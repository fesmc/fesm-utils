module fesmdata
    ! Helpers for reading FesmData v2 files (fesmdata/<Domain>/<GRID>/),
    ! shared by the regions and topodata modules and the model hosts:
    !   - flag tables: the codes and names of a categorical variable, from its
    !     CF attributes flag_values and flag_meanings
    !   - the grid of a file: its global attribute grid_name, and the grid
    !     itself from the cdo description grid_<GRID>.txt next to the file
    !     (as FesmData writes it), else in maps/
    !   - placing a field of a file on a target grid: remapped from the grid
    !     of the file if its name differs, else read as it is (a file
    !     without grid_name must then be on the target grid)
    !   - reading a 2D field, remapped onto a target grid if needed
    !   - {domain}, {grid_name} and {key} placeholders in paths
    !   - units of a variable (its units attribute): the factor of a water
    !     flux to other units, the reference temperature of a temperature

    use fesm_log, only : fesm_error
    use fesm_units, only : units_conv

    use precision
    use constants,   only : mv
    use ncio
    use nml,         only : nml_replace
    use coordinates, only : grid_class
    use grid_cdo,    only : grid_cdo_read_desc
    use mapping,     only : map_class, map_field, map_init

    implicit none

    integer, parameter :: len_name = 128

    type flag_table_class
        integer,                 allocatable :: codes(:)
        character(len=len_name), allocatable :: names(:)
    end type

    interface fesmdata_read_field
        module procedure fesmdata_read_field_int, fesmdata_read_field_wp
    end interface

    private
    public :: len_name
    public :: flag_table_class, flag_table_read, flag_table_merge, flag_table_from_values
    public :: flag_codes, flag_name, flag_names_joined
    public :: fesmdata_grid_name, fesmdata_grid_read, fesmdata_grid_desc_exists
    public :: fesmdata_map_init, fesmdata_read_field
    public :: fesmdata_parse_path
    public :: fesmdata_flux_conv, fesmdata_temp_ref

contains

    ! ===== Flag tables ========================================================

    subroutine flag_table_read(tab, fname, varname)
        ! Codes and names of a variable from its flag_values and flag_meanings.

        implicit none

        type(flag_table_class), intent(OUT) :: tab
        character(len=*),       intent(IN)  :: fname, varname

        character(len=:), allocatable :: meanings
        integer :: n, nm, k, i0, i1

        if (.not. nc_exists_attr(fname, varname, "flag_values")) then
            call fesm_error("flag_table_read", "variable without flag_values.", &
                "file = "//trim(fname)//new_line("a")//"variable = "//trim(varname))
        end if

        n = nc_size_attr(fname, varname, "flag_values")
        allocate(tab%codes(n), tab%names(n))
        call nc_read_attr(fname, varname, "flag_values", tab%codes)

        nm = nc_size_attr(fname, varname, "flag_meanings")
        allocate(character(len=nm) :: meanings)
        call nc_read_attr(fname, varname, "flag_meanings", meanings)

        ! Space-separated names, one per code
        i0 = 1
        do k = 1, n
            do while (i0 .le. nm)
                if (meanings(i0:i0) .ne. " ") exit
                i0 = i0 + 1
            end do
            if (i0 .gt. nm) exit
            i1 = index(meanings(i0:), " ")
            if (i1 .eq. 0) then
                i1 = nm
            else
                i1 = i0 + i1 - 2
            end if
            tab%names(k) = meanings(i0:i1)
            i0 = i1 + 2
        end do

        if (k .le. n .or. len_trim(meanings(min(i0,nm+1):)) .gt. 0) then
            call fesm_error("flag_table_read", "flag_values and flag_meanings differ in length.", &
                "file = "//trim(fname)//new_line("a")//"variable = "//trim(varname))
        end if

        return

    end subroutine flag_table_read

    subroutine flag_table_from_values(tab, field)
        ! A table of the positive values of a field, named by their values
        ! (for fields without flag attributes).

        implicit none

        type(flag_table_class), intent(OUT) :: tab
        integer,                intent(IN)  :: field(:,:)

        integer :: c

        allocate(tab%codes(0), tab%names(0))
        if (.not. any(field .gt. 0)) return

        c = minval(field, mask=field .gt. 0)
        do
            tab%codes = [tab%codes, c]
            tab%names = [tab%names, int_to_name(c)]
            if (.not. any(field .gt. c)) exit
            c = minval(field, mask=field .gt. c)
        end do

        return

    end subroutine flag_table_from_values

    pure function int_to_name(c) result(name)
        integer, intent(IN) :: c
        character(len=len_name) :: name
        write(name, "(i0)") c
    end function int_to_name

    subroutine flag_table_merge(tab, tab_add)
        ! Add the entries of tab_add with codes not yet in tab, keeping tab
        ! sorted by code.

        implicit none

        type(flag_table_class), intent(INOUT) :: tab
        type(flag_table_class), intent(IN)    :: tab_add

        integer,                 allocatable :: codes(:)
        character(len=len_name), allocatable :: names(:)
        logical, allocatable :: new(:)
        integer :: k, j, n

        allocate(new(size(tab_add%codes)))
        do k = 1, size(tab_add%codes)
            new(k) = .not. any(tab%codes .eq. tab_add%codes(k))
        end do

        codes = [tab%codes, pack(tab_add%codes, new)]
        names = [tab%names, pack(tab_add%names, new)]

        ! Insertion sort by code
        n = size(codes)
        do k = 2, n
            j = k
            do while (j .gt. 1)
                if (codes(j-1) .le. codes(j)) exit
                codes(j-1:j) = codes([j,j-1])
                names(j-1:j) = names([j,j-1])
                j = j - 1
            end do
        end do

        call move_alloc(codes, tab%codes)
        call move_alloc(names, tab%names)

        return

    end subroutine flag_table_merge

    function flag_codes(tab, name) result(codes)
        ! All codes with the name (ignoring case; spaces read as underscores).
        ! A name may have several codes, e.g. "Africa" in both hemispheres.

        implicit none

        type(flag_table_class), intent(IN) :: tab
        character(len=*),       intent(IN) :: name
        integer, allocatable :: codes(:)

        character(len=len_name) :: key
        logical, allocatable :: match(:)
        integer :: k

        key = normalize_name(name)

        allocate(match(size(tab%codes)))
        do k = 1, size(tab%codes)
            match(k) = (normalize_name(tab%names(k)) .eq. key)
        end do
        codes = pack(tab%codes, match)

        return

    end function flag_codes

    function flag_name(tab, code) result(name)
        ! Name of a code ("" if the code is not in the table).

        implicit none

        type(flag_table_class), intent(IN) :: tab
        integer,                intent(IN) :: code
        character(len=len_name) :: name

        integer :: k

        name = ""
        do k = 1, size(tab%codes)
            if (tab%codes(k) .eq. code) then
                name = tab%names(k)
                exit
            end if
        end do

        return

    end function flag_name

    function flag_names_joined(tab) result(str)
        ! All names, space-separated (for messages).

        implicit none

        type(flag_table_class), intent(IN) :: tab
        character(len=:), allocatable :: str

        integer :: k

        str = ""
        do k = 1, size(tab%names)
            str = str//trim(tab%names(k))//" "
        end do

    end function flag_names_joined

    pure function normalize_name(name) result(key)
        ! Lower case, spaces between words as underscores.

        implicit none

        character(len=*), intent(IN) :: name
        character(len=len_name) :: key

        integer :: k, c

        key = adjustl(name)
        do k = 1, len_trim(key)
            c = ichar(key(k:k))
            if (c .ge. ichar("A") .and. c .le. ichar("Z")) key(k:k) = char(c + 32)
            if (key(k:k) .eq. " ") key(k:k) = "_"
        end do

    end function normalize_name

    ! ===== Units of a variable ================================================

    function fesmdata_flux_conv(fname, varname, units_out) result(fac)
        ! Factor converting a water flux from the units of the variable (its
        ! units attribute) to units_out: var_out = fac*var. The mass is
        ! kg m-2 or mm (water equivalent, the same), the time s, d or yr
        ! (1 yr = 365.25 d, as in FesmData), e.g. "kg m-2 yr-1", "mm d-1",
        ! "mm w.e. a-1", "mm/d". Stops on missing or unknown units. The
        ! factor is that of units_conv (fesm_units).

        implicit none

        character(len=*), intent(IN) :: fname
        character(len=*), intent(IN) :: varname
        character(len=*), intent(IN) :: units_out
        real(wp) :: fac

        character(len=256) :: units
        real(dp) :: scale, offset, s2, o2
        integer  :: stat, stat2

        units = var_units(fname, varname, "fesmdata_flux_conv")

        ! A water flux on both sides (the canonical class of units_conv)
        call units_conv(units, units_out, scale, offset, stat=stat)
        call units_conv(units_out, "kg m-2 yr-1", s2, o2, stat=stat2)
        if (stat .ne. 0 .or. stat2 .ne. 0) then
            call fesm_error("fesmdata_flux_conv", "water-flux units not recognized.", &
                "file      = "//trim(fname)//new_line("a")// &
                "varname   = "//trim(varname)//new_line("a")// &
                "units     = "//trim(units)//new_line("a")// &
                "units_out = "//trim(units_out))
        end if

        fac = real(scale, wp)

    end function fesmdata_flux_conv

    function fesmdata_temp_ref(fname, varname) result(t_ref)
        ! Reference temperature of the units of a temperature variable (its
        ! units attribute), to add to get kelvin: 0 for K, 273.15 for degC.
        ! Stops on missing or unknown units. The offset of units_conv to K.

        implicit none

        character(len=*), intent(IN) :: fname
        character(len=*), intent(IN) :: varname
        real(wp) :: t_ref

        character(len=256) :: units
        real(dp) :: scale, offset
        integer  :: stat

        units = var_units(fname, varname, "fesmdata_temp_ref")

        call units_conv(units, "K", scale, offset, stat=stat)
        if (stat .ne. 0 .or. scale .ne. 1.0_dp) then
            call fesm_error("fesmdata_temp_ref", "temperature units not recognized.", &
                "file    = "//trim(fname)//new_line("a")// &
                "varname = "//trim(varname)//new_line("a")// &
                "units   = "//trim(units))
        end if

        t_ref = real(offset, wp)

    end function fesmdata_temp_ref

    function var_units(fname, varname, proc) result(units)
        ! The units attribute of a variable; stops without one.

        implicit none

        character(len=*), intent(IN) :: fname, varname, proc
        character(len=256) :: units

        if (.not. nc_exists_attr(fname, varname, "units")) then
            call fesm_error(proc, "variable without a units attribute.", &
                "file    = "//trim(fname)//new_line("a")// &
                "varname = "//trim(varname))
        end if
        call nc_read_attr(fname, varname, "units", units)

    end function var_units

    ! ===== Grid of a file =====================================================

    function fesmdata_grid_name(fname) result(grid_name)
        ! The grid of a FesmData file (its global attribute grid_name).

        implicit none

        character(len=*), intent(IN) :: fname
        character(len=256) :: grid_name

        if (.not. nc_exists_attr(fname, "grid_name")) then
            call fesm_error("fesmdata_grid_name", "file without the global attribute grid_name.", &
                "file = "//trim(fname))
        end if
        call nc_read_attr(fname, "grid_name", grid_name)

    end function fesmdata_grid_name

    subroutine fesmdata_grid_read(grid, fname, grid_name)
        ! The grid of a file (grid_name, by default its global attribute
        ! grid_name), from its cdo description grid_<grid_name>.txt next to
        ! the file, else in maps/.

        implicit none

        type(grid_class), intent(OUT) :: grid
        character(len=*), intent(IN)  :: fname
        character(len=*), intent(IN), optional :: grid_name

        character(len=256)  :: gname
        character(len=1024) :: fldr
        logical :: found

        if (present(grid_name)) then
            gname = grid_name
        else
            gname = fesmdata_grid_name(fname)
        end if

        found = grid_desc_folder(fname, gname, fldr)
        if (.not. found) then
            call fesm_error("fesmdata_grid_read", "no description of the grid of the file.", &
                "file      = "//trim(fname)//new_line("a")// &
                "grid_name = "//trim(gname)//new_line("a")// &
                "looked for grid_"//trim(gname)//".txt next to the file and in maps/")
        end if

        call grid_cdo_read_desc(grid, trim(gname), trim(fldr))

        return

    end subroutine fesmdata_grid_read

    function fesmdata_grid_desc_exists(fname, grid_name) result(found)
        ! A cdo description of the grid grid_name exists next to the file
        ! (grid_<grid_name>.txt, as FesmData writes it), else in maps/.

        implicit none

        character(len=*), intent(IN) :: fname, grid_name
        logical :: found

        character(len=1024) :: fldr

        found = grid_desc_folder(fname, grid_name, fldr)

    end function fesmdata_grid_desc_exists

    function grid_desc_folder(fname, grid_name, fldr) result(found)
        ! The folder of the cdo description grid_<grid_name>.txt: that of the
        ! file, else maps/.

        implicit none

        character(len=*), intent(IN)  :: fname, grid_name
        character(len=*), intent(OUT) :: fldr
        logical :: found

        integer :: q

        q = index(fname, "/", back=.true.)
        fldr = "."
        if (q .gt. 0) fldr = fname(1:q-1)

        inquire(file=trim(fldr)//"/grid_"//trim(grid_name)//".txt", exist=found)
        if (.not. found) then
            fldr = "maps"
            inquire(file=trim(fldr)//"/grid_"//trim(grid_name)//".txt", exist=found)
        end if

    end function grid_desc_folder

    ! ===== Fields on a target grid ============================================

    subroutine fesmdata_map_init(map, remap, nx, ny, fname, varname, method, grid, grid_name, grid_src)
        ! How a field of a file gets onto the target grid: with grid, and a
        ! grid of the file (grid_name, else its global attribute grid_name)
        ! of another name, the field is remapped (remap, with map by method);
        ! else it is read as it is, and must have the size of grid. nx, ny:
        ! the size of the field as read. grid_src: the grid of the file ("":
        ! unknown).

        implicit none

        type(map_class),  intent(INOUT) :: map
        logical,          intent(OUT)   :: remap
        integer,          intent(OUT)   :: nx, ny
        character(len=*), intent(IN)    :: fname, varname, method
        type(grid_class), intent(IN), optional  :: grid
        character(len=*), intent(IN), optional  :: grid_name
        character(len=*), intent(OUT), optional :: grid_src

        character(len=256) :: gsrc
        type(grid_class) :: grid_file
        integer, allocatable :: dims(:)

        gsrc = ""
        if (present(grid_name)) gsrc = grid_name
        if (len_trim(gsrc) .eq. 0 .and. nc_exists_attr(fname, "grid_name")) gsrc = fesmdata_grid_name(fname)
        if (present(grid_src)) grid_src = gsrc

        remap = .false.
        if (present(grid)) remap = (len_trim(gsrc) .gt. 0 .and. trim(gsrc) .ne. trim(grid%name))

        if (remap) then
            call fesmdata_grid_read(grid_file, fname, gsrc)
            call map_init(map, grid_file, grid, method=trim(method), fldr="maps")
            nx = grid%G%nx
            ny = grid%G%ny
        else
            call nc_dims(fname, varname, dims=dims)
            nx = dims(1)
            ny = dims(2)
            if (present(grid)) then
                if (nx .ne. grid%G%nx .or. ny .ne. grid%G%ny) then
                    if (len_trim(gsrc) .eq. 0) gsrc = "(none: the file has no grid_name)"
                    call fesm_error("fesmdata_map_init", "the field is not on the target grid.", &
                        "file        = "//trim(fname)//new_line("a")// &
                        "variable    = "//trim(varname)//new_line("a")// &
                        "grid        = "//trim(gsrc)//new_line("a")// &
                        "target grid = "//trim(grid%name))
                end if
            end if
        end if

        return

    end subroutine fesmdata_map_init

    ! ===== Reading fields =====================================================

    subroutine fesmdata_read_field_int(fname, varname, var, nx, ny, remap, map, none)
        ! Read an integer 2D field into var(nx,ny). With remap, the field is
        ! mapped from the grid of the file with map, and target cells without
        ! a source value get none, as do negative values (e.g. fill values):
        ! these are a class of their own, not missing values that the map
        ! would fill from neighbouring cells.

        implicit none

        character(len=*),     intent(IN)  :: fname, varname
        integer, allocatable, intent(OUT) :: var(:,:)
        integer,              intent(IN)  :: nx, ny
        logical,              intent(IN)  :: remap
        type(map_class),      intent(IN)  :: map
        integer,              intent(IN)  :: none

        integer, allocatable :: src(:,:), dims(:)
        logical, allocatable :: filled(:,:)

        allocate(var(nx,ny))

        if (remap) then
            call nc_dims(fname, varname, dims=dims)
            allocate(src(dims(1),dims(2)))
            allocate(filled(nx,ny))
            call nc_read(fname, varname, src)
            where (src .lt. 0) src = none
            ! (none-1 does not occur: no source value is missing)
            call map_field(map, varname, src, var, missing_value=none-1, mask2=filled, reset=.true.)
            where (.not. filled) var = none
        else
            call nc_read(fname, varname, var)
            where (var .lt. 0) var = none
        end if

        return

    end subroutine fesmdata_read_field_int

    subroutine fesmdata_read_field_wp(fname, varname, var, nx, ny, remap, map)
        ! Read a real 2D field into var(nx,ny), missing values as mv. With
        ! remap, the field is mapped from the grid of the file with map.

        implicit none

        character(len=*),      intent(IN)  :: fname, varname
        real(wp), allocatable, intent(OUT) :: var(:,:)
        integer,               intent(IN)  :: nx, ny
        logical,               intent(IN)  :: remap
        type(map_class),       intent(IN)  :: map

        real(wp), allocatable :: src(:,:)
        integer,  allocatable :: dims(:)

        allocate(var(nx,ny))

        if (remap) then
            call nc_dims(fname, varname, dims=dims)
            allocate(src(dims(1),dims(2)))
            call nc_read(fname, varname, src, missing_value=mv)
            call map_field(map, varname, src, var, missing_value=mv, reset=.true.)
        else
            call nc_read(fname, varname, var, missing_value=mv)
        end if

        return

    end subroutine fesmdata_read_field_wp

    ! ===== Paths ==============================================================

    subroutine fesmdata_parse_path(path, domain, grid_name, subs)
        ! Substitute {domain}, {grid_name} and any extra {key}->value pairs
        ! (subs(k,1) = key without braces, subs(k,2) = value).

        implicit none

        character(len=*), intent(INOUT) :: path
        character(len=*), intent(IN), optional :: domain, grid_name
        character(len=*), intent(IN), optional :: subs(:,:)

        integer :: k

        if (present(domain))    call nml_replace(path, "{domain}",    trim(domain))
        if (present(grid_name)) call nml_replace(path, "{grid_name}", trim(grid_name))

        if (present(subs)) then
            do k = 1, size(subs,1)
                call nml_replace(path, "{"//trim(subs(k,1))//"}", trim(subs(k,2)))
            end do
        end if

        return

    end subroutine fesmdata_parse_path

end module fesmdata
