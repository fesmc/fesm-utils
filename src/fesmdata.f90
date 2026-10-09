module fesmdata
    ! Helpers for reading FesmData v2 files ($ICE_DATA/v2/<Domain>/<GRID>/),
    ! shared by the regions and topodata modules:
    !   - flag tables: the codes and names of a categorical variable, from its
    !     CF attributes flag_values and flag_meanings
    !   - the grid of a file: its global attribute grid_name, and the grid
    !     itself from the cdo description grid_<GRID>.txt next to the file
    !     (as FesmData writes it), else in maps/
    !   - reading a 2D field, remapped onto a target grid if needed
    !   - {domain}, {grid_name} and {key} placeholders in paths

    use, intrinsic :: iso_fortran_env, only : error_unit

    use precision
    use constants,   only : mv
    use ncio
    use nml,         only : nml_replace
    use coordinates, only : grid_class
    use grid_cdo,    only : grid_cdo_read_desc
    use mapping,     only : map_class, map_field

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
    public :: flag_table_class, flag_table_read, flag_table_merge
    public :: flag_codes, flag_name, flag_names_joined
    public :: fesmdata_grid_name, fesmdata_grid_read
    public :: fesmdata_read_field
    public :: fesmdata_parse_path

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
            call fesmdata_error("flag_table_read", "variable without flag_values.", &
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
            call fesmdata_error("flag_table_read", "flag_values and flag_meanings differ in length.", &
                "file = "//trim(fname)//new_line("a")//"variable = "//trim(varname))
        end if

        return

    end subroutine flag_table_read

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

    ! ===== Grid of a file =====================================================

    function fesmdata_grid_name(fname) result(grid_name)
        ! The grid of a FesmData file (its global attribute grid_name).

        implicit none

        character(len=*), intent(IN) :: fname
        character(len=256) :: grid_name

        if (.not. nc_exists_attr(fname, "grid_name")) then
            call fesmdata_error("fesmdata_grid_name", "file without the global attribute grid_name.", &
                "file = "//trim(fname))
        end if
        call nc_read_attr(fname, "grid_name", grid_name)

    end function fesmdata_grid_name

    subroutine fesmdata_grid_read(grid, fname)
        ! The grid of a FesmData file, from its cdo description
        ! grid_<grid_name>.txt next to the file, else in maps/.

        implicit none

        type(grid_class), intent(OUT) :: grid
        character(len=*), intent(IN)  :: fname

        character(len=256)  :: grid_name
        character(len=1024) :: fldr
        logical :: found
        integer :: q

        grid_name = fesmdata_grid_name(fname)

        q = index(fname, "/", back=.true.)
        fldr = "."
        if (q .gt. 0) fldr = fname(1:q-1)

        inquire(file=trim(fldr)//"/grid_"//trim(grid_name)//".txt", exist=found)
        if (.not. found) then
            fldr = "maps"
            inquire(file=trim(fldr)//"/grid_"//trim(grid_name)//".txt", exist=found)
        end if
        if (.not. found) then
            call fesmdata_error("fesmdata_grid_read", "no description of the grid of the file.", &
                "file      = "//trim(fname)//new_line("a")// &
                "grid_name = "//trim(grid_name)//new_line("a")// &
                "looked for grid_"//trim(grid_name)//".txt next to the file and in maps/")
        end if

        call grid_cdo_read_desc(grid, trim(grid_name), trim(fldr))

        return

    end subroutine fesmdata_grid_read

    ! ===== Reading fields =====================================================

    subroutine fesmdata_read_field_int(fname, varname, var, nx, ny, remap, map, none)
        ! Read an integer 2D field into var(nx,ny). With remap, the field is
        ! mapped from the grid of the file with map, and target cells without
        ! a source value get none.

        implicit none

        character(len=*),     intent(IN)  :: fname, varname
        integer, allocatable, intent(OUT) :: var(:,:)
        integer,              intent(IN)  :: nx, ny
        logical,              intent(IN)  :: remap
        type(map_class),      intent(IN)  :: map
        integer,              intent(IN)  :: none

        integer, allocatable :: src(:,:)
        logical, allocatable :: filled(:,:)

        allocate(var(nx,ny))

        if (remap) then
            allocate(src(nc_size(fname,"xc"),nc_size(fname,"yc")))
            allocate(filled(nx,ny))
            call nc_read(fname, varname, src)
            call map_field(map, varname, src, var, mask2=filled, reset=.true.)
            where (.not. filled) var = none
        else
            call nc_read(fname, varname, var)
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

        allocate(var(nx,ny))

        if (remap) then
            allocate(src(nc_size(fname,"xc"),nc_size(fname,"yc")))
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

    subroutine fesmdata_error(proc, msg, detail)
        ! Abort with a framed message on error_unit.

        implicit none

        character(len=*), intent(IN)           :: proc
        character(len=*), intent(IN)           :: msg
        character(len=*), intent(IN), optional :: detail

        integer :: p0, p1

        write(error_unit,"(a)") ""
        write(error_unit,"(a)") "fesmdata:: error in "//trim(proc)
        write(error_unit,"(a)") "    "//trim(msg)

        if (present(detail)) then
            p0 = 1
            do
                p1 = index(detail(p0:), new_line("a"))
                if (p1 == 0) then
                    write(error_unit,"(a)") "    "//trim(detail(p0:))
                    exit
                end if
                write(error_unit,"(a)") "    "//trim(detail(p0:p0+p1-2))
                p0 = p0 + p1
            end do
        end if

        write(error_unit,"(a)") "  stopped by fesmdata."
        error stop 1

    end subroutine fesmdata_error

end module fesmdata
