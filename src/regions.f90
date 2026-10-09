module regions
    ! Regions, zones and basins of a FesmData v2 grid, and masks derived from
    ! them by selection expressions.
    !
    ! Files (FesmData/Regions, $ICE_DATA/v2/<Domain>/<GRID>/):
    !   <GRID>_REGIONS.nc       region_1, region_2, region_3, zone, dist_shelfbreak
    !   <GRID>_BASINS-<set>.nc  basin, basin_group (optional), basin_mask
    ! The names of the codes come from the CF attributes flag_values and
    ! flag_meanings of each variable, which list only the codes present on
    ! the grid.
    !
    ! Region codes have two decimal digits per level below the first: the
    ! path "1.3.1" is the code 10301. region_n holds the code of the deepest
    ! region at or above level n, so region_3 places every cell at all levels
    ! and in_region(region_3, code) selects a region with its subregions.
    !
    ! Selection expressions (regions_select, and the named masks of the
    ! namelist) combine terms field:value,value,... where a comma is an OR:
    !   region:Greenland,1.5        names, paths or codes, at any level
    !   zone:land,continental_shelf names or codes of the zone
    !   <set>:11,12                 basins of a loaded basin set (names or ids)
    !   <set>.group:1               basin groups of the set
    !   <set>.mask:1                original extent of the basins of the set
    ! A term prefixed by ~ is negated, & joins terms (AND) and | joins
    ! clauses (OR, binding weaker than &). "all" and "none" select every cell
    ! and no cell. Names are matched ignoring case, with spaces read as
    ! underscores. Example: "region:Greenland & ~zone:open_ocean | Zwally2012:11"
    !
    ! With a target grid whose name differs from the grid of the files, all
    ! fields are remapped by nearest neighbour onto it. The grid of the files
    ! (their global attribute grid_name) is read from its cdo description
    ! grid_<name>.txt next to the files, else in maps/.

    use, intrinsic :: iso_fortran_env, only : error_unit

    use precision
    use constants,   only : mv
    use ncio
    use nml
    use coordinates, only : grid_class
    use grid_cdo,    only : grid_cdo_read_desc
    use mapping,     only : map_class, map_init, map_field

    implicit none

    integer, parameter :: len_name = 128
    integer, parameter :: len_expr = 1000
    integer, parameter :: n_sets_max  = 20
    integer, parameter :: n_masks_max = 50

    ! Zone of a remapped target cell without a source value
    integer, parameter :: zone_undefined = -1

    type flag_table_class
        integer,                 allocatable :: codes(:)
        character(len=len_name), allocatable :: names(:)
    end type

    type basin_set_class
        character(len=56)    :: name                ! e.g. "Zwally2012"
        character(len=1024)  :: filename
        logical              :: with_group
        integer, allocatable :: basin(:,:)
        integer, allocatable :: basin_group(:,:)
        integer, allocatable :: basin_mask(:,:)
        type(flag_table_class) :: tab_basin, tab_group
    end type

    type region_mask_class
        character(len=56)       :: name
        character(len=len_expr) :: expr
        logical, allocatable    :: mask(:,:)
    end type

    type regions_param_class
        character(len=1024) :: path_regions         ! <GRID>_REGIONS.nc
        character(len=1024) :: path_basins          ! template with {set}
        character(len=56),       allocatable :: basin_sets(:)
        character(len=56),       allocatable :: mask_names(:)
        character(len=len_expr), allocatable :: mask_exprs(:)

        ! Internal parameters
        character(len=256) :: grid_src              ! grid of the files
        character(len=256) :: grid_tgt              ! grid of the fields
    end type

    type regions_class
        type(regions_param_class) :: par
        integer :: nx, ny

        integer,  allocatable :: region_1(:,:), region_2(:,:), region_3(:,:)
        integer,  allocatable :: zone(:,:)
        real(wp), allocatable :: dist_shelfbreak(:,:)

        ! Region names: the tables of region_1..3 together
        type(flag_table_class) :: tab_region, tab_zone

        type(basin_set_class),   allocatable :: basins(:)
        type(region_mask_class), allocatable :: masks(:)

        ! Online remapping from the grid of the files
        logical         :: remap = .false.
        type(map_class) :: map
    end type

    private
    public :: flag_table_class, basin_set_class, region_mask_class
    public :: regions_param_class, regions_class
    public :: regions_init_nml, regions_init_arg, regions_end
    public :: regions_select, regions_mask
    public :: flag_codes, flag_name
    public :: region_code, region_level, region_ancestor, region_path, in_region
    public :: zone_undefined

contains

    ! ===== Initialization =====================================================

    subroutine regions_init_nml(reg, filename, group, domain, grid_name, subs, grid, verbose)
        ! Load regions, basin sets and named masks as given by a namelist group:
        !   path_regions = "ice_data/v2/{domain}/{grid_name}/{grid_name}_REGIONS.nc"
        !   path_basins  = "ice_data/v2/{domain}/{grid_name}/{grid_name}_BASINS-{set}.nc"
        !   basin_sets   = "Zwally2012"          (optional; none by default)
        !   masks        = "grl_shelf"           (optional; none by default)
        !   mask_grl_shelf = "region:Greenland & zone:continental_shelf"
        ! {grid_name} is the grid of the files. grid: the target grid, onto
        ! which the fields are remapped if it differs from the grid of the files.

        implicit none

        type(regions_class), intent(INOUT) :: reg
        character(len=*),    intent(IN)    :: filename
        character(len=*),    intent(IN)    :: group
        character(len=*),    intent(IN), optional :: domain
        character(len=*),    intent(IN), optional :: grid_name
        character(len=*),    intent(IN), optional :: subs(:,:)   ! extra {key}->value path substitutions
        type(grid_class),    intent(IN), optional :: grid
        logical,             intent(IN), optional :: verbose

        call regions_par_load(reg%par, filename, group, domain, grid_name, subs, verbose)

        call regions_init_data(reg, grid)

        return

    end subroutine regions_init_nml

    subroutine regions_init_arg(reg, path_regions, path_basins, basin_sets, mask_names, mask_exprs, grid)
        ! Load regions, basin sets and named masks from explicit arguments
        ! (see regions_init_nml). path_basins is a template with {set}.

        implicit none

        type(regions_class), intent(INOUT) :: reg
        character(len=*),    intent(IN)    :: path_regions
        character(len=*),    intent(IN), optional :: path_basins
        character(len=*),    intent(IN), optional :: basin_sets(:)
        character(len=*),    intent(IN), optional :: mask_names(:)
        character(len=*),    intent(IN), optional :: mask_exprs(:)
        type(grid_class),    intent(IN), optional :: grid

        integer :: n_sets, n_masks

        n_sets = 0
        if (present(basin_sets)) n_sets = size(basin_sets)
        n_masks = 0
        if (present(mask_names)) n_masks = size(mask_names)

        if (n_sets .gt. 0 .and. .not. present(path_basins)) then
            call regions_error("regions_init_arg", "basin_sets are given, but no path_basins.")
        end if
        if (n_masks .gt. 0) then
            if (.not. present(mask_exprs)) then
                call regions_error("regions_init_arg", "mask_names are given, but no mask_exprs.")
            end if
            if (size(mask_exprs) .ne. n_masks) then
                call regions_error("regions_init_arg", "mask_names and mask_exprs differ in length.")
            end if
        end if

        reg%par%path_regions = trim(path_regions)
        reg%par%path_basins  = ""
        if (present(path_basins)) reg%par%path_basins = trim(path_basins)

        allocate(reg%par%basin_sets(n_sets))
        if (n_sets .gt. 0) reg%par%basin_sets = basin_sets

        allocate(reg%par%mask_names(n_masks), reg%par%mask_exprs(n_masks))
        if (n_masks .gt. 0) then
            reg%par%mask_names = mask_names
            reg%par%mask_exprs = mask_exprs
        end if

        call regions_init_data(reg, grid)

        return

    end subroutine regions_init_arg

    subroutine regions_init_data(reg, grid)
        ! Read the files given by reg%par (remapping onto grid if needed) and
        ! evaluate the named masks.

        implicit none

        type(regions_class), intent(INOUT) :: reg
        type(grid_class),    intent(IN), optional :: grid

        character(len=1024) :: fname
        type(flag_table_class) :: tab
        integer :: k

        fname = reg%par%path_regions

        if (.not. nc_exists_var(fname, "region_1")) then
            call regions_error("regions_init_data", &
                "not a FesmData v2 regions file (no variable region_1).", &
                "file = "//trim(fname)//new_line("a")// &
                "FesmData v1 regions files (one variable mask) are not supported.")
        end if

        ! Grid of the files, and the map onto the target grid if it differs
        call nc_read_attr(fname, "grid_name", reg%par%grid_src)
        reg%par%grid_tgt = reg%par%grid_src
        reg%remap = .false.
        if (present(grid)) then
            reg%par%grid_tgt = grid%name
            if (trim(grid%name) .ne. trim(reg%par%grid_src)) then
                call regions_remap_init(reg, fname, grid)
            end if
        end if

        if (reg%remap) then
            reg%nx = grid%G%nx
            reg%ny = grid%G%ny
        else
            reg%nx = nc_size(fname, "xc")
            reg%ny = nc_size(fname, "yc")
            if (present(grid)) then
                if (grid%G%nx .ne. reg%nx .or. grid%G%ny .ne. reg%ny) then
                    call regions_error("regions_init_data", &
                        "the target grid has the name of the grid of the file, but another size.", &
                        "file = "//trim(fname)//new_line("a")//"grid = "//trim(grid%name))
                end if
            end if
        end if

        ! Regions and zones
        call read_int(reg, fname, "region_1", 0, reg%region_1)
        call read_int(reg, fname, "region_2", 0, reg%region_2)
        call read_int(reg, fname, "region_3", 0, reg%region_3)
        call read_int(reg, fname, "zone", zone_undefined, reg%zone)
        call read_real(reg, fname, "dist_shelfbreak", reg%dist_shelfbreak)

        call flag_table_read(reg%tab_region, fname, "region_1")
        do k = 2, 3
            call flag_table_read(tab, fname, "region_"//char(ichar("0")+k))
            call flag_table_merge(reg%tab_region, tab)
        end do
        call flag_table_read(reg%tab_zone, fname, "zone")

        ! Basin sets
        if (allocated(reg%basins)) deallocate(reg%basins)
        allocate(reg%basins(size(reg%par%basin_sets)))
        do k = 1, size(reg%basins)
            call basin_set_load(reg, reg%basins(k), reg%par%basin_sets(k))
        end do

        ! Named masks
        if (allocated(reg%masks)) deallocate(reg%masks)
        allocate(reg%masks(size(reg%par%mask_names)))
        do k = 1, size(reg%masks)
            reg%masks(k)%name = reg%par%mask_names(k)
            reg%masks(k)%expr = reg%par%mask_exprs(k)
            allocate(reg%masks(k)%mask(reg%nx,reg%ny))
            reg%masks(k)%mask = regions_select(reg, reg%masks(k)%expr)
        end do

        return

    end subroutine regions_init_data

    subroutine basin_set_load(reg, bs, name)

        implicit none

        type(regions_class),   intent(INOUT) :: reg
        type(basin_set_class), intent(INOUT) :: bs
        character(len=*),      intent(IN)    :: name

        character(len=256) :: grid_name

        bs%name     = trim(name)
        bs%filename = reg%par%path_basins
        call nml_replace(bs%filename, "{set}", trim(name))

        ! All files must be on the same grid (one map)
        call nc_read_attr(bs%filename, "grid_name", grid_name)
        if (trim(grid_name) .ne. trim(reg%par%grid_src)) then
            call regions_error("basin_set_load", &
                "the basins file is on another grid than the regions file.", &
                "file      = "//trim(bs%filename)//new_line("a")// &
                "grid_name = "//trim(grid_name)//new_line("a")// &
                "regions   = "//trim(reg%par%grid_src))
        end if

        call read_int(reg, bs%filename, "basin",      0, bs%basin)
        call read_int(reg, bs%filename, "basin_mask", 0, bs%basin_mask)
        call flag_table_read(bs%tab_basin, bs%filename, "basin")

        bs%with_group = nc_exists_var(bs%filename, "basin_group")
        if (bs%with_group) then
            call read_int(reg, bs%filename, "basin_group", 0, bs%basin_group)
            call flag_table_read(bs%tab_group, bs%filename, "basin_group")
        end if

        return

    end subroutine basin_set_load

    subroutine regions_end(reg)

        implicit none

        type(regions_class), intent(INOUT) :: reg

        type(regions_class) :: reg0

        reg = reg0

        return

    end subroutine regions_end

    subroutine regions_par_load(par, filename, group, domain, grid_name, subs, verbose)

        implicit none

        type(regions_param_class), intent(OUT) :: par
        character(len=*), intent(IN) :: filename
        character(len=*), intent(IN) :: group
        character(len=*), intent(IN), optional :: domain
        character(len=*), intent(IN), optional :: grid_name
        character(len=*), intent(IN), optional :: subs(:,:)
        logical,          intent(IN), optional :: verbose

        character(len=56) :: sets(n_sets_max), names(n_masks_max)
        logical :: print_summary
        integer :: k, n

        print_summary = .true.
        if (present(verbose)) print_summary = verbose

        call nml_read(filename, group, "path_regions", par%path_regions)
        call parse_path(par%path_regions, domain, grid_name, subs)

        par%path_basins = ""
        sets = ""
        if (nml_has_param(filename, group, "basin_sets")) then
            call nml_read(filename, group, "basin_sets",  sets)
            call nml_read(filename, group, "path_basins", par%path_basins)
            call parse_path(par%path_basins, domain, grid_name, subs)
        end if
        n = count(len_trim(sets) .gt. 0)
        allocate(par%basin_sets(n))
        par%basin_sets = pack(sets, len_trim(sets) .gt. 0)

        names = ""
        if (nml_has_param(filename, group, "masks")) then
            call nml_read(filename, group, "masks", names)
        end if
        n = count(len_trim(names) .gt. 0)
        allocate(par%mask_names(n), par%mask_exprs(n))
        par%mask_names = pack(names, len_trim(names) .gt. 0)
        do k = 1, n
            call nml_read(filename, group, "mask_"//trim(par%mask_names(k)), par%mask_exprs(k))
        end do

        if (print_summary) then
            write(*,*) "Loading: ", trim(filename), ":: ", trim(group)
            write(*,*) "path_regions  = ", trim(par%path_regions)
            if (size(par%basin_sets) .gt. 0) then
                write(*,*) "path_basins   = ", trim(par%path_basins)
                write(*,*) "basin_sets    = ", (trim(par%basin_sets(k))//" ", k=1,size(par%basin_sets))
            end if
            do k = 1, size(par%mask_names)
                write(*,*) "mask_"//trim(par%mask_names(k))//" = ", trim(par%mask_exprs(k))
            end do
        end if

        return

    end subroutine regions_par_load

    subroutine parse_path(path, domain, grid_name, subs)
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

    end subroutine parse_path

    ! ===== Reading and remapping ==============================================

    subroutine regions_remap_init(reg, fname, grid)
        ! Nearest-neighbour map from the grid of the files onto grid. The grid
        ! of the files is read from grid_<name>.txt next to the file, else in
        ! maps/; the map is cached in maps/.

        implicit none

        type(regions_class), intent(INOUT) :: reg
        character(len=*),    intent(IN)    :: fname
        type(grid_class),    intent(IN)    :: grid

        type(grid_class)    :: grid_src
        character(len=1024) :: fldr
        logical :: found
        integer :: q

        q = index(fname, "/", back=.true.)
        fldr = "."
        if (q .gt. 0) fldr = fname(1:q-1)

        inquire(file=trim(fldr)//"/grid_"//trim(reg%par%grid_src)//".txt", exist=found)
        if (.not. found) then
            fldr = "maps"
            inquire(file=trim(fldr)//"/grid_"//trim(reg%par%grid_src)//".txt", exist=found)
        end if
        if (.not. found) then
            call regions_error("regions_remap_init", &
                "no description of the grid of the files.", &
                "grid_src = "//trim(reg%par%grid_src)//new_line("a")// &
                "looked for grid_"//trim(reg%par%grid_src)//".txt next to "//trim(fname)// &
                " and in maps/")
        end if

        call grid_cdo_read_desc(grid_src, trim(reg%par%grid_src), trim(fldr))
        call map_init(reg%map, grid_src, grid, method="nn", fldr="maps")
        reg%remap = .true.

        return

    end subroutine regions_remap_init

    subroutine read_int(reg, fname, varname, none, var)
        ! Read an integer field, remapped onto the target grid if needed.
        ! none: the value of target cells without a source value.

        implicit none

        type(regions_class),  intent(IN)  :: reg
        character(len=*),     intent(IN)  :: fname, varname
        integer,              intent(IN)  :: none
        integer, allocatable, intent(OUT) :: var(:,:)

        integer, allocatable :: src(:,:)
        logical, allocatable :: filled(:,:)

        allocate(var(reg%nx,reg%ny))

        if (reg%remap) then
            allocate(src(nc_size(fname,"xc"),nc_size(fname,"yc")))
            allocate(filled(reg%nx,reg%ny))
            call nc_read(fname, varname, src)
            call map_field(reg%map, varname, src, var, mask2=filled, reset=.true.)
            where (.not. filled) var = none
        else
            call nc_read(fname, varname, var)
        end if

        return

    end subroutine read_int

    subroutine read_real(reg, fname, varname, var)
        ! Read a real field, remapped onto the target grid if needed (missing
        ! values are mv).

        implicit none

        type(regions_class),   intent(IN)  :: reg
        character(len=*),      intent(IN)  :: fname, varname
        real(wp), allocatable, intent(OUT) :: var(:,:)

        real(wp), allocatable :: src(:,:)

        allocate(var(reg%nx,reg%ny))

        if (reg%remap) then
            allocate(src(nc_size(fname,"xc"),nc_size(fname,"yc")))
            call nc_read(fname, varname, src, missing_value=mv)
            call map_field(reg%map, varname, src, var, missing_value=mv, reset=.true.)
        else
            call nc_read(fname, varname, var, missing_value=mv)
        end if

        return

    end subroutine read_real

    ! ===== Flag tables ========================================================

    subroutine flag_table_read(tab, fname, varname)
        ! Codes and names of a variable from its flag_values and flag_meanings.

        implicit none

        type(flag_table_class), intent(OUT) :: tab
        character(len=*),       intent(IN)  :: fname, varname

        character(len=:), allocatable :: meanings
        integer :: n, nm, k, i0, i1

        if (.not. nc_exists_attr(fname, varname, "flag_values")) then
            call regions_error("flag_table_read", "variable without flag_values.", &
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
            call regions_error("flag_table_read", "flag_values and flag_meanings differ in length.", &
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

    ! ===== Selection ==========================================================

    function regions_mask(reg, name) result(mask)
        ! A named mask of the namelist (or of regions_init_arg).

        implicit none

        type(regions_class), intent(IN) :: reg
        character(len=*),    intent(IN) :: name
        logical :: mask(reg%nx,reg%ny)

        integer :: k

        do k = 1, size(reg%masks)
            if (trim(reg%masks(k)%name) .eq. trim(name)) then
                mask = reg%masks(k)%mask
                return
            end if
        end do

        call regions_error("regions_mask", "no mask with this name.", "name = "//trim(name))

    end function regions_mask

    function regions_select(reg, expr) result(mask)
        ! Cells selected by an expression (see the module header).

        implicit none

        type(regions_class), intent(IN) :: reg
        character(len=*),    intent(IN) :: expr
        logical :: mask(reg%nx,reg%ny)

        logical :: clause(reg%nx,reg%ny)
        integer :: i0, i1, j0, j1, n

        n = len_trim(expr)
        if (n .eq. 0) then
            call regions_error("regions_select", "empty expression.")
        end if

        mask = .false.

        ! OR over the clauses (|), AND over the terms of a clause (&)
        i0 = 1
        do while (i0 .le. n+1)
            i1 = index(expr(i0:n), "|")
            i1 = merge(n, i0+i1-2, i1 .eq. 0)

            clause = .true.
            j0 = i0
            do while (j0 .le. i1+1)
                j1 = index(expr(j0:i1), "&")
                j1 = merge(i1, j0+j1-2, j1 .eq. 0)
                clause = clause .and. select_term(reg, expr(j0:j1), expr)
                j0 = j1 + 2
            end do

            mask = mask .or. clause
            i0 = i1 + 2
        end do

        return

    end function regions_select

    function select_term(reg, term, expr) result(mask)
        ! Cells selected by one term [~]field:value,value,... (or all, none).

        implicit none

        type(regions_class), intent(IN) :: reg
        character(len=*),    intent(IN) :: term
        character(len=*),    intent(IN) :: expr       ! for error messages
        logical :: mask(reg%nx,reg%ny)

        character(len=len_expr) :: t, field, values
        type(flag_table_class) :: tab_mask
        integer, allocatable :: codes(:)
        logical :: negate
        integer :: q, k, ks

        t = adjustl(term)
        negate = (t(1:1) .eq. "~")
        if (negate) t = adjustl(t(2:))

        if (len_trim(t) .eq. 0) then
            call regions_error("regions_select", "empty term.", "expression = "//trim(expr))
        end if

        q = index(t, ":")
        if (q .eq. 0) then
            select case(lower(t))
                case("all")
                    mask = .true.
                case("none")
                    mask = .false.
                case default
                    call regions_error("regions_select", &
                        "a term must be field:values, all or none.", &
                        "term       = "//trim(t)//new_line("a")//"expression = "//trim(expr))
            end select
        else
            field  = lower(adjustl(t(1:q-1)))
            values = t(q+1:)

            select case(trim(field))

                case("region")
                    codes = resolve_values(reg%tab_region, values, .true., "region", expr)
                    mask = .false.
                    do k = 1, size(codes)
                        mask = mask .or. in_region(reg%region_3, codes(k))
                    end do

                case("zone")
                    codes = resolve_values(reg%tab_zone, values, .false., "zone", expr)
                    mask = .false.
                    do k = 1, size(codes)
                        mask = mask .or. (reg%zone .eq. codes(k))
                    end do

                case default
                    ! <set>, <set>.group or <set>.mask
                    q  = index(field, ".")
                    ks = find_set(reg, field(1:merge(len_trim(field), q-1, q .eq. 0)))
                    if (ks .eq. 0) then
                        call regions_error("regions_select", &
                            "unknown field (region, zone or a loaded basin set).", &
                            "field      = "//trim(field)//new_line("a")//"expression = "//trim(expr))
                    end if

                    associate(bs => reg%basins(ks))
                    mask = .false.
                    if (q .eq. 0) then
                        codes = resolve_values(bs%tab_basin, values, .false., field, expr)
                        do k = 1, size(codes)
                            mask = mask .or. (bs%basin .eq. codes(k))
                        end do
                    else if (trim(field(q+1:)) .eq. "group") then
                        if (.not. bs%with_group) then
                            call regions_error("regions_select", "the basin set has no basin_group.", &
                                "set        = "//trim(bs%name)//new_line("a")//"expression = "//trim(expr))
                        end if
                        codes = resolve_values(bs%tab_group, values, .false., field, expr)
                        do k = 1, size(codes)
                            mask = mask .or. (bs%basin_group .eq. codes(k))
                        end do
                    else if (trim(field(q+1:)) .eq. "mask") then
                        tab_mask%codes = [0, 1]
                        tab_mask%names = [character(len=len_name) :: "0", "1"]
                        codes = resolve_values(tab_mask, values, .false., field, expr)
                        do k = 1, size(codes)
                            mask = mask .or. (bs%basin_mask .eq. codes(k))
                        end do
                    else
                        call regions_error("regions_select", &
                            "a basin set field is <set>, <set>.group or <set>.mask.", &
                            "field      = "//trim(field)//new_line("a")//"expression = "//trim(expr))
                    end if
                    end associate

            end select
        end if

        if (negate) mask = .not. mask

        return

    end function select_term

    function resolve_values(tab, values, is_region, field, expr) result(codes)
        ! Codes of a comma-separated list of names, codes and (for regions)
        ! paths such as 1.3.1.

        implicit none

        type(flag_table_class), intent(IN) :: tab
        character(len=*),       intent(IN) :: values
        logical,                intent(IN) :: is_region
        character(len=*),       intent(IN) :: field, expr
        integer, allocatable :: codes(:)

        character(len=len_name) :: v
        integer, allocatable :: c(:)
        integer :: i0, i1, n, code, ios

        allocate(codes(0))

        n = len_trim(values)
        i0 = 1
        do while (i0 .le. n+1)
            i1 = index(values(i0:n), ",")
            i1 = merge(n, i0+i1-2, i1 .eq. 0)
            v = adjustl(values(i0:i1))
            i0 = i1 + 2

            if (len_trim(v) .eq. 0) then
                call regions_error("regions_select", "empty value.", &
                    "field      = "//trim(field)//new_line("a")//"expression = "//trim(expr))
            end if

            if (verify(trim(v), "0123456789") .eq. 0) then
                read(v, *, iostat=ios) code
                codes = [codes, code]
            else if (is_region .and. verify(trim(v), "0123456789.") .eq. 0) then
                code = region_code(v)
                if (code .le. 0) then
                    call regions_error("regions_select", "invalid region path.", &
                        "path       = "//trim(v)//new_line("a")//"expression = "//trim(expr))
                end if
                codes = [codes, code]
            else
                c = flag_codes(tab, v)
                if (size(c) .eq. 0) then
                    call regions_error("regions_select", "unknown name (not on this grid).", &
                        "field      = "//trim(field)//new_line("a")// &
                        "name       = "//trim(v)//new_line("a")// &
                        "expression = "//trim(expr)//new_line("a")// &
                        "names      = "//join_names(tab))
                end if
                codes = [codes, c]
            end if
        end do

        return

    end function resolve_values

    integer function find_set(reg, name)

        implicit none

        type(regions_class), intent(IN) :: reg
        character(len=*),    intent(IN) :: name

        integer :: k

        find_set = 0
        do k = 1, size(reg%basins)
            if (lower(reg%basins(k)%name) .eq. lower(name)) then
                find_set = k
                exit
            end if
        end do

    end function find_set

    ! ===== Region codes (as FesmUtils.jl region_codes.jl) =====================

    pure function region_code(path) result(code)
        ! Code of a region path ("1.3.1" => 10301); 0 if the path is invalid
        ! (each part must be in 1..99).

        implicit none

        character(len=*), intent(IN) :: path
        integer :: code

        integer :: i0, i1, n, part, ios

        code = 0
        n = len_trim(path)
        if (n .eq. 0) return

        i0 = 1
        do while (i0 .le. n+1)
            i1 = index(path(i0:n), ".")
            i1 = merge(n, i0+i1-2, i1 .eq. 0)
            if (i1 .lt. i0) then
                code = 0
                return
            end if
            read(path(i0:i1), *, iostat=ios) part
            if (ios .ne. 0 .or. part .lt. 1 .or. part .gt. 99) then
                code = 0
                return
            end if
            code = code*100 + part
            i0 = i1 + 2
        end do

    end function region_code

    elemental function region_level(code) result(level)
        ! Level of a region code (1 for 1..99, 2 for 101..9999, ...); 0 for
        ! code < 1 (no region).

        implicit none

        integer, intent(IN) :: code
        integer :: level

        integer :: c

        level = 0
        if (code .lt. 1) return

        level = 1
        c = code
        do while (c .ge. 100)
            c = c / 100
            level = level + 1
        end do

    end function region_level

    elemental function region_ancestor(code, level) result(anc)
        ! Code of the region at level that contains region code (code itself
        ! at its own level); 0 if there is none.

        implicit none

        integer, intent(IN) :: code, level
        integer :: anc

        integer :: n

        anc = 0
        n = region_level(code)
        if (level .lt. 1 .or. level .gt. n) return

        anc = code / 100**(n-level)

    end function region_ancestor

    pure function region_path(code) result(path)
        ! Path of a region code (10301 => "1.3.1"); "" if the code is invalid.

        implicit none

        integer, intent(IN) :: code
        character(len=32) :: path

        character(len=4) :: part
        integer :: k, n, c, p

        path = ""
        n = region_level(code)
        c = code
        do k = 1, n
            p = mod(c, 100)
            if (p .eq. 0) then
                path = ""
                return
            end if
            write(part, "(i0)") p
            if (k .eq. 1) then
                path = trim(part)
            else
                path = trim(part)//"."//trim(path)
            end if
            c = c / 100
        end do

    end function region_path

    elemental function in_region(codes, code) result(inside)
        ! Whether a cell with region code codes lies in region code or one of
        ! its subregions (cells with code 0 are in no region).

        implicit none

        integer, intent(IN) :: codes, code
        logical :: inside

        integer :: level

        level  = region_level(code)
        inside = (codes .gt. 0 .and. level .gt. 0)
        if (inside) inside = (region_ancestor(codes, level) .eq. code)

    end function in_region

    ! ===== Helpers ============================================================

    pure function lower(str) result(res)

        implicit none

        character(len=*), intent(IN) :: str
        character(len=len_trim(str)) :: res

        integer :: k, c

        res = str
        do k = 1, len(res)
            c = ichar(res(k:k))
            if (c .ge. ichar("A") .and. c .le. ichar("Z")) res(k:k) = char(c + 32)
        end do

    end function lower

    pure function normalize_name(name) result(key)
        ! Lower case, spaces between words as underscores.

        implicit none

        character(len=*), intent(IN) :: name
        character(len=len_name) :: key

        integer :: k

        key = lower(adjustl(name))
        do k = 1, len_trim(key)
            if (key(k:k) .eq. " ") key(k:k) = "_"
        end do

    end function normalize_name

    function join_names(tab) result(str)

        implicit none

        type(flag_table_class), intent(IN) :: tab
        character(len=:), allocatable :: str

        integer :: k

        str = ""
        do k = 1, size(tab%names)
            str = str//trim(tab%names(k))//" "
        end do

    end function join_names

    subroutine regions_error(proc, msg, detail)
        ! Abort with a framed message on error_unit.

        implicit none

        character(len=*), intent(IN)           :: proc
        character(len=*), intent(IN)           :: msg
        character(len=*), intent(IN), optional :: detail

        integer :: p0, p1

        write(error_unit,"(a)") ""
        write(error_unit,"(a)") "regions:: error in "//trim(proc)
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

        write(error_unit,"(a)") "  stopped by regions."
        error stop 1

    end subroutine regions_error

end module regions
