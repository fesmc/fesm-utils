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
    ! A basin set may also be a custom file (path_basins_<set>, var_basins_<set>
    ! in the namelist): a 2D integer field of basin ids on the grid of the
    ! regions file; without flag attributes its names are the ids, and values
    ! below 1 (or missing) are no basin.
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
    use ncio
    use nml
    use coordinates, only : grid_class
    use mapping,     only : map_class, map_init
    use fesmdata

    implicit none

    integer, parameter :: len_expr = 1000
    integer, parameter :: n_sets_max  = 20
    integer, parameter :: n_masks_max = 50

    ! Zone of a remapped target cell without a source value
    integer, parameter :: zone_undefined = -1

    type basin_set_class
        character(len=56)    :: name                ! e.g. "Zwally2012"
        character(len=1024)  :: filename
        character(len=56)    :: varname             ! "basin" for FesmData files
        logical              :: with_group, with_mask
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
        character(len=1024),     allocatable :: basin_paths(:)  ! per set (template or custom)
        character(len=56),       allocatable :: basin_vars(:)   ! per set ("basin" or custom)
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
    public :: regions_select, regions_mask, regions_basin_ids
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

    subroutine regions_init_arg(reg, path_regions, path_basins, basin_sets, mask_names, mask_exprs, grid, &
                                basin_paths, basin_vars)
        ! Load regions, basin sets and named masks from explicit arguments
        ! (see regions_init_nml). path_basins is a template with {set};
        ! basin_paths and basin_vars give the file and variable of each set
        ! instead (custom sets; "" = from the template, "basin").

        implicit none

        type(regions_class), intent(INOUT) :: reg
        character(len=*),    intent(IN)    :: path_regions
        character(len=*),    intent(IN), optional :: path_basins
        character(len=*),    intent(IN), optional :: basin_sets(:)
        character(len=*),    intent(IN), optional :: mask_names(:)
        character(len=*),    intent(IN), optional :: mask_exprs(:)
        type(grid_class),    intent(IN), optional :: grid
        character(len=*),    intent(IN), optional :: basin_paths(:)
        character(len=*),    intent(IN), optional :: basin_vars(:)

        integer :: n_sets, n_masks, k

        n_sets = 0
        if (present(basin_sets)) n_sets = size(basin_sets)
        n_masks = 0
        if (present(mask_names)) n_masks = size(mask_names)

        if (n_sets .gt. 0 .and. .not. (present(path_basins) .or. present(basin_paths))) then
            call regions_error("regions_init_arg", "basin_sets are given, but no path_basins or basin_paths.")
        end if
        if (present(basin_paths)) then
            if (size(basin_paths) .ne. n_sets) &
                call regions_error("regions_init_arg", "basin_sets and basin_paths differ in length.")
        end if
        if (present(basin_vars)) then
            if (size(basin_vars) .ne. n_sets) &
                call regions_error("regions_init_arg", "basin_sets and basin_vars differ in length.")
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

        allocate(reg%par%basin_sets(n_sets), reg%par%basin_paths(n_sets), reg%par%basin_vars(n_sets))
        reg%par%basin_paths = ""
        reg%par%basin_vars  = "basin"
        if (n_sets .gt. 0) reg%par%basin_sets = basin_sets
        if (present(basin_paths)) reg%par%basin_paths = basin_paths
        if (present(basin_vars))  reg%par%basin_vars  = basin_vars
        do k = 1, n_sets
            if (len_trim(reg%par%basin_paths(k)) .eq. 0) then
                reg%par%basin_paths(k) = reg%par%path_basins
                call nml_replace(reg%par%basin_paths(k), "{set}", trim(reg%par%basin_sets(k)))
            end if
            if (len_trim(reg%par%basin_vars(k)) .eq. 0) reg%par%basin_vars(k) = "basin"
        end do

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
        reg%par%grid_src = fesmdata_grid_name(fname)
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
        call fesmdata_read_field(fname, "region_1", reg%region_1, reg%nx, reg%ny, reg%remap, reg%map, 0)
        call fesmdata_read_field(fname, "region_2", reg%region_2, reg%nx, reg%ny, reg%remap, reg%map, 0)
        call fesmdata_read_field(fname, "region_3", reg%region_3, reg%nx, reg%ny, reg%remap, reg%map, 0)
        call fesmdata_read_field(fname, "zone", reg%zone, reg%nx, reg%ny, reg%remap, reg%map, zone_undefined)
        call fesmdata_read_field(fname, "dist_shelfbreak", reg%dist_shelfbreak, reg%nx, reg%ny, reg%remap, reg%map)

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
            call basin_set_load(reg, reg%basins(k), reg%par%basin_sets(k), &
                                reg%par%basin_paths(k), reg%par%basin_vars(k))
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

    subroutine basin_set_load(reg, bs, name, filename, varname)
        ! A basin set: a FesmData basins file (varname = "basin", with
        ! basin_group and basin_mask if present) or a custom field of basin ids.

        implicit none

        type(regions_class),   intent(INOUT) :: reg
        type(basin_set_class), intent(INOUT) :: bs
        character(len=*),      intent(IN)    :: name, filename, varname

        character(len=256)   :: grid_name
        character(len=64), allocatable :: dim_names(:)
        integer, allocatable :: dims(:)
        logical :: fesm

        bs%name     = trim(name)
        bs%filename = trim(filename)
        bs%varname  = trim(varname)
        fesm = (trim(varname) .eq. "basin")

        ! All files must be on the grid of the regions file (one map): by the
        ! name of the grid if the file has one, else by the size of the field
        if (nc_exists_attr(bs%filename, "grid_name")) then
            grid_name = fesmdata_grid_name(bs%filename)
            if (trim(grid_name) .ne. trim(reg%par%grid_src)) then
                call regions_error("basin_set_load", &
                    "the basins file is on another grid than the regions file.", &
                    "file      = "//trim(bs%filename)//new_line("a")// &
                    "grid_name = "//trim(grid_name)//new_line("a")// &
                    "regions   = "//trim(reg%par%grid_src))
            end if
        else
            call nc_dims(bs%filename, bs%varname, dim_names, dims)
            if (size(dims) .lt. 2) then
                call regions_error("basin_set_load", "the basins are not a 2D field.", &
                    "file = "//trim(bs%filename)//new_line("a")//"variable = "//trim(bs%varname))
            end if
            if (dims(1) .ne. nc_size(reg%par%path_regions, "xc") .or. &
                dims(2) .ne. nc_size(reg%par%path_regions, "yc")) then
                call regions_error("basin_set_load", &
                    "the basins file (without grid_name) differs in size from the regions file.", &
                    "file = "//trim(bs%filename)//new_line("a")//"variable = "//trim(bs%varname))
            end if
        end if

        call fesmdata_read_field(bs%filename, bs%varname, bs%basin, reg%nx, reg%ny, reg%remap, reg%map, 0)
        where (bs%basin .lt. 1) bs%basin = 0

        if (nc_exists_attr(bs%filename, bs%varname, "flag_values")) then
            call flag_table_read(bs%tab_basin, bs%filename, bs%varname)
        else
            call flag_table_from_values(bs%tab_basin, bs%basin)
        end if

        bs%with_mask = fesm .and. nc_exists_var(bs%filename, "basin_mask")
        if (bs%with_mask) then
            call fesmdata_read_field(bs%filename, "basin_mask", bs%basin_mask, reg%nx, reg%ny, reg%remap, reg%map, 0)
        end if

        bs%with_group = fesm .and. nc_exists_var(bs%filename, "basin_group")
        if (bs%with_group) then
            call fesmdata_read_field(bs%filename, "basin_group", bs%basin_group, reg%nx, reg%ny, reg%remap, reg%map, 0)
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
        call fesmdata_parse_path(par%path_regions, domain, grid_name, subs)

        ! Basin sets: from the template path_basins, or custom (path_basins_<set>,
        ! var_basins_<set>)
        par%path_basins = ""
        sets = ""
        if (nml_has_param(filename, group, "basin_sets")) then
            call nml_read(filename, group, "basin_sets",  sets)
        end if
        if (nml_has_param(filename, group, "path_basins")) then
            call nml_read(filename, group, "path_basins", par%path_basins)
            call fesmdata_parse_path(par%path_basins, domain, grid_name, subs)
        end if
        n = count(len_trim(sets) .gt. 0)
        allocate(par%basin_sets(n), par%basin_paths(n), par%basin_vars(n))
        par%basin_sets = pack(sets, len_trim(sets) .gt. 0)
        do k = 1, n
            if (nml_has_param(filename, group, "path_basins_"//trim(par%basin_sets(k)))) then
                call nml_read(filename, group, "path_basins_"//trim(par%basin_sets(k)), par%basin_paths(k))
                call fesmdata_parse_path(par%basin_paths(k), domain, grid_name, subs)
            else if (len_trim(par%path_basins) .gt. 0) then
                par%basin_paths(k) = par%path_basins
                call nml_replace(par%basin_paths(k), "{set}", trim(par%basin_sets(k)))
            else
                call regions_error("regions_par_load", "a basin set without a file.", &
                    "set = "//trim(par%basin_sets(k))//new_line("a")// &
                    "give path_basins (template with {set}) or path_basins_"//trim(par%basin_sets(k)))
            end if
            par%basin_vars(k) = "basin"
            if (nml_has_param(filename, group, "var_basins_"//trim(par%basin_sets(k)))) then
                call nml_read(filename, group, "var_basins_"//trim(par%basin_sets(k)), par%basin_vars(k))
            end if
        end do

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
            do k = 1, size(par%basin_sets)
                write(*,*) "basins "//trim(par%basin_sets(k))//" = ", trim(par%basin_paths(k)), &
                           " :: ", trim(par%basin_vars(k))
            end do
            do k = 1, size(par%mask_names)
                write(*,*) "mask_"//trim(par%mask_names(k))//" = ", trim(par%mask_exprs(k))
            end do
        end if

        return

    end subroutine regions_par_load

    ! ===== Reading and remapping ==============================================

    subroutine regions_remap_init(reg, fname, grid)
        ! Nearest-neighbour map from the grid of the files onto grid (cached
        ! in maps/).

        implicit none

        type(regions_class), intent(INOUT) :: reg
        character(len=*),    intent(IN)    :: fname
        type(grid_class),    intent(IN)    :: grid

        type(grid_class) :: grid_src

        call fesmdata_grid_read(grid_src, fname)
        call map_init(reg%map, grid_src, grid, method="nn", fldr="maps")
        reg%remap = .true.

        return

    end subroutine regions_remap_init

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

    function regions_basin_ids(reg, spec, extent) result(ids)
        ! Basin ids of a loaded basin set: spec = "<set>" (basin) or
        ! "<set>.group" (basin_group); 0 = no basin. extent: the original
        ! extent of the basins (basin_mask), or where ids > 0 without one.

        implicit none

        type(regions_class), intent(IN) :: reg
        character(len=*),    intent(IN) :: spec
        logical, optional,   intent(OUT) :: extent(:,:)
        integer :: ids(reg%nx,reg%ny)

        integer :: q, ks

        q  = index(spec, ".")
        ks = find_set(reg, spec(1:merge(len_trim(spec), q-1, q .eq. 0)))
        if (ks .eq. 0) then
            call regions_error("regions_basin_ids", "not a loaded basin set.", "spec = "//trim(spec))
        end if

        if (q .eq. 0) then
            ids = reg%basins(ks)%basin
        else if (lower(spec(q+1:)) .eq. "group") then
            if (.not. reg%basins(ks)%with_group) then
                call regions_error("regions_basin_ids", "the basin set has no basin_group.", "spec = "//trim(spec))
            end if
            ids = reg%basins(ks)%basin_group
        else
            call regions_error("regions_basin_ids", "spec is <set> or <set>.group.", "spec = "//trim(spec))
        end if

        if (present(extent)) then
            if (reg%basins(ks)%with_mask) then
                extent = (reg%basins(ks)%basin_mask .eq. 1)
            else
                extent = (ids .gt. 0)
            end if
        end if

    end function regions_basin_ids

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
                        if (.not. bs%with_mask) then
                            call regions_error("regions_select", "the basin set has no basin_mask.", &
                                "set        = "//trim(bs%name)//new_line("a")//"expression = "//trim(expr))
                        end if
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
                        "names      = "//flag_names_joined(tab))
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
