module regions
    ! Named categorical layers on a grid (regions, zones, basins, ... from any
    ! source), and masks derived from them by selection expressions.
    !
    ! A layer is a 2D integer field with a table of names of its codes.
    ! Layers come from any file (regions_load_layer; namelist keys layers,
    ! path_<layer>, ...), from the program (regions_add_layer), or from the
    ! files of FesmData v2 ($ICE_DATA/v2/<Domain>/<GRID>/):
    !   <GRID>_REGIONS.nc       layers region (hierarchical, from region_1..3)
    !                           and zone (regions_load_fesmdata)
    !   <GRID>_BASINS-<set>.nc  layers <set> (basin), <set>.group (basin_group)
    !                           and <set>.mask (basin_mask) (regions_load_basins)
    ! The names of the codes are given (codes, names), else the CF attributes
    ! flag_values and flag_meanings of the variable, else the codes
    ! themselves. Negative values (e.g. fill values) and cells without a
    ! source value after remapping are regions_undefined.
    !
    ! Hierarchical layers (region) have codes of two decimal digits per level
    ! below the first: the path "1.3.1" is the code 10301, and a cell holds
    ! the code of its deepest region, so in_region(values, code) selects a
    ! region with its subregions.
    !
    ! Selection expressions (regions_select, and the named masks) combine
    ! terms layer:value,value,... where a comma is an OR:
    !   region:Greenland,1.5        names, paths or codes, at any level
    !   zone:land,continental_shelf names or codes
    !   Zwally2012:11,12            any layer, e.g. a basin set
    !   Zwally2012.group:1
    ! A term prefixed by ~ is negated, & joins terms (AND) and | joins
    ! clauses (OR, binding weaker than &). "all" and "none" select every cell
    ! and no cell. Names are matched ignoring case, with spaces read as
    ! underscores. Example: "region:Greenland & ~zone:open_ocean | Zwally2012:11"
    !
    ! With a target grid, a file whose grid (its global attribute grid_name,
    ! or grid_<layer> in the namelist) has another name is remapped onto it by
    ! nearest neighbour; the grid is read from grid_<name>.txt next to the
    ! file, else in maps/. A file without a grid name must be on the target
    ! grid.

    use, intrinsic :: iso_fortran_env, only : error_unit

    use precision
    use ncio
    use nml
    use coordinates, only : grid_class
    use mapping,     only : map_class
    use fesmdata

    implicit none

    integer, parameter :: len_expr     = 1000
    integer, parameter :: len_layer    = 56
    integer, parameter :: n_layers_max = 50
    integer, parameter :: n_masks_max  = 50
    integer, parameter :: n_codes_max  = 500

    ! Code of an undefined cell
    integer, parameter :: regions_undefined = -1

    ! Names that are not layers (namelist keys path_regions and path_basins,
    ! and the terms all and none)
    character(len=8), parameter :: names_reserved(4) = &
        [character(len=8) :: "regions", "basins", "all", "none"]

    type region_layer_class
        character(len=len_layer) :: name
        logical                  :: hier = .false.  ! hierarchical codes
        integer, allocatable     :: values(:,:)
        type(flag_table_class)   :: tab
    end type

    type region_mask_class
        character(len=56)       :: name
        character(len=len_expr) :: expr
        logical, allocatable    :: mask(:,:)
    end type

    type regions_class
        integer :: nx = 0, ny = 0
        logical          :: with_grid = .false.
        type(grid_class) :: grid                    ! target grid (if with_grid)

        type(region_layer_class), allocatable :: layers(:)
        type(region_mask_class),  allocatable :: masks(:)
    end type

    private
    public :: flag_table_class, region_layer_class, region_mask_class, regions_class
    public :: regions_init, regions_init_nml, regions_end
    public :: regions_load_layer, regions_load_fesmdata, regions_load_basins
    public :: regions_add_layer, regions_add_mask, regions_write
    public :: regions_find, regions_select, regions_mask, regions_basin_ids
    public :: flag_codes, flag_name
    public :: region_code, region_level, region_ancestor, region_path, in_region
    public :: regions_undefined

contains

    ! ===== Initialization =====================================================

    subroutine regions_init(reg, grid, nx, ny)
        ! An empty set of layers on grid (or of size nx, ny; else of the size
        ! of the first layer). "all", "none", "None" and "domain" work on it.

        implicit none

        type(regions_class), intent(INOUT) :: reg
        type(grid_class),    intent(IN), optional :: grid
        integer,             intent(IN), optional :: nx, ny

        call regions_end(reg)

        allocate(reg%layers(0), reg%masks(0))

        if (present(grid)) then
            reg%with_grid = .true.
            reg%grid = grid
            reg%nx = grid%G%nx
            reg%ny = grid%G%ny
        else if (present(nx) .and. present(ny)) then
            reg%nx = nx
            reg%ny = ny
        end if

        return

    end subroutine regions_init

    subroutine regions_init_nml(reg, filename, group, domain, grid_name, subs, grid, verbose)
        ! Load layers and named masks as given by a namelist group (all keys
        ! optional; without any, an empty set on grid):
        !   path_regions = "ice_data/v2/{domain}/{grid_name}/{grid_name}_REGIONS.nc"
        !   path_basins  = "ice_data/v2/{domain}/{grid_name}/{grid_name}_BASINS-{set}.nc"
        !   basin_sets   = "Zwally2012"
        !   layers       = "glaciers"            layers of any file:
        !   path_glaciers  = "glaciers.nc"
        !   var_glaciers   = "id"                (default: the layer name)
        !   codes_glaciers = 1 2                 (default: flag_values, else the codes)
        !   names_glaciers = "Rhone" "Aletsch"   (default: flag_meanings, else the codes)
        !   hier_glaciers  = False               (default: the attribute hierarchical)
        !   grid_glaciers  = "ALPS-1KM"          (default: the attribute grid_name)
        !   masks        = "grl_shelf"
        !   mask_grl_shelf = "region:Greenland & zone:continental_shelf"
        ! {grid_name} in the paths: grid_name. grid: the target grid, onto
        ! which the layers are remapped from the grids of their files.

        implicit none

        type(regions_class), intent(INOUT) :: reg
        character(len=*),    intent(IN)    :: filename
        character(len=*),    intent(IN)    :: group
        character(len=*),    intent(IN), optional :: domain
        character(len=*),    intent(IN), optional :: grid_name
        character(len=*),    intent(IN), optional :: subs(:,:)   ! extra {key}->value path substitutions
        type(grid_class),    intent(IN), optional :: grid
        logical,             intent(IN), optional :: verbose

        character(len=1024)      :: path, path_basins
        character(len=len_layer) :: names(n_layers_max), var, gname
        character(len=len_name)  :: code_names(n_codes_max)
        character(len=len_expr)  :: expr
        integer :: codes(n_codes_max)
        logical :: print_summary, hier
        integer :: k, n, nc

        print_summary = .true.
        if (present(verbose)) print_summary = verbose

        call regions_init(reg, grid)

        if (print_summary) write(*,*) "Loading: ", trim(filename), ":: ", trim(group)

        ! FesmData regions and zones
        if (nml_has_param(filename, group, "path_regions")) then
            call nml_read(filename, group, "path_regions", path)
            call fesmdata_parse_path(path, domain, grid_name, subs)
            if (print_summary) write(*,*) "path_regions = ", trim(path)
            call regions_load_fesmdata(reg, path)
        end if

        ! FesmData basin sets
        names = ""
        if (nml_has_param(filename, group, "basin_sets")) call nml_read(filename, group, "basin_sets", names)
        n = count(len_trim(names) .gt. 0)
        if (n .gt. 0) then
            if (.not. nml_has_param(filename, group, "path_basins")) then
                call regions_error("regions_init_nml", "basin_sets without path_basins.", "group = "//trim(group))
            end if
            call nml_read(filename, group, "path_basins", path_basins)
            call fesmdata_parse_path(path_basins, domain, grid_name, subs)
            names(1:n) = pack(names, len_trim(names) .gt. 0)
            do k = 1, n
                path = path_basins
                call nml_replace(path, "{set}", trim(names(k)))
                if (print_summary) write(*,*) "basins "//trim(names(k))//" = ", trim(path)
                call regions_load_basins(reg, names(k), path)
            end do
        end if

        ! Layers of any file
        names = ""
        if (nml_has_param(filename, group, "layers")) call nml_read(filename, group, "layers", names)
        n = count(len_trim(names) .gt. 0)
        names(1:n) = pack(names, len_trim(names) .gt. 0)
        do k = 1, n
            associate(l => names(k))
                call nml_read(filename, group, "path_"//trim(l), path)
                call fesmdata_parse_path(path, domain, grid_name, subs)

                var = l
                if (nml_has_param(filename, group, "var_"//trim(l))) &
                    call nml_read(filename, group, "var_"//trim(l), var)
                gname = ""
                if (nml_has_param(filename, group, "grid_"//trim(l))) &
                    call nml_read(filename, group, "grid_"//trim(l), gname)

                if (print_summary) write(*,*) "layer "//trim(l)//" = ", trim(path), " :: ", trim(var)

                nc = -1
                if (nml_has_param(filename, group, "codes_"//trim(l))) then
                    codes = huge(0)
                    code_names = ""
                    call nml_read(filename, group, "codes_"//trim(l), codes)
                    call nml_read(filename, group, "names_"//trim(l), code_names)
                    nc = count(codes .ne. huge(0))
                    if (count(len_trim(code_names) .gt. 0) .ne. nc) then
                        call regions_error("regions_init_nml", "codes_ and names_ of a layer differ in length.", &
                            "group = "//trim(group)//new_line("a")//"layer = "//trim(l))
                    end if
                end if

                if (nml_has_param(filename, group, "hier_"//trim(l))) then
                    call nml_read(filename, group, "hier_"//trim(l), hier)
                    if (nc .ge. 0) then
                        call regions_load_layer(reg, l, path, var, hier=hier, codes=codes(1:nc), &
                                                names=code_names(1:nc), grid_name=gname)
                    else
                        call regions_load_layer(reg, l, path, var, hier=hier, grid_name=gname)
                    end if
                else
                    if (nc .ge. 0) then
                        call regions_load_layer(reg, l, path, var, codes=codes(1:nc), &
                                                names=code_names(1:nc), grid_name=gname)
                    else
                        call regions_load_layer(reg, l, path, var, grid_name=gname)
                    end if
                end if
            end associate
        end do

        ! Named masks
        names = ""
        if (nml_has_param(filename, group, "masks")) call nml_read(filename, group, "masks", names)
        n = count(len_trim(names) .gt. 0)
        names(1:n) = pack(names, len_trim(names) .gt. 0)
        do k = 1, n
            call nml_read(filename, group, "mask_"//trim(names(k)), expr)
            if (print_summary) write(*,*) "mask_"//trim(names(k))//" = ", trim(expr)
            call regions_add_mask(reg, names(k), expr)
        end do

        return

    end subroutine regions_init_nml

    subroutine regions_end(reg)

        implicit none

        type(regions_class), intent(INOUT) :: reg

        type(regions_class) :: reg0

        reg = reg0

        return

    end subroutine regions_end

    ! ===== Layers =============================================================

    subroutine regions_load_layer(reg, name, path, var, hier, codes, names, grid_name)
        ! A layer from the integer variable var (default: name) of any file,
        ! on the target grid (see the module header). codes, names: the names
        ! of the codes. hier: hierarchical codes (default: the attribute
        ! hierarchical of the variable). grid_name: the grid of a file
        ! without the attribute grid_name.

        implicit none

        type(regions_class), intent(INOUT) :: reg
        character(len=*),    intent(IN)    :: name, path
        character(len=*),    intent(IN), optional :: var
        logical,             intent(IN), optional :: hier
        integer,             intent(IN), optional :: codes(:)
        character(len=*),    intent(IN), optional :: names(:)
        character(len=*),    intent(IN), optional :: grid_name

        call check_layer_name(name, "regions_load_layer")
        call load_layer(reg, name, path, var, hier, codes, names, grid_name)

        return

    end subroutine regions_load_layer

    subroutine regions_load_fesmdata(reg, path)
        ! The layers region (hierarchical; values of region_3, names of
        ! region_1..3) and zone of a FesmData v2 regions file.

        implicit none

        type(regions_class), intent(INOUT) :: reg
        character(len=*),    intent(IN)    :: path

        type(flag_table_class) :: tab
        integer :: k

        if (.not. nc_exists_var(path, "region_1")) then
            call regions_error("regions_load_fesmdata", &
                "not a FesmData v2 regions file (no variable region_1).", &
                "file = "//trim(path)//new_line("a")// &
                "FesmData v1 regions files (one variable mask) are not supported.")
        end if

        call load_layer(reg, "region", path, "region_3", hier=.true.)
        associate(tab_region => reg%layers(size(reg%layers))%tab)
            do k = 1, 2
                call flag_table_read(tab, path, "region_"//char(ichar("0")+k))
                call flag_table_merge(tab_region, tab)
            end do
        end associate

        call load_layer(reg, "zone", path, "zone")

        return

    end subroutine regions_load_fesmdata

    subroutine regions_load_basins(reg, set, path)
        ! The layers <set> (basin), <set>.group (basin_group) and <set>.mask
        ! (basin_mask) of a FesmData v2 basins file (the last two if present).

        implicit none

        type(regions_class), intent(INOUT) :: reg
        character(len=*),    intent(IN)    :: set, path

        call check_layer_name(set, "regions_load_basins")

        call load_layer(reg, set, path, "basin")
        if (nc_exists_var(path, "basin_group")) call load_layer(reg, trim(set)//".group", path, "basin_group")
        if (nc_exists_var(path, "basin_mask"))  call load_layer(reg, trim(set)//".mask",  path, "basin_mask")

        return

    end subroutine regions_load_basins

    subroutine regions_add_layer(reg, name, field, codes, names, hier)
        ! A layer from a field of the program (on the grid of the set).

        implicit none

        type(regions_class), intent(INOUT) :: reg
        character(len=*),    intent(IN)    :: name
        integer,             intent(IN)    :: field(:,:)
        integer,             intent(IN), optional :: codes(:)
        character(len=*),    intent(IN), optional :: names(:)
        logical,             intent(IN), optional :: hier

        type(region_layer_class) :: layer

        call check_layer_name(name, "regions_add_layer")

        layer%name   = name
        layer%values = field
        where (layer%values .lt. 0) layer%values = regions_undefined
        if (present(hier)) layer%hier = hier
        call layer_table(layer, codes, names)

        call add_layer(reg, layer)

        return

    end subroutine regions_add_layer

    subroutine load_layer(reg, name, path, var, hier, codes, names, grid_name)
        ! Read a layer (see regions_load_layer), without checking its name.

        implicit none

        type(regions_class), intent(INOUT) :: reg
        character(len=*),    intent(IN)    :: name, path
        character(len=*),    intent(IN), optional :: var
        logical,             intent(IN), optional :: hier
        integer,             intent(IN), optional :: codes(:)
        character(len=*),    intent(IN), optional :: names(:)
        character(len=*),    intent(IN), optional :: grid_name

        type(region_layer_class) :: layer
        character(len=len_layer) :: varname
        type(map_class) :: map
        logical :: remap
        integer :: nx, ny, ival

        varname = name
        if (present(var)) varname = var

        if (.not. nc_exists_var(path, varname)) then
            call regions_error("regions_load_layer", "variable not in the file.", &
                "layer    = "//trim(name)//new_line("a")// &
                "file     = "//trim(path)//new_line("a")//"variable = "//trim(varname))
        end if

        if (reg%with_grid) then
            call fesmdata_map_init(map, remap, nx, ny, path, varname, "nn", grid=reg%grid, grid_name=grid_name)
        else
            call fesmdata_map_init(map, remap, nx, ny, path, varname, "nn", grid_name=grid_name)
        end if

        layer%name = name
        call fesmdata_read_field(path, varname, layer%values, nx, ny, remap, map, regions_undefined)
        where (layer%values .lt. 0) layer%values = regions_undefined

        if (present(hier)) then
            layer%hier = hier
        else if (nc_exists_attr(path, varname, "hierarchical")) then
            call nc_read_attr(path, varname, "hierarchical", ival)
            layer%hier = (ival .eq. 1)
        end if

        if (.not. present(codes) .and. nc_exists_attr(path, varname, "flag_values")) then
            call flag_table_read(layer%tab, path, varname)
        else
            call layer_table(layer, codes, names)
        end if

        call add_layer(reg, layer)

        return

    end subroutine load_layer

    subroutine layer_table(layer, codes, names)
        ! The table of a layer: codes and names, else its positive values
        ! (named by themselves).

        implicit none

        type(region_layer_class), intent(INOUT) :: layer
        integer,          intent(IN), optional  :: codes(:)
        character(len=*), intent(IN), optional  :: names(:)

        if (present(codes) .neqv. present(names)) then
            call regions_error("regions", "codes and names of a layer go together.", "layer = "//trim(layer%name))
        end if

        if (present(codes)) then
            if (size(codes) .ne. size(names)) then
                call regions_error("regions", "codes and names of a layer differ in length.", &
                    "layer = "//trim(layer%name))
            end if
            layer%tab%codes = codes
            allocate(layer%tab%names(size(names)))
            layer%tab%names = names
        else
            call flag_table_from_values(layer%tab, layer%values)
        end if

        return

    end subroutine layer_table

    subroutine add_layer(reg, layer)
        ! Append a layer (of the size of the set; the first sets it).

        implicit none

        type(regions_class),      intent(INOUT) :: reg
        type(region_layer_class), intent(IN)    :: layer

        if (.not. allocated(reg%layers)) call regions_init(reg)

        if (regions_find(reg, layer%name) .gt. 0) then
            call regions_error("regions", "a layer of this name exists already.", "layer = "//trim(layer%name))
        end if

        if (reg%nx .eq. 0 .and. reg%ny .eq. 0) then
            reg%nx = size(layer%values,1)
            reg%ny = size(layer%values,2)
        else if (size(layer%values,1) .ne. reg%nx .or. size(layer%values,2) .ne. reg%ny) then
            call regions_error("regions", "the layer differs in size from the set.", "layer = "//trim(layer%name))
        end if

        reg%layers = [reg%layers, layer]

        return

    end subroutine add_layer

    subroutine check_layer_name(name, proc)
        ! Names of layers: not reserved, and without "." (kept for the
        ! layers <set>.group and <set>.mask of basin sets).

        implicit none

        character(len=*), intent(IN) :: name, proc

        if (len_trim(name) .eq. 0 .or. index(name, ".") .gt. 0 .or. &
            any(names_reserved .eq. lower(name))) then
            call regions_error(proc, "invalid layer name (empty, with '.', or reserved: regions basins all none).", &
                "name = "//trim(name))
        end if

    end subroutine check_layer_name

    integer function regions_find(reg, name)
        ! Index of the layer name (ignoring case), 0 if there is none.

        implicit none

        type(regions_class), intent(IN) :: reg
        character(len=*),    intent(IN) :: name

        integer :: k

        regions_find = 0
        if (.not. allocated(reg%layers)) return
        do k = 1, size(reg%layers)
            if (lower(reg%layers(k)%name) .eq. lower(adjustl(name))) then
                regions_find = k
                exit
            end if
        end do

    end function regions_find

    subroutine regions_write(reg, filename, grid, layers)
        ! Write layers (default: all) with their tables (flag_values,
        ! flag_meanings) and the attribute hierarchical, so that they load
        ! again as they are. With grid, a new file on it (dimensions xc, yc,
        ! global attribute grid_name); else into the existing file filename
        ! (dimensions xc, yc).

        implicit none

        type(regions_class), intent(IN) :: reg
        character(len=*),    intent(IN) :: filename
        type(grid_class),    intent(IN), optional :: grid
        character(len=*),    intent(IN), optional :: layers(:)

        character(len=:), allocatable :: meanings
        integer :: k, j

        if (present(grid)) then
            call nc_create(filename)
            call nc_write_dim(filename, "xc", x=grid%G%x, units="km")
            call nc_write_dim(filename, "yc", x=grid%G%y, units="km")
            call nc_write_attr(filename, "grid_name", trim(grid%name))
        end if

        do k = 1, size(reg%layers)
            associate(l => reg%layers(k))
                if (present(layers)) then
                    if (.not. any([(lower(layers(j)) .eq. lower(l%name), j=1,size(layers))])) cycle
                end if

                call nc_write(filename, trim(l%name), l%values, dim1="xc", dim2="yc", units="1")
                if (size(l%tab%codes) .gt. 0) then
                    meanings = trim(l%tab%names(1))
                    do j = 2, size(l%tab%names)
                        meanings = meanings//" "//trim(l%tab%names(j))
                    end do
                    call nc_write_attr(filename, trim(l%name), "flag_values", l%tab%codes)
                    call nc_write_attr(filename, trim(l%name), "flag_meanings", meanings)
                end if
                if (l%hier) call nc_write_attr(filename, trim(l%name), "hierarchical", 1)
            end associate
        end do

        return

    end subroutine regions_write

    ! ===== Masks ==============================================================

    subroutine regions_add_mask(reg, name, expr)
        ! A named mask from a selection expression.

        implicit none

        type(regions_class), intent(INOUT) :: reg
        character(len=*),    intent(IN)    :: name, expr

        type(region_mask_class) :: m

        if (.not. allocated(reg%masks)) allocate(reg%masks(0))

        m%name = name
        m%expr = expr
        m%mask = regions_select(reg, expr)
        reg%masks = [reg%masks, m]

        return

    end subroutine regions_add_mask

    function regions_mask(reg, name) result(mask)
        ! A named mask (of the namelist, or of regions_add_mask).

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
        ! Basin ids of a layer (spec = its name, e.g. "Zwally2012" or
        ! "Zwally2012.group"); 0 = no basin. Also "None" (no basins, 0
        ! everywhere) and "domain" (the whole domain is one basin, 1).
        ! extent: the original extent of the basins (the layer <set>.mask of
        ! a basin set), else where ids > 0.

        implicit none

        type(regions_class), intent(IN) :: reg
        character(len=*),    intent(IN) :: spec
        logical, optional,   intent(OUT) :: extent(:,:)
        integer :: ids(reg%nx,reg%ny)

        integer :: q, k, km

        select case(trim(spec))
            case("None")
                ids = 0
                if (present(extent)) extent = .false.
                return
            case("domain")
                ids = 1
                if (present(extent)) extent = .true.
                return
        end select

        k = regions_find(reg, spec)
        if (k .eq. 0) then
            call regions_error("regions_basin_ids", "no layer of this name (or None, domain).", &
                "spec   = "//trim(spec)//new_line("a")//"layers = "//layer_names(reg))
        end if
        ids = max(reg%layers(k)%values, 0)

        if (present(extent)) then
            q  = index(spec, ".")
            km = regions_find(reg, spec(1:merge(len_trim(spec), q-1, q .eq. 0))//".mask")
            if (km .gt. 0) then
                extent = (reg%layers(km)%values .eq. 1)
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
        ! Cells selected by one term [~]layer:value,value,... (or all, none).

        implicit none

        type(regions_class), intent(IN) :: reg
        character(len=*),    intent(IN) :: term
        character(len=*),    intent(IN) :: expr       ! for error messages
        logical :: mask(reg%nx,reg%ny)

        character(len=len_expr) :: t
        integer, allocatable :: codes(:)
        logical :: negate
        integer :: q, k, kl

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
                        "a term must be layer:values, all or none.", &
                        "term       = "//trim(t)//new_line("a")//"expression = "//trim(expr))
            end select
        else
            kl = regions_find(reg, t(1:q-1))
            if (kl .eq. 0) then
                call regions_error("regions_select", "unknown layer.", &
                    "layer      = "//trim(adjustl(t(1:q-1)))//new_line("a")// &
                    "expression = "//trim(expr)//new_line("a")//"layers     = "//layer_names(reg))
            end if

            associate(l => reg%layers(kl))
                codes = resolve_values(l%tab, t(q+1:), l%hier, l%name, expr)
                mask = .false.
                do k = 1, size(codes)
                    if (l%hier) then
                        mask = mask .or. in_region(l%values, codes(k))
                    else
                        mask = mask .or. (l%values .eq. codes(k))
                    end if
                end do
            end associate
        end if

        if (negate) mask = .not. mask

        return

    end function select_term

    function resolve_values(tab, values, hier, layer, expr) result(codes)
        ! Codes of a comma-separated list of names, codes and (for
        ! hierarchical layers) paths such as 1.3.1.

        implicit none

        type(flag_table_class), intent(IN) :: tab
        character(len=*),       intent(IN) :: values
        logical,                intent(IN) :: hier
        character(len=*),       intent(IN) :: layer, expr
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
                    "layer      = "//trim(layer)//new_line("a")//"expression = "//trim(expr))
            end if

            if (verify(trim(v), "0123456789") .eq. 0) then
                read(v, *, iostat=ios) code
                codes = [codes, code]
            else if (hier .and. verify(trim(v), "0123456789.") .eq. 0) then
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
                        "layer      = "//trim(layer)//new_line("a")// &
                        "name       = "//trim(v)//new_line("a")// &
                        "expression = "//trim(expr)//new_line("a")// &
                        "names      = "//flag_names_joined(tab))
                end if
                codes = [codes, c]
            end if
        end do

        return

    end function resolve_values

    function layer_names(reg) result(str)
        ! All names of layers, space-separated (for messages).

        implicit none

        type(regions_class), intent(IN) :: reg
        character(len=:), allocatable :: str

        integer :: k

        str = ""
        if (.not. allocated(reg%layers)) return
        do k = 1, size(reg%layers)
            str = str//trim(reg%layers(k)%name)//" "
        end do

    end function layer_names

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
