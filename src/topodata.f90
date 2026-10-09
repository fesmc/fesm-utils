module topodata
    ! Topography of a FesmData v2 grid (<GRID>_TOPO-<product>.nc, from
    ! FesmData/Topo), read as it is in the file. Only the listed variables
    ! are loaded; the others stay unallocated.
    !
    !   z_bed, z_srf, H_ice, z_bed_sd    elevations and thickness (m)
    !   f_ocn, f_land, f_grnd, f_flt     area fractions of ice-free ocean,
    !                                    ice-free land, grounded and floating ice
    !   mask                             dominant surface type (topo_mask_*)
    !   src_id                           source with the largest weight
    !
    ! With a target grid whose name differs from the grid of the file, the
    ! real fields are remapped onto it with par%remap (default "con",
    ! conservative), mask and src_id by nearest neighbour. Missing values are
    ! mv; mask and src_id of target cells without a source value are -1.
    !
    ! Other files (e.g. ISMIP7) may be read too: par%names gives the variable
    ! of the file for each of par%vars. A file without the global attribute
    ! grid_name must be on the target grid (not remapped).

    use, intrinsic :: iso_fortran_env, only : error_unit

    use precision
    use ncio
    use nml
    use coordinates, only : grid_class
    use mapping,     only : map_class, map_init
    use fesmdata

    implicit none

    ! Classes of mask (FesmData/Topo)
    integer, parameter :: topo_mask_ocean = 0
    integer, parameter :: topo_mask_land  = 1
    integer, parameter :: topo_mask_grnd  = 2
    integer, parameter :: topo_mask_flt   = 3

    integer, parameter :: n_vars_max = 10
    character(len=8), parameter :: vars_all(n_vars_max) = &
        [character(len=8) :: "z_bed", "z_srf", "H_ice", "z_bed_sd", &
                             "f_ocn", "f_land", "f_grnd", "f_flt", "mask", "src_id"]
    character(len=8), parameter :: vars_default(4) = &
        [character(len=8) :: "z_bed", "z_srf", "H_ice", "z_bed_sd"]

    type topodata_param_class
        character(len=1024) :: path                 ! <GRID>_TOPO-<product>.nc
        character(len=56), allocatable :: vars(:)   ! variables to load
        character(len=56), allocatable :: names(:)  ! their names in the file
        character(len=32)  :: remap                 ! method for the real fields

        ! Internal parameters
        character(len=256) :: grid_src              ! grid of the file ("": none)
        character(len=256) :: grid_tgt              ! grid of the fields
        character(len=256) :: product               ! global attributes of the file
        character(len=512) :: sources
    end type

    type topodata_class
        type(topodata_param_class) :: par
        integer :: nx, ny

        real(wp), allocatable :: z_bed(:,:), z_srf(:,:), H_ice(:,:), z_bed_sd(:,:)
        real(wp), allocatable :: f_ocn(:,:), f_land(:,:), f_grnd(:,:), f_flt(:,:)
        integer,  allocatable :: mask(:,:), src_id(:,:)
        type(flag_table_class) :: tab_mask, tab_src

        ! Online remapping from the grid of the file: map for the real
        ! fields, map_nn for mask and src_id (unused if remap = "nn")
        logical         :: remap = .false.
        type(map_class) :: map, map_nn
    end type

    private
    public :: topo_mask_ocean, topo_mask_land, topo_mask_grnd, topo_mask_flt
    public :: topodata_param_class, topodata_class
    public :: topodata_init_nml, topodata_init_arg, topodata_end

contains

    subroutine topodata_init_nml(td, filename, group, domain, grid_name, subs, grid, verbose)
        ! Load the topography as given by a namelist group:
        !   path  = "ice_data/v2/{domain}/{grid_name}/{grid_name}_TOPO-BedMachine-v6.nc"
        !   vars  = "z_bed" "z_srf" "H_ice" "z_bed_sd"   (optional; this default)
        !   names = "bed" "surface" "thickness" ""       (optional; "" = as vars)
        !   remap = "con"                                (optional; this default)
        ! {grid_name} is the grid of the file. grid: the target grid, onto which
        ! the fields are remapped if it differs from the grid of the file.

        implicit none

        type(topodata_class), intent(INOUT) :: td
        character(len=*),     intent(IN)    :: filename
        character(len=*),     intent(IN)    :: group
        character(len=*),     intent(IN), optional :: domain
        character(len=*),     intent(IN), optional :: grid_name
        character(len=*),     intent(IN), optional :: subs(:,:)   ! extra {key}->value path substitutions
        type(grid_class),     intent(IN), optional :: grid
        logical,              intent(IN), optional :: verbose

        call topodata_par_load(td%par, filename, group, domain, grid_name, subs, verbose)

        call topodata_init_data(td, grid)

        return

    end subroutine topodata_init_nml

    subroutine topodata_init_arg(td, path, vars, names, remap, grid)
        ! Load the topography from explicit arguments (see topodata_init_nml).

        implicit none

        type(topodata_class), intent(INOUT) :: td
        character(len=*),     intent(IN)    :: path
        character(len=*),     intent(IN), optional :: vars(:)
        character(len=*),     intent(IN), optional :: names(:)
        character(len=*),     intent(IN), optional :: remap
        type(grid_class),     intent(IN), optional :: grid

        td%par%path = trim(path)

        if (present(vars)) then
            allocate(td%par%vars(size(vars)))
            td%par%vars = vars
        else
            allocate(td%par%vars(size(vars_default)))
            td%par%vars = vars_default
        end if

        allocate(td%par%names(size(td%par%vars)))
        td%par%names = ""
        if (present(names)) then
            if (size(names) .ne. size(td%par%vars)) then
                call topodata_error("topodata_init_arg", "vars and names differ in length.")
            end if
            td%par%names = names
        end if
        where (len_trim(td%par%names) .eq. 0) td%par%names = td%par%vars

        td%par%remap = "con"
        if (present(remap)) td%par%remap = trim(remap)

        call topodata_init_data(td, grid)

        return

    end subroutine topodata_init_arg

    subroutine topodata_init_data(td, grid)
        ! Read the variables given by td%par (remapping onto grid if needed).

        implicit none

        type(topodata_class), intent(INOUT) :: td
        type(grid_class),     intent(IN), optional :: grid

        character(len=1024) :: fname
        type(grid_class) :: grid_src
        logical :: with_int
        integer :: k
        integer, allocatable :: dims(:)

        fname = td%par%path

        do k = 1, size(td%par%vars)
            if (.not. any(vars_all .eq. td%par%vars(k))) then
                call topodata_error("topodata_init_data", "unknown variable.", &
                    "variable  = "//trim(td%par%vars(k))//new_line("a")// &
                    "variables = z_bed z_srf H_ice z_bed_sd f_ocn f_land f_grnd f_flt mask src_id")
            end if
            if (.not. nc_exists_var(fname, td%par%names(k))) then
                call topodata_error("topodata_init_data", "variable not in the file.", &
                    "file     = "//trim(fname)//new_line("a")//"variable = "//trim(td%par%names(k)))
            end if
        end do

        td%par%grid_src = ""
        if (nc_exists_attr(fname, "grid_name")) td%par%grid_src = fesmdata_grid_name(fname)
        td%par%product  = ""
        td%par%sources  = ""
        if (nc_exists_attr(fname, "product")) call nc_read_attr(fname, "product", td%par%product)
        if (nc_exists_attr(fname, "sources")) call nc_read_attr(fname, "sources", td%par%sources)

        ! Maps onto the target grid if it differs from the grid of the file
        td%par%grid_tgt = td%par%grid_src
        td%remap = .false.
        if (present(grid)) then
            td%par%grid_tgt = grid%name
            if (len_trim(td%par%grid_src) .gt. 0 .and. &
                trim(grid%name) .ne. trim(td%par%grid_src)) then
                with_int = any(td%par%vars .eq. "mask") .or. any(td%par%vars .eq. "src_id")
                call fesmdata_grid_read(grid_src, fname)
                call map_init(td%map, grid_src, grid, method=trim(td%par%remap), fldr="maps")
                if (with_int .and. trim(td%par%remap) .ne. "nn") then
                    call map_init(td%map_nn, grid_src, grid, method="nn", fldr="maps")
                end if
                td%remap = .true.
            end if
        end if

        if (td%remap) then
            td%nx = grid%G%nx
            td%ny = grid%G%ny
        else
            call nc_dims(fname, td%par%names(1), dims=dims)
            td%nx = dims(1)
            td%ny = dims(2)
            if (present(grid)) then
                if (grid%G%nx .ne. td%nx .or. grid%G%ny .ne. td%ny) then
                    if (len_trim(td%par%grid_src) .gt. 0) then
                        call topodata_error("topodata_init_data", &
                            "the target grid has the name of the grid of the file, but another size.", &
                            "file = "//trim(fname)//new_line("a")//"grid = "//trim(grid%name))
                    else
                        call topodata_error("topodata_init_data", &
                            "the file (without grid_name) differs in size from the target grid.", &
                            "file = "//trim(fname)//new_line("a")//"grid = "//trim(grid%name))
                    end if
                end if
            end if
        end if

        do k = 1, size(td%par%vars)
            associate(name => td%par%names(k))
                select case(trim(td%par%vars(k)))
                    case("z_bed")
                        call fesmdata_read_field(fname, trim(name), td%z_bed,    td%nx, td%ny, td%remap, td%map)
                    case("z_srf")
                        call fesmdata_read_field(fname, trim(name), td%z_srf,    td%nx, td%ny, td%remap, td%map)
                    case("H_ice")
                        call fesmdata_read_field(fname, trim(name), td%H_ice,    td%nx, td%ny, td%remap, td%map)
                    case("z_bed_sd")
                        call fesmdata_read_field(fname, trim(name), td%z_bed_sd, td%nx, td%ny, td%remap, td%map)
                    case("f_ocn")
                        call fesmdata_read_field(fname, trim(name), td%f_ocn,    td%nx, td%ny, td%remap, td%map)
                    case("f_land")
                        call fesmdata_read_field(fname, trim(name), td%f_land,   td%nx, td%ny, td%remap, td%map)
                    case("f_grnd")
                        call fesmdata_read_field(fname, trim(name), td%f_grnd,   td%nx, td%ny, td%remap, td%map)
                    case("f_flt")
                        call fesmdata_read_field(fname, trim(name), td%f_flt,    td%nx, td%ny, td%remap, td%map)
                    case("mask")
                        call read_field_nn(td, fname, trim(name), td%mask)
                        call flag_table_read(td%tab_mask, fname, trim(name))
                    case("src_id")
                        call read_field_nn(td, fname, trim(name), td%src_id)
                        call flag_table_read(td%tab_src, fname, trim(name))
                end select
            end associate
        end do

        return

    end subroutine topodata_init_data

    subroutine read_field_nn(td, fname, varname, var)
        ! An integer field, remapped by nearest neighbour if needed.

        implicit none

        type(topodata_class), intent(IN)  :: td
        character(len=*),     intent(IN)  :: fname, varname
        integer, allocatable, intent(OUT) :: var(:,:)

        if (trim(td%par%remap) .eq. "nn") then
            call fesmdata_read_field(fname, varname, var, td%nx, td%ny, td%remap, td%map,    -1)
        else
            call fesmdata_read_field(fname, varname, var, td%nx, td%ny, td%remap, td%map_nn, -1)
        end if

        return

    end subroutine read_field_nn

    subroutine topodata_end(td)

        implicit none

        type(topodata_class), intent(INOUT) :: td

        type(topodata_class) :: td0

        td = td0

        return

    end subroutine topodata_end

    subroutine topodata_par_load(par, filename, group, domain, grid_name, subs, verbose)

        implicit none

        type(topodata_param_class), intent(OUT) :: par
        character(len=*), intent(IN) :: filename
        character(len=*), intent(IN) :: group
        character(len=*), intent(IN), optional :: domain
        character(len=*), intent(IN), optional :: grid_name
        character(len=*), intent(IN), optional :: subs(:,:)
        logical,          intent(IN), optional :: verbose

        character(len=56) :: vars(n_vars_max), names(n_vars_max)
        logical :: print_summary
        integer :: k

        print_summary = .true.
        if (present(verbose)) print_summary = verbose

        call nml_read(filename, group, "path", par%path)
        call fesmdata_parse_path(par%path, domain, grid_name, subs)

        vars = ""
        vars(1:size(vars_default)) = vars_default
        if (nml_has_param(filename, group, "vars")) then
            vars = ""
            call nml_read(filename, group, "vars", vars)
        end if
        allocate(par%vars(count(len_trim(vars) .gt. 0)))
        par%vars = pack(vars, len_trim(vars) .gt. 0)

        ! Names in the file, by position in vars ("" = as vars)
        names = ""
        if (nml_has_param(filename, group, "names")) then
            call nml_read(filename, group, "names", names)
        end if
        if (any(len_trim(vars) .eq. 0 .and. len_trim(names) .gt. 0)) then
            call topodata_error("topodata_par_load", "more names than vars.", &
                "file = "//trim(filename)//new_line("a")//"group = "//trim(group))
        end if
        allocate(par%names(size(par%vars)))
        par%names = pack(names, len_trim(vars) .gt. 0)
        where (len_trim(par%names) .eq. 0) par%names = par%vars

        par%remap = "con"
        if (nml_has_param(filename, group, "remap")) then
            call nml_read(filename, group, "remap", par%remap)
        end if

        if (print_summary) then
            write(*,*) "Loading: ", trim(filename), ":: ", trim(group)
            write(*,*) "path  = ", trim(par%path)
            write(*,*) "vars  = ", (trim(par%vars(k))//" ", k=1,size(par%vars))
            write(*,*) "names = ", (trim(par%names(k))//" ", k=1,size(par%names))
            write(*,*) "remap = ", trim(par%remap)
        end if

        return

    end subroutine topodata_par_load

    subroutine topodata_error(proc, msg, detail)
        ! Abort with a framed message on error_unit.

        implicit none

        character(len=*), intent(IN)           :: proc
        character(len=*), intent(IN)           :: msg
        character(len=*), intent(IN), optional :: detail

        integer :: p0, p1

        write(error_unit,"(a)") ""
        write(error_unit,"(a)") "topodata:: error in "//trim(proc)
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

        write(error_unit,"(a)") "  stopped by topodata."
        error stop 1

    end subroutine topodata_error

end module topodata
