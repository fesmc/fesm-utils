module esm

    use, intrinsic :: iso_fortran_env, only : input_unit, output_unit, error_unit

    use nml
    use ncio
    use varslice

    implicit none

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)

    ! Working precision
    integer,  parameter :: wp = sp

    ! Default missing value
    real(wp), parameter :: mv = -9999.0_wp

    type esm_state_class

        character(len=64)    :: name       ! label, e.g. "lgm", "pd", "hist", "proj"

        ! Segment-mode timing metadata
        real(wp)             :: t_start    ! start of valid time window [yr]
        real(wp)             :: t_end      ! end   of valid time window [yr]

        ! Interp-mode metadata
        real(wp)             :: anchor     ! alpha value at which this state is "pure"

        ! Atmospheric fields
        type(varslice_class) :: ts          ! surface temperature
        type(varslice_class) :: pr          ! precipitation

        ! Oceanic fields
        type(varslice_class) :: to          ! ocean temperature
        type(varslice_class) :: so          ! ocean salinity

        ! Surface / topographic fields
        type(varslice_class) :: zs          ! surface elevation

    end type esm_state_class

    type esm_forcing_class

        ! Experiment metadata
        character(len=256)   :: gcm
        character(len=256)   :: scenario
        character(len=256)   :: experiment
        character(len=256)   :: domain
        character(len=256)   :: grid_name
        character(len=256)   :: ctrl_run_type
        real(wp)             :: lapse(2)
        real(wp)             :: beta_p
        real(wp)             :: f_ocn
        real(wp)             :: f_polar
        real(wp)             :: dT_lim
        character(len=256)   :: grid_src

        ! Update mode
        ! "segment" : find the state whose [t_start,t_end] contains current time
        ! "interp"  : blend the two states bracketing a scalar alpha index
        character(len=16)    :: forcing_mode

        integer                             :: n_clim
        character(len=56), allocatable      :: clim_names
        type(esm_state_class), allocatable  :: clim(:)   ! size n_clim

        ! General fields
        type(varslice_class) :: basins

    end type esm_forcing_class

    type esm_class
        ! Diagnostic / output fields

        integer :: nx
        integer :: ny 
        integer :: nm

        ! Atmosphere (monthly)
        real(wp), allocatable :: ts(:,:,:)        ! near-surface temperature [K]
        real(wp), allocatable :: pr(:,:,:)        ! precipitation [mm/yr]

        ! Anomalies
        real(wp), allocatable :: dts(:,:,:)       ! temperature anomaly [K]
        real(wp), allocatable :: dpr(:,:,:)       ! precipitation relative anomaly [-]
        real(wp), allocatable :: dts_var(:,:,:)   ! temperature variability anomaly [K]
        real(wp), allocatable :: dpr_var(:,:,:)   ! precipitation variability anomaly [-]

        ! Ocean
        real(wp), allocatable :: dto(:,:)         ! ocean temperature anomaly [K]
        real(wp), allocatable :: dso(:,:)         ! ocean salinity anomaly [psu]
        real(wp), allocatable :: dto_var(:,:)
        real(wp), allocatable :: dso_var(:,:)

        ! Mean fields
        real(wp), allocatable :: ts_sum(:,:)
        real(wp), allocatable :: ts_ann(:,:)
        real(wp), allocatable :: pr_ann(:,:)

        ! Store all internal data related to building the current forcing here
        type(esm_forcing_class) :: f

    end type

    ! Ancillary types
    type esm_ice_var_class
        character(len=56)  :: name
        character(len=128) :: long_name
        character(len=12)  :: var_type
        character(len=128) :: standard_name
        character(len=128) :: units_in
        character(len=128) :: units_out
        real(wp) :: unit_scale
        real(wp) :: unit_offset
    end type

    type esm_experiment_class
        character(len=56)   :: expname
        character(len=56)   :: group
        character(len=56)   :: model
        character(len=256)  :: experiment
        character(len=256)  :: file_suffix
    end type

    type esm_ice_class
        type(esm_ice_var_class), allocatable :: vars(:)
    end type

    private

    public :: esm_state_class
    public :: esm_forcing_class
    public :: esm_class

    public :: esm_experiment_class
    public :: esm_ice_class

contains


    ! =========================================================================
    ! varslice wrapper routines
    ! =========================================================================
    subroutine varslice_init_nml_esm(vs, filename, group, domain, grid_name, &
                                     gcm, scenario, verbose)

        implicit none

        type(varslice_class), intent(INOUT) :: vs
        character(len=*),     intent(IN)    :: filename
        character(len=*),     intent(IN)    :: group
        character(len=*),     intent(IN), optional :: domain
        character(len=*),     intent(IN), optional :: grid_name
        character(len=*),     intent(IN)    :: gcm
        character(len=*),     intent(IN)    :: scenario
        logical,              intent(IN), optional :: verbose

        call varslice_par_load_esm(vs%par, filename, group, domain, grid_name, &
                                   gcm, scenario, verbose)
        call varslice_init_data(vs)

        return

    end subroutine varslice_init_nml_esm

    subroutine varslice_par_load_esm(par, filename, group, domain, grid_name, &
                                     gcm, scenario, verbose)

        implicit none

        type(varslice_param_class), intent(OUT) :: par
        character(len=*), intent(IN) :: filename
        character(len=*), intent(IN) :: group
        character(len=*), intent(IN), optional :: domain
        character(len=*), intent(IN), optional :: grid_name
        character(len=*), intent(IN) :: gcm
        character(len=*), intent(IN) :: scenario
        logical, optional :: verbose

        logical :: init_pars, print_summary
        integer :: i

        character(len=56) :: meta(3)
        real(wp) :: scaling(2)

        init_pars     = .FALSE.
        print_summary = .TRUE.
        if (present(verbose)) print_summary = verbose

        call nml_read(filename, group, "filename",    par%filename,     init=init_pars)
        call nml_read(filename, group, "meta",        meta,             init=init_pars)
        call nml_read(filename, group, "scaling",     scaling,          init=init_pars)
        call nml_read(filename, group, "time_info",   par%time_par,     init=init_pars)

        par%name = meta(1)
        par%units_in = meta(2)
        par%units_out = meta(3)

        par%unit_scale = scaling(1)
        par%unit_offset = scaling(2)

        if (par%time_par(3) .eq. 0.0) then
            par%with_time = .FALSE.
        end if

        call parse_path(par%filename, domain, grid_name)
        call parse_path_esm(par%filename, gcm, scenario)
        call get_matching_files(par%filenames, par%filename)

        if (par%time_par(3) .eq. 0.0) par%time_par(2) = par%time_par(1)

        if (par%time_par(4) .gt. 1.0) then
            par%with_time_sub = .TRUE.
        else
            par%with_time_sub = .FALSE.
        end if

        if (print_summary) then
            write(*,*) "Loading: ", trim(filename), ":: ", trim(group)
            write(*,*) "filename      = ", trim(par%filename)
            if (size(par%filenames,1) .gt. 1) then
                write(*,*) "filenames     = "
                do i = 1, size(par%filenames,1)
                    write(*,*) "                  ", trim(par%filenames(i))
                end do
            end if
            write(*,*) "name          = ", trim(par%name)
            write(*,*) "units_in      = ", trim(par%units_in)
            write(*,*) "units_out     = ", trim(par%units_out)
            write(*,*) "unit_scale    = ", par%unit_scale
            write(*,*) "unit_offset   = ", par%unit_offset
            write(*,*) "with_time     = ", par%with_time
            write(*,*) "with_time_sub = ", par%with_time_sub
            if (par%with_time) write(*,*) "time_par      = ", par%time_par
        end if

        return

    end subroutine varslice_par_load_esm

    subroutine parse_path_esm(path, gcm, scenario)

        implicit none

        character(len=*), intent(INOUT) :: path
        character(len=*), intent(IN)    :: gcm, scenario

        call nml_replace(path, "{gcm}",      trim(gcm))
        call nml_replace(path, "{scenario}", trim(scenario))

        return

    end subroutine parse_path_esm

    subroutine parse_path(path, domain, grid_name)

        implicit none

        character(len=*), intent(INOUT) :: path
        character(len=*), intent(IN)    :: domain, grid_name

        call nml_replace(path, "{domain}",    trim(domain))
        call nml_replace(path, "{grid_name}", trim(grid_name))

        return

    end subroutine parse_path

    subroutine parse_clim_string(input_str, name, var_names, num_vars)
        
        implicit none
        
        character(len=32), intent(OUT) :: name
        character(len=10), intent(OUT) :: var_names(:)
        integer,           intent(OUT) :: num_vars
        character(len=*),  intent(IN)  :: input_str
        
        ! Internal variables
        integer :: bracket_open, bracket_close
        integer :: i, start_pos, comma_pos
        character(len=len(input_str)) :: var_part
        
        ! Initialize
        name = ""
        var_names = ""
        num_vars = 0
        
        bracket_open  = index(input_str, "[")
        bracket_close = index(input_str, "]")
        
        if (bracket_open == 0) then
            ! Case 1: No brackets - Use default 5 variables
            name = trim(adjustl(input_str))
            var_names(1:5) = ["ts", "pr", "zs", "to", "so"]
            num_vars = 5
        else
            ! Case 2: Brackets found - Parse name and list
            name = trim(adjustl(input_str(1:bracket_open-1)))
            
            ! Extract the string inside the brackets
            var_part = input_str(bracket_open+1 : bracket_close-1)
            
            start_pos = 1
            do i = 1, size(var_names)
                comma_pos = index(var_part(start_pos:), ",")
                
                if (comma_pos == 0) then
                    ! Last (or only) variable in the list
                    var_names(i) = trim(adjustl(var_part(start_pos:)))
                    num_vars = i
                    exit
                else
                    ! Extract variable before the comma
                    var_names(i) = trim(adjustl(var_part(start_pos : start_pos + comma_pos - 2)))
                    start_pos = start_pos + comma_pos
                end if
            end do
        end if
        
    end subroutine parse_clim_string

    ! =========================================================================
    ! Allocation helpers
    ! =========================================================================
    subroutine esm_allocate(esm, nx, ny)

        implicit none

        type(esm_class), intent(INOUT) :: esm
        integer, intent(IN) :: nx, ny

        call esm_deallocate(esm)

        esm%nx = nx
        esm%ny = ny
        esm%nm = 12

        allocate(esm%ts(nx,ny,12))
        allocate(esm%pr(nx,ny,12))
        allocate(esm%dts(nx,ny,12))
        allocate(esm%dpr(nx,ny,12))
        allocate(esm%dto(nx,ny))
        allocate(esm%dso(nx,ny))
        allocate(esm%dts_var(nx,ny,12))
        allocate(esm%dpr_var(nx,ny,12))
        allocate(esm%dto_var(nx,ny))
        allocate(esm%dso_var(nx,ny))
        allocate(esm%ts_ann(nx,ny))
        allocate(esm%ts_sum(nx,ny))
        allocate(esm%pr_ann(nx,ny))

        return

    end subroutine esm_allocate

    subroutine esm_deallocate(esm)

        implicit none

        type(esm_class), intent(INOUT) :: esm

        if (allocated(esm%ts))     deallocate(esm%ts)
        if (allocated(esm%pr))      deallocate(esm%pr)
        if (allocated(esm%dts))     deallocate(esm%dts)
        if (allocated(esm%dpr))     deallocate(esm%dpr)
        if (allocated(esm%dto))     deallocate(esm%dto)
        if (allocated(esm%dso))     deallocate(esm%dso)
        if (allocated(esm%dts_var)) deallocate(esm%dts_var)
        if (allocated(esm%dpr_var)) deallocate(esm%dpr_var)
        if (allocated(esm%dto_var)) deallocate(esm%dto_var)
        if (allocated(esm%dso_var)) deallocate(esm%dso_var)
        if (allocated(esm%ts_sum)) deallocate(esm%ts_sum)
        if (allocated(esm%ts_ann)) deallocate(esm%ts_ann)
        if (allocated(esm%pr_ann))  deallocate(esm%pr_ann)

        return

    end subroutine esm_deallocate

    ! ===== ROUTINES from marine_shelf ========
    ! These were brought in here and marine_shelf dependency removed, 
    ! but eventually may need to be organized better (ajr, 2026-03-02)

    subroutine interp_shelf(out2D,in3d,H_ice,z_bed,f_grnd,z_sl,depth, &
        interp_depth, interp_method, depth_const, depth_range )
        ! Calculate 2D fields from 3D ocean fields representative

        implicit none

        real(wp), intent(INOUT) :: out2d(:,:)
        real(wp), intent(IN) :: in3d(:,:,:)
        real(wp), intent(IN) :: H_ice(:,:)
        real(wp), intent(IN) :: z_bed(:,:)
        real(wp), intent(IN) :: f_grnd(:,:)
        real(wp), intent(IN) :: z_sl(:,:)
        real(wp), intent(IN) :: depth(:)
        character(len=*), intent(IN) :: interp_depth
        character(len=*), intent(IN) :: interp_method
        real(wp), intent(IN) :: depth_const
        real(wp), intent(IN), optional :: depth_range(2)

        ! Local variables
        integer :: i, j, nx, ny, nz
        real(wp), allocatable :: depth_shlf
        real(wp), allocatable :: wt_shlf(:)

        real(wp), parameter :: rho_ice_sw = 917.0 / 1028.0         ! ajr, should define as a physical constant better, or pass in as an argument

        nx = size(H_ice,1)
        ny = size(H_ice,2)
        nz = size(depth,1)

        ! Allocate objects
        allocate(wt_shlf(nz))
        wt_shlf = 0.0
        out2d   = 0.0

        ! Loop over domain and update variables at each point (vertical interpolation)
        do j = 1, ny
        do i = 1, nx

            ! 1. Calculate the depth of the current shelf base

            select case(trim(interp_depth))

                case("shlf")
                ! Assign the depth to the shelf depth

                if(H_ice(i,j) .gt. 0.0 .and. f_grnd(i,j) .lt. 1.0) then
                ! Floating ice, depth == z_ice_base
                depth_shlf = H_ice(i,j)*rho_ice_sw
                else if(H_ice(i,j) .gt. 0.0 .and. f_grnd(i,j) .eq. 1.0) then
                ! Grounded ice, depth == H_ocn = z_sl-z_bed
                depth_shlf = z_sl(i,j) - z_bed(i,j)
                else
                ! Open ocean, depth == constant value, eg 2000 m.
                depth_shlf = depth_const
                end if

                case("bed")
                ! Assign the depth corresponding to the bedrock

                depth_shlf = z_sl(i,j) - z_bed(i,j)

                case("const")

                depth_shlf = depth_const

                case DEFAULT

                write(*,*) "interp_shelf:: Error: interp_depth method not recognized."
                write(*,*) "interp_depth = ", trim(interp_depth)

            end select

            ! 2. Calculate weighting function for vertical depths ===========================

            select case(trim(interp_method))

                case("mean")
                ! Equal weighting of layers within a specified depth range

                call calc_shelf_variable_mean(wt_shlf,depth,depth_range)

                case("layer")
                ! All weight given to the nearest layer to depth of shelf

                call calc_shelf_variable_layer(wt_shlf,depth,depth_shlf)

                case("interp")
                ! Interpolation weights from the two nearest layers to depth of shelf

                call  calc_shelf_variable_depth(wt_shlf,depth,depth_shlf)

                case DEFAULT
                write(*,*) "interp_shelf:: error: interp_method not recognized: ", interp_method
                write(*,*) "Must be one of [mean, layer, interp]"
                stop

            end select

            ! Normalize weighting function
            if (sum(wt_shlf) .gt. 0.0_wp) then
                wt_shlf = wt_shlf / sum(wt_shlf)
            else
                write(*,*) "marshelf_interp_shelf:: Error: weighting should be > 0."
                stop
            end if

            ! 3. Calculate water properties at depths of interest ============================
            out2d(i,j)  = sum(in3d(i,j,:) * wt_shlf)

        end do
        end do 

        return

    end subroutine interp_shelf

    subroutine calc_shelf_variable_mean(wt_shlf,depth,depth_range)
        ! Calculate average variable value for a given range of depths
        ! at a specific point (x,y). 

        implicit none 
        
        real(wp), intent(OUT)   :: wt_shlf(:)
        real(wp), intent(IN)    :: depth(:)
        real(wp), intent(IN), optional :: depth_range(:)

        ! Local variables 
        integer :: k0, k1 
        
        ! depth_range was set to optional, to facilitate flexibility in another routine.
        ! Error tracking is performed here to make sure that depth_range is actually provided as an argument.
        if (.not. present(depth_range)) then
            write(*,*) "calc_shelf_temperature_mean:: depth_range is needed, but has not been provided as an argument."
            stop
        end if

        ! Note: this requires that k1 > k0, and it should be
        ! weighted by the thickness of the layers (to do!)
        ! Note: depth is z-coordinate, ie positive below sea level  
        k0 = minloc(abs(depth-depth_range(1)),dim=1)
        k1 = minloc(abs(depth-depth_range(2)),dim=1)

        if (k1 < k0) then 
            write(*,*) "calc_shelf_temperature_mean:: error in depth_range calculation."
            write(*,*) "depth_min, depth_max: ", depth_range 
            write(*,*) "indices(k0,k1): ", k0, k1
            write(*,*) "depths(k0,k1):  ", depth(k0), depth(k1) 
            stop 
        end if 

        ! Get index weights to produce mean water temperature for these depths
        wt_shlf        = 0.0 
        wt_shlf(k0:k1) = 1.0 

        return 

    end subroutine calc_shelf_variable_mean

    subroutine calc_shelf_variable_layer(wt_shlf,depth,depth_shlf)
        ! Calculates the water temperature at the depth of the ice shelf
        ! It assigns the temperature of the nearest layer
        
        implicit none

        real(wp), intent(OUT)   :: wt_shlf(:)
        real(wp), intent(IN)    :: depth(:)
        real(wp), intent(IN)    :: depth_shlf

        ! Local variables
        integer :: k0, i, j

        ! Determine layer closest to target depth 
        k0 = minloc(abs(depth-depth_shlf),dim=1)

        ! Assign weighting function 
        wt_shlf     = 0.0 
        wt_shlf(k0) = 1.0 

        return

    end subroutine calc_shelf_variable_layer

    subroutine calc_shelf_variable_depth(wt_shlf,depth,depth_shlf)
        ! Calculates the water temperature from linear interpolation of vertical profile
        
        implicit none

        real(wp), intent(OUT)   :: wt_shlf(:)
        real(wp), intent(IN)    :: depth(:)
        real(wp), intent(IN)    :: depth_shlf

        ! Local variables 
        integer  :: j, n 
        real(wp) :: alpha 

        n    = size(depth) 

        do j = 1, n 
            if (depth(j) .ge. depth_shlf) exit 
        end do

        if (j .eq. 1) then 
            wt_shlf    = 0.0 
            wt_shlf(1) = 1.0 
        else if (j .eq. n+1) then 
            wt_shlf    = 0.0 
            wt_shlf(n) = 1.0 
        else 
            alpha = (depth_shlf - depth(j-1)) / (depth(j) - depth(j-1))
            wt_shlf      = 0.0 
            wt_shlf(j-1) = 1.0-alpha 
            wt_shlf(j)   = alpha
        end if  

        return

    end subroutine calc_shelf_variable_depth

end module esm
