module varslice

    use, intrinsic :: iso_fortran_env, only : input_unit, output_unit, error_unit
    
    use precision
    use constants, only: mv, TOL
    use ncio
    use nml
    use mapping, only : map_class, map_field

    implicit none 

    logical, parameter :: verbose = .FALSE.

    ! Values with magnitude >= var_lim_crazy are treated as unflagged
    ! fill/garbage and reset to the missing value. This catches NetCDF
    ! default fill values (~9.97e36) that were not tagged via missing_value.
    ! Note: legitimate data with magnitude above this limit would also be
    ! discarded — raise it if a variable can genuinely exceed 1e10.
    real(wp), parameter :: var_lim_crazy = 1.0e10_wp

    type varslice_param_class

        character(len=1024) :: filename
        character(len=1024), allocatable :: filenames(:)
        character(len=56)   :: name 
        character(len=56)   :: units_in
        character(len=56)   :: units_out
        real(wp) :: unit_scale 
        real(wp) :: unit_offset
        logical  :: with_time 
        logical  :: with_time_sub

        ! Internal parameters
        integer  :: ndim 
        real(dp) :: time_par(4) 

    end type

    type varslice_class 

        type(varslice_param_class) :: par 

        ! Parameters defined during update call
        real(dp)          :: time_range(2)
        character(len=56) :: slice_method 
        integer           :: range_rep 
        
        ! Variable information
        integer,  allocatable :: dim(:)  
        real(wp), allocatable :: x(:) 
        real(wp), allocatable :: y(:)
        real(wp), allocatable :: z(:)  
        real(dp), allocatable :: time(:)
        real(dp), allocatable :: time_sub(:)
        integer, allocatable  :: idx(:)
        integer, allocatable  :: nt_files(:)   ! per-file time length (cached at init)

        real(wp), allocatable :: var(:,:,:,:)

    contains

        ! Convenience accessors returning a natural-rank copy of the data,
        ! indexed by the trailing time/level position (default 1). The classic
        ! vs%var(:,:,k,1) indexing continues to work unchanged.
        !   v1(k) => var(:,k,1,1)   v2(k) => var(:,:,k,1)   v3(k) => var(:,:,:,k)
        procedure :: v1 => varslice_v1
        procedure :: v2 => varslice_v2
        procedure :: v3 => varslice_v3

    end type

    interface axis_init
        module procedure axis_init_sp
        module procedure axis_init_dp
    end interface

    ! Render scalars / 1-D arrays into the `detail` string of varslice_error.
    interface to_str
        module procedure to_str_int, to_str_sp, to_str_dp, to_str_log
        module procedure to_str_int1, to_str_sp1, to_str_dp1
    end interface

    private
    public :: varslice_param_class
    public :: varslice_class
    public :: varslice_update
    public :: varslice_init_nml 
    public :: varslice_init_arg
    public :: varslice_init_data
    public :: varslice_end 

    public :: varslice_map_to_grid

    public :: get_matching_files
    
    public :: print_var_range
    
contains

    subroutine varslice_map_to_grid(vs_tgt,vs_src,map,mask_out,stat,method,reset,missing_value, &
                                    mask_pack,fill_method,filt_method,filt_par,verbose)
        ! Given a input vs object on a source grid src,
        ! map the variable to a target grid tgt and store
        ! in new vs object.

        implicit none

        type(varslice_class),  intent(INOUT) :: vs_tgt
        type(varslice_class),  intent(IN)    :: vs_src
        type(map_class),       intent(IN)    :: map
        integer,                intent(OUT), optional :: mask_out(:,:)   ! Mask showing where interpolation was done
        character(len=*),       intent(IN),  optional :: stat            ! Aggregation: mean (default), count, stdev
        character(len=*),       intent(IN),  optional :: method          ! Interpolation-kernel override (distance maps)
        logical,                intent(IN),  optional :: reset           ! Reset var_tgt initially to missing_value?
        real(wp),               intent(IN),  optional :: missing_value   ! Points not included in mapping
        logical,                intent(IN),  optional :: mask_pack(:,:)  ! Mask for where to interpolate
        character(len=*),       intent(IN),  optional :: fill_method     ! Method to fill in remaining missing values
        character(len=*),       intent(IN),  optional :: filt_method     ! Method to use for filtering
        real(wp),               intent(IN),  optional :: filt_par(:)     ! gaussian=[sigma,dx]; poisson=[tol]
        logical,                intent(IN),  optional :: verbose         ! Print information
        
        ! Local variables
        integer :: nx, ny, nz, nt
        integer :: i, j, k, t
        integer :: ndim
        real(wp) :: miss                        ! missing value used for mapping
        logical, allocatable :: mask_int(:,:)   ! per-slice interpolation mask (logical form of mask_out)

        ! Missing value: caller override, else the package default
        miss = mv
        if (present(missing_value)) miss = missing_value

        ! Consistency check: is this at least a 2D field?
        ndim = size(vs_src%dim,1)

        if (ndim .lt. 2 .or. (vs_src%par%with_time .and. ndim .lt. 3)) then
            call varslice_error("varslice_map_to_grid", &
                "mapping can only be done for fields with two spatial dimensions.", &
                "name = "//trim(vs_src%par%name)//new_line("a")// &
                "ndim = "//to_str(ndim)//new_line("a")// &
                "dim  = "//to_str(vs_src%dim))
        end if

        ! Determine size of target var
        nx = map%G%nx
        ny = map%G%ny
        nz = size(vs_src%var,3)         ! z should keep the same dimension
        nt = size(vs_src%var,4)         ! time should keep the same dimension
        
        ! Initialize meta information by copying entire vs object
        vs_tgt = vs_src

        ! Re-allocate target varslice object fields as needed
        deallocate(vs_tgt%x)
        deallocate(vs_tgt%y)
        deallocate(vs_tgt%var)
        
        allocate(vs_tgt%x(nx))
        allocate(vs_tgt%y(ny))
        allocate(vs_tgt%var(nx,ny,nz,nt))

        ! Per-slice interpolation mask (logical); converted to mask_out below
        allocate(mask_int(nx,ny))

        ! Loop over dimensions and remap variable
        do t = 1, nt
        do k = 1, nz
            call map_field(map,trim(vs_tgt%par%name),vs_src%var(:,:,k,t),vs_tgt%var(:,:,k,t), &
                            stat=stat,method=method,missing_value=miss,reset=reset,mask_pack=mask_pack, &
                            fill_method=fill_method,filt_method=filt_method,filt_par=filt_par, &
                            verbose=verbose,mask2=mask_int)
        end do
        end do

        ! Report where interpolation was done (last slice), as an integer 0/1 mask
        if (present(mask_out)) mask_out = merge(1, 0, mask_int)

        ! Done! target field should now be available.

        return

    end subroutine varslice_map_to_grid

    subroutine varslice_update(vs,time,method,fill,rep,print_summary)
        ! Routine to update transient climate forcing to match 
        ! current `time`. 

        ! time = a specific time or a range of times

        ! method = ["exact","range","interp","extrap","range_min","range_mean","range_max"]
        ! difference time-slicing methods applied to the data read from the file.
        ! Check routine `get_indices()` to see how data is loaded for each slice_method.
        ! ["exact","range"] produce data identical to the input file, where
        ! "exact" loads the data from the file for the time index matching the desired time
        ! exactly - if it is not available, missing values are returned.
        ! "range" loads all the data within a certain time range.
        ! ["interp","extrap"] both return a time slice by interpolating to the
        ! desired time from the two closest bracketing timesteps available in the file. "extrap"
        ! allows for setting the time slice equal to the first or last timestep available, if
        ! the desired time is out of bounds, while "interp" returns missing values in this case. 
        ! "range_*" methods return one time slice with the method applied to the data within the 
        ! range given by `time`. 

        ! rep: frequency to apply slice_method over time. If rep=1, then calculation (mean/sd/etc)
        ! will be applied to each time index, returning a field with no time dimension. 
        ! If rep=12, then calculation will be applied to every 12th index, resulting in 
        ! 12 values along dimension.

        ! fill: method to handle missing values. By default, no treatment and missing values
        ! are included in returned fields. TODO

        implicit none 

        type(varslice_class), target, intent(INOUT) :: vs
        real(wp),         optional, intent(IN)    :: time(:)        ! [yr] Current time, or time range
        character(len=*), optional, intent(IN)    :: method         ! slice_method (only if with_time==True)
        character(len=*), optional, intent(IN)    :: fill           ! none, min, max, mean (how to fill in missing values)
        integer,          optional, intent(IN)    :: rep            ! Only if with_time==True
        logical,          optional, intent(IN)    :: print_summary  ! Print summary of updated variable

        ! Local variables 
        integer :: k, k0, k1, nt
        integer :: nt_tot, nt_rep, nt_major, nt_out
        integer :: n1, n2, i, j, l    
        type(varslice_param_class) :: par 
        logical  :: with_time
        logical  :: with_time_sub
        real(wp) :: time_range(2) 
        real(wp), allocatable :: time_wt(:)
        character(len=56) :: slice_method
        character(len=56) :: fill_method
        character(len=56) :: vec_method 
        integer  :: range_rep
        integer,  allocatable :: kk(:)
        real(wp), allocatable, target :: var(:,:,:,:)
        integer  :: s, nsp
        real(wp), pointer, contiguous :: src2d(:,:), dst2d(:,:)

        ! Define shortcuts
        par = vs%par 
        with_time = par%with_time 

        if (present(rep)) then
            range_rep     = rep
        else
            range_rep     = 1 
        end if

        if (range_rep .gt. 1) then
            with_time_sub = .TRUE.
        else
            with_time_sub = .FALSE.
        end if

        if (with_time) then 

            if (.not. present(time)) then
                call varslice_error("varslice_update", &
                    "current time or time range must be given as an argument (1D array).")
            end if

            ! Consistency check
            if (size(time,1) .eq. 2) then
                if (time(2) .lt. time(1)) then
                    call varslice_error("varslice_update", &
                        "time(2) should be >= time(1).", &
                        "time = "//to_str(time))
                end if
            end if

        end if 

        slice_method = "exact"
        if (present(method)) slice_method = trim(method)

        fill_method = "none"
        if (present(fill)) fill_method = trim(fill) 

        if (trim(fill_method) .ne. "none") then
            call varslice_error("varslice_update", &
                "fill methods have not yet been implemented. Set fill_method='none' for now.", &
                "fill_method = "//trim(fill_method))
        end if
        
        if (present(time)) then

            if (size(time) .eq. 2) then 
                time_range = time
            else 
                time_range(1:2) = time(1) 
            end if 

        else 

            time_range = vs%time_range 

        end if 

        
        if ( with_time .and. trim(slice_method) .eq. trim(vs%slice_method) &
                .and. range_rep .eq. vs%range_rep &
                .and. vs%time_range(1) .eq. time_range(1) &
                .and. vs%time_range(2) .eq. time_range(2) ) then 

            ! Do nothing, the varslice object is already up to date 
            ! fo the current time and method. 

        else 

            ! Set parameters in varslice object 
            vs%slice_method = trim(slice_method)
            vs%range_rep    = range_rep 
            vs%time_range   = time_range 

            ! 2. Read variable and convert units as needed

            ! Make sure var variable is deallocated and ready to be modified
            if (allocated(vs%var)) deallocate(vs%var)
            
            if (.not. with_time) then 
                ! Handle cases that do not have time dimension (simpler)

                select case(par%ndim)

                    case(1)

                        ! Allocate var to the right size 
                        allocate(vs%var(vs%dim(1),1,1,1))

                        ! 1D variable
                        call nc_read(par%filename,par%name,vs%var(:,1,1,1),missing_value=mv)

                    case(2)

                        ! Allocate var to the right size 
                        allocate(vs%var(vs%dim(1),vs%dim(2),1,1))

                        ! 2D variable 
                        call nc_read(par%filename,par%name,vs%var(:,:,1,1),missing_value=mv)

                    case(3)

                        ! Allocate var to the right size 
                        allocate(vs%var(vs%dim(1),vs%dim(2),vs%dim(3),1))

                        ! 3D variable 
                        call nc_read(par%filename,par%name,vs%var(:,:,:,1),missing_value=mv)

                    case DEFAULT 

                        call varslice_error("varslice_update", &
                            "ndim >= 4 with no time dimension is not allowed.", &
                            "ndim = "//to_str(par%ndim))

                end select

            else 
                ! Cases with a time dimension (more complicated)

                if (trim(slice_method) .eq. "interp" .or. &
                    trim(slice_method) .eq. "extrap") then 
                    ! Update time range for interp/extrap methods 

                    ! Additional consistency check 
                    if (size(time,1) .ne. 1) then
                        call varslice_error("varslice_update", &
                            "to use slice_method=['interp','extrap'], only one time should be "// &
                            "provided as an argument.", &
                            "time = "//to_str(time))
                    end if

                end if

                ! Determine indices of data to load 

                call get_indices(vs%idx,vs%time,vs%time_range,slice_method,with_time_sub)
                
                k0 = minval(vs%idx)
                k1 = maxval(vs%idx)

                !write(*,*) "idx: ", k0, k1
                !if (k0 .gt. 0) write(*,*) vs%time(k0)
                !if (k1 .gt. 0) write(*,*) vs%time(k1)

                if (k0 .gt. 0 .and. k1 .gt. 0) then 
                    ! Dimension range is available for loading, proceed 

                    ! Get size of time dimension needed for loading
                    nt_tot = k1-k0+1

                    ! Get size of time dimensions for major axis and sub axis
                    nt_rep   = vs%range_rep
                    nt_major = max(nt_tot / nt_rep, 1)

                    if (nt_major .ne. int(real(nt_tot)/real(nt_rep))) then
                        call varslice_error("varslice_update", &
                            "number of major time axis points does not match number of total "// &
                            "points divided by number of sub-time points.", &
                            "nt_rep   = "//to_str(nt_rep)//new_line("a")// &
                            "nt_tot   = "//to_str(nt_tot)//new_line("a")// &
                            "nt_major = "//to_str(nt_major)//new_line("a")// &
                            "int(real(nt_tot)/real(nt_rep)) = "//to_str(int(real(nt_tot)/real(nt_rep))))
                    end if

                    if (verbose) then
                        write(*,*) "nt_tot:   ", nt_tot
                        write(*,*) "nt_rep:   ", nt_rep
                        write(*,*) "nt_major: ", nt_major
                        ! write(*,*) "time: ", vs%time_range, k0, k1 
                        ! write(*,*) "      ", vs%time(k0), vs%time(k1)
                    end if

                    select case(par%ndim)

                        case(1)

                            ! Allocate local var to the right size 
                            allocate(var(nt_tot,1,1,1))

                            ! 0D (point) variable plus time dimension
                            call nc_read_multifile(par%filenames,par%name,var,missing_value=mv, &
                                    start=[k0],count=[nt_tot],nt_files=vs%nt_files)

                        case(2)

                            ! Allocate local var to the right size 
                            allocate(var(vs%dim(1),nt_tot,1,1))

                            ! 1D variable plus time dimension
                            call nc_read_multifile(par%filenames,par%name,var,missing_value=mv, &
                                    start=[1,k0],count=[vs%dim(1),nt_tot],nt_files=vs%nt_files)

                        case(3)
        
                            ! Allocate local var to the right size 
                            allocate(var(vs%dim(1),vs%dim(2),nt_tot,1))

                            ! 2D variable plus time dimension
                            call nc_read_multifile(par%filenames,par%name,var,missing_value=mv, &
                                    start=[1,1,k0],count=[vs%dim(1),vs%dim(2),nt_tot],nt_files=vs%nt_files)

                        case(4)

                            ! Allocate local var to the right size 
                            allocate(var(vs%dim(1),vs%dim(2),vs%dim(3),nt_tot))

                            ! 3D variable plus time dimension
                            call nc_read_multifile(par%filenames,par%name,var,missing_value=mv, &
                                    start=[1,1,1,k0],count=[vs%dim(1),vs%dim(2),vs%dim(3),nt_tot],nt_files=vs%nt_files)

                        case DEFAULT 

                            call varslice_error("varslice_update", &
                                "ndim > 4 with time dimension not allowed.", &
                                "ndim = "//to_str(par%ndim))

                    end select

                    ! At this point, the local var variable has been defined 
                    ! with data from the file for the appropriate time indices 

                    ! Next, we need to allocate the vs%var variable to the 
                    ! appropriate size and perform any calculations on the time 
                    ! indices of the local var variable as needed. 

                    ! Handle special case: if only one time is available 
                    ! for interp/extrap methods, then change method 
                    ! to exact 
                    if ( (trim(vs%slice_method) .eq. "interp" .or. & 
                          trim(vs%slice_method) .eq. "extrap") .and. &
                          nt_major .eq. 1) then 
                        ! Same time is given for upper and lower bound

                        slice_method = "exact"
                    end if

                    ! Now, allocate the vs%var variable to the right size

                    select case(trim(slice_method)) 

                        case("exact","range")
                            ! Allocate vs%var to the same size as var 
                            ! and store all values 

                            
                            if (allocated(vs%var)) deallocate(vs%var)
                            allocate(vs%var(size(var,1),size(var,2),size(var,3),size(var,4)))

                            ! Store data in vs%var 
                            vs%var = var 

                        case("interp","extrap","range_mean","range_sd","range_min","range_max","range_sum")
                            ! var should have two major time dimensions to interpolate
                            ! between. Allocate vs%var to appropriate size and 
                            ! perform interpolation 

                            ! Note: slice_method='interp' and 'extrap' use the same method, since 
                            ! the indices determine interpolation weights (ie, 
                            ! for slice_method='interp', if the time of interest lies
                            ! outside of the bounds of the data, then k0=k1=-1 and 
                            ! output data are set to missing values) 

                            select case(trim(slice_method))

                                case("interp","extrap")
                                    
                                    vec_method = "mean"     ! interp methods use the (weighted) mean

                                case("range_mean","range_sd","range_min","range_max","range_sum")

                                    ! Define 'vec_method'
                                    n1 = index(slice_method,"_")
                                    n2 = len_trim(slice_method)
                                    vec_method = slice_method(n1+1:n2)

                            end select

                            ! Size of dimension out is the size of the 
                            ! repitition desired. Ie, range_rep = 1 means 
                            ! to apply mean/sd/etc to all values along dimension
                            ! with the result of calculating 1 value. 
                            ! range_rep = 12 means apply calculation to every 12th 
                            ! value, resulting in 12 values along dimension.

                            nt_out = vs%range_rep

                            select case(trim(slice_method))

                                case("interp","extrap")
                                    
                                    if (nt_tot .ne. 2*nt_out) then
                                        ! Remember that if nt_tot==nt_out, this means that nt_major=1,
                                        ! and so only one time slice was available. So the method was
                                        ! changed to 'exact'.
                                        call varslice_error("varslice_update", &
                                            "something went wrong during interpolation. Exactly 2 major "// &
                                            "time slices should be available to interpolate between. Check!", &
                                            "nt_tot   = "//to_str(nt_tot)//new_line("a")// &
                                            "rep      = "//to_str(nt_out)//new_line("a")// &
                                            "2*rep    = "//to_str(2*nt_out)//new_line("a")// &
                                            "time     = "//to_str(time)//new_line("a")// &
                                            "k0, k1   = "//to_str([k0,k1])//new_line("a")// &
                                            "times    = "//to_str(vs%time(k0:k1)))
                                    end if

                                    ! Calculate time weighting between two extremes
                                    allocate(time_wt(2))
                                    if (with_time_sub) then
                                        time_wt(2) = real(floor(time(1))-floor(vs%time(k0))) / real(floor(vs%time(k1)) - floor(vs%time(k0)))
                                        time_wt(1) = 1.0_wp - time_wt(2) 
                                    else
                                        time_wt(2) = (time(1)-vs%time(k0)) / (vs%time(k1) - vs%time(k0))
                                        time_wt(1) = 1.0_wp - time_wt(2) 
                                    end if

                                case("range_mean","range_sd","range_min","range_max","range_sum")

                                    ! Equal weights to all values on major axis
                                    allocate(time_wt(nt_major))
                                    time_wt = 1.0

                            end select

                            if (minval(time_wt) .lt. 0.0_wp .or. maxval(time_wt) .gt. 1.0_wp) then
                                call varslice_error("varslice_update", &
                                    "interpolation weights are incorrect.", &
                                    "time_wt  = "//to_str(time_wt)//new_line("a")// &
                                    "time     = "//to_str(time)//new_line("a")// &
                                    "time(k0) = "//to_str(vs%time(k0))//new_line("a")// &
                                    "time(k1) = "//to_str(vs%time(k1)))
                            end if
                            
                            ! Make sure that var has at least as many values as we expect 
                            if (nt_out .gt. nt_tot) then
                                call varslice_error("varslice_update", &
                                    "the specified time range does not provide enough data points to "// &
                                    "be consistent with the specified value of range_rep "// &
                                    "(range_rep must be <= nt).", &
                                    "time_range      = "//to_str(vs%time_range)//new_line("a")// &
                                    "nt (time_range) = "//to_str(nt_tot)//new_line("a")// &
                                    "range_rep       = "//to_str(vs%range_rep))
                            end if

                            ! Allocate output at the natural rank (spatial dims
                            ! from the loaded slab, nt_out in the time slot).
                            call alloc_var(vs%var, shape(var), par%ndim, nt_out)

                            ! Collapse every rank to a single (space, time) view:
                            ! nsp is the flattened spatial extent, and the reduction
                            ! (mean/sd/interp/...) is one loop regardless of ndim.
                            nsp = size(var) / nt_tot
                            src2d(1:nsp,1:nt_tot) => var
                            dst2d(1:nsp,1:nt_out) => vs%var

                            do k = 1, nt_out

                                ! Get indices for current repitition
                                call get_rep_indices(kk,i0=k,i1=nt_tot,nrep=vs%range_rep)

                                do s = 1, nsp
                                    ! Calculate the vector value desired (mean,sd,...)
                                    call calc_vec_value(dst2d(s,k),src2d(s,kk),vec_method,mv,wt=time_wt)
                                end do

                            end do


                    end select

                else 
                    ! Dimension range was not found, set variable to missing values 

                    select case(par%ndim)
                        case(1)
                            allocate(vs%var(vs%dim(1),1,1,1))
                        case(2)
                            allocate(vs%var(vs%dim(1),vs%dim(2),1,1))
                        case(3)
                            allocate(vs%var(vs%dim(1),vs%dim(2),vs%dim(3),1))
                    end select

                    vs%var = mv 

                end if

            end if 

            ! Make sure crazy values have been set to missing (for safety)
            where (abs(vs%var) .ge. var_lim_crazy) vs%var = mv

            ! Apply scaling 
            where (vs%var .ne. mv) 
                vs%var = vs%var*par%unit_scale + par%unit_offset
            end where 
            
        end if 

        if (present(print_summary)) then
            if (print_summary) then
                call print_var_range(vs%var, trim(slice_method), mv) 
            end if
        end if

        return 

    end subroutine varslice_update

    subroutine get_indices(idx, x, xrange, slice_method, with_sub)
        ! Get indices in range x0 <= x <= x1
        ! This way all decimal values within x0 and x1 are included.
        
        implicit none
        
        integer,  allocatable, intent(INOUT) :: idx(:)          ! Output indices
        real(dp),              intent(IN)    :: x(:)            ! Time array in years
        real(dp),              intent(IN)    :: xrange(:)       ! [Start, End] x (inclusive) or [x_interp]
        character(len=*),      intent(IN)    :: slice_method    ! method = "exact", "interp", "extrap"
        logical,               intent(IN)    :: with_sub        ! Use fractional time unit too

        ! Local variables
        integer :: i, n, nidx
        integer,  allocatable :: ii(:)
        real(dp), allocatable :: xmain(:)
        real(dp), allocatable :: dist(:)
        real(dp) :: dist_min_lo, dist_min_hi
        real(dp) :: x0, x1

        n     = size(x)

        ! Size work arrays to the time axis (no fixed upper bound), and
        ! keep them in double precision to match the incoming time axis x.
        allocate(ii(n))
        allocate(xmain(n))
        allocate(dist(n))

        xmain = 0
        nidx  = 0
        ii    = 0

        if (trim(slice_method) .eq. "interp" .or. trim(slice_method) .eq. "extrap") then
            x0 = xrange(1)
            x1 = xrange(1)
        else
            x0 = xrange(1)
            x1 = xrange(2) 
        end if

        if (with_sub) then
            if ( abs(x0-floor(x0)) .gt. TOL .or. abs(x1-floor(x1)) .gt. TOL) then
                call varslice_error("get_indices", &
                    "when time sub-axis is used, the time range should be specified by whole numbers.", &
                    "xrange = "//to_str(xrange))
            end if
        end if

        if (with_sub) then
            xmain(1:n) = floor(x)
        else
            xmain(1:n) = x
        end if

        select case(trim(slice_method))

            case("exact","range","range_mean","range_sd","range_min","range_max","range_sum")

                if ( x0 .ge. xmain(1)-TOL .and. x1 .le. xmain(n)+TOL) then
                    ! All values within the range are available, proceed
                    
                    do i = 1, n
                        if (xmain(i) >= x0-TOL .and. xmain(i) <= x1+TOL) ii(i) = i
                    end do
                    
                end if 

            case("interp","extrap")

                if (x0 .ge. xmain(1)-TOL .and. x0 .le. xmain(n)+TOL) then
                    ! Interpolation is possible, proceed
                    
                    dist = x0 - xmain(1:n)
                    dist_min_lo = maxval(dist(1:n),mask=dist(1:n).lt.0.0)
                    dist_min_hi = minval(dist(1:n),mask=dist(1:n).ge.0.0)

                    do i = 1, n
                        if (dist(i) .eq. dist_min_lo .or. dist(i) .eq. dist_min_hi) then
                            ii(i) = i
                        end if
                    end do
                
                else if (trim(slice_method) .eq. "extrap") then
                    ! Interp value is above or below bounds, and extrapolation is desired

                    if (x0 .lt. xmain(1)-TOL) then
                        dist = x0 - xmain(1:n)
                        dist_min_lo = maxval(dist(1:n),mask=dist(1:n).lt.0.0)
                        do i = 1, n
                            if (dist(i) .eq. dist_min_lo) then
                                ii(i) = i
                            end if
                        end do
                    else if (x0 .gt. xmain(n)+TOL) then
                        dist = x0 - xmain(1:n)
                        dist_min_hi = minval(dist(1:n),mask=dist(1:n).ge.0.0)
                        do i = 1, n
                            if (dist(i) .eq. dist_min_hi) then
                                ii(i) = i
                            end if
                        end do
                    end if

                end if

            case DEFAULT
                ! Guard against a slice_method that has no index rule here.
                ! Without this, an unhandled method silently yields no
                ! indices and the whole field is returned as missing.
                call varslice_error("get_indices", &
                    "slice_method not recognized.", &
                    "slice_method = "//trim(slice_method))

        end select

        nidx = count(ii .gt. 0)

        if (nidx .gt. 0) then
            if (allocated(idx)) deallocate(idx)
            allocate(idx(nidx)) 
            idx = pack(ii, ii.gt.0)
        else
            if (allocated(idx)) deallocate(idx)
            allocate(idx(1))
            idx(1) = -1
        end if

        if (verbose) then
            write(*,*) "get_indices: ", x0, x1, slice_method, with_sub
            !write(*,*) "x: ", xmain(1:n)
            write(*,*) "nidx: ", nidx
            !write(*,*) "idx: ", pack(ii,ii.gt.0)
            write(*,*) "idx: ", idx
            write(*,*) "x:   ", x(idx)
        end if

        return

    end subroutine get_indices

    subroutine get_rep_indices(ii,i0,i1,nrep)
        ! Given starting value i0 and final value i1, 
        ! and number of values to skip nrep, generate 
        ! a vector of indices ii. 
        ! eg, i0 = 1, i1 = 36, nrep = 12
        ! => ii = [1,13,25]
        ! eg, i0 = 2, i1 = 36, nrep = 12
        ! => ii = [2,14,26]
        
        implicit none 

        integer, allocatable, intent(OUT) :: ii(:) 
        integer, intent(IN) :: i0 
        integer, intent(IN) :: i1 
        integer, intent(IN) :: nrep 

        ! Local variables
        integer :: i, ni, ntot
        integer, allocatable :: jj(:)

        ni = i1-i0+1

        allocate(jj(ni))
        jj = 0

        do i = 1, ni
            jj(i) = i0 + nrep*(i-1) 
            if (jj(i) .gt. i1) then 
                jj(i) = 0
                exit 
            end if 
        end do 

        ntot = count(jj .gt. 0)
        if (allocated(ii)) deallocate(ii)
        allocate(ii(ntot)) 

        ii = jj(1:ntot) 

        return 

    end subroutine get_rep_indices

    subroutine calc_vec_value(val,var,method,mv,wt)

        implicit none 

        real(wp),         intent(OUT) :: val 
        real(wp),         intent(IN)  :: var(:) 
        character(len=*), intent(IN)  :: method 
        real(wp),         intent(IN)  :: mv 
        real(wp),         intent(IN), optional :: wt(:) 

        ! Local variables 
        integer  :: ntot 
        real(wp) :: mean, variance 
        real(wp), allocatable :: wt_now(:) 

        ntot = count(var .ne. mv)
        
        if (ntot .gt. 0) then 
            ! Values available for calculations 

            allocate(wt_now(size(var)))
            if (present(wt)) then 
                wt_now = wt

                ! Safety check
                if (size(wt,1) .ne. size(var,1)) then
                    call varslice_error("calc_vec_value", &
                        "wt vector must be the same length as the var vector.", &
                        "size(wt)  = "//to_str(size(wt,1))//new_line("a")// &
                        "size(var) = "//to_str(size(var,1)))
                end if
            else 
                wt_now = 1.0_wp
            end if 

            where(var .eq. mv) wt_now = 0.0_wp 

            select case(trim(method))

                case("sum","mean")
                    ! Calculate a weighted sum/mean with the right weights for each case 

                    ! Normalize weights to sum to one
                    if (trim(method) .eq. "mean") wt_now = wt_now / sum(wt_now)
                    
                    val = sum(wt_now*var)

                case("sd")
                    ! Calculate the weighted mean and then weighted standard deviation 

                    ! Normalize weights to sum to one
                    wt_now = wt_now / sum(wt_now)

                    if (ntot .ge. 2) then

                        mean     = sum(wt_now*var)
                        variance = real(ntot/(ntot-1),wp) * sum( wt_now*(var-mean)**2 )
                        val      = sqrt(variance)

                    else 

                        val = mv 

                    end if 

                case("min") 

                    val = minval(var,mask=var.ne.mv)

                case("max")

                    val = maxval(var,mask=var.ne.mv)

                case DEFAULT 

                    call varslice_error("calc_vec_value", &
                        "method not recognized.", &
                        "method = "//trim(method))

            end select

        else 
            ! No values available in vector, set to missing value 

            val = mv 

        end if 

        return 

    end subroutine calc_vec_value

    subroutine varslice_init_nml(vs,filename,group,domain,grid_name,verbose)
        ! Routine to load information related to a given 
        ! transient variable, so that it can be processed properly.

        implicit none 

        type(varslice_class),   intent(INOUT) :: vs
        character(len=*),       intent(IN)    :: filename
        character(len=*),       intent(IN)    :: group
        character(len=*),       intent(IN), optional :: domain
        character(len=*),       intent(IN), optional :: grid_name
        logical,                intent(IN), optional :: verbose 
        ! Local variables 
        
        ! First load parameters from nml file 
        call varslice_par_load(vs%par,filename,group,domain,grid_name,verbose)

        ! Perform remaining init operations 
        call varslice_init_data(vs) 

        return 

    end subroutine varslice_init_nml

    subroutine varslice_init_arg(vs,filename,name,units_in,units_out,scale,offset, &
                                                            with_time,time_par)
        ! Routine to load information related to a given 
        ! transient variable, so that it can be processed properly.

        implicit none 

        type(varslice_class), intent(INOUT) :: vs
        character(len=*),      intent(IN)   :: filename
        character(len=*),      intent(IN)   :: name
        character(len=*),      intent(IN)   :: units_in 
        character(len=*),      intent(IN)   :: units_out 
        real(wp), optional,    intent(IN)   :: scale 
        real(wp), optional,    intent(IN)   :: offset 
        logical,  optional,    intent(IN)   :: with_time 
        real(wp), optional,    intent(IN)   :: time_par(:) 

        ! Local variables 
        
        ! Define parameters from subroutine arguments 
        vs%par%filename  = trim(filename) 
        vs%par%name      = trim(name) 
        vs%par%units_in  = trim(units_in) 
        vs%par%units_out = trim(units_out) 
        
        vs%par%unit_scale = 1.0_wp 
        if (present(scale)) vs%par%unit_scale = scale 

        vs%par%unit_offset = 0.0_wp 
        if (present(offset)) vs%par%unit_offset = offset 
        
        vs%par%with_time = .TRUE. 
        if (present(with_time)) vs%par%with_time = with_time 

        vs%par%time_par = [0.0,0.0,0.0,0.0]
        if (present(time_par)) vs%par%time_par(1:size(time_par)) = time_par

        ! Resolve file list and derive internal time parameters
        ! (mirrors the nml path; without this par%filenames and
        !  par%with_time_sub would be left unset)
        call varslice_par_finalize(vs%par)

        ! Perform remaining init operations
        call varslice_init_data(vs)

        return

    end subroutine varslice_init_arg

    subroutine varslice_init_data(vs)

        implicit none 

        type(varslice_class), intent(INOUT) :: vs 

        ! Local variables  
        character(len=12), allocatable :: dim_names(:) 
        logical :: with_time 
        logical :: with_time_sub
        real(dp) :: dt 
        integer  :: i, k 
        integer  :: nt
        integer, allocatable :: dim_now(:)
        character(len=1024)  :: filename_dims
        character(len=:), allocatable :: fnames

        ! Local shortcut
        with_time       = vs%par%with_time 
        with_time_sub   = vs%par%with_time_sub 

        ! Store the first filename locally, to be used for getting variable dimensions (except time)
        filename_dims = vs%par%filenames(1)

        ! First make sure all data objects are deallocated 
        call varslice_end(vs)

        ! Get information from netcdf file - first file in case their are multiple sources 
        call nc_dims(filename_dims,vs%par%name,dim_names,vs%dim)
        vs%par%ndim = size(vs%dim,1)

        ! Determine total time dimension size and cache the per-file time
        ! lengths. nt_files is reused by nc_read_multifile so it need not
        ! re-query the files on every update.
        if (with_time) then
            if (allocated(vs%nt_files)) deallocate(vs%nt_files)
            allocate(vs%nt_files(size(vs%par%filenames,1)))
            do i = 1, size(vs%par%filenames,1)
                call nc_dims(vs%par%filenames(i),vs%par%name,dim_names,dim_now)
                vs%nt_files(i) = dim_now(vs%par%ndim)
            end do
            vs%dim(vs%par%ndim) = sum(vs%nt_files)
        end if

! ======== TO DO =============
! In the case, of using nc_read_interp, data in arrays will have shape nx,ny
! of target grid, not necessarily of input data file. Adjust dims here
! based on target grid definition.

        !call nc_read()


! ============================


        if (with_time) then

            if (with_time_sub) then
                ! Need to generate a sub-annual axis too

                dt = 1.0 / vs%par%time_par(4)

                call axis_init(vs%time, &
                                x0=vs%par%time_par(1) + dt/2.0, &
                                x1=vs%par%time_par(2) + 1.0 - dt/2.0, &
                                dx=dt)

                ! Store number of sub-value (1:time_par(4)), so 1:12 for monthly values, repeating for all of time.

                if (allocated(vs%time_sub)) deallocate(vs%time_sub)
                allocate(vs%time_sub(size(vs%time)))

                do k = 1, size(vs%time)
                    vs%time_sub(k) = ceiling( (vs%time(k)-floor(vs%time(k)))/dt )
                    !write(*,*) vs%time(k), vs%time_sub(k)
                end do

                ! write(*,*) "Testing axis generation: "
                ! write(*,*) 
                ! do k = 1, size(vs%time)
                !     write(*,*) vs%time(k), floor(vs%time(k)), ceiling( (vs%time(k)-floor(vs%time(k)))/dt )
                ! end do
                ! write(*,*) 
                ! write(*,*) "nt = ", size(vs%time)
                ! stop 

                vs%par%with_time_sub = .TRUE. 

            else
                ! No special sub-annual axis needed.
                ! Note: time_par(4) must be greater than 1 to work
                ! as a fractional amount of time between major units.

                ! Initialize time vector from user parameters 
                call axis_init(vs%time,x0=vs%par%time_par(1), &
                                    x1=vs%par%time_par(2), &
                                    dx=vs%par%time_par(3))

            end if

            ! Check to make sure time vector matches netcdf file length 
            if (size(vs%time,1) .ne. vs%dim(vs%par%ndim)) then
                fnames = ""
                if (size(vs%par%filenames,1) .gt. 1) then
                    fnames = "filenames    ="
                    do i = 1, size(vs%par%filenames,1)
                        fnames = fnames//new_line("a")//"               "//trim(vs%par%filenames(i))
                    end do
                    fnames = fnames//new_line("a")
                end if
                call varslice_error("varslice_init_data", &
                    "generated time coordinate does not match the length of the time "// &
                    "dimension in the netcdf file.", &
                    "time_par    = "//to_str(vs%par%time_par)//new_line("a")// &
                    "size(time)  = "//to_str(size(vs%time,1))//new_line("a")// &
                    "nt (netcdf) = "//to_str(vs%dim(vs%par%ndim))//new_line("a")// &
                    "filename    = "//trim(vs%par%filename)//new_line("a")// &
                    fnames// &
                    "time        = "//to_str(vs%time))
            end if

        end if 

        ! Allocate coordinate variables and load coordinates too. The number
        ! of spatial axes is ndim, minus one when the last dimension is time.
        select case(vs%par%ndim)

            case(1)
                ! Point (0D space), or a 1D-in-space variable with no time.
                ! No spatial axis coordinates to load here.

            case(2)
                if (with_time) then
                    ! 1D space + time
                    call load_or_generate_axis(filename_dims,dim_names(1),vs%dim(1),vs%x)
                else
                    ! 2D space
                    call load_or_generate_axis(filename_dims,dim_names(1),vs%dim(1),vs%x)
                    call load_or_generate_axis(filename_dims,dim_names(2),vs%dim(2),vs%y)
                end if

            case(3)
                if (with_time) then
                    ! 2D space + time
                    call load_or_generate_axis(filename_dims,dim_names(1),vs%dim(1),vs%x)
                    call load_or_generate_axis(filename_dims,dim_names(2),vs%dim(2),vs%y)
                else
                    ! 3D space
                    call load_or_generate_axis(filename_dims,dim_names(1),vs%dim(1),vs%x)
                    call load_or_generate_axis(filename_dims,dim_names(2),vs%dim(2),vs%y)
                    call load_or_generate_axis(filename_dims,dim_names(3),vs%dim(3),vs%z)
                end if

            case(4)
                if (with_time) then
                    ! 3D space + time
                    call load_or_generate_axis(filename_dims,dim_names(1),vs%dim(1),vs%x)
                    call load_or_generate_axis(filename_dims,dim_names(2),vs%dim(2),vs%y)
                    call load_or_generate_axis(filename_dims,dim_names(3),vs%dim(3),vs%z)
                else
                    call varslice_error("varslice_init_data", &
                        "4D array without time dimension is not yet supported.")
                end if

            case DEFAULT
                call varslice_error("varslice_init_data", &
                    "ndim > 4 not allowed.", &
                    "ndim = "//to_str(vs%par%ndim))

        end select

        ! Make sure all axis variables exist
        if (.not. allocated(vs%x)) then 
            allocate(vs%x(1))
            vs%x = 0.0_wp 
        end if 
        if (.not. allocated(vs%y)) then 
            allocate(vs%y(1))
            vs%y = 0.0_wp 
        end if 
        if (.not. allocated(vs%z)) then 
            allocate(vs%z(1))
            vs%z = 0.0_wp 
        end if 
        if (.not. allocated(vs%time)) then 
            allocate(vs%time(1))
            vs%time = 0.0_wp 
        end if 
        if (.not. allocated(vs%time_sub)) then 
            allocate(vs%time_sub(1))
            vs%time_sub = 0.0_wp 
        end if 
        if (.not. allocated(vs%idx)) then 
            allocate(vs%idx(1))
            vs%idx = 0 
        end if 

        return  

    end subroutine varslice_init_data

    subroutine alloc_var(var, src_shape, ndim, nout)
        ! (Re)allocate the 4D storage array so that the first ndim-1 axes are
        ! the spatial extents (taken from src_shape, the shape of the loaded
        ! slab), the ndim-th axis has length nout (the time-like output), and
        ! any remaining trailing axes are singleton. This keeps a single code
        ! path for every rank instead of a per-ndim allocate.

        implicit none

        real(wp), allocatable, intent(INOUT) :: var(:,:,:,:)
        integer,               intent(IN)    :: src_shape(4)
        integer,               intent(IN)    :: ndim
        integer,               intent(IN)    :: nout

        integer :: shp(4)

        shp = 1
        if (ndim .ge. 2) shp(1:ndim-1) = src_shape(1:ndim-1)
        shp(ndim) = nout

        if (allocated(var)) deallocate(var)
        allocate(var(shp(1),shp(2),shp(3),shp(4)))

        return

    end subroutine alloc_var

    function varslice_v1(vs, k) result(v)
        ! 1D (vector) copy at trailing index k: var(:,k,1,1)
        class(varslice_class), intent(IN) :: vs
        integer, optional,     intent(IN) :: k
        real(wp), allocatable :: v(:)
        integer :: kk
        kk = 1
        if (present(k)) kk = k
        v = vs%var(:,kk,1,1)
    end function varslice_v1

    function varslice_v2(vs, k) result(v)
        ! 2D (field) copy at trailing index k: var(:,:,k,1)
        class(varslice_class), intent(IN) :: vs
        integer, optional,     intent(IN) :: k
        real(wp), allocatable :: v(:,:)
        integer :: kk
        kk = 1
        if (present(k)) kk = k
        v = vs%var(:,:,kk,1)
    end function varslice_v2

    function varslice_v3(vs, k) result(v)
        ! 3D (volume) copy at trailing index k: var(:,:,:,k)
        class(varslice_class), intent(IN) :: vs
        integer, optional,     intent(IN) :: k
        real(wp), allocatable :: v(:,:,:)
        integer :: kk
        kk = 1
        if (present(k)) kk = k
        v = vs%var(:,:,:,kk)
    end function varslice_v3

    subroutine load_or_generate_axis(filename, varname, n, x)
        ! Allocate axis x to length n and populate it: read the coordinate
        ! variable from file if it exists, otherwise generate a 1..n index axis.

        implicit none

        character(len=*),      intent(IN)    :: filename
        character(len=*),      intent(IN)    :: varname
        integer,               intent(IN)    :: n
        real(wp), allocatable, intent(INOUT) :: x(:)

        if (allocated(x)) deallocate(x)
        allocate(x(n))

        if (nc_exists_var(filename, varname)) then
            call nc_read(filename, varname, x)
        else
            call axis_init(x, nx=n)
        end if

        return

    end subroutine load_or_generate_axis

    subroutine varslice_end(vs)
        ! Deallocate all variables

        implicit none 

        type(varslice_class), intent(INOUT) :: vs 

        if (allocated(vs%dim))          deallocate(vs%dim)
        if (allocated(vs%x))            deallocate(vs%x)
        if (allocated(vs%y))            deallocate(vs%y)
        if (allocated(vs%z))            deallocate(vs%z)
        if (allocated(vs%time))         deallocate(vs%time)
        if (allocated(vs%time_sub))     deallocate(vs%time_sub)
        if (allocated(vs%idx))          deallocate(vs%idx)
        if (allocated(vs%nt_files))     deallocate(vs%nt_files)
        if (allocated(vs%var))          deallocate(vs%var)
        
        return 

    end subroutine varslice_end

    subroutine varslice_par_load(par,filename,group,domain,grid_name,verbose)

        type(varslice_param_class), intent(OUT) :: par 
        character(len=*), intent(IN) :: filename
        character(len=*), intent(IN) :: group
        character(len=*), intent(IN), optional :: domain
        character(len=*), intent(IN), optional :: grid_name   
        logical, optional :: verbose 

        ! Local variables
        logical  :: print_summary
        integer  :: i
        character(len=56) :: units(2)     ! [units_in, units_out]
        real(wp)          :: scaling(2)   ! [unit_scale, unit_offset]
        real(dp)          :: time(5)      ! [active, x0, x1, dx, sub]

        print_summary = .TRUE.
        if (present(verbose)) print_summary = verbose

        ! Read the condensed namelist keys. Each value bundles related
        ! quantities so the group stays short but fully explicit:
        !   units   = "<in>" "<out>"
        !   scaling = <unit_scale> <unit_offset>
        !   time    = <active> <x0> <x1> <dx> <sub>
        ! `active` (1.0/0.0) flags whether the variable is time-varying, or
        ! whether x0..x1 merely document the period the (static) data represents.
        ! `sub` is optional (trailing) and defaults to 0; every other value is
        ! required, keeping the provenance visible in the file.

        ! Pre-set defaults so an omitted trailing element falls back sensibly
        ! (the keys themselves remain mandatory).
        units   = ""
        scaling = [1.0_wp, 0.0_wp]
        time    = 0.0_dp

        call nml_read(filename,group,"filename", par%filename)
        call nml_read(filename,group,"name",     par%name)
        call nml_read(filename,group,"units",    units)
        call nml_read(filename,group,"scaling",  scaling)
        call nml_read(filename,group,"time",     time)

        ! Unpack into the internal parameter fields
        par%units_in    = units(1)
        par%units_out   = units(2)
        par%unit_scale  = scaling(1)
        par%unit_offset = scaling(2)
        par%with_time   = (time(1) .ge. 0.5_dp)
        par%time_par    = time(2:5)

        ! Parse filename as needed
        if (present(domain) .and. present(grid_name)) then
            call parse_path(par%filename,domain,grid_name)
        end if

        ! Resolve file list and derive internal time parameters
        call varslice_par_finalize(par)

        ! Summary 
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
            if (par%with_time) then
                write(*,*) "time_par    = ", par%time_par
            end if
        end if 

        return

    end subroutine varslice_par_load

    subroutine varslice_par_finalize(par)
        ! Resolve the (possibly wildcard) filename into the sorted list of
        ! matching files and derive the internal time parameters. Shared by
        ! both the nml- and argument-based initialization paths so they stay
        ! consistent.

        implicit none

        type(varslice_param_class), intent(INOUT) :: par

        ! See if multiple files are available (also allocates par%filenames)
        call get_matching_files(par%filenames, par%filename)

        ! For a time-varying field with no step given, collapse to a single
        ! time slice (x1=x0). Skip this for a static field, where x0..x1 are
        ! only documenting the represented period and must be preserved.
        if (par%with_time .and. par%time_par(3) .eq. 0.0) par%time_par(2) = par%time_par(1)

        if (par%time_par(4) .gt. 1.0) then
            par%with_time_sub = .TRUE.
        else
            par%with_time_sub = .FALSE.
        end if

        return

    end subroutine varslice_par_finalize

    subroutine parse_path(path,domain,grid_name)

        implicit none 

        character(len=*), intent(INOUT) :: path 
        character(len=*), intent(IN)    :: domain, grid_name 

        call nml_replace(path,"{domain}",   trim(domain))
        call nml_replace(path,"{grid_name}",trim(grid_name))
        
        return 

    end subroutine parse_path
    
    subroutine axis_init_sp(x,x0,x1,dx,nx)

        implicit none 

        real(sp), allocatable, intent(OUT) :: x(:)
        real(sp), optional :: x0
        real(sp), optional :: x1
        real(sp), optional :: dx
        integer,  optional :: nx 

        ! Local variables 
        integer :: i  
        real(sp) :: x0_now
        real(sp) :: x1_now
        real(sp) :: dx_now
        integer  :: nx_now 
        real(sp) :: nx_check

        dx_now = 1.0_sp 
        if (present(dx)) dx_now = dx 

        x0_now = 0.0_sp 
        if (present(x0)) x0_now = x0 
        
        ! Note: if x1 is present, nx is ignored 
        if (present(x1)) then 
            x1_now = x1 
        else if (present(nx)) then 
            x1_now = x0_now + (nx-1)*dx_now 
        else
            call varslice_error("axis_init", &
                "either x1 or nx must be present.")
        end if

        if (allocated(x)) deallocate(x)

        nx_now   = nint((x1_now-x0_now)/dx_now) + 1
        nx_check = (x1_now-x0_now)/dx_now + 1

        if ( abs(nx_now - nx_check) .gt. TOL ) then
            ! Make sure nx is a round number.
            call varslice_error("axis_init", &
                "desired axis bounds [x0,x1] do not divide evenly with dx.", &
                "x0 = "//to_str(x0_now)//new_line("a")// &
                "x1 = "//to_str(x1_now)//new_line("a")// &
                "dx = "//to_str(dx_now)//new_line("a")// &
                "nx = "//to_str((x1_now-x0_now)/dx_now + 1))
        end if
        
        allocate(x(nx_now))

        do i = 1, nx_now 
            x(i) = x0_now + dx_now*real(i-1,wp)
        end do 

        return

    end subroutine axis_init_sp
    
    subroutine axis_init_dp(x,x0,x1,dx,nx)

        implicit none 

        real(dp), allocatable, intent(OUT) :: x(:)
        real(dp), optional :: x0
        real(dp), optional :: x1
        real(dp), optional :: dx
        integer,  optional :: nx 

        ! Local variables 
        integer :: i  
        real(dp) :: x0_now
        real(dp) :: x1_now
        real(dp) :: dx_now
        integer  :: nx_now 
        real(dp) :: nx_check

        dx_now = 1.0_dp 
        if (present(dx)) dx_now = dx 

        x0_now = 0.0_dp 
        if (present(x0)) x0_now = x0 
        
        ! Note: if x1 is present, nx is ignored 
        if (present(x1)) then 
            x1_now = x1 
        else if (present(nx)) then 
            x1_now = x0_now + (nx-1)*dx_now 
        else
            call varslice_error("axis_init", &
                "either x1 or nx must be present.")
        end if

        if (allocated(x)) deallocate(x)

        nx_now   = nint((x1_now-x0_now)/dx_now) + 1
        nx_check = (x1_now-x0_now)/dx_now + 1

        if ( abs(nx_now - nx_check) .gt. TOL ) then
            ! Make sure nx is a round number.
            call varslice_error("axis_init", &
                "desired axis bounds [x0,x1] do not divide evenly with dx.", &
                "x0 = "//to_str(x0_now)//new_line("a")// &
                "x1 = "//to_str(x1_now)//new_line("a")// &
                "dx = "//to_str(dx_now)//new_line("a")// &
                "nx = "//to_str((x1_now-x0_now)/dx_now + 1))
        end if
        
        allocate(x(nx_now))

        do i = 1, nx_now 
            x(i) = x0_now + dx_now*real(i-1,dp)
        end do 

        return

    end subroutine axis_init_dp

    subroutine print_var_range(var,name,missing_value,time)

        implicit none 

        real(wp),         intent(IN) :: var(:,:,:,:) 
        character(len=*), intent(IN) :: name
        real(wp),         intent(IN) :: missing_value 
        real(wp), intent(IN), optional :: time 

        ! Local variables 
        real(wp) :: vmin, vmax 

        if (count(var.ne.missing_value) .gt. 0) then 
            vmin = minval(var,mask=var.ne.missing_value)
            vmax = maxval(var,mask=var.ne.missing_value)
        else 
            vmin = missing_value 
            vmax = missing_value 
        end if 

        if (present(time)) then 
            write(*,"(f10.1,2x,a16,a3,2f14.3)") time, trim(name), ": ", vmin, vmax
        else 
            write(*,"(10x,2x,a16,a3,2f14.3)") trim(name), ": ", vmin, vmax
        end if 

        return 

    end subroutine print_var_range

    subroutine get_matching_files(filenames, pattern)

        implicit none

        character(len=*), allocatable, intent(OUT) :: filenames(:)
        character(len=*), intent(in) :: pattern
        
        ! Local variables
        character(len=1024) :: command
        character(len=256) :: temp_filename
        integer :: i, ios, unit
        integer :: num_files

        integer, parameter :: max_num_files_allowed = 10000

        ! Temporary file to store file names. Make the name unique per
        ! process (PID) so concurrent runs sharing a working directory
        ! (e.g. batch experiments) do not clobber each other's list.
        write(temp_filename,"(a,i0,a)") "filelist_", getpid(), ".tmp"

        ! Create command to list files sorted alphabetically
        command = "ls -1 " // trim(pattern) // " | sort > " // trim(temp_filename)
        call execute_command_line(trim(command))

        ! Open the temporary file
        open(newunit=unit, file=trim(temp_filename), status="old", action="read", iostat=ios)
        
        ! Make sure at least one file was found
        if (ios /= 0) then
            call varslice_error("get_matching_files", &
                "temporary file could not be opened.", &
                "temp_filename = "//trim(temp_filename)//new_line("a")// &
                "iostat        = "//to_str(ios))
        end if

        ! First, count the number of files
        num_files = 0
        
        do i = 1, max_num_files_allowed
            read(unit, '(a)', iostat=ios)
            if (ios /= 0) exit
            num_files = num_files + 1
        end do
        rewind(unit)

        ! Make sure at least one file was found
        if (num_files .eq. 0) then
            call varslice_error("get_matching_files", &
                "filename(s) not found.", &
                "filename = "//trim(pattern))
        end if

        ! Allocate the output array
        if (allocated(filenames)) deallocate(filenames)
        allocate(filenames(num_files))

        ! Read the filenames into the array
        do i = 1, num_files
            read(unit, '(a)', iostat=ios) filenames(i)
        end do
        close(unit)

        ! Remove the temporary file
        call execute_command_line("rm " // trim(temp_filename))

        return

    end subroutine get_matching_files


    ! === Extensions to nc_read to handle multiple filenames ====

    ! Should this eventually be incorporated into ncio itself??
    ! It might be hard to do, since it would have to propagate through all routines,
    ! not just nc4_read_internal, but also nc_dims, etc...
    ! Now routine benefits from generic var[4D], but perhaps harder to collapse to 1D reading.

    subroutine nc_read_multifile(filenames,name,var,start,count,nt_files,missing_value)

        implicit none

        character (len=*),  intent(IN)      :: filenames(:)
        character (len=*),  intent(IN)      :: name
        real(wp),           intent(INOUT)   :: var(:,:,:,:)
        integer,            intent(IN)      :: start(:)
        integer,            intent(IN)      :: count(:)
        integer,            intent(IN)      :: nt_files(:)   ! per-file time length (cached)
        real(wp),           intent(IN), optional :: missing_value

        ! Local variables
        integer :: i, k0, nk, j0, nj, t0, nt, num_files
        integer :: ndim
        character(len=:), allocatable :: fnames

        ! Get number of dimensions we are working with
        ndim = size(start,1)

        num_files = size(filenames)

        ! Consistency check
        if (sum(nt_files) .lt. count(ndim)) then
            fnames = "per-file nt:"
            do i = 1, num_files
                fnames = fnames//new_line("a")// &
                    "  "//trim(filenames(i))//"  "//to_str(nt_files(i))
            end do
            call varslice_error("nc_read_multifile", &
                "number of time axis values read in is not sufficient to cover count.", &
                "count        = "//to_str(count)//new_line("a")// &
                "sum(nt_files) = "//to_str(sum(nt_files))//new_line("a")// &
                fnames)
        end if

        ! To do: figure out where index k0 begins within nt_files.
        k0 = start(size(start,1))
        nk = 0
        j0 = 0
        nt = 0
        do i = 1, num_files
            nk = nk + nt_files(i)       ! count maximum available including this file
            !nk = min(nk,count(ndim))       ! Limit maximum to total count in case it is less
            !write(*,*) "k0: ", k0, nk
            if (k0 .le. nk) then 
                j0 = k0 - j0            ! start index for current file
                nj = nk - k0 + 1        ! count available from current file
                nj = min(nj,count(ndim)-nt)  ! limit to values still needed (not total)
                nt = nt + nj            ! store total to be loaded so far
                t0 = nt - nj + 1        ! start index in output array
                !write(*,*) "j: ", j0, nj, k0, nk, t0, nt

                select case(ndim)

                    case(1)
                        call nc_read(filenames(i),name,var(t0:nt,1,1,1),missing_value=missing_value, &
                                                start=[j0],count=[nj])
                    case(2)
                        call nc_read(filenames(i),name,var(:,t0:nt,1,1),missing_value=missing_value, &
                                                start=[start(1),j0],count=[count(1),nj])
                    case(3)
                        call nc_read(filenames(i),name,var(:,:,t0:nt,1),missing_value=missing_value, &
                                                start=[start(1),start(2),j0],count=[count(1),count(2),nj])
                    case(4)
                        call nc_read(filenames(i),name,var(:,:,:,t0:nt),missing_value=missing_value, &
                                                start=[start(1),start(2),start(3),j0],count=[count(1),count(2),count(3),nj])
                end select

                k0 = k0 + nj            ! start index for whole dimension over all files
            end if
            j0 = nk                 ! reset j0 index to end of current total

            if (nt .ge. count(ndim)) exit
        end do

        if (nt .ne. count(ndim)) then
            call varslice_error("nc_read_multifile", &
                "number of time axis values read in does not match the expected total.", &
                "count = "//to_str(count)//new_line("a")// &
                "nk    = "//to_str(nk)//new_line("a")// &
                "nt    = "//to_str(nt))
        end if

        return

    end subroutine nc_read_multifile

    ! ===== Error reporting (varslice-local) =================================
    ! Self-contained helpers so every varslice failure aborts with a uniform,
    ! context-rich message on error_unit. Deliberately private and
    ! dependency-free (no shared error module) to avoid coupling varslice to
    ! the other utils modules. `varslice_error` frames the failing routine and
    ! message (plus an optional multi-line `detail`) and halts via `error stop`
    ! (nonzero exit, and a backtrace if built with -fbacktrace). `to_str`
    ! renders scalars / 1-D arrays into the `detail` string; build detail with
    ! several values by joining them with new_line("a").

    subroutine varslice_error(proc, msg, detail)
        ! Abort with a framed, informative message on error_unit.
        implicit none
        character(len=*), intent(in)           :: proc    ! failing routine, e.g. "varslice_update"
        character(len=*), intent(in)           :: msg     ! what went wrong
        character(len=*), intent(in), optional :: detail  ! extra lines (join with new_line("a"))

        integer :: p0, p1

        write(error_unit,"(a)") ""
        write(error_unit,"(a)") "varslice:: error in "//trim(proc)
        write(error_unit,"(a)") "    "//trim(msg)

        if (present(detail)) then
            ! Emit `detail` line by line, indenting each, splitting on newlines.
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

        write(error_unit,"(a)") "  stopped by varslice."
        flush(error_unit)
        error stop 1

    end subroutine varslice_error

    function to_str_int(v) result(s)
        implicit none
        integer, intent(in) :: v
        character(len=:), allocatable :: s
        character(len=32) :: buf
        write(buf,"(i0)") v
        s = trim(adjustl(buf))
    end function to_str_int

    function to_str_sp(v) result(s)
        implicit none
        real(sp), intent(in) :: v
        character(len=:), allocatable :: s
        character(len=64) :: buf
        write(buf,"(g0)") v
        s = trim(adjustl(buf))
    end function to_str_sp

    function to_str_dp(v) result(s)
        implicit none
        real(dp), intent(in) :: v
        character(len=:), allocatable :: s
        character(len=64) :: buf
        write(buf,"(g0)") v
        s = trim(adjustl(buf))
    end function to_str_dp

    function to_str_log(v) result(s)
        implicit none
        logical, intent(in) :: v
        character(len=:), allocatable :: s
        if (v) then
            s = "T"
        else
            s = "F"
        end if
    end function to_str_log

    function to_str_int1(v) result(s)
        implicit none
        integer, intent(in) :: v(:)
        character(len=:), allocatable :: s
        integer :: i
        s = ""
        do i = 1, size(v)
            s = trim(s)//" "//to_str_int(v(i))
        end do
        s = adjustl(s)
    end function to_str_int1

    function to_str_sp1(v) result(s)
        implicit none
        real(sp), intent(in) :: v(:)
        character(len=:), allocatable :: s
        integer :: i
        s = ""
        do i = 1, size(v)
            s = trim(s)//" "//to_str_sp(v(i))
        end do
        s = adjustl(s)
    end function to_str_sp1

    function to_str_dp1(v) result(s)
        implicit none
        real(dp), intent(in) :: v(:)
        character(len=:), allocatable :: s
        integer :: i
        s = ""
        do i = 1, size(v)
            s = trim(s)//" "//to_str_dp(v(i))
        end do
        s = adjustl(s)
    end function to_str_dp1

end module varslice
