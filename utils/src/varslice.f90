module varslice

    use, intrinsic :: iso_fortran_env, only : input_unit, output_unit, error_unit
    
    use precision
    use ncio 
    use nml 
    use mapping_scrip

    implicit none 

    logical, parameter :: verbose = .FALSE.

    type varslice_param_class

        character(len=1024) :: filename
        character(len=56)   :: name 
        character(len=56)   :: units_in
        character(len=56)   :: units_out
        real(wp) :: unit_scale 
        real(wp) :: unit_offset
        logical  :: with_time 
        logical  :: with_time_sub

        ! Internal parameters
        integer  :: ndim 
        real(wp) :: time_par(4) 

    end type

    type varslice_class 

        type(varslice_param_class) :: par 

        ! Parameters defined during update call
        real(wp)          :: time_range(2)
        character(len=56) :: slice_method 
        integer           :: range_rep 
        
        ! Variable information
        integer,  allocatable :: dim(:)  
        real(wp), allocatable :: x(:) 
        real(wp), allocatable :: y(:)
        real(wp), allocatable :: z(:)  
        real(wp), allocatable :: time(:)
        real(wp), allocatable :: time_sub(:)
        integer, allocatable  :: idx(:)
        
        real(wp), allocatable :: var(:,:,:,:)
                
    end type 

    private 
    public :: varslice_param_class
    public :: varslice_class
    public :: varslice_update
    public :: varslice_init_nml 
    public :: varslice_init_arg
    public :: varslice_init_data
    public :: varslice_end 

    public :: print_var_range
    
contains

    subroutine varslice_map_to_grid(vs_tgt,vs_src,mps)

        implicit none

        type(varslice_class),  intent(INOUT) :: vs_tgt
        type(varslice_class),  intent(IN)    :: vs_src
        type(map_scrip_class), intent(IN)    :: mps

        ! Local variables
        integer :: nx, ny 


        ! Determine size of target grid
        nx = mps%dst_grid_dims(1)
        ny = mps%dst_grid_dims(2)

        ! Initialize meta information
        vs_tgt%par = vs_src%par

        ! Re-allocate target varslice object


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

        type(varslice_class),       intent(INOUT) :: vs
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
        real(wp), allocatable :: var(:,:,:,:) 

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
                write(*,*) "varslice_update:: Error: current time or time range &
                            &must be given as an argument (1D array)."
                stop 
            end if 

            ! Consistency check 
            if (size(time,1) .eq. 2) then 
                if (time(2) .lt. time(1)) then 
                    write(*,*) "varslice_update:: Error: time(2) should be >= time(1)."
                    write(*,*) "time = ", time
                    stop 
                end if
            end if 

        end if 

        slice_method = "exact"
        if (present(method)) slice_method = trim(method)

        fill_method = "none"
        if (present(fill)) fill_method = trim(fill) 

        if (trim(fill_method) .ne. "none") then 
            write(error_unit,*) "Error: varslice: fill methods have not yet been implemented. &
            &Set fill_method='none' for now."
            stop
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

            if (.not. with_time) then 
                ! Handle cases that do not have time dimension (simpler)

                select case(par%ndim)

                    case(1)

                        ! 1D variable
                        call nc_read(par%filename,par%name,vs%var(:,1,1,1),missing_value=mv)

                    case(2)

                        ! 2D variable 
                        call nc_read(par%filename,par%name,vs%var(:,:,1,1),missing_value=mv)

                    case(3)

                        ! 3D variable 
                        call nc_read(par%filename,par%name,vs%var(:,:,:,1),missing_value=mv)

                    case DEFAULT 

                        write(*,*) "varslice_update:: ndim >= 4 with no time dimension is not allowed."
                        write(*,*) "ndim = ", par%ndim 
                        stop 

                end select

            else 
                ! Cases with a time dimension (more complicated)

                if (trim(slice_method) .eq. "interp" .or. &
                    trim(slice_method) .eq. "extrap") then 
                    ! Update time range for interp/extrap methods 

                    ! Additional consistency check 
                    if (size(time,1) .ne. 1) then 
                        write(*,*) "varslice_update:: Error: to use slice_method=['interp','extrap'], &
                        &only one time should be provided as an argument."
                        write(*,*) "time = ", time 
                        stop 
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
                        write(error_unit,*) "varslice_update:: Error: number of major time axis points &
                        &does not match number of total points divided by number of sub-time points."
                        write(error_unit,*) "nt_rep   = ", nt_rep
                        write(error_unit,*) "nt_tot   = ", nt_tot
                        write(error_unit,*) "nt_major = ", nt_major
                        write(error_unit,*) "int(real(nt_tot)/real(nt_rep)) = ", int(real(nt_tot)/real(nt_rep))
                        stop
                    end if

                    if (verbose) then
                        write(*,*) "nt_tot:   ", nt_tot
                        write(*,*) "nt_rep:   ", nt_rep
                        write(*,*) "nt_major: ", nt_major
                        ! write(*,*) "time: ", vs%time_range, k0, k1 
                        ! write(*,*) "      ", vs%time(k0), vs%time(k1)
                    end if
                        

                    if (allocated(var)) deallocate(var) 

                    select case(par%ndim)

                        case(1)

                            ! Allocate local var to the right size 
                            allocate(var(nt_tot,1,1,1))

                            ! 0D (point) variable plus time dimension 
                            call nc_read(par%filename,par%name,var,missing_value=mv, &
                                    start=[k0],count=[nt_tot])

                        case(2)

                            ! Allocate local var to the right size 
                            allocate(var(vs%dim(1),nt_tot,1,1))

                            ! 1D variable plus time dimension 
                            call nc_read(par%filename,par%name,var,missing_value=mv, &
                                    start=[1,k0],count=[vs%dim(1),nt_tot])

                        case(3)
        
                            ! Allocate local var to the right size 
                            allocate(var(vs%dim(1),vs%dim(2),nt_tot,1))

                            ! 2D variable plus time dimension 
                            call nc_read(par%filename,par%name,var,missing_value=mv, &
                                    start=[1,1,k0],count=[vs%dim(1),vs%dim(2),nt_tot])

                        case(4)

                            ! Allocate local var to the right size 
                            allocate(var(vs%dim(1),vs%dim(2),vs%dim(3),nt_tot))

                            ! 3D variable plus time dimension 
                            call nc_read(par%filename,par%name,var,missing_value=mv, &
                                    start=[1,1,1,k0],count=[vs%dim(1),vs%dim(2),vs%dim(3),nt_tot]) 

                        case DEFAULT 

                            write(*,*) "varslice_update:: ndim > 4 with time dimension not allowed."
                            write(*,*) "ndim = ", par%ndim 
                            stop 

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

                            if (size(vs%var,1) .eq. size(var,1) .and. &
                                size(vs%var,2) .eq. size(var,2) .and. &
                                size(vs%var,3) .eq. size(var,3) .and. &
                                size(vs%var,4) .eq. size(var,4) ) then 

                                ! No allocation needed size is the same 

                            else 

                                deallocate(vs%var)
                                allocate(vs%var(size(vs%var,1),size(vs%var,2),size(vs%var,3),size(vs%var,4)))

                            end if 

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
                                        write(*,*) "varslice_update:: Error: something went wrong during &
                                        &interpolation. Exactly 2 major time slices should be available to interpolate &
                                        &between. Check!"
                                        write(*,*) "nt_tot   = ", nt_tot
                                        write(*,*) "rep      = ", nt_out 
                                        write(*,*) "2*rep    = ", 2*nt_out 
                                        write(*,*) "time     = ", time 
                                        write(*,*) "indices k0, k1: ", k0, k1 
                                        write(*,*) "times : ", vs%time(k0:k1) 
                                        stop 
                                        ! Remember that if nt_tot==nt_out, this means that nt_major=1,
                                        ! and so only one time slice was available. So the method was
                                        ! changed to 'exact'.
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
                                write(*,*) "varslice_update:: Error: interpolation weights are incorrect."
                                write(*,*) "time_wt  = ", time_wt 
                                write(*,*) "time     = ", time 
                                write(*,*) "time(k0) = ", vs%time(k0)
                                write(*,*) "time(k1) = ", vs%time(k1) 
                                stop
                            end if
                            
                            ! Make sure that var has at least as many values as we expect 
                            if (nt_out .gt. nt_tot) then 
                                write(*,*) "varslice_update:: Error: the specified time range &
                                    & does not provide enough data points to be consistent with &
                                    & the specified value of range_rep."
                                write(*,*) "time_range      = ", vs%time_range 
                                write(*,*) "nt (time_range) = ", nt_tot 
                                write(*,*) "range_rep       = ", vs%range_rep 
                                write(*,*) "range_rep must be <= nt."
                                stop 
                            end if 

                            deallocate(vs%var)

                            select case(par%ndim)

                                case(1)
                                    allocate(vs%var(nt_out,1,1,1))

                                    ! Calculate each slice 
                                    do k = 1, nt_out 

                                        ! Get indices for current repitition
                                        call get_rep_indices(kk,i0=k,i1=nt_tot,nrep=vs%range_rep)

                                        ! Calculate the vector value desired (mean,sd,...)
                                        call calc_vec_value(vs%var(k,1,1,1),var(kk,1,1,1),vec_method,mv,wt=time_wt)

                                    end do 

                                case(2)
                                    allocate(vs%var(size(var,1),nt_out,1,1))

                                    ! Calculate each slice 
                                    do k = 1, nt_out 

                                        ! Get indices for current repitition
                                        call get_rep_indices(kk,i0=k,i1=nt_tot,nrep=vs%range_rep)

                                        do i = 1, size(vs%var,1)
                                            ! Calculate the vector value desired (mean,sd,...)
                                            call calc_vec_value(vs%var(i,k,1,1),var(i,kk,1,1),vec_method,mv,wt=time_wt)
                                        end do 

                                    end do 
                                    
                                case(3)
                                    allocate(vs%var(size(var,1),size(var,2),nt_out,1))
                                
                                    ! Calculate each slice 
                                    do k = 1, nt_out 

                                        ! Get indices for current repitition
                                        call get_rep_indices(kk,i0=k,i1=nt_tot,nrep=vs%range_rep)

                                        do j = 1, size(vs%var,2)
                                        do i = 1, size(vs%var,1)
                                            ! Calculate the vector value desired (mean,sd,...)
                                            call calc_vec_value(vs%var(i,j,k,1),var(i,j,kk,1),vec_method,mv,wt=time_wt)
                                        end do
                                        end do 
                                        
                                    end do 
                                    
                                case(4)
                                    allocate(vs%var(size(var,1),size(var,2),size(var,3),nt_out))

                                    ! Calculate each slice 
                                    do k = 1, nt_out 

                                        ! Get indices for current repitition
                                        call get_rep_indices(kk,i0=k,i1=nt_tot,nrep=vs%range_rep)

                                        do l = 1, size(vs%var,3)
                                        do j = 1, size(vs%var,2)
                                        do i = 1, size(vs%var,1)
                                            ! Calculate the vector value desired (mean,sd,...)
                                            call calc_vec_value(vs%var(i,j,l,k),var(i,j,l,kk),vec_method,mv,wt=time_wt)
                                        end do
                                        end do
                                        end do 
                                        
                                    end do 
                                    
                            end select


                    end select

                else 
                    ! Dimension range was not found, set variable to missing values 

                    vs%var = mv 

                end if

            end if 

            ! Make sure crazy values have been set to missing (for safety)
            where (abs(vs%var) .ge. 1e10) vs%var = mv 

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
        real(wp),              intent(IN)    :: x(:)            ! Time array in years
        real(wp),              intent(IN)    :: xrange(:)       ! [Start, End] x (inclusive) or [x_interp]
        character(len=*),      intent(IN)    :: slice_method    ! method = "exact", "interp", "extrap"
        logical,               intent(IN)    :: with_sub        ! Use fractional time unit too

        ! Local variables
        integer :: i, n, nidx, i1, ni
        integer :: ii(10000)
        real(wp) :: xmain(10000)
        real(wp) :: dist(10000)
        real(wp) :: dist_min_lo, dist_min_hi
        real(wp) :: x0, x1

        n     = size(x)
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
                write(error_unit,*) "get_indices:: Error: when time sub-axis is used, then the time range &
                & should be specified by whole numbers."
                write(error_unit,*) "xrange: ", xrange
                stop
            end if
        end if

        if (with_sub) then
            xmain(1:n) = floor(x)
        else
            xmain(1:n) = x
        end if

        select case(trim(slice_method))

            case("exact","range","range_mean","range_sd","range_min","range_max")

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
        integer :: jj(10000)

        ni = i1-i0+1 

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
                    write(error_unit,*) "calc_vec_value:: Error: wt vector must be the same length as the var vector."
                    write(error_unit,*) "size(wt):  ", size(wt,1)
                    write(error_unit,*) "size(var): ", size(var,1)
                    stop
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

                    write(*,*) "calc_vec_value:: Error: method not recognized."
                    write(*,*) "method = ", trim(method) 
                    stop 

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
        real(wp) :: dt 
        integer  :: k 

        ! Local shortcut
        with_time       = vs%par%with_time 
        with_time_sub   = vs%par%with_time_sub 

        ! First make sure all data objects are deallocated 
        call varslice_end(vs)

        ! Get information from netcdf file 
        call nc_dims(vs%par%filename,vs%par%name,dim_names,vs%dim)
        vs%par%ndim = size(vs%dim,1)



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
                write(*,*) "varslice_init_data:: Error: generated time coordinate &
                &does not match the length of the time dimension in the netcdf file."
                write(*,*) "time_par:    ", vs%par%time_par 
                write(*,*) "size(time):  ", size(vs%time,1)
                write(*,*) "nt (netcdf): ", vs%dim(vs%par%ndim)
                write(*,*) "filename:    ", trim(vs%par%filename)
                write(*,*) 
                write(*,*) "time = ", vs%time
                stop 
            end if 

        end if 

        ! Allocate coordinate and data variables,
        ! load coordinates too 
        select case(vs%par%ndim)

            case(1)
                if (with_time) then 
                    allocate(vs%var(1,1,1,1))
                else 
                    allocate(vs%var(vs%dim(1),1,1,1))
                end if 
            case(2)
                if (with_time) then 
                    allocate(vs%x(vs%dim(1)))
                    allocate(vs%var(vs%dim(1),1,1,1))

                    if (nc_exists_var(vs%par%filename,dim_names(1))) then 
                        call nc_read(vs%par%filename,dim_names(1),vs%x)
                    else
                        call axis_init(vs%x,nx=vs%dim(1))
                    end if
                else 
                    allocate(vs%x(vs%dim(1)))
                    allocate(vs%y(vs%dim(2)))
                    allocate(vs%var(vs%dim(1),vs%dim(2),1,1))

                    if (nc_exists_var(vs%par%filename,dim_names(1))) then 
                        call nc_read(vs%par%filename,dim_names(1),vs%x)
                    else
                        call axis_init(vs%x,nx=vs%dim(1))
                    end if
                    if (nc_exists_var(vs%par%filename,dim_names(2))) then 
                        call nc_read(vs%par%filename,dim_names(2),vs%y)
                    else
                        call axis_init(vs%y,nx=vs%dim(2))
                    end if
                    
                end if 

            case(3)
                if (with_time) then
                    allocate(vs%x(vs%dim(1)))
                    allocate(vs%y(vs%dim(2)))
                    allocate(vs%var(vs%dim(1),vs%dim(2),1,1))

                    if (nc_exists_var(vs%par%filename,dim_names(1))) then 
                        call nc_read(vs%par%filename,dim_names(1),vs%x)
                    else
                        call axis_init(vs%x,nx=vs%dim(1))
                    end if
                    if (nc_exists_var(vs%par%filename,dim_names(2))) then 
                        call nc_read(vs%par%filename,dim_names(2),vs%y)
                    else
                        call axis_init(vs%y,nx=vs%dim(2))
                    end if
                    
                else 
                    allocate(vs%x(vs%dim(1)))
                    allocate(vs%y(vs%dim(2)))
                    allocate(vs%z(vs%dim(3)))
                    allocate(vs%var(vs%dim(1),vs%dim(2),vs%dim(3),1))

                    if (nc_exists_var(vs%par%filename,dim_names(1))) then 
                        call nc_read(vs%par%filename,dim_names(1),vs%x)
                    else
                        call axis_init(vs%x,nx=vs%dim(1))
                    end if
                    if (nc_exists_var(vs%par%filename,dim_names(2))) then 
                        call nc_read(vs%par%filename,dim_names(2),vs%y)
                    else
                        call axis_init(vs%y,nx=vs%dim(2))
                    end if
                    if (nc_exists_var(vs%par%filename,dim_names(3))) then 
                        call nc_read(vs%par%filename,dim_names(3),vs%z)
                    else
                        call axis_init(vs%z,nx=vs%dim(3))
                    end if
                    
                end if
                  
            case(4)
                if (with_time) then 
                    allocate(vs%x(vs%dim(1)))
                    allocate(vs%y(vs%dim(2)))
                    allocate(vs%z(vs%dim(3)))
                    allocate(vs%var(vs%dim(1),vs%dim(2),vs%dim(3),1))

                    if (nc_exists_var(vs%par%filename,dim_names(1))) then 
                        call nc_read(vs%par%filename,dim_names(1),vs%x)
                    else
                        call axis_init(vs%x,nx=vs%dim(1))
                    end if
                    if (nc_exists_var(vs%par%filename,dim_names(2))) then 
                        call nc_read(vs%par%filename,dim_names(2),vs%y)
                    else
                        call axis_init(vs%y,nx=vs%dim(2))
                    end if
                    if (nc_exists_var(vs%par%filename,dim_names(3))) then 
                        call nc_read(vs%par%filename,dim_names(3),vs%z)
                    else
                        call axis_init(vs%z,nx=vs%dim(3))
                    end if
                    
                else 
                    write(*,*) "varslice_init_data:: 4D array without time dimension is not yet supported."
                    stop 
                end if

                    
            case DEFAULT 
                write(*,*) "varslice_init_data:: ndim > 4 not allowed."
                write(*,*) "ndim = ", vs%par%ndim 
                stop 

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
        logical  :: init_pars
        logical  :: print_summary 

        init_pars = .FALSE.

        print_summary = .TRUE. 
        if (present(verbose)) print_summary = verbose 

        call nml_read(filename,group,"filename",       par%filename,     init=init_pars)
        call nml_read(filename,group,"name",           par%name,         init=init_pars)
        call nml_read(filename,group,"units_in",       par%units_in,     init=init_pars)
        call nml_read(filename,group,"units_out",      par%units_out,    init=init_pars)
        call nml_read(filename,group,"unit_scale",     par%unit_scale,   init=init_pars)   
        call nml_read(filename,group,"unit_offset",    par%unit_offset,  init=init_pars)   
        call nml_read(filename,group,"with_time",      par%with_time,    init=init_pars)   
        call nml_read(filename,group,"time_par",       par%time_par,     init=init_pars)   
        
        ! Parse filename as needed
        if (present(domain) .and. present(grid_name)) then
            call parse_path(par%filename,domain,grid_name)
        end if 

        ! Make sure time parameters are consistent time_par=[x0,x1,dx]
        if (par%time_par(3) .eq. 0.0) par%time_par(2) = par%time_par(1) 

        if (par%time_par(4) .gt. 1.0) then
            par%with_time_sub = .TRUE.
        else
            par%with_time_sub = .FALSE.
        end if

        ! Summary 
        if (print_summary) then  
            write(*,*) "Loading: ", trim(filename), ":: ", trim(group)
            write(*,*) "filename      = ", trim(par%filename)
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

    subroutine parse_path(path,domain,grid_name)

        implicit none 

        character(len=*), intent(INOUT) :: path 
        character(len=*), intent(IN)    :: domain, grid_name 

        call nml_replace(path,"{domain}",   trim(domain))
        call nml_replace(path,"{grid_name}",trim(grid_name))
        
        return 

    end subroutine parse_path
    
    subroutine axis_init(x,x0,x1,dx,nx)

        implicit none 

        real(wp), allocatable, intent(OUT) :: x(:)
        real(wp), optional :: x0
        real(wp), optional :: x1
        real(wp), optional :: dx
        integer,  optional :: nx 

        ! Local variables 
        integer :: i  
        real(wp) :: x0_now
        real(wp) :: x1_now
        real(wp) :: dx_now
        integer  :: nx_now 
        real(wp) :: nx_check

        dx_now = 1.0_wp 
        if (present(dx)) dx_now = dx 

        x0_now = 0.0_wp 
        if (present(x0)) x0_now = x0 
        
        ! Note: if x1 is present, nx is ignored 
        if (present(x1)) then 
            x1_now = x1 
        else if (present(nx)) then 
            x1_now = x0_now + (nx-1)*dx_now 
        else 
            write(*,*) "axis_init:: Error: either x1 or nx must be present."
            stop 
        end if 

        if (allocated(x)) deallocate(x)

        nx_now   = nint((x1_now-x0_now)/dx_now) + 1
        nx_check = (x1_now-x0_now)/dx_now + 1

        if ( abs(nx_now - nx_check) .gt. TOL ) then
            ! Make sure nx is a round number.
            write(error_unit,*) "axis_init:: Error: desired axis bounds [x0,x1] do &
            & not divide evenly with dx."
            write(error_unit,*) "  x0: ", x0_now
            write(error_unit,*) "  x1: ", x1_now
            write(error_unit,*) "  dx: ", dx_now
            write(error_unit,*) "  nx: ", ((x1_now-x0_now)/dx_now + 1)
            stop
        end if
        
        allocate(x(nx_now))

        do i = 1, nx_now 
            x(i) = x0_now + dx_now*real(i-1,wp)
        end do 

        return

    end subroutine axis_init

    subroutine print_var_range(var,name,mv,time)

        implicit none 

        real(wp),         intent(IN) :: var(:,:,:,:) 
        character(len=*), intent(IN) :: name
        real(wp),         intent(IN) :: mv 
        real(wp), intent(IN), optional :: time 

        ! Local variables 
        real(wp) :: vmin, vmax 

        if (count(var.ne.mv) .gt. 0) then 
            vmin = minval(var,mask=var.ne.mv)
            vmax = maxval(var,mask=var.ne.mv)
        else 
            vmin = mv 
            vmax = mv 
        end if 

        if (present(time)) then 
            write(*,"(f10.1,2x,a10,a3,2f14.3)") time, trim(name), ": ", vmin, vmax
        else 
            write(*,"(10x,2x,a10,a3,2f14.3)") trim(name), ": ", vmin, vmax
        end if 

        return 

    end subroutine print_var_range

end module varslice
