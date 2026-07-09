
module tsgen 

    use nml 
    use ncio 

    use precision, only: wp
    use constants, only: mv, pi, TOL
    implicit none

    ! Minimum forcing-rate magnitude [f/yr]: small but strictly nonzero, so the
    ! PI-controller methods never divide by (or take the log of) zero.
    real(wp), parameter :: DF_DT_MIN = 1e-9_wp

    type tsgen_par_class
        ! Configuration (read from namelist)
        character(len=56) :: label
        character(len=56) :: method
        logical  :: with_kill
        real(wp) :: dt_init
        real(wp) :: dt_ramp
        real(wp) :: dt_conv
        real(wp) :: dt_ave
        real(wp) :: df_sign
        real(wp) :: eps
        real(wp) :: df_dt_max
        real(wp) :: sigma
        real(wp) :: f_min
        real(wp) :: f_max
        real(wp) :: f_conv
    end type

    type tsgen_class

        type(tsgen_par_class) :: par

        ! History buffers over the time-averaging window
        real(wp), allocatable :: time(:)
        real(wp), allocatable :: var(:)
        real(wp), allocatable :: dv_dt(:)

        ! PI-controller history (n:n-2)
        real(wp) :: pi_df(3)
        real(wp) :: pi_eta(3)

        ! Runtime state
        real(wp) :: dt
        real(wp) :: dv_dt_ave
        real(wp) :: df_dt
        real(wp) :: f_now
        real(wp) :: f_mean_now
        real(wp) :: eta_now
        logical  :: kill
        real(wp) :: time_init
        logical  :: after_step      ! ramp-time-step: on the return leg of the triangle
    end type

    private
    public :: tsgen_class
    public :: tsgen_init 
    public :: tsgen_update

contains

    subroutine tsgen_init(ts,filename,time,label)

        type(tsgen_class), intent(INOUT) :: ts 
        character(len=*),   intent(IN)    :: filename 
        real(wp),           intent(IN)    :: time 
        character(len=*),   intent(IN), optional :: label 
        
        integer :: ntot 
        character(len=56) :: par_label 

        par_label = "tsgen"
        if (present(label)) par_label = trim(par_label)//"_"//trim(label)

        ! Load parameters
        call nml_read(filename,trim(par_label),"method",      ts%par%method)
        call nml_read(filename,trim(par_label),"with_kill",   ts%par%with_kill)
        call nml_read(filename,trim(par_label),"dt_ave",      ts%par%dt_ave)
        call nml_read(filename,trim(par_label),"dt_init",     ts%par%dt_init)
        call nml_read(filename,trim(par_label),"dt_ramp",     ts%par%dt_ramp)
        call nml_read(filename,trim(par_label),"dt_conv",      ts%par%dt_conv)
        call nml_read(filename,trim(par_label),"df_sign",     ts%par%df_sign)
        call nml_read(filename,trim(par_label),"eps",         ts%par%eps)
        call nml_read(filename,trim(par_label),"df_dt_max",   ts%par%df_dt_max)
        call nml_read(filename,trim(par_label),"sigma",       ts%par%sigma)
        call nml_read(filename,trim(par_label),"f_min",       ts%par%f_min)
        call nml_read(filename,trim(par_label),"f_max",       ts%par%f_max)
        call nml_read(filename,trim(par_label),"f_conv",      ts%par%f_conv)
        
        ! Make sure sign is only +1/-1 
        ts%par%df_sign = sign(1.0_wp,ts%par%df_sign)
        
        ! Consistency check 
        if (ts%par%df_sign .eq. -1.0_wp .and. trim(ts%par%method) .eq. "sin") then 
            write(*,*) "tsgen:: Error: method='sin' can only be used with &
            &df_sign=1.0 (non-negative). Please set df_sign=1.0 or use &
            &a different method." 
            stop
        end if 

        ! Set after_step false to start, since the first step is the initial ramp
        ts%after_step = .FALSE.

        ! Define label for this tsgen object
        ts%par%label = "tsgen"
        if (present(label)) ts%par%label = trim(ts%par%label)//"_"//trim(label)

        ! (Re)initialize history vectors to a large value
        ! to store many timesteps.
        ntot = 2000
        if (allocated(ts%time))  deallocate(ts%time)
        if (allocated(ts%var))   deallocate(ts%var)
        if (allocated(ts%dv_dt)) deallocate(ts%dv_dt)
        allocate(ts%time(ntot))
        allocate(ts%var(ntot))
        allocate(ts%dv_dt(ntot))

        ! Initialize variable values
        ts%time  = MV
        ts%var   = MV
        ts%dv_dt = MV

        ts%dv_dt_ave = 0.0_wp
        ts%df_dt     = 0.0_wp

        ts%pi_df  = DF_DT_MIN
        ts%pi_eta = ts%par%eps

        ! Override above choice for sin forcing 
        if (trim(ts%par%method) .eq. "sin") then 

            ! Calculate initial forcing value based on sin function
            call calc_sin_now(ts%f_mean_now,0.0_wp,ts%par%dt_ramp, &
                                ts%par%f_min,ts%par%f_max,x_offset=0.0_wp)

        else 
            ! Determine initial forcing value from parameters 

            if (ts%par%df_sign .gt. 0.0_wp) then 
                ts%f_mean_now = ts%par%f_min 
            else 
                ts%f_mean_now = ts%par%f_max 
            end if 

        end if 

        ! Set noise to zero for now
        ts%eta_now = 0.0_wp
        ts%f_now = ts%f_mean_now
        
        ! Set kill switch to false to start 
        ts%kill = .FALSE. 

        ! Make sure kill is not active for some methods
        select case(trim(ts%par%method))
            case("sin")
                ts%par%with_kill = .FALSE. 
        end select

        ! Store initial simulation time for reference (for ramp method)
        ts%time_init = time 

        return 

    end subroutine tsgen_init

  
    subroutine tsgen_update(ts,time,var,dv_dt)
        ! Generate the transient forcing value f_now for the current time,
        ! given the model response variable var (and optionally its rate dv_dt).

        type(tsgen_class), intent(INOUT) :: ts
        real(wp),           intent(IN)    :: time
        real(wp),           intent(IN)    :: var 
        real(wp), optional, intent(IN)    :: dv_dt

        ! Local variables 
        real(wp) :: time_elapsed
        real(wp) :: dv_dt_now
        real(wp) :: f_scale 
        integer  :: ntot, kmin, kmax, nk
        real(wp) :: dt_tot 
        real(wp) :: dvdt_fac 
        real(wp) :: pi_df_now

        ! For periodic forcing
        real(wp) :: f_tmp

        ! Since dv_dt is typically calculated over an averaging period,
        ! assume second-order PI controller parameters are needed. 
        integer, parameter :: pi_order = 2
        
        ! Get size of ts vectors 
        ntot = size(ts%time,1) 

        ! Get current timestep
        if (ts%time(ntot) .ne. MV) then
            ts%dt = time - ts%time(ntot)
        else
            ts%dt = 0.0_wp
        end if

        ! Get current derivative
        if (present(dv_dt)) then 
            dv_dt_now = dv_dt
        else if (ts%dt .gt. 0.0_wp) then
            dv_dt_now = (var-ts%var(ntot))/ts%dt
        else
            dv_dt_now = 0.0_wp
        end if

        ! Remove oldest point from beginning and add current one to the end
        ts%time  = eoshift(ts%time,  1,boundary=time)
        ts%var   = eoshift(ts%var,   1,boundary=var)
        ts%dv_dt = eoshift(ts%dv_dt, 1,boundary=dv_dt_now)
        

        ! Determine range of indices of times within our time-averaging window. 
        ! ajr: `findloc` only available for gfotran9 and above:
        ! kmin = findloc(ts%time .ge. time - ts%par%dt_ave,value=.TRUE., &
        !                                              dim=1,mask=ts%time.ne.MV)
        kmin = minloc(ts%time,dim=1, &
                mask=(ts%time .ge. time - ts%par%dt_ave) .and. ts%time.ne.MV)
        
        kmax = ntot

        ! Determine currently available time window
        ! Note: do not use kmin here, in case time step does not match dt_ave,
        ! rather, use all available times in the vector to see if enough 
        ! time has passed. 
        dt_tot = ts%time(kmax) - minval(ts%time,mask=ts%time.ne.MV)

        ! Get currently elapsed time overall
        time_elapsed = time - ts%time_init 

         
        ! Calculate average derivative over time steps of interest
        
        ! Get current number of averaging points 
        nk = kmax - kmin + 1 
        

        ! Calculate the current average value of the derivative
        ts%dv_dt_ave = sum(ts%dv_dt(kmin:kmax)) / real(nk,wp)
    

        if ( time_elapsed .le. ts%par%dt_init) then 
            ! Initialization time period, no transient methods applied

            ts%df_dt = 0.0_wp 

        else 
            ! Initialization time has passed, apply transient methods


            ! Determine the magnitude of rate of change (without sign)
            ! depending on method to be used.
            select case(trim(ts%par%method))

                case("const") 
                    ! Apply a constant rate of change, independent of dv_dt.
                    ! Use the df_dt_max parameter as a constant.

                    ts%df_dt = ts%par%df_dt_max

                case("ramp-time","ramp-time-step")
                    ! Ramp up to the constant rate of change for the first N years. 
                    ! Then maintain a constant anomaly (independent of dv_dt). 

                        if (trim(ts%par%method) .eq. "ramp-time-step" .and. &
                                (.not. ts%after_step) ) then
                            ! When first extreme value is reached, 
                            ! new extreme value should be imposed and 
                            ! df_sign possibly reversed.

                            if ( (ts%par%df_sign .lt. 0.0 .and.&
                                    ts%f_mean_now .le. ts%par%f_min) ) then
                                ! Ramp-down complete, start second step

                                if (ts%par%f_conv .le. ts%par%f_min) then 
                                    ! Rate still going down or constant in second step.
                                    ts%par%f_max   = ts%par%f_min
                                    ts%par%f_min   = ts%par%f_conv
                                else
                                    ! Rate changes sign in second step.
                                    ts%par%f_max   = ts%par%f_conv
                                    ts%par%df_sign = -ts%par%df_sign
                                end if 

                                ! Update duration of current step
                                ts%par%dt_ramp = ts%par%dt_conv 
                                
                                ! Set switch to know we are on the return part of the triangle
                                ts%after_step = .TRUE.

                            else if ( (ts%par%df_sign .gt. 0.0 .and.&
                                    ts%f_mean_now .ge. ts%par%f_max) ) then  
                            ! Ramp-up complete, switch directions

                                if (ts%par%f_conv .ge. ts%par%f_max) then 
                                    ! Rate still going up or constant in second step.
                                    ts%par%f_min = ts%par%f_max
                                    ts%par%f_max = ts%par%f_conv
                                else
                                    ! Rate changes sign in second step.
                                    ts%par%f_min   = ts%par%f_conv
                                    ts%par%df_sign = -ts%par%df_sign
                                end if 

                                ! Update duration of current step
                                ts%par%dt_ramp = ts%par%dt_conv 
                                
                                ! Set switch to know we are on the return part of the triangle
                                ts%after_step = .TRUE.

                            end if
                            
                        end if

                        if ( (ts%par%df_sign .lt. 0.0 .and.&
                                    ts%f_mean_now .le. ts%par%f_min) .or. &
                             (ts%par%df_sign .gt. 0.0 .and.&
                                    ts%f_mean_now .ge. ts%par%f_max) ) then
                            ! Ramp-up complete, no more forcing change

                            ts%df_dt = 0.0_wp

                        else
                            ! Linear rate of change from f_max to f_min (or vice versa) over
                            ! the time of interest dt_ramp.

                            ts%df_dt = abs(ts%par%f_max-ts%par%f_min)/ts%par%dt_ramp

                        end if

                case("ramp-slope")
                    ! Ramp up to the constant rate of change for the first N years. 
                    ! Then maintain a constant anomaly (independent of dv_dt). 

                    
                        if ( (ts%par%df_sign .lt. 0.0 .and.&
                                    ts%f_mean_now .le. ts%par%f_min) .or. &
                             (ts%par%df_sign .gt. 0.0 .and.&
                                    ts%f_mean_now .ge. ts%par%f_max) ) then  
                            ! Ramp-up complete, no more forcing change 

                            ts%df_dt = 0.0_wp 

                        else 
                            ! Linear rate of change from f_max to f_min (or vice versa) over 
                            ! the time of interest dt_ramp. 

                            ts%df_dt = ts%par%df_dt_max

                        end if 

                case("sin")
                    ! Apply periodic forcing with a given period, amplitude and x/y offset

                    ! Calculate expected current forcing value based on time elapsed
                    call calc_sin_now(f_tmp,time_elapsed,ts%par%dt_ramp, &
                                        ts%par%f_min,ts%par%f_max,x_offset=0.0_wp)

                    ! Get rate of change to apply later 
                    if (ts%dt .ne. 0.0_wp) then 
                        ts%df_dt = (f_tmp-ts%f_mean_now)/ts%dt
                    else 
                        ts%df_dt = 0.0_wp 
                    end if 
                    
                case("exp")

                    if (dt_tot .lt. ts%par%dt_ave) then 
                        ! Not enough time has passed, maintain constant forcing 
                        ! (to avoid affects of noisy derivatives)

                        ts%df_dt = 0.0_wp 

                    else 
                        ! Calculate the current forcing rate, df_dt
                        ! BASED ON EXPONENTIAL (sharp transition, tuneable)
                        ! Returns scalar in range [0-1], 0.6 at dv_dt==eps
                        ! Note: apply limit to dvdt_fac of a maximum value of 10, so 
                        ! that exp function doesn't explode (exp(-10)=>0.0)
                        dvdt_fac = min(abs(ts%dv_dt_ave)/ts%par%eps,10.0_wp)
                        f_scale  = exp(-dvdt_fac)

                        ! Get forcing rate of change magnitude
                        ts%df_dt = ( DF_DT_MIN + f_scale*(ts%par%df_dt_max-DF_DT_MIN) )

                    end if 

                case("PI42","H312b","H312PID","H321PID","PID1")

                    if (dt_tot .lt. ts%par%dt_ave) then 
                        ! Not enough time has passed, maintain constant forcing 
                        ! (to avoid affects of noisy derivatives)

                        ts%df_dt = 0.0_wp 

                    else 
                        ! Calculate the current forcing rate, df_dt

                        ! Calculate adaptive dfdt value using proportional-integral (PI) methods
                        call set_adaptive_timestep_pc(pi_df_now,ts%pi_df,ts%pi_eta,ts%par%eps, &
                                            DF_DT_MIN,ts%par%df_dt_max,pi_order,ts%par%method)

                        ! Remove oldest point from the end and insert latest point in beginning
                        ts%pi_eta = eoshift(ts%pi_eta,-1,boundary=abs(ts%dv_dt_ave))
                        ts%pi_df  = eoshift(ts%pi_df, -1,boundary=abs(pi_df_now))

                        ! Apply limits to eta so that algorithm works well. 
                        ! pi_eta should be greater than zero
                        ts%pi_eta(1) = max(ts%pi_eta(1),1e-3_wp)

                        ! Get forcing rate of change magnitude in [f/yr]
                        ts%df_dt = ts%pi_df(1) 

                    end if 

            end select 

        end if 

        ! Apply sign of change
        ts%df_dt = ts%par%df_sign*ts%df_dt

        ! Avoid underflow errors
        if (abs(ts%df_dt) .lt. TOL) ts%df_dt = 0.0_wp


        ! ajr: Note: for now, keep diagnosing df_dt even when limits have been reached.
        ! df_dt will not be applied beyond forcing limits though.
        ! To actually set df_dt to zero too, uncomment following lines:

        ! Set df_dt to zero if desired forcing limits have been reached 
        !if (ts%f_mean_now .le. ts%par%f_min) ts%df_dt = 0.0
        !if (ts%f_mean_now .ge. ts%par%f_max) ts%df_dt = 0.0

        if (ts%dt .gt. 0.0_wp) then 
            ! Update f_now, etc. if time step is non-zero. 

            ! Update the mean forcing value 
            ts%f_mean_now = ts%f_mean_now + (ts%df_dt*ts%dt) 

            ! Ensure f_min/f_max bounds are not exceeded
            if (ts%f_mean_now .lt. ts%par%f_min) ts%f_mean_now = ts%par%f_min 
            if (ts%f_mean_now .gt. ts%par%f_max) ts%f_mean_now = ts%par%f_max 

            ! If desired, generate some noise 
            if (ts%par%sigma .gt. 0.0) then 
                call gen_random_normal(ts%eta_now,mu=0.0_wp,sigma=ts%par%sigma)
            else 
                ts%eta_now = 0.0_wp 
            end if 

        end if 
        
        ! Update the real forcing value 
        ts%f_now = ts%f_mean_now + ts%eta_now 

        ! Check if kill should be activated 
        if (ts%par%with_kill .and. &
            abs(ts%dv_dt_ave) .lt. ts%par%eps) then 

            if (ts%par%df_sign .gt. 0.0 .and. &
                ts%f_mean_now .ge. ts%par%f_max) then 
                ts%kill = .TRUE.
            end if 

            if (ts%par%df_sign .lt. 0.0 .and. &
                ts%f_mean_now .le. ts%par%f_min) then 
                ts%kill = .TRUE.
            end if 

        end if 

        return

    end subroutine tsgen_update
            
    subroutine set_adaptive_timestep_pc(dt_new,dt,eta,eps,dtmin,dtmax,pc_k,controller)
        ! Calculate the timestep following algorithm for 
        ! a general predictor-corrector (pc) method.
        ! Implemented followig Cheng et al (2017, GMD)

        implicit none 

        real(wp), intent(OUT) :: dt_new               ! [yr]   Timestep (n+1)
        real(wp), intent(IN)  :: dt(:)                ! [yr]   Timesteps (n:n-2)
        real(wp), intent(IN)  :: eta(:)               ! [X/yr] Maximum truncation error (n:n-2)
        real(wp), intent(IN)  :: eps                  ! [--]   Tolerance value (eg, eps=1e-4)
        real(wp), intent(IN)  :: dtmin                ! [yr]   Minimum allowed timestep, must be > 0
        real(wp), intent(IN)  :: dtmax                ! [yr]   Maximum allowed timestep
        integer,    intent(IN)  :: pc_k                 ! pc_k gives the order of the timestepping scheme (pc_k=2 for FE-SBE, pc_k=3 for AB-SAM)
        character(len=*), intent(IN) :: controller      ! Adaptive controller to use [PI42, H312b, H312PID]

        ! Local variables
        real(wp) :: dt_n, dt_nm1, dt_nm2          ! [yr]   Timesteps (n:n-2)
        real(wp) :: eta_n, eta_nm1, eta_nm2       ! [X/yr] Maximum truncation error (n:n-2)
        real(wp) :: rho_n, rho_nm1, rho_nm2
        real(wp) :: rhohat_n
        real(wp) :: k_i
        real(wp) :: k_p, k_d 

        ! Smoothing parameter; Söderlind and Wang (2006) method, Eq. 10
        ! Values on the order of [0.7,2.0] are reasonable. Higher kappa slows variation in dt
        real(wp), parameter :: kappa = 2.0_wp 
        
        ! Step 1: Save information needed for adapative controller algorithms 

        ! Save dt from several timesteps (potentially more available)
        dt_n    = max(dt(1),dtmin) 
        dt_nm1  = max(dt(2),dtmin) 
        dt_nm2  = max(dt(3),dtmin)

        ! Save eta from several timesteps (potentially more available)
        eta_n   = eta(1)
        eta_nm1 = eta(2)
        eta_nm2 = eta(3)

        ! Calculate rho from several timesteps 
        rho_nm1 = (dt_n   / dt_nm1) 
        rho_nm2 = (dt_nm1 / dt_nm2) 

        ! Step 2: calculate scaling for the next timestep (dt,n+1)
        select case(trim(controller))

            case("PI42")
                ! Söderlind and Wang, 2006; Cheng et al., 2017
                ! Deeper discussion in Söderlind, 2002. Note for example, 
                ! that Söderlind (2002) recommends:
                ! k*k_i >= 0.3 
                ! k*k_p >= 0.2 
                ! k*k_i + k*k_p <= 0.7 
                ! However, the default values of Söderlind and Wang (2006) and Cheng et al (2017)
                ! are outside of these bounds. 

                ! Default parameter values 
                k_i = 2.0_wp / (pc_k*5.0_wp)
                k_p = 1.0_wp / (pc_k*5.0_wp)
                
                ! Improved parameter values (reduced oscillations)
!                 k_i = 4.0_wp / (pc_k*10.0_wp)
!                 k_p = 3.0_wp / (pc_k*10.0_wp)

                ! Experimental parameter values (minimal oscillations, does not access largest timesteps)
!                 k_i = 0.5_wp / (pc_k*10.0_wp)
!                 k_p = 6.5_wp / (pc_k*10.0_wp)

                ! Default parameter values
                rho_n = calc_pi_rho_pi42(eta_n,eta_nm1,rho_nm1,eps,k_i,k_p,alpha_2=0.0_wp)

            case("H312b") 
                ! Söderlind (2003) H312b, Eq. 31+ (unlabeled) 
                
                rho_n = calc_pi_rho_H312b(eta_n,eta_nm1,eta_nm2,rho_nm1,rho_nm2,eps,k=real(pc_k,wp),b=8.0_wp)

            case("H312PID") 
                ! Söderlind (2003) H312PD, Eq. 38
                ! Note: Suggested k_i =(2/9)*1/pc_k, but lower value gives more stable solution

                !k_i = (2.0_wp/9.0_wp)*1.0_wp/real(pc_k,wp)
                k_i = 0.08_wp/real(pc_k,wp)

                rho_n = calc_pi_rho_H312PID(eta_n,eta_nm1,eta_nm2,eps,k_i)

            case("H321PID")

                k_i = 0.1  / real(pc_k,wp)
                k_p = 0.45 / real(pc_k,wp) 

                rho_n = calc_pi_rho_H321PID(eta_n,eta_nm1,eta_nm2,dt_n,dt_nm1,eps,k_i,k_p)
                
            case("PID1")

                k_i = 0.175 
                k_p = 0.075
                k_d = 0.01 

                rho_n = calc_pi_rho_PID1(eta_n,eta_nm1,eta_nm2,eps,k_i,k_p,k_d)
                
            case DEFAULT 

                write(*,*) "set_adaptive_timestep_pc:: Error: controller not recognized."
                write(*,*) "controller = ", trim(controller) 
                stop 

        end select 

        ! Scale rho_n for smoothness 
        rhohat_n = rho_n
        !rhohat_n = min(rho_n,1.1)
        !rhohat_n = 1.0_wp + kappa * atan((rho_n-1.0_wp)/kappa) ! Söderlind and Wang, 2006, Eq. 10
        
        ! Step 3: calculate the next time timestep (dt,n+1)
        dt_new = rhohat_n * dt_n

        ! Ensure timestep is also within parameter limits 
        dt_new = max(dtmin,dt_new)  ! dt >= dtmin
        dt_new = min(dtmax,dt_new)  ! dt <= dtmax

        return 

    end subroutine set_adaptive_timestep_pc

    function calc_pi_rho_pi42(eta_n,eta_nm1,rho_nm1,eps,k_i,k_p,alpha_2) result(rho_n)

        implicit none 

        real(wp), intent(IN) :: eta_n 
        real(wp), intent(IN) :: eta_nm1 
        real(wp), intent(IN) :: rho_nm1 
        real(wp), intent(IN) :: eps 
        real(wp), intent(IN) :: k_i 
        real(wp), intent(IN) :: k_p
        real(wp), intent(IN) :: alpha_2 
        real(wp) :: rho_n 

        ! Söderlind and Wang, 2006; Cheng et al., 2017
        ! Original formulation: Söderlind, 2002, Eq. 3.12:
        rho_n   = (eps/eta_n)**(k_i+k_p) * (eps/eta_nm1)**(-k_p) * rho_nm1**(-alpha_2)

        return 

    end function calc_pi_rho_pi42

    function calc_pi_rho_H312b(eta_n,eta_nm1,eta_nm2,rho_nm1,rho_nm2,eps,k,b) result(rho_n)

        implicit none 

        real(wp), intent(IN) :: eta_n 
        real(wp), intent(IN) :: eta_nm1 
        real(wp), intent(IN) :: eta_nm2 
        real(wp), intent(IN) :: rho_nm1
        real(wp), intent(IN) :: rho_nm2
        real(wp), intent(IN) :: eps 
        real(wp), intent(IN) :: k 
        real(wp), intent(IN) :: b 
        real(wp) :: rho_n 

        ! Local variables 
        real(wp) :: beta_1, beta_2, beta_3 
        real(wp) :: alpha_2, alpha_3 

        beta_1  =  1.0_wp / (k*b)
        beta_2  =  2.0_wp / (k*b)
        beta_3  =  1.0_wp / (k*b)
        alpha_2 = -3.0_wp / b 
        alpha_3 = -1.0_wp / b 

        ! Söderlind (2003) H312b, Eq. 31+ (unlabeled) 
        rho_n   = (eps/eta_n)**beta_1 * (eps/eta_nm1)**beta_2 * (eps/eta_nm2)**beta_3 &
                            * rho_nm1**alpha_2 * rho_nm2**alpha_3 

        return 

    end function calc_pi_rho_H312b

    function calc_pi_rho_H312PID(eta_n,eta_nm1,eta_nm2,eps,k_i) result(rho_n)

        implicit none 

        real(wp), intent(IN) :: eta_n 
        real(wp), intent(IN) :: eta_nm1 
        real(wp), intent(IN) :: eta_nm2 
        real(wp), intent(IN) :: eps 
        real(wp), intent(IN) :: k_i  
        real(wp) :: rho_n 

        ! Local variables 
        real(wp) :: k_i_1, k_i_2, k_i_3

        k_i_1   = k_i / 4.0_wp 
        k_i_2   = k_i / 2.0_wp 
        k_i_3   = k_i / 4.0_wp 

        ! Söderlind (2003) H312PID, Eq. 38
        rho_n   = (eps/eta_n)**k_i_1 * (eps/eta_nm1)**k_i_2 * (eps/eta_nm2)**k_i_3

        return 

    end function calc_pi_rho_H312PID

    function calc_pi_rho_H321PID(eta_n,eta_nm1,eta_nm2,dt_n,dt_nm1,eps,k_i,k_p) result(rho_n)

        implicit none 

        real(wp), intent(IN) :: eta_n 
        real(wp), intent(IN) :: eta_nm1 
        real(wp), intent(IN) :: eta_nm2 
        real(wp), intent(IN) :: dt_n 
        real(wp), intent(IN) :: dt_nm1 
        real(wp), intent(IN) :: eps 
        real(wp), intent(IN) :: k_i  
        real(wp), intent(IN) :: k_p
        real(wp) :: rho_n 

        ! Local variables 
        real(wp) :: k_i_1, k_i_2, k_i_3

        k_i_1   =   0.75_wp*k_i + 0.50_wp*k_p 
        k_i_2   =   0.50_wp*k_i 
        k_i_3   = -(0.25_wp*k_i + 0.50_wp*k_p)

        ! Söderlind (2003) H321PID, Eq. 42
        rho_n   = (eps/eta_n)**k_i_1 * (eps/eta_nm1)**k_i_2 * (eps/eta_nm2)**k_i_3 * (dt_n / dt_nm1)

        return 

    end function calc_pi_rho_H321PID

    function calc_pi_rho_PID1(eta_n,eta_nm1,eta_nm2,eps,k_i,k_p,k_d) result(rho_n)

        implicit none 

        real(wp), intent(IN) :: eta_n 
        real(wp), intent(IN) :: eta_nm1 
        real(wp), intent(IN) :: eta_nm2 
        real(wp), intent(IN) :: eps 
        real(wp), intent(IN) :: k_i  
        real(wp), intent(IN) :: k_p
        real(wp), intent(IN) :: k_d
        real(wp) :: rho_n 

        ! https://www.mathematik.uni-dortmund.de/~kuzmin/cfdintro/lecture8.pdf
        ! Page 20 (theoretical basis unclear/unknown)
        rho_n   = (eps/eta_n)**k_i * (eta_nm1/eta_n)**k_p * (eta_nm1**2/(eta_n*eta_nm2))**k_d

        return 

    end function calc_pi_rho_PID1


    subroutine calc_sin_now(f_now,time_now,p,f_min,f_max,x_offset)

        implicit none 

        real(wp), intent(OUT) :: f_now 
        real(wp), intent(IN)  :: time_now 
        real(wp), intent(IN)  :: p
        real(wp), intent(IN)  :: f_min
        real(wp), intent(IN)  :: f_max 
        real(wp), intent(IN)  :: x_offset

        ! Local variables 
        real(wp) :: amp 
        real(wp) :: y_offset 

        ! Get sinusoidal parameters
        amp      = 0.5_wp * (f_max - f_min)
        y_offset = 0.5_wp * (f_max + f_min)

        ! Calculate expected current forcing value based on time elapsed
        f_now = amp*sin(2.0_wp*pi*(time_now+x_offset)/p) + y_offset

        return 

    end subroutine calc_sin_now


    subroutine gen_random_normal(ynrm,mu,sigma)
        ! Calculate a random number from a normal distribution 
        ! following the Box-Mueller algorithm 
        ! https://en.wikipedia.org/wiki/Normal_distribution#Generating_values_from_normal_distribution 

        implicit none 

        real(wp), intent(OUT) :: ynrm
        real(wp), intent(IN)  :: mu 
        real(wp), intent(IN)  :: sigma 

        ! Local variables
        real(wp) :: yuni(2)

        ! Get 2 numbers from uniform distribution between 0 and 1
        call random_number(yuni)
        
        ! Convert to normal distribution using the Box-Mueller algorithm
        ynrm = mu + sigma * sqrt(-2.0_wp*log(yuni(1))) * cos(2.0_wp*pi*yuni(2))

        return 

    end subroutine gen_random_normal

end module tsgen
