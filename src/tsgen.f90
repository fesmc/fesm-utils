
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
        ! Configuration (read from namelist, immutable after init)
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

        ! Derived at init: .TRUE. for response-driven (feedback) methods
        logical  :: feedback
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
        real(wp) :: time_prev       ! time seen at the previous update (MV if none)
        real(wp) :: dv_dt_ave
        real(wp) :: df_dt
        real(wp) :: f_now
        real(wp) :: f_mean_now
        real(wp) :: eta_now
        logical  :: kill
        real(wp) :: time_init
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

        ! Classify the method: response-driven (feedback) vs time-driven
        select case(trim(ts%par%method))
            case("exp","PI42","H312b","H312PID","H321PID","PID1")
                ts%par%feedback = .TRUE.
            case default
                ts%par%feedback = .FALSE.
        end select

        ! Kill is not meaningful for periodic forcing
        if (trim(ts%par%method) .eq. "sin") ts%par%with_kill = .FALSE.

        ! Define label for this tsgen object
        ts%par%label = "tsgen"
        if (present(label)) ts%par%label = trim(ts%par%label)//"_"//trim(label)

        ! The history buffers are only needed for feedback control and/or the
        ! kill switch (both use the windowed response derivative). Allocate them
        ! only in that case; time-driven forcing is stateless/analytic.
        if (allocated(ts%time))  deallocate(ts%time)
        if (allocated(ts%var))   deallocate(ts%var)
        if (allocated(ts%dv_dt)) deallocate(ts%dv_dt)
        if (ts%par%feedback .or. ts%par%with_kill) then
            ntot = 2000
            allocate(ts%time(ntot))
            allocate(ts%var(ntot))
            allocate(ts%dv_dt(ntot))
            ts%time  = MV
            ts%var   = MV
            ts%dv_dt = MV
        end if

        ts%dv_dt_ave = 0.0_wp
        ts%df_dt     = 0.0_wp

        ts%pi_df  = DF_DT_MIN
        ts%pi_eta = ts%par%eps

        ! Initial forcing value (analytic series at elapsed time zero)
        ts%f_mean_now = eval_open(ts%par,0.0_wp)
        ts%eta_now    = 0.0_wp
        ts%f_now      = ts%f_mean_now

        ! Kill switch off to start
        ts%kill = .FALSE.

        ! Store initial simulation time; no previous time seen yet
        ts%time_init = time
        ts%time_prev = MV

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
        real(wp) :: f_mean_prev

        ! Current timestep, from the time seen at the previous update
        if (ts%time_prev .ne. MV) then
            ts%dt = time - ts%time_prev
        else
            ts%dt = 0.0_wp
        end if

        ! Update the history buffers / windowed derivative when they are in use
        ! (feedback control and/or the kill switch need dv_dt_ave)
        if (ts%par%feedback .or. ts%par%with_kill) call update_history(ts,time,var,dv_dt)

        ! Elapsed time since initialization
        time_elapsed = time - ts%time_init

        ! Remember previous mean forcing (for the diagnostic rate df_dt)
        f_mean_prev = ts%f_mean_now

        if (ts%par%feedback) then
            ! Response-driven (feedback) control: modulate the forcing rate from
            ! the averaged response derivative, then integrate to update f_mean.

            if (time_elapsed .le. ts%par%dt_init) then
                ! Initialization period: no forcing change yet
                ts%df_dt = 0.0_wp
            else
                call calc_df_dt_feedback(ts)

                ! Apply sign of change and suppress underflow
                ts%df_dt = ts%par%df_sign * ts%df_dt
                if (abs(ts%df_dt) .lt. TOL) ts%df_dt = 0.0_wp
            end if

            if (ts%dt .gt. 0.0_wp) then
                ! Integrate and keep within [f_min,f_max]
                ts%f_mean_now = ts%f_mean_now + ts%df_dt*ts%dt
                ts%f_mean_now = min(max(ts%f_mean_now,ts%par%f_min),ts%par%f_max)
            end if

        else
            ! Time-driven forcing: evaluate the analytic series directly.
            ! The reported rate df_dt is the finite difference over the step.

            ts%f_mean_now = eval_open(ts%par,time_elapsed)

            if (ts%dt .gt. 0.0_wp) then
                ts%df_dt = (ts%f_mean_now - f_mean_prev)/ts%dt
            else
                ts%df_dt = 0.0_wp
            end if
        end if

        ! Optional noise on top of the mean forcing (only on a non-zero step)
        if (ts%dt .gt. 0.0_wp) then
            if (ts%par%sigma .gt. 0.0_wp) then
                call gen_random_normal(ts%eta_now,mu=0.0_wp,sigma=ts%par%sigma)
            else
                ts%eta_now = 0.0_wp
            end if
        end if

        ! Real forcing value = mean + noise
        ts%f_now = ts%f_mean_now + ts%eta_now

        ! Kill switch: stop once the response has equilibrated at a bound
        if (ts%par%with_kill .and. abs(ts%dv_dt_ave) .lt. ts%par%eps) then
            if (ts%par%df_sign .gt. 0.0_wp .and. ts%f_mean_now .ge. ts%par%f_max) ts%kill = .TRUE.
            if (ts%par%df_sign .lt. 0.0_wp .and. ts%f_mean_now .le. ts%par%f_min) ts%kill = .TRUE.
        end if

        ! Remember this time for the next call
        ts%time_prev = time

        return

    end subroutine tsgen_update

    subroutine update_history(ts,time,var,dv_dt)
        ! Push the current (time,var,dv_dt) onto the history buffers and refresh
        ! the windowed average response derivative dv_dt_ave.

        type(tsgen_class), intent(INOUT) :: ts
        real(wp),           intent(IN)    :: time
        real(wp),           intent(IN)    :: var
        real(wp), optional, intent(IN)    :: dv_dt

        ! Local variables
        integer  :: ntot, kmin, kmax, nk
        real(wp) :: dv_dt_now

        ntot = size(ts%time,1)

        ! Current response derivative
        if (present(dv_dt)) then
            dv_dt_now = dv_dt
        else if (ts%dt .gt. 0.0_wp) then
            dv_dt_now = (var-ts%var(ntot))/ts%dt
        else
            dv_dt_now = 0.0_wp
        end if

        ! Drop the oldest point and append the current one
        ts%time  = eoshift(ts%time,  1,boundary=time)
        ts%var   = eoshift(ts%var,   1,boundary=var)
        ts%dv_dt = eoshift(ts%dv_dt, 1,boundary=dv_dt_now)

        ! Indices within the averaging window [time-dt_ave, time].
        ! (findloc needs gfortran>=9, so use minloc with a mask instead.)
        kmin = minloc(ts%time,dim=1, &
                mask=(ts%time .ge. time - ts%par%dt_ave) .and. ts%time.ne.MV)
        kmax = ntot
        nk   = kmax - kmin + 1

        ! Windowed average of the response derivative
        ts%dv_dt_ave = sum(ts%dv_dt(kmin:kmax)) / real(nk,wp)

        return

    end subroutine update_history

    function eval_open(par,time_elapsed) result(f_mean)
        ! Analytic (time-driven) forcing value as a pure function of elapsed time.
        ! Shared by tsgen_update (time-driven methods) and tsgen_tabulate.

        type(tsgen_par_class), intent(IN) :: par
        real(wp),              intent(IN) :: time_elapsed
        real(wp) :: f_mean

        ! Local variables
        real(wp) :: tau, rate, f_start, f_ext

        ! Time since the end of the initialization period
        tau = max(0.0_wp, time_elapsed - par%dt_init)

        ! Forcing value at the start of the ramp (depends on ramp direction)
        f_start = merge(par%f_min, par%f_max, par%df_sign .gt. 0.0_wp)

        select case(trim(par%method))

            case("const","ramp-slope")
                ! Constant-rate ramp (df_dt_max), clamped to [f_min,f_max]
                f_mean = f_start + par%df_sign*par%df_dt_max*tau
                f_mean = min(max(f_mean,par%f_min),par%f_max)

            case("ramp-time")
                ! Linear ramp from f_start to the far bound over dt_ramp
                rate   = abs(par%f_max-par%f_min)/par%dt_ramp
                f_mean = f_start + par%df_sign*rate*tau
                f_mean = min(max(f_mean,par%f_min),par%f_max)

            case("ramp-time-step")
                ! Triangle waveform: f_start -> first extreme (over dt_ramp)
                !                             -> f_conv       (over dt_conv) -> hold
                f_ext = merge(par%f_max, par%f_min, par%df_sign .gt. 0.0_wp)
                if (tau .le. par%dt_ramp) then
                    f_mean = f_start + (f_ext-f_start)*(tau/par%dt_ramp)
                else if (tau .le. par%dt_ramp+par%dt_conv) then
                    f_mean = f_ext + (par%f_conv-f_ext)*((tau-par%dt_ramp)/par%dt_conv)
                else
                    f_mean = par%f_conv
                end if

            case("sin")
                ! Periodic forcing; held flat during the initialization period
                if (time_elapsed .le. par%dt_init) then
                    call calc_sin_now(f_mean,0.0_wp,par%dt_ramp, &
                                        par%f_min,par%f_max,x_offset=0.0_wp)
                else
                    call calc_sin_now(f_mean,time_elapsed,par%dt_ramp, &
                                        par%f_min,par%f_max,x_offset=0.0_wp)
                end if

            case default
                ! Feedback methods have no analytic form: report the start value
                f_mean = f_start

        end select

        return

    end function eval_open

    subroutine calc_df_dt_feedback(ts)
        ! Response-driven forcing-rate magnitude (unsigned), stored in ts%df_dt.
        ! Uses the windowed response derivative dv_dt_ave and therefore requires
        ! the history buffer, which is allocated for feedback methods.

        type(tsgen_class), intent(INOUT) :: ts

        ! Local variables
        real(wp) :: dt_tot, dvdt_fac, f_scale, pi_df_now

        ! dv_dt is a windowed average, so use 2nd-order PI controller parameters
        integer, parameter :: pi_order = 2

        ! Time span currently covered by the history buffer
        dt_tot = ts%time(size(ts%time,1)) - minval(ts%time,mask=ts%time.ne.MV)

        if (dt_tot .lt. ts%par%dt_ave) then
            ! Not enough time has passed; hold the forcing to avoid reacting
            ! to noisy derivatives.
            ts%df_dt = 0.0_wp
            return
        end if

        select case(trim(ts%par%method))

            case("exp")
                ! Exponential response (sharp, tuneable): scale in [0,1], ~0.6 at
                ! dv_dt==eps. Cap the exponent so exp() cannot overflow.
                dvdt_fac = min(abs(ts%dv_dt_ave)/ts%par%eps,10.0_wp)
                f_scale  = exp(-dvdt_fac)
                ts%df_dt = DF_DT_MIN + f_scale*(ts%par%df_dt_max-DF_DT_MIN)

            case default    ! PI42, H312b, H312PID, H321PID, PID1
                ! Adaptive rate via proportional-integral(-derivative) controller
                call calc_forcing_rate_pc(pi_df_now,ts%pi_df,ts%pi_eta,ts%par%eps, &
                                    DF_DT_MIN,ts%par%df_dt_max,pi_order,ts%par%method)

                ! Push newest error/rate to the front of the controller history
                ts%pi_eta = eoshift(ts%pi_eta,-1,boundary=abs(ts%dv_dt_ave))
                ts%pi_df  = eoshift(ts%pi_df, -1,boundary=abs(pi_df_now))

                ! Keep the newest error strictly positive for stability
                ts%pi_eta(1) = max(ts%pi_eta(1),1e-3_wp)

                ts%df_dt = ts%pi_df(1)

        end select

        return

    end subroutine calc_df_dt_feedback

    subroutine calc_forcing_rate_pc(df_new,df,eta,eps,df_min,df_max,pc_k,controller)
        ! Adaptive forcing-rate update following the general predictor-corrector
        ! (pc) controller algorithms of Cheng et al. (2017, GMD). The forcing
        ! increment df plays the role of the "timestep" in the original scheme.

        implicit none

        real(wp), intent(OUT) :: df_new               ! [f/yr] Forcing rate (n+1)
        real(wp), intent(IN)  :: df(:)                ! [f/yr] Forcing rates (n:n-2)
        real(wp), intent(IN)  :: eta(:)               ! [X/yr] Error signal, |dv_dt_ave| (n:n-2)
        real(wp), intent(IN)  :: eps                  ! [--]   Tolerance value (eg, eps=1e-4)
        real(wp), intent(IN)  :: df_min               ! [f/yr] Minimum allowed rate, must be > 0
        real(wp), intent(IN)  :: df_max               ! [f/yr] Maximum allowed rate
        integer,    intent(IN)  :: pc_k                 ! pc_k gives the order of the timestepping scheme (pc_k=2 for FE-SBE, pc_k=3 for AB-SAM)
        character(len=*), intent(IN) :: controller      ! Adaptive controller to use [PI42, H312b, H312PID]

        ! Local variables
        real(wp) :: df_n, df_nm1, df_nm2          ! [f/yr] Forcing rates (n:n-2)
        real(wp) :: eta_n, eta_nm1, eta_nm2       ! [X/yr] Error signal (n:n-2)
        real(wp) :: rho_n, rho_nm1, rho_nm2
        real(wp) :: rhohat_n
        real(wp) :: k_i
        real(wp) :: k_p, k_d

        ! Smoothing parameter; Söderlind and Wang (2006) method, Eq. 10
        ! Values on the order of [0.7,2.0] are reasonable. Higher kappa slows variation in df
        real(wp), parameter :: kappa = 2.0_wp

        ! Step 1: Save information needed for adapative controller algorithms

        ! Save df from several steps (potentially more available)
        df_n    = max(df(1),df_min)
        df_nm1  = max(df(2),df_min)
        df_nm2  = max(df(3),df_min)

        ! Save eta from several steps (potentially more available)
        eta_n   = eta(1)
        eta_nm1 = eta(2)
        eta_nm2 = eta(3)

        ! Calculate rho from several steps
        rho_nm1 = (df_n   / df_nm1)
        rho_nm2 = (df_nm1 / df_nm2)

        ! Step 2: calculate scaling for the next rate (df,n+1)
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

                rho_n = calc_pi_rho_H321PID(eta_n,eta_nm1,eta_nm2,df_n,df_nm1,eps,k_i,k_p)
                
            case("PID1")

                k_i = 0.175 
                k_p = 0.075
                k_d = 0.01 

                rho_n = calc_pi_rho_PID1(eta_n,eta_nm1,eta_nm2,eps,k_i,k_p,k_d)
                
            case DEFAULT 

                write(*,*) "calc_forcing_rate_pc:: Error: controller not recognized."
                write(*,*) "controller = ", trim(controller)
                stop

        end select

        ! Scale rho_n for smoothness
        rhohat_n = rho_n
        !rhohat_n = min(rho_n,1.1)
        !rhohat_n = 1.0_wp + kappa * atan((rho_n-1.0_wp)/kappa) ! Söderlind and Wang, 2006, Eq. 10

        ! Step 3: calculate the next forcing rate (df,n+1)
        df_new = rhohat_n * df_n

        ! Ensure rate is also within parameter limits
        df_new = max(df_min,df_new)  ! df >= df_min
        df_new = min(df_max,df_new)  ! df <= df_max

        return

    end subroutine calc_forcing_rate_pc

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
