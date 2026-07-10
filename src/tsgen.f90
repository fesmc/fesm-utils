
module tsgen

    use nml
    use ncio
    use series

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
        real(wp) :: tol
        real(wp) :: df_dt_max
        real(wp) :: sigma
        real(wp) :: f_min
        real(wp) :: f_max
        real(wp) :: f_conv

        ! Tabulated-series method ("series"): file + (netCDF) variable names
        character(len=512) :: series_file
        character(len=56)  :: series_var
        character(len=56)  :: series_time_name

        ! Derived at init: .TRUE. for response-driven (feedback) methods
        logical  :: feedback
    end type

    type tsgen_vec_class
        ! Per-channel outputs, allocated only when nc > 1 (a multi-channel
        ! "series" file, e.g. 12 monthly values). For nc == 1 the scalar
        ! outputs on tsgen_class are used and these stay unallocated.
        real(wp), allocatable :: f_now(:)
        real(wp), allocatable :: f_mean(:)
        real(wp), allocatable :: eps(:)
    end type

    type tsgen_class

        type(tsgen_par_class) :: par

        ! Loaded forcing series (only for method="series")
        type(series_class) :: ser

        ! History buffers over the time-averaging window
        real(wp), allocatable :: time(:)
        real(wp), allocatable :: var(:)
        real(wp), allocatable :: dv_dt(:)

        ! PI-controller history (n:n-2)
        real(wp) :: pi_df(3)
        real(wp) :: pi_eta(3)

        ! Runtime state
        integer  :: nc              ! number of output channels (1 except a 2-D series)
        real(wp) :: dt
        real(wp) :: time_prev       ! time seen at the previous update (MV if none)
        real(wp) :: dv_dt_ave
        real(wp) :: df_dt
        real(wp) :: f_now           ! scalar forcing (channel mean when nc>1)
        real(wp) :: f_mean_now
        real(wp) :: eps             ! realized noise (channel mean when nc>1)
        logical  :: kill
        real(wp) :: time_init

        ! Per-channel outputs (allocated only when nc>1)
        type(tsgen_vec_class) :: vec
    end type

    private
    public :: tsgen_class
    public :: tsgen_init
    public :: tsgen_update
    public :: tsgen_tabulate
    public :: tsgen_write
    public :: tsgen_write_step
    public :: tsgen_restart_write
    public :: tsgen_restart_read

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

        ! The method determines which parameters are relevant, so read it first.
        call nml_read(filename,trim(par_label),"method",      ts%par%method)

        ! Noise amplitude applies to every method (it is added on top of the
        ! mean forcing), so it is always read.
        call nml_read(filename,trim(par_label),"sigma",       ts%par%sigma)

        ! Classify the method: response-driven (feedback) vs time-driven
        select case(trim(ts%par%method))
            case("exp","PI42","H312b","H312PID","H321PID","PID1")
                ts%par%feedback = .TRUE.
            case default
                ts%par%feedback = .FALSE.
        end select

        if (trim(ts%par%method) .eq. "series") then
            ! Tabulated series: read only the file/variable parameters, load it,
            ! and derive the channel count. Ramp/feedback parameters are unused.
            call tsgen_par_defaults(ts%par)
            ts%par%with_kill = .FALSE.

            call nml_read(filename,trim(par_label),"series_file",ts%par%series_file)

            if (series_is_netcdf(ts%par%series_file)) then
                call nml_read(filename,trim(par_label),"series_var",      ts%par%series_var)
                call nml_read(filename,trim(par_label),"series_time_name",ts%par%series_time_name)
                call series_load(ts%ser,ts%par%series_file, &
                                    varname=ts%par%series_var, &
                                    time_name=ts%par%series_time_name, &
                                    sigma_name=trim(ts%par%series_var)//"_sd")
            else
                call series_load(ts%ser,ts%par%series_file)
            end if

            ts%nc = ts%ser%nc

        else
            ! Time-driven (analytic) or response-driven (feedback) scalar forcing
            call nml_read(filename,trim(par_label),"with_kill",   ts%par%with_kill)
            call nml_read(filename,trim(par_label),"dt_ave",      ts%par%dt_ave)
            call nml_read(filename,trim(par_label),"dt_init",     ts%par%dt_init)
            call nml_read(filename,trim(par_label),"dt_ramp",     ts%par%dt_ramp)
            call nml_read(filename,trim(par_label),"dt_conv",     ts%par%dt_conv)
            call nml_read(filename,trim(par_label),"df_sign",     ts%par%df_sign)
            call nml_read(filename,trim(par_label),"tol",         ts%par%tol)
            call nml_read(filename,trim(par_label),"df_dt_max",   ts%par%df_dt_max)
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

            ! Kill is not meaningful for periodic forcing
            if (trim(ts%par%method) .eq. "sin") ts%par%with_kill = .FALSE.

            ts%nc = 1
        end if

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

        ! Allocate per-channel output buffers for a multi-channel series
        if (allocated(ts%vec%f_now))  deallocate(ts%vec%f_now)
        if (allocated(ts%vec%f_mean)) deallocate(ts%vec%f_mean)
        if (allocated(ts%vec%eps))    deallocate(ts%vec%eps)
        if (ts%nc .gt. 1) then
            allocate(ts%vec%f_now(ts%nc))
            allocate(ts%vec%f_mean(ts%nc))
            allocate(ts%vec%eps(ts%nc))
            ts%vec%eps = 0.0_wp
        end if

        ts%dv_dt_ave = 0.0_wp
        ts%df_dt     = 0.0_wp

        ts%pi_df  = DF_DT_MIN
        ts%pi_eta = ts%par%tol

        ! Initial forcing value
        if (trim(ts%par%method) .eq. "series") then
            ! Series interpolates on the absolute time axis stored in the file,
            ! so sample it at the initial simulation time.
            call eval_series(ts,time,noise=.FALSE.)
        else
            ! Analytic series at elapsed time zero
            ts%f_mean_now = eval_open(ts%par,0.0_wp)
            ts%eps        = 0.0_wp
            ts%f_now      = ts%f_mean_now
        end if

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

        ! Remember previous mean forcing (for the diagnostic rate df_dt)
        f_mean_prev = ts%f_mean_now

        if (trim(ts%par%method) .eq. "series") then
            ! Tabulated series: interpolate on the file's own (absolute) axis.
            call eval_series(ts,time,noise=.TRUE.)

            if (ts%dt .gt. 0.0_wp) then
                ts%df_dt = (ts%f_mean_now - f_mean_prev)/ts%dt
            else
                ts%df_dt = 0.0_wp
            end if

            ts%time_prev = time
            return
        end if

        ! Update the history buffers / windowed derivative when they are in use
        ! (feedback control and/or the kill switch need dv_dt_ave)
        if (ts%par%feedback .or. ts%par%with_kill) call update_history(ts,time,var,dv_dt)

        ! Elapsed time since initialization
        time_elapsed = time - ts%time_init

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
                call gen_random_normal(ts%eps,mu=0.0_wp,sigma=ts%par%sigma)
            else
                ts%eps = 0.0_wp
            end if
        end if

        ! Real forcing value = mean + noise
        ts%f_now = ts%f_mean_now + ts%eps

        ! Kill switch: stop once the response has equilibrated at a bound
        if (ts%par%with_kill .and. abs(ts%dv_dt_ave) .lt. ts%par%tol) then
            if (ts%par%df_sign .gt. 0.0_wp .and. ts%f_mean_now .ge. ts%par%f_max) ts%kill = .TRUE.
            if (ts%par%df_sign .lt. 0.0_wp .and. ts%f_mean_now .le. ts%par%f_min) ts%kill = .TRUE.
        end if

        ! Remember this time for the next call
        ts%time_prev = time

        return

    end subroutine tsgen_update

    subroutine eval_series(ts,time,noise)
        ! Evaluate the tabulated forcing series at the (absolute) time `time`,
        ! filling the scalar outputs (channel mean when nc>1) and, for a
        ! multi-channel series, the per-channel `vec` buffers.

        type(tsgen_class), intent(INOUT) :: ts
        real(wp),           intent(IN)    :: time
        logical,            intent(IN)    :: noise

        ! Local automatic arrays sized by the channel count
        real(wp) :: vmean(ts%nc)
        real(wp) :: vsig(ts%nc)
        real(wp) :: vnoise(ts%nc)
        integer  :: c

        vmean = series_interp(ts%ser,time)

        ! Per-channel noise (only on a non-zero step): use the series' own
        ! standard deviation if present, otherwise fall back to par%sigma.
        vnoise = 0.0_wp
        if (noise .and. ts%dt .gt. 0.0_wp) then
            vsig = series_interp_sig(ts%ser,time)
            do c = 1, ts%nc
                if (vsig(c) .le. 0.0_wp) vsig(c) = ts%par%sigma
                if (vsig(c) .gt. 0.0_wp) call gen_random_normal(vnoise(c),mu=0.0_wp,sigma=vsig(c))
            end do
        end if

        ! Scalar outputs (channel mean for a multi-channel series)
        ts%f_mean_now = sum(vmean)  / real(ts%nc,wp)
        ts%eps        = sum(vnoise) / real(ts%nc,wp)
        ts%f_now      = ts%f_mean_now + ts%eps

        ! Per-channel outputs
        if (ts%nc .gt. 1) then
            ts%vec%f_mean = vmean
            ts%vec%eps    = vnoise
            ts%vec%f_now  = vmean + vnoise
        end if

        return

    end subroutine eval_series

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
                ! dv_dt==tol. Cap the exponent so exp() cannot overflow.
                dvdt_fac = min(abs(ts%dv_dt_ave)/ts%par%tol,10.0_wp)
                f_scale  = exp(-dvdt_fac)
                ts%df_dt = DF_DT_MIN + f_scale*(ts%par%df_dt_max-DF_DT_MIN)

            case default    ! PI42, H312b, H312PID, H321PID, PID1
                ! Adaptive rate via proportional-integral(-derivative) controller
                call calc_forcing_rate_pc(pi_df_now,ts%pi_df,ts%pi_eta,ts%par%tol, &
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

    subroutine tsgen_tabulate(ts,time,f,df_dt,verbose)
        ! Sample the (time-driven) forcing series over a caller-supplied time
        ! axis, for validation / diagnostics / output. Valid for the analytic
        ! and tabulated-series methods (the latter reported as the channel mean);
        ! feedback methods have no closed form and must be stepped online via
        ! tsgen_update.

        type(tsgen_class),  intent(IN)  :: ts
        real(wp),           intent(IN)  :: time(:)      ! [yr] simulation-time axis
        real(wp),           intent(OUT) :: f(:)         ! [f]  mean forcing at each time
        real(wp), optional, intent(OUT) :: df_dt(:)     ! [f/yr] finite-difference rate
        logical,  optional, intent(IN)  :: verbose      ! print a short summary

        ! Local variables
        integer  :: i, n
        logical  :: verb
        logical  :: is_series

        verb = .FALSE.
        if (present(verbose)) verb = verbose

        if (ts%par%feedback) then
            write(*,*) "tsgen_tabulate:: Error: method '"//trim(ts%par%method)// &
                "' is response-driven and has no analytic series. Step it with tsgen_update."
            stop
        end if

        is_series = (trim(ts%par%method) .eq. "series")

        n = size(time,1)

        do i = 1, n
            if (is_series) then
                ! Series carries its own absolute time axis (no init offset);
                ! report the channel mean to match the scalar f_now semantics.
                f(i) = sum(series_interp(ts%ser,time(i))) / real(ts%nc,wp)
            else
                f(i) = eval_open(ts%par, time(i) - ts%time_init)
            end if
        end do

        if (present(df_dt)) then
            ! Central finite differences (one-sided at the ends)
            if (n .gt. 1) then
                df_dt(1) = (f(2)-f(1)) / (time(2)-time(1))
                do i = 2, n-1
                    df_dt(i) = (f(i+1)-f(i-1)) / (time(i+1)-time(i-1))
                end do
                df_dt(n) = (f(n)-f(n-1)) / (time(n)-time(n-1))
            else
                df_dt = 0.0_wp
            end if
        end if

        if (verb) then
            write(*,"(a)")        "tsgen_tabulate:: "//trim(ts%par%label)// &
                                    " (method="//trim(ts%par%method)//")"
            write(*,"(a,i0)")     "  n points       = ", n
            write(*,"(a,2g14.6)") "  time range     = ", time(1), time(n)
            write(*,"(a,3g14.6)") "  f min/mean/max = ", minval(f), sum(f)/real(n,wp), maxval(f)
        end if

        return

    end subroutine tsgen_tabulate

    subroutine tsgen_write(filename,ts,time,f,df_dt)
        ! Write a tabulated forcing series (from tsgen_tabulate) to a netCDF file.

        character(len=*),   intent(IN) :: filename
        type(tsgen_class),  intent(IN) :: ts
        real(wp),           intent(IN) :: time(:)
        real(wp),           intent(IN) :: f(:)
        real(wp), optional, intent(IN) :: df_dt(:)

        call nc_create(filename)
        call nc_write_attr(filename,"tsgen_label", trim(ts%par%label))
        call nc_write_attr(filename,"tsgen_method",trim(ts%par%method))

        call nc_write_dim(filename,"time",x=time,units="yr")
        call nc_write(filename,"f",f,dim1="time", &
                        long_name="Transient forcing value")
        if (present(df_dt)) then
            call nc_write(filename,"df_dt",df_dt,dim1="time", &
                        long_name="Forcing rate of change",units="1/yr")
        end if

        return

    end subroutine tsgen_write

    subroutine tsgen_write_step(filename,ts,time,label)
        ! Append the current forcing state (f_now, df_dt, dv_dt_ave) to an
        ! existing 1D timeseries file at the time index for `time`. The file must
        ! already carry a "time" dimension (the caller's own timeseries file); the
        ! variable names are prefixed by the tsgen label (or `label`, if given) so
        ! several tsgen objects can share one file without clashing.

        character(len=*),  intent(IN) :: filename
        type(tsgen_class), intent(IN) :: ts
        real(wp),          intent(IN) :: time
        character(len=*),  intent(IN), optional :: label

        integer            :: ncid, n
        character(len=64)  :: pre

        pre = trim(ts%par%label)
        if (present(label)) pre = trim(label)

        call nc_open(filename, ncid, writable=.TRUE.)
        n = nc_time_index(filename, "time", time, ncid)

        call nc_write(filename,"time",time,dim1="time",start=[n],count=[1],ncid=ncid)
        call nc_write(filename,trim(pre)//"_f",ts%f_now,dim1="time", &
                        start=[n],count=[1],ncid=ncid, &
                        long_name="Transient forcing value")
        call nc_write(filename,trim(pre)//"_df_dt",ts%df_dt,dim1="time", &
                        start=[n],count=[1],ncid=ncid,units="1/yr", &
                        long_name="Forcing rate of change")
        call nc_write(filename,trim(pre)//"_dv_dt",ts%dv_dt_ave,dim1="time", &
                        start=[n],count=[1],ncid=ncid, &
                        long_name="Windowed response rate of change")

        call nc_close(ncid)

        return

    end subroutine tsgen_write_step

    subroutine tsgen_restart_write(filename,ts,time)
        ! Write the prognostic state of a tsgen object to a self-contained restart
        ! file. Time-driven methods (ramp/sin/const) resume exactly from time_init;
        ! feedback methods (exp, PI/PID) additionally need the integrated forcing
        ! (f_mean_now), the PI-controller history (pi_df/pi_eta) and the response
        ! history buffers, which are written when allocated. "series" carries no
        ! prognostic state (it interpolates the file on the absolute axis).

        character(len=*),  intent(IN) :: filename
        type(tsgen_class), intent(IN) :: ts
        real(wp),          intent(IN) :: time

        integer  :: ncid
        real(wp) :: kill_r

        kill_r = 0.0_wp
        if (ts%kill) kill_r = 1.0_wp

        call nc_create(filename)
        call nc_write_attr(filename,"tsgen_label", trim(ts%par%label))
        call nc_write_attr(filename,"tsgen_method",trim(ts%par%method))

        call nc_write_dim(filename,"time",x=time,dx=1.0_wp,nx=1,units="year",unlimited=.TRUE.)
        call nc_write_dim(filename,"pt3", x=1,dx=1,nx=3)

        call nc_open(filename, ncid, writable=.TRUE.)

        call nc_write(filename,"time",       time,           dim1="time",start=[1],count=[1],ncid=ncid)
        call nc_write(filename,"time_init",  ts%time_init,   dim1="time",start=[1],count=[1],ncid=ncid)
        call nc_write(filename,"time_prev",  ts%time_prev,   dim1="time",start=[1],count=[1],ncid=ncid)
        call nc_write(filename,"f_mean_now", ts%f_mean_now,  dim1="time",start=[1],count=[1],ncid=ncid)
        call nc_write(filename,"f_now",      ts%f_now,       dim1="time",start=[1],count=[1],ncid=ncid)
        call nc_write(filename,"eps",        ts%eps,         dim1="time",start=[1],count=[1],ncid=ncid)
        call nc_write(filename,"df_dt",      ts%df_dt,       dim1="time",start=[1],count=[1],ncid=ncid)
        call nc_write(filename,"dv_dt_ave",  ts%dv_dt_ave,   dim1="time",start=[1],count=[1],ncid=ncid)
        call nc_write(filename,"kill",       kill_r,         dim1="time",start=[1],count=[1],ncid=ncid)

        call nc_write(filename,"pi_df",  ts%pi_df,  dim1="pt3",ncid=ncid)
        call nc_write(filename,"pi_eta", ts%pi_eta, dim1="pt3",ncid=ncid)

        ! Response history buffers (feedback / kill switch only)
        if (allocated(ts%time)) then
            call nc_write_dim(filename,"hist",x=1,dx=1,nx=size(ts%time),ncid=ncid)
            call nc_write(filename,"hist_time", ts%time, dim1="hist",ncid=ncid)
            call nc_write(filename,"hist_var",  ts%var,  dim1="hist",ncid=ncid)
            call nc_write(filename,"hist_dv_dt",ts%dv_dt,dim1="hist",ncid=ncid)
        end if

        call nc_close(ncid)

        write(*,*) "tsgen_restart_write:: wrote "//trim(filename)

        return

    end subroutine tsgen_restart_write

    subroutine tsgen_restart_read(ts,filename)
        ! Restore tsgen state written by tsgen_restart_write. Must be called AFTER
        ! tsgen_init (which reads the namelist, classifies the method and allocates
        ! the history buffers); this routine overwrites the runtime state so the
        ! series continues from where it stopped. Missing variables are tolerated
        ! (older restart files, or a method that never wrote buffers).

        type(tsgen_class), intent(INOUT) :: ts
        character(len=*),  intent(IN)    :: filename

        real(wp) :: kill_r

        if (.not. nc_exists_var(filename,"f_mean_now")) then
            write(*,*) "tsgen_restart_read:: WARNING: no tsgen state in "// &
                trim(filename)//"; keeping cold-start values."
            return
        end if

        call nc_read(filename,"time_init",  ts%time_init,  start=[1],count=[1])
        call nc_read(filename,"time_prev",  ts%time_prev,  start=[1],count=[1])
        call nc_read(filename,"f_mean_now", ts%f_mean_now, start=[1],count=[1])
        call nc_read(filename,"f_now",      ts%f_now,      start=[1],count=[1])
        call nc_read(filename,"eps",        ts%eps,        start=[1],count=[1])
        call nc_read(filename,"df_dt",      ts%df_dt,      start=[1],count=[1])
        call nc_read(filename,"dv_dt_ave",  ts%dv_dt_ave,  start=[1],count=[1])

        call nc_read(filename,"kill",       kill_r,        start=[1],count=[1])
        ts%kill = (kill_r .gt. 0.5_wp)

        if (nc_exists_var(filename,"pi_df"))  call nc_read(filename,"pi_df", ts%pi_df)
        if (nc_exists_var(filename,"pi_eta")) call nc_read(filename,"pi_eta",ts%pi_eta)

        if (allocated(ts%time) .and. nc_exists_var(filename,"hist_time")) then
            call nc_read(filename,"hist_time", ts%time, start=[1],count=[size(ts%time)])
            call nc_read(filename,"hist_var",  ts%var,  start=[1],count=[size(ts%var)])
            call nc_read(filename,"hist_dv_dt",ts%dv_dt,start=[1],count=[size(ts%dv_dt)])
        end if

        write(*,*) "tsgen_restart_read:: restored f_now = ", ts%f_now, " from "//trim(filename)

        return

    end subroutine tsgen_restart_read

    subroutine tsgen_par_defaults(par)
        ! Neutral values for the ramp/feedback parameters, used by methods (e.g.
        ! "series") that do not read them, so downstream code never sees
        ! uninitialized fields.

        type(tsgen_par_class), intent(INOUT) :: par

        par%with_kill = .FALSE.
        par%dt_init   = 0.0_wp
        par%dt_ramp   = 1.0_wp
        par%dt_conv   = 1.0_wp
        par%dt_ave    = 1.0_wp
        par%df_sign   = 1.0_wp
        par%tol       = 1.0_wp
        par%df_dt_max = 0.0_wp
        par%f_min     = 0.0_wp
        par%f_max     = 0.0_wp
        par%f_conv    = 0.0_wp

        return

    end subroutine tsgen_par_defaults

    logical function series_is_netcdf(filename)
        ! .TRUE. when the filename extension denotes a netCDF file.

        character(len=*), intent(IN) :: filename
        integer :: idot

        series_is_netcdf = .FALSE.
        idot = index(filename,".",back=.TRUE.)
        if (idot .gt. 0) then
            select case(trim(adjustl(filename(idot+1:))))
                case("nc","nc4","cdf","netcdf")
                    series_is_netcdf = .TRUE.
            end select
        end if

        return

    end function series_is_netcdf

    subroutine calc_forcing_rate_pc(df_new,df,eta,tol,df_min,df_max,pc_k,controller)
        ! Adaptive forcing-rate update following the general predictor-corrector
        ! (pc) controller algorithms of Cheng et al. (2017, GMD). The forcing
        ! increment df plays the role of the "timestep" in the original scheme.

        implicit none

        real(wp), intent(OUT) :: df_new               ! [f/yr] Forcing rate (n+1)
        real(wp), intent(IN)  :: df(:)                ! [f/yr] Forcing rates (n:n-2)
        real(wp), intent(IN)  :: eta(:)               ! [X/yr] Error signal, |dv_dt_ave| (n:n-2)
        real(wp), intent(IN)  :: tol                  ! [--]   Tolerance value (eg, tol=1e-4)
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
                rho_n = calc_pi_rho_pi42(eta_n,eta_nm1,rho_nm1,tol,k_i,k_p,alpha_2=0.0_wp)

            case("H312b")
                ! Söderlind (2003) H312b, Eq. 31+ (unlabeled)

                rho_n = calc_pi_rho_H312b(eta_n,eta_nm1,eta_nm2,rho_nm1,rho_nm2,tol,k=real(pc_k,wp),b=8.0_wp)

            case("H312PID")
                ! Söderlind (2003) H312PD, Eq. 38
                ! Note: Suggested k_i =(2/9)*1/pc_k, but lower value gives more stable solution

                !k_i = (2.0_wp/9.0_wp)*1.0_wp/real(pc_k,wp)
                k_i = 0.08_wp/real(pc_k,wp)

                rho_n = calc_pi_rho_H312PID(eta_n,eta_nm1,eta_nm2,tol,k_i)

            case("H321PID")

                k_i = 0.1  / real(pc_k,wp)
                k_p = 0.45 / real(pc_k,wp)

                rho_n = calc_pi_rho_H321PID(eta_n,eta_nm1,eta_nm2,df_n,df_nm1,tol,k_i,k_p)

            case("PID1")

                k_i = 0.175
                k_p = 0.075
                k_d = 0.01

                rho_n = calc_pi_rho_PID1(eta_n,eta_nm1,eta_nm2,tol,k_i,k_p,k_d)

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

    function calc_pi_rho_pi42(eta_n,eta_nm1,rho_nm1,tol,k_i,k_p,alpha_2) result(rho_n)

        implicit none

        real(wp), intent(IN) :: eta_n
        real(wp), intent(IN) :: eta_nm1
        real(wp), intent(IN) :: rho_nm1
        real(wp), intent(IN) :: tol
        real(wp), intent(IN) :: k_i
        real(wp), intent(IN) :: k_p
        real(wp), intent(IN) :: alpha_2
        real(wp) :: rho_n

        ! Söderlind and Wang, 2006; Cheng et al., 2017
        ! Original formulation: Söderlind, 2002, Eq. 3.12:
        rho_n   = (tol/eta_n)**(k_i+k_p) * (tol/eta_nm1)**(-k_p) * rho_nm1**(-alpha_2)

        return

    end function calc_pi_rho_pi42

    function calc_pi_rho_H312b(eta_n,eta_nm1,eta_nm2,rho_nm1,rho_nm2,tol,k,b) result(rho_n)

        implicit none

        real(wp), intent(IN) :: eta_n
        real(wp), intent(IN) :: eta_nm1
        real(wp), intent(IN) :: eta_nm2
        real(wp), intent(IN) :: rho_nm1
        real(wp), intent(IN) :: rho_nm2
        real(wp), intent(IN) :: tol
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
        rho_n   = (tol/eta_n)**beta_1 * (tol/eta_nm1)**beta_2 * (tol/eta_nm2)**beta_3 &
                            * rho_nm1**alpha_2 * rho_nm2**alpha_3

        return

    end function calc_pi_rho_H312b

    function calc_pi_rho_H312PID(eta_n,eta_nm1,eta_nm2,tol,k_i) result(rho_n)

        implicit none

        real(wp), intent(IN) :: eta_n
        real(wp), intent(IN) :: eta_nm1
        real(wp), intent(IN) :: eta_nm2
        real(wp), intent(IN) :: tol
        real(wp), intent(IN) :: k_i
        real(wp) :: rho_n

        ! Local variables
        real(wp) :: k_i_1, k_i_2, k_i_3

        k_i_1   = k_i / 4.0_wp
        k_i_2   = k_i / 2.0_wp
        k_i_3   = k_i / 4.0_wp

        ! Söderlind (2003) H312PID, Eq. 38
        rho_n   = (tol/eta_n)**k_i_1 * (tol/eta_nm1)**k_i_2 * (tol/eta_nm2)**k_i_3

        return

    end function calc_pi_rho_H312PID

    function calc_pi_rho_H321PID(eta_n,eta_nm1,eta_nm2,dt_n,dt_nm1,tol,k_i,k_p) result(rho_n)

        implicit none

        real(wp), intent(IN) :: eta_n
        real(wp), intent(IN) :: eta_nm1
        real(wp), intent(IN) :: eta_nm2
        real(wp), intent(IN) :: dt_n
        real(wp), intent(IN) :: dt_nm1
        real(wp), intent(IN) :: tol
        real(wp), intent(IN) :: k_i
        real(wp), intent(IN) :: k_p
        real(wp) :: rho_n

        ! Local variables
        real(wp) :: k_i_1, k_i_2, k_i_3

        k_i_1   =   0.75_wp*k_i + 0.50_wp*k_p
        k_i_2   =   0.50_wp*k_i
        k_i_3   = -(0.25_wp*k_i + 0.50_wp*k_p)

        ! Söderlind (2003) H321PID, Eq. 42
        rho_n   = (tol/eta_n)**k_i_1 * (tol/eta_nm1)**k_i_2 * (tol/eta_nm2)**k_i_3 * (dt_n / dt_nm1)

        return

    end function calc_pi_rho_H321PID

    function calc_pi_rho_PID1(eta_n,eta_nm1,eta_nm2,tol,k_i,k_p,k_d) result(rho_n)

        implicit none

        real(wp), intent(IN) :: eta_n
        real(wp), intent(IN) :: eta_nm1
        real(wp), intent(IN) :: eta_nm2
        real(wp), intent(IN) :: tol
        real(wp), intent(IN) :: k_i
        real(wp), intent(IN) :: k_p
        real(wp), intent(IN) :: k_d
        real(wp) :: rho_n

        ! https://www.mathematik.uni-dortmund.de/~kuzmin/cfdintro/lecture8.pdf
        ! Page 20 (theoretical basis unclear/unknown)
        rho_n   = (tol/eta_n)**k_i * (eta_nm1/eta_n)**k_p * (eta_nm1**2/(eta_n*eta_nm2))**k_d

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
