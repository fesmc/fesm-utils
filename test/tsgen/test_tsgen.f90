program test_tsgen
    ! Exercise the tsgen (time-series generator) module:
    !   * time-driven methods (ramp-time, ramp-time-step, sin) via the analytic
    !     tabulation path, checked against known values;
    !   * a feedback method (PI42) stepped online, checked to stay within bounds;
    !   * the kill switch (const forcing that saturates at f_max).
    !   * defaults: keys missing from the par group fall back to a template group
    !     in a defaults file (defaults_file/defaults_group).
    !   * the window modifier (win_*), analytic and series methods;
    !   * the RNG state in the restart: noise after a restart as without it.
    ! The namelist is written by the test itself, so no external files are needed.

    use precision, only: dp, wp
    use tsgen

    implicit none

    character(len=*), parameter :: nml_file = "test_tsgen_tmp.nml"
    character(len=*), parameter :: ser_file = "test_tsgen_series.dat"
    character(len=*), parameter :: def_file = "test_tsgen_defaults_tmp.nml"
    character(len=*), parameter :: win_file = "test_tsgen_win_series.dat"
    character(len=*), parameter :: rst_file = "test_tsgen_restart_tmp.nc"
    real(wp),         parameter :: tol = real(1e-4,wp)
    integer :: nfail

    nfail = 0

    write(*,*) "==================================================="
    write(*,*) " tsgen test suite"
    write(*,*) "==================================================="
    write(*,"(a)") ""

    call write_namelist()

    call test_ramp_time()
    call test_ramp_time_step()
    call test_sin()
    call test_series_method()
    call test_feedback_bounds()
    call test_kill()
    call test_defaults()
    call test_window()
    call test_rng_restart()

    write(*,"(a)") ""
    if (nfail == 0) then
        write(*,"(a)") " All tsgen tests PASSED."
    else
        write(*,"(a,i0,a)") " tsgen tests FAILED: ", nfail, " case(s)."
    end if

    ! Clean up the temporary files
    call remove_file(nml_file)
    call remove_file(ser_file)
    call remove_file(def_file)
    call remove_file(win_file)
    call remove_file(rst_file)

    if (nfail > 0) stop 1

contains

    subroutine test_series_method()
        ! Tabulated ascii series [time value]: forcing must equal the clamped
        ! linear interpolation of the file at the update time (no noise).
        type(tsgen_class) :: ts
        integer :: u

        open(newunit=u, file=ser_file, status="replace", action="write")
        write(u,"(a)") "   0.0    0.0"
        write(u,"(a)") "1000.0   10.0"
        write(u,"(a)") "2000.0   30.0"
        close(u)

        call tsgen_init(ts, nml_file, 0.0_dp, group="tsgen_ser")

        call tsgen_update(ts, 500.0_dp, 0.0_wp)
        call check("series     f(500)  == 5 ", ts%f_now,  5.0_wp)
        call tsgen_update(ts, 1500.0_dp, 0.0_wp)
        call check("series     f(1500) == 20", ts%f_now, 20.0_wp)
        call tsgen_update(ts, 3000.0_dp, 0.0_wp)
        call check("series     clamp(3000) == 30", ts%f_now, 30.0_wp)
    end subroutine test_series_method

    subroutine test_ramp_time()
        ! Linear ramp 0 -> 10 over dt_ramp=1000, then held at 10.
        type(tsgen_class) :: ts
        real(dp) :: time(31)
        real(wp) :: f(31)
        integer  :: i

        do i = 1, 31
            time(i) = real(i-1,dp)*100.0_dp
        end do
        call tsgen_init(ts, nml_file, 0.0_dp, group="tsgen_ramp")
        call tsgen_tabulate(ts, time, f)

        call check("ramp-time  f(0)   == 0 ", f(1),   0.0_wp)
        call check("ramp-time  f(500) == 5 ", f(6),   5.0_wp)
        call check("ramp-time  f(1000)== 10", f(11), 10.0_wp)
        call check("ramp-time  f(2000)== 10", f(21), 10.0_wp)
    end subroutine test_ramp_time

    subroutine test_ramp_time_step()
        ! Triangle: 0 -> 10 (over dt_ramp) -> f_conv=5 (over dt_conv) -> held.
        type(tsgen_class) :: ts
        real(dp) :: time(31)
        real(wp) :: f(31)
        integer  :: i

        do i = 1, 31
            time(i) = real(i-1,dp)*100.0_dp
        end do
        call tsgen_init(ts, nml_file, 0.0_dp, group="tsgen_step")
        call tsgen_tabulate(ts, time, f)

        call check("ramp-step  f(1000)== 10 ", f(11), 10.0_wp)
        call check("ramp-step  f(1500)== 7.5", f(16),  7.5_wp)
        call check("ramp-step  f(2000)== 5  ", f(21),  5.0_wp)
        call check("ramp-step  f(3000)== 5  ", f(31),  5.0_wp)
    end subroutine test_ramp_time_step

    subroutine test_sin()
        ! Sinusoid within [0,10], period 2000: mid at t=0, peak at 500, mid at 1000.
        type(tsgen_class) :: ts
        real(dp) :: time(21)
        real(wp) :: f(21)
        integer  :: i

        do i = 1, 21
            time(i) = real(i-1,dp)*100.0_dp
        end do
        call tsgen_init(ts, nml_file, 0.0_dp, group="tsgen_sin")
        call tsgen_tabulate(ts, time, f)

        call check("sin        f(0)   == 5 ", f(1),   5.0_wp)
        call check("sin        f(500) == 10", f(6),  10.0_wp)
        call check("sin        f(1000)== 5 ", f(11),  5.0_wp)
        call check("sin        f(1500)== 0 ", f(16),  0.0_wp)
    end subroutine test_sin

    subroutine test_feedback_bounds()
        ! PI42 controller stepped online with a decaying response derivative.
        ! Forcing must stay within [f_min,f_max] and remain finite throughout.
        type(tsgen_class) :: ts
        real(dp) :: tt, dt
        real(wp) :: v, dvdt
        integer  :: i
        logical  :: ok

        call tsgen_init(ts, nml_file, 0.0_dp, group="tsgen_pi")
        tt = 0.0_dp
        dt = 50.0_dp
        v  = 100.0_wp
        ok = .TRUE.
        do i = 1, 300
            dvdt = -0.5_wp*exp(-tt/1000.0_dp)
            v = v + dvdt*dt
            call tsgen_update(ts, tt, v, dv_dt=dvdt)
            if (ts%f_now .lt. ts%par%f_min - tol) ok = .FALSE.
            if (ts%f_now .gt. ts%par%f_max + tol) ok = .FALSE.
            if (ts%f_now /= ts%f_now)             ok = .FALSE.   ! NaN guard
            tt = tt + dt
        end do
        call check_logical("PI42       stays within [f_min,f_max]", ok)
    end subroutine test_feedback_bounds

    subroutine test_kill()
        ! Constant forcing that saturates at f_max with an equilibrated response
        ! (dv_dt=0 < eps): the kill switch must activate once f reaches f_max.
        type(tsgen_class) :: ts
        real(dp) :: tt, dt
        integer  :: i

        call tsgen_init(ts, nml_file, 0.0_dp, group="tsgen_kill")
        tt = 0.0_dp
        dt = 50.0_dp
        do i = 1, 40
            call tsgen_update(ts, tt, 0.0_wp, dv_dt=0.0_wp)
            tt = tt + dt
        end do
        call check_logical("const+kill activates at saturation ", ts%kill)
    end subroutine test_kill

    subroutine test_defaults()
        ! Template group "tsgen_tmpl" in def_file holds every key (ramp 0 -> 10
        ! over dt_ramp=1000). Par group "tsgen_user" sets only method and f_max=20,
        ! so dt_ramp etc. come from the template and f_max is overlaid: f(500) == 10.
        type(tsgen_class) :: ts
        real(dp) :: time(21)
        real(wp) :: f(21)
        integer  :: i, u

        open(newunit=u, file=def_file, status="replace", action="write")
        call put(u, "tmpl", "const", ".FALSE.", 100.0_wp, 1000.0_wp, 1000.0_wp, 5.0_wp, 1.0e-4_wp, 0.01_wp, 1.0_wp, &
                 win=.TRUE.)
        close(u)

        do i = 1, 21
            time(i) = real(i-1,dp)*100.0_dp
        end do
        call tsgen_init(ts, nml_file, 0.0_dp, group="tsgen_user", &
                        defaults_file=def_file, defaults_group="tsgen_tmpl")
        call tsgen_tabulate(ts, time, f)

        call check("defaults   f_max overlaid == 20", ts%par%f_max, 20.0_wp)
        call check("defaults   f(500)  == 10       ", f(6),  10.0_wp)
        call check("defaults   f(1000) == 20       ", f(11), 20.0_wp)
    end subroutine test_defaults

    subroutine test_window()
        ! Constant value 2 with win_ref = 1, win_fac = 0.5, win_above in the
        ! window -12 ka < t < -1 ka: 1 + 0.5*(2 - 1) = 1.5 inside, 2 outside.
        ! A taper of 1 kyr gives half the effect 500 yr inside the edges (1.75).
        ! A value below win_ref is untouched with win_above (scaled without).
        type(tsgen_class) :: ts
        real(dp) :: time(5)
        real(wp) :: f(5)
        integer  :: u

        time = [-15e3_dp, -12e3_dp, -6e3_dp, -1e3_dp, 0.0_dp]

        call tsgen_init(ts, nml_file, -20e3_dp, group="tsgen_win")
        call check("window     f_now(init) outside == 2", ts%f_now, 2.0_wp)
        call tsgen_tabulate(ts, time, f)
        call check("window     f(-15ka) == 2  ", f(1), 2.0_wp)
        call check("window     f(-12ka) == 2  ", f(2), 2.0_wp)
        call check("window     f(-6ka)  == 1.5", f(3), 1.5_wp)
        call check("window     f(-1ka)  == 2  ", f(4), 2.0_wp)
        call check("window     f(0)     == 2  ", f(5), 2.0_wp)
        call tsgen_update(ts, -6e3_dp, 0.0_wp)
        call check("window     f_now(-6ka) == 1.5", ts%f_now, 1.5_wp)
        call check("window     f_mean(-6ka) == 2 ", ts%f_mean_now, 2.0_wp)

        time = [-12e3_dp, -11.5e3_dp, -11e3_dp, -1.5e3_dp, -1e3_dp]
        call tsgen_init(ts, nml_file, -20e3_dp, group="tsgen_win_taper")
        call tsgen_tabulate(ts, time, f)
        call check("taper      f(-12ka)   == 2   ", f(1), 2.0_wp)
        call check("taper      f(-11.5ka) == 1.75", f(2), 1.75_wp)
        call check("taper      f(-11ka)   == 1.5 ", f(3), 1.5_wp)
        call check("taper      f(-1.5ka)  == 1.75", f(4), 1.75_wp)
        call check("taper      f(-1ka)    == 2   ", f(5), 2.0_wp)

        call tsgen_init(ts, nml_file, -20e3_dp, group="tsgen_win_below")
        call tsgen_update(ts, -6e3_dp, 0.0_wp)
        call check("win_above  below ref untouched == 0.5", ts%f_now, 0.5_wp)

        call tsgen_init(ts, nml_file, -20e3_dp, group="tsgen_win_both")
        call tsgen_update(ts, -6e3_dp, 0.0_wp)
        call check("no above   below ref scaled == 0.75 ", ts%f_now, 0.75_wp)

        ! Series path (eval_series): constant series of 2
        open(newunit=u, file=win_file, status="replace", action="write")
        write(u,"(a)") "-20000.0   2.0"
        write(u,"(a)") "     0.0   2.0"
        close(u)
        call tsgen_init(ts, nml_file, -20e3_dp, group="tsgen_win_ser")
        call tsgen_update(ts, -15e3_dp, 0.0_wp)
        call check("win series f_now(-15ka) == 2  ", ts%f_now, 2.0_wp)
        call tsgen_update(ts, -6e3_dp, 0.0_wp)
        call check("win series f_now(-6ka)  == 1.5", ts%f_now, 1.5_wp)
    end subroutine test_window

    subroutine test_rng_restart()
        ! Noisy ramp (sigma > 0) run to t1 = 100, restart written, continued to
        ! t2 = 200. A fresh object restarted from the file (after other draws
        ! from the generator) must give the identical noise from t1 to t2.
        type(tsgen_class) :: ts
        real(wp) :: f_a(10), f_b(10), eps_a(10), eps_b(10), junk(7)
        integer  :: i

        call tsgen_init(ts, nml_file, 0.0_dp, group="tsgen_noise")
        do i = 1, 10
            call tsgen_update(ts, real(i,dp)*10.0_dp, 0.0_wp)
        end do
        call tsgen_restart_write(rst_file, ts, 100.0_dp)
        do i = 1, 10
            call tsgen_update(ts, 100.0_dp + real(i,dp)*10.0_dp, 0.0_wp)
            f_a(i)   = ts%f_now
            eps_a(i) = ts%eps
        end do

        call tsgen_init(ts, nml_file, 100.0_dp, group="tsgen_noise")
        call random_number(junk)            ! move the generator on
        call tsgen_restart_read(ts, rst_file)
        do i = 1, 10
            call tsgen_update(ts, 100.0_dp + real(i,dp)*10.0_dp, 0.0_wp)
            f_b(i)   = ts%f_now
            eps_b(i) = ts%eps
        end do

        call check_logical("rng        noise drawn (eps /= 0)        ", any(eps_a .ne. 0.0_wp))
        call check_logical("rng        restart: identical noise t1-t2", all(eps_a .eq. eps_b))
        call check_logical("rng        restart: identical f_now t1-t2", all(f_a .eq. f_b))
    end subroutine test_rng_restart

    ! ---- helpers ----

    subroutine check(name, got, expected)
        character(len=*), intent(in) :: name
        real(wp),         intent(in) :: got, expected
        logical :: ok
        character(len=4) :: tag
        ok = abs(got-expected) .le. tol
        tag = "PASS"
        if (.not. ok) then
            tag = "FAIL"
            nfail = nfail + 1
        end if
        write(*,"(a,a,a,f10.4,a,f10.4)") "  ["//tag//"] ", name, "  got=", got, " exp=", expected
    end subroutine check

    subroutine check_logical(name, ok)
        character(len=*), intent(in) :: name
        logical,          intent(in) :: ok
        character(len=4) :: tag
        tag = "PASS"
        if (.not. ok) then
            tag = "FAIL"
            nfail = nfail + 1
        end if
        write(*,"(a,a)") "  ["//tag//"] ", name
    end subroutine check_logical

    subroutine write_namelist()
        integer :: u
        open(newunit=u, file=nml_file, status="replace", action="write")
        call put(u, "ramp", "ramp-time",      ".FALSE.", 100.0_wp, 1000.0_wp, 1000.0_wp, 5.0_wp, 1.0e-4_wp, 0.01_wp, 1.0_wp)
        call put(u, "step", "ramp-time-step", ".FALSE.", 100.0_wp, 1000.0_wp, 1000.0_wp, 5.0_wp, 1.0e-4_wp, 0.01_wp, 1.0_wp)
        call put(u, "sin",  "sin",            ".FALSE.", 100.0_wp, 2000.0_wp, 1000.0_wp, 5.0_wp, 1.0e-4_wp, 0.01_wp, 1.0_wp)
        call put(u, "pi",   "PI42",           ".TRUE.",  200.0_wp, 1000.0_wp, 1000.0_wp, 5.0_wp, 1.0e-2_wp, 0.05_wp, 1.0_wp)
        call put(u, "kill", "const",          ".TRUE.",  100.0_wp, 1000.0_wp, 1000.0_wp, 5.0_wp, 1.0e-4_wp, 0.02_wp, 1.0_wp)

        ! Partial group for test_defaults: other keys come from the template.
        write(u,"(a)") "&tsgen_user"
        write(u,"(a)") "    method = ""ramp-time"""
        write(u,"(a)") "    f_max  = 20.0"
        write(u,"(a)") "/"
        write(u,"(a)") ""

        ! Window modifier on a constant value (const, f_min = f_max)
        call put_win(u, "win",       2.0_wp, "0.0",    ".TRUE.")
        call put_win(u, "win_taper", 2.0_wp, "1000.0", ".TRUE.")
        call put_win(u, "win_below", 0.5_wp, "0.0",    ".TRUE.")
        call put_win(u, "win_both",  0.5_wp, "0.0",    ".FALSE.")

        write(u,"(a)") "&tsgen_win_ser"
        write(u,"(a)") "    method      = ""series"""
        write(u,"(a)") "    sigma       = 0.0"
        write(u,"(a)") "    series_file = """//win_file//""""
        write(u,"(a)") "    win_t0      = -12e3"
        write(u,"(a)") "    win_t1      = -1e3"
        write(u,"(a)") "    win_fac     = 0.5"
        write(u,"(a)") "    win_ref     = 1.0"
        write(u,"(a)") "    win_above   = .TRUE."
        write(u,"(a)") "/"
        write(u,"(a)") ""

        ! Noisy ramp for the RNG restart test
        write(u,"(a)") "&tsgen_noise"
        write(u,"(a)") "    method    = ""ramp-time"""
        write(u,"(a)") "    with_kill = .FALSE."
        write(u,"(a)") "    dt_ave    = 100.0"
        write(u,"(a)") "    dt_init   = 0.0"
        write(u,"(a)") "    dt_ramp   = 1000.0"
        write(u,"(a)") "    dt_conv   = 1000.0"
        write(u,"(a)") "    df_sign   = 1.0"
        write(u,"(a)") "    tol       = 1e-4"
        write(u,"(a)") "    df_dt_max = 0.01"
        write(u,"(a)") "    sigma     = 0.5"
        write(u,"(a)") "    f_min     = 0.0"
        write(u,"(a)") "    f_max     = 10.0"
        write(u,"(a)") "    f_conv    = 5.0"
        write(u,"(a)") "/"
        write(u,"(a)") ""

        ! Tabulated-series method: only method/sigma/series_file are consulted.
        write(u,"(a)") "&tsgen_ser"
        write(u,"(a)") "    method      = ""series"""
        write(u,"(a)") "    sigma       = 0.0"
        write(u,"(a)") "    series_file = """//ser_file//""""
        write(u,"(a)") "/"
        write(u,"(a)") ""

        close(u)
    end subroutine write_namelist

    subroutine put(u, label, method, kill, dt_ave, dt_ramp, dt_conv, f_conv, tolv, df_dt_max, df_sign, win)
        ! win: also write the window keys at their defaults (a template group
        ! must declare every key tsgen_init reads)
        integer,          intent(in) :: u
        character(len=*), intent(in) :: label, method, kill
        real(wp),         intent(in) :: dt_ave, dt_ramp, dt_conv, f_conv, tolv, df_dt_max, df_sign
        logical, optional, intent(in) :: win
        write(u,"(a)")        "&tsgen_"//trim(label)
        write(u,"(a)")        "    method    = """//trim(method)//""""
        write(u,"(a)")        "    with_kill = "//trim(kill)
        write(u,"(a,g0)")     "    dt_ave    = ", dt_ave
        write(u,"(a)")        "    dt_init   = 0.0"
        write(u,"(a,g0)")     "    dt_ramp   = ", dt_ramp
        write(u,"(a,g0)")     "    dt_conv   = ", dt_conv
        write(u,"(a,g0)")     "    df_sign   = ", df_sign
        write(u,"(a,g0)")     "    tol       = ", tolv
        write(u,"(a,g0)")     "    df_dt_max = ", df_dt_max
        write(u,"(a)")        "    sigma     = 0.0"
        write(u,"(a)")        "    f_min     = 0.0"
        write(u,"(a)")        "    f_max     = 10.0"
        write(u,"(a,g0)")     "    f_conv    = ", f_conv
        if (present(win)) then
            if (win) then
                write(u,"(a)") "    win_t0    = 0.0"
                write(u,"(a)") "    win_t1    = 0.0"
                write(u,"(a)") "    win_taper = 0.0"
                write(u,"(a)") "    win_fac   = 1.0"
                write(u,"(a)") "    win_ref   = 0.0"
                write(u,"(a)") "    win_above = .FALSE."
            end if
        end if
        write(u,"(a)")        "/"
        write(u,"(a)")        ""
    end subroutine put

    subroutine put_win(u, label, fval, taper, above)
        ! const method held at fval (f_min = f_max = fval, df_dt_max = 0), with
        ! the window -12 ka < t < -1 ka, win_ref = 1, win_fac = 0.5.
        integer,          intent(in) :: u
        character(len=*), intent(in) :: label, taper, above
        real(wp),         intent(in) :: fval
        write(u,"(a)")        "&tsgen_"//trim(label)
        write(u,"(a)")        "    method    = ""const"""
        write(u,"(a)")        "    with_kill = .FALSE."
        write(u,"(a)")        "    dt_ave    = 100.0"
        write(u,"(a)")        "    dt_init   = 0.0"
        write(u,"(a)")        "    dt_ramp   = 1000.0"
        write(u,"(a)")        "    dt_conv   = 1000.0"
        write(u,"(a)")        "    df_sign   = 1.0"
        write(u,"(a)")        "    tol       = 1e-4"
        write(u,"(a)")        "    df_dt_max = 0.0"
        write(u,"(a)")        "    sigma     = 0.0"
        write(u,"(a,g0)")     "    f_min     = ", fval
        write(u,"(a,g0)")     "    f_max     = ", fval
        write(u,"(a,g0)")     "    f_conv    = ", fval
        write(u,"(a)")        "    win_t0    = -12e3"
        write(u,"(a)")        "    win_t1    = -1e3"
        write(u,"(a)")        "    win_taper = "//trim(taper)
        write(u,"(a)")        "    win_fac   = 0.5"
        write(u,"(a)")        "    win_ref   = 1.0"
        write(u,"(a)")        "    win_above = "//trim(above)
        write(u,"(a)")        "/"
        write(u,"(a)")        ""
    end subroutine put_win

    subroutine remove_file(fname)
        character(len=*), intent(in) :: fname
        integer :: u, ios
        open(newunit=u, file=fname, status="old", iostat=ios)
        if (ios == 0) close(u, status="delete")
    end subroutine remove_file

end program test_tsgen
