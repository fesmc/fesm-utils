program test_tsgen
    ! Exercise the tsgen (time-series generator) module:
    !   * time-driven methods (ramp-time, ramp-time-step, sin) via the analytic
    !     tabulation path, checked against known values;
    !   * a feedback method (PI42) stepped online, checked to stay within bounds;
    !   * the kill switch (const forcing that saturates at f_max).
    ! The namelist is written by the test itself, so no external files are needed.

    use precision, only: wp
    use tsgen

    implicit none

    character(len=*), parameter :: nml_file = "test_tsgen_tmp.nml"
    character(len=*), parameter :: ser_file = "test_tsgen_series.dat"
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

    write(*,"(a)") ""
    if (nfail == 0) then
        write(*,"(a)") " All tsgen tests PASSED."
    else
        write(*,"(a,i0,a)") " tsgen tests FAILED: ", nfail, " case(s)."
    end if

    ! Clean up the temporary files
    call remove_file(nml_file)
    call remove_file(ser_file)

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

        call tsgen_init(ts, nml_file, 0.0_wp, label="ser")

        call tsgen_update(ts, 500.0_wp, 0.0_wp)
        call check("series     f(500)  == 5 ", ts%f_now,  5.0_wp)
        call tsgen_update(ts, 1500.0_wp, 0.0_wp)
        call check("series     f(1500) == 20", ts%f_now, 20.0_wp)
        call tsgen_update(ts, 3000.0_wp, 0.0_wp)
        call check("series     clamp(3000) == 30", ts%f_now, 30.0_wp)
    end subroutine test_series_method

    subroutine test_ramp_time()
        ! Linear ramp 0 -> 10 over dt_ramp=1000, then held at 10.
        type(tsgen_class) :: ts
        real(wp) :: time(31), f(31)
        integer  :: i

        do i = 1, 31
            time(i) = real(i-1,wp)*100.0_wp
        end do
        call tsgen_init(ts, nml_file, 0.0_wp, label="ramp")
        call tsgen_tabulate(ts, time, f)

        call check("ramp-time  f(0)   == 0 ", f(1),   0.0_wp)
        call check("ramp-time  f(500) == 5 ", f(6),   5.0_wp)
        call check("ramp-time  f(1000)== 10", f(11), 10.0_wp)
        call check("ramp-time  f(2000)== 10", f(21), 10.0_wp)
    end subroutine test_ramp_time

    subroutine test_ramp_time_step()
        ! Triangle: 0 -> 10 (over dt_ramp) -> f_conv=5 (over dt_conv) -> held.
        type(tsgen_class) :: ts
        real(wp) :: time(31), f(31)
        integer  :: i

        do i = 1, 31
            time(i) = real(i-1,wp)*100.0_wp
        end do
        call tsgen_init(ts, nml_file, 0.0_wp, label="step")
        call tsgen_tabulate(ts, time, f)

        call check("ramp-step  f(1000)== 10 ", f(11), 10.0_wp)
        call check("ramp-step  f(1500)== 7.5", f(16),  7.5_wp)
        call check("ramp-step  f(2000)== 5  ", f(21),  5.0_wp)
        call check("ramp-step  f(3000)== 5  ", f(31),  5.0_wp)
    end subroutine test_ramp_time_step

    subroutine test_sin()
        ! Sinusoid within [0,10], period 2000: mid at t=0, peak at 500, mid at 1000.
        type(tsgen_class) :: ts
        real(wp) :: time(21), f(21)
        integer  :: i

        do i = 1, 21
            time(i) = real(i-1,wp)*100.0_wp
        end do
        call tsgen_init(ts, nml_file, 0.0_wp, label="sin")
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
        real(wp) :: tt, dt, v, dvdt
        integer  :: i
        logical  :: ok

        call tsgen_init(ts, nml_file, 0.0_wp, label="pi")
        tt = 0.0_wp
        dt = 50.0_wp
        v  = 100.0_wp
        ok = .TRUE.
        do i = 1, 300
            dvdt = -0.5_wp*exp(-tt/1000.0_wp)
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
        real(wp) :: tt, dt
        integer  :: i

        call tsgen_init(ts, nml_file, 0.0_wp, label="kill")
        tt = 0.0_wp
        dt = 50.0_wp
        do i = 1, 40
            call tsgen_update(ts, tt, 0.0_wp, dv_dt=0.0_wp)
            tt = tt + dt
        end do
        call check_logical("const+kill activates at saturation ", ts%kill)
    end subroutine test_kill

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

        ! Tabulated-series method: only method/sigma/series_file are consulted.
        write(u,"(a)") "&tsgen_ser"
        write(u,"(a)") "    method      = ""series"""
        write(u,"(a)") "    sigma       = 0.0"
        write(u,"(a)") "    series_file = """//ser_file//""""
        write(u,"(a)") "/"
        write(u,"(a)") ""

        close(u)
    end subroutine write_namelist

    subroutine put(u, label, method, kill, dt_ave, dt_ramp, dt_conv, f_conv, tolv, df_dt_max, df_sign)
        integer,          intent(in) :: u
        character(len=*), intent(in) :: label, method, kill
        real(wp),         intent(in) :: dt_ave, dt_ramp, dt_conv, f_conv, tolv, df_dt_max, df_sign
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
        write(u,"(a)")        "/"
        write(u,"(a)")        ""
    end subroutine put

    subroutine remove_file(fname)
        character(len=*), intent(in) :: fname
        integer :: u, ios
        open(newunit=u, file=fname, status="old", iostat=ios)
        if (ios == 0) close(u, status="delete")
    end subroutine remove_file

end program test_tsgen
