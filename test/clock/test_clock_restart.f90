program test_clock_restart
    ! fesm_clock restart bundles (clock + cadence state, ncio dp scalars):
    !   * written at 100 by a run that started at 0;
    !   * continuation (time_init = 100): continued, n_origin and n_last restored;
    !   * branch (time_init = 50): n_origin and n_last reset to time_init;
    !   * a bundle written at N = 1800 read at N = 36000 (converted exactly);
    !   * a bundle at N = 7 holding 1/7 yr: not convertible to N = 36000, rejected.
    ! The layout of a bundle (clock and cadence variables and units, the
    ! attribute clock_method, the cadence appended to the clock's file) is
    ! checked with ncio.
    ! Writes logs/clock/restart_*.nc (run from the repo root).

    use precision,  only : dp
    use ncio,       only : nc_exists_var, nc_read_attr, nc_read
    use fesm_clock

    implicit none

    character(len=*), parameter :: fa = "logs/clock/restart_a.nc"
    character(len=*), parameter :: fb = "logs/clock/restart_b1800.nc"
    character(len=*), parameter :: fc = "logs/clock/restart_c7.nc"
    type(clock_t)   :: clk
    type(cadence_t) :: c_geo
    integer :: nfail, stat
    character(len=32) :: str
    real(dp) :: x

    nfail = 0
    call execute_command_line("mkdir -p logs/clock && rm -f "//fa//" "//fb//" "//fc)

    ! ----- write at 100 ----------------------------------------------------------
    call run_and_write(0.0_dp, 100.0_dp, 1.0_dp, 36000, fa)

    ! ----- layout of the bundle ----------------------------------------------------
    call check("bundle holds the clock and the cadence variables", &
               nc_exists_var(fa, "clock_ticks_per_year") .and. nc_exists_var(fa, "clock_n") .and. &
               nc_exists_var(fa, "clock_n_origin") .and. nc_exists_var(fa, "clock_dtt") .and. &
               nc_exists_var(fa, "cadence_geo_n_last") .and. nc_exists_var(fa, "cadence_geo_dt") .and. &
               nc_exists_var(fa, "cadence_geo_phase"))
    str = ""
    call nc_read_attr(fa, "clock_n", "units", str)
    call check("clock_n: units = ticks", trim(str) == "ticks")
    str = ""
    call nc_read_attr(fa, "clock_ticks_per_year", "long_name", str)
    call check("clock_ticks_per_year: long_name", trim(str) == "ticks per year")
    str = ""
    call nc_read_attr(fa, "clock_method", str)
    call check("global clock_method = rel", trim(str) == "rel")
    call nc_read(fa, "clock_ticks_per_year", x)
    call check("clock_ticks_per_year = 36000", x == 36000.0_dp)
    call nc_read(fa, "clock_n", x)
    call check("clock_n = 100 yr in ticks", x == 100.0_dp*36000.0_dp)
    call nc_read(fa, "cadence_geo_dt", x)
    call check("cadence_geo_dt = 10 yr in ticks", x == 10.0_dp*36000.0_dp)

    ! ----- continuation ------------------------------------------------------------
    call clock_init(clk, 100.0_dp, 200.0_dp, 0.5_dp, "365_day", "tropical")
    call cadence_init(c_geo, clk, "geo", dt=10.0_dp)
    call clock_restart_read(clk, fa)
    call cadence_restart_read(c_geo, clk, fa)
    call clock_start(clk)
    call check("continuation: continued = .true.", clk%continued)
    call check("continuation: elapsed origin 0 restored (time_elapsed = 100)", clock_time_elapsed(clk) == 100.0_dp)
    call check("continuation: n_last restored (t_last = 100)", cadence_t_last(c_geo, clk) == 100.0_dp)
    call check("continuation: geo not due before the first advance", .not. cadence_due(c_geo, clk))
    call clock_advance(clk)
    call check("continuation: time_elapsed runs on (100.5)", clock_time_elapsed(clk) == 100.5_dp)

    ! ----- branch ------------------------------------------------------------------
    call clock_init(clk, 50.0_dp, 150.0_dp, 1.0_dp, "365_day", "tropical")
    call cadence_init(c_geo, clk, "geo", dt=10.0_dp)
    call clock_restart_read(clk, fa)
    call cadence_restart_read(c_geo, clk, fa)
    call clock_start(clk)
    call check("branch: continued = .false.", .not. clk%continued)
    call check("branch: elapsed origin reset (time_elapsed = 0)", clock_time_elapsed(clk) == 0.0_dp)
    call check("branch: n_last = time_init (t_last = 50)", cadence_t_last(c_geo, clk) == 50.0_dp)

    ! ----- N = 1800 -> 36000 ----------------------------------------------------------
    call run_and_write(0.0_dp, 100.0_dp, 0.5_dp, 1800, fb)
    call clock_init(clk, 100.0_dp, 200.0_dp, 1.0_dp, "365_day", "tropical")
    call cadence_init(c_geo, clk, "geo", dt=10.0_dp)
    call clock_restart_read(clk, fb, stat=stat)
    call check("N = 1800 bundle read at N = 36000", stat == 0)
    call cadence_restart_read(c_geo, clk, fb, stat=stat)
    call check("N = 1800 cadence state read at N = 36000", stat == 0)
    call check("... continuation, time_elapsed = 100", clk%continued .and. clock_time_elapsed(clk) == 100.0_dp)
    call check("... n_last = 100", cadence_t_last(c_geo, clk) == 100.0_dp)

    ! ----- N = 7 bundle at 1/7 yr: not convertible -----------------------------------
    call clock_init(clk, 0.0_dp, 1.0_dp, 1.0_dp/7.0_dp, "none", "365d", ticks_per_year=7)
    call clock_start(clk)
    call clock_advance(clk)
    call clock_restart_write(clk, fc)
    call clock_init(clk, 0.0_dp, 10.0_dp, 1.0_dp, "365_day", "tropical")
    call clock_restart_read(clk, fc, stat=stat)
    call check("N = 7 bundle at 1/7 yr rejected at N = 36000", stat == 1)

    call clock_init(clk, 0.0_dp, 10.0_dp, 1.0_dp, "365_day", "tropical")
    call clock_restart_read(clk, "logs/clock/does_not_exist.nc", stat=stat)
    call check("missing bundle rejected", stat == 1)

    call summary("test_clock_restart")

contains

    subroutine run_and_write(ti, te, dtt, tpy, filename)
        ! Run ti -> te with a 10-yr cadence; write the bundle at te.
        real(dp),         intent(in) :: ti, te, dtt
        integer,          intent(in) :: tpy
        character(len=*), intent(in) :: filename
        type(clock_t)    :: clk
        type(cadence_t)  :: c_geo
        type(schedule_t) :: s_rst

        call clock_init(clk, ti, te, dtt, "365_day", "tropical", ticks_per_year=tpy)
        call cadence_init(c_geo, clk, "geo", dt=10.0_dp)
        call schedule_init(s_rst, clk, "restart", times=[te], restart=.true.)
        call clock_start(clk)
        do
            if (schedule_due(s_rst, clk)) then
                call clock_restart_write(clk, filename)
                call cadence_restart_write(c_geo, filename)
            end if
            if (clock_finished(clk)) exit
            call clock_advance(clk)
            if (cadence_due(c_geo, clk)) call cadence_mark(c_geo, clk)
        end do
    end subroutine run_and_write

    subroutine check(name, ok)
        character(len=*), intent(in) :: name
        logical,          intent(in) :: ok
        if (ok) then
            write(*,"(a)") "  PASS  "//name
        else
            write(*,"(a)") "  FAIL  "//name
            nfail = nfail + 1
        end if
    end subroutine check

    subroutine summary(prog)
        character(len=*), intent(in) :: prog
        write(*,"(a)") ""
        if (nfail == 0) then
            write(*,"(a)") " All "//prog//" tests PASSED."
        else
            write(*,"(a,i0,a)") " "//prog//" tests FAILED: ", nfail, " case(s)."
            error stop 1
        end if
    end subroutine summary

end program test_clock_restart
