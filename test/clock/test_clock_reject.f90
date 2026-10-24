program test_clock_reject
    ! fesm_clock: every invalid setup is rejected at init/start (via the
    ! optional stat; the messages with the value and the nearest valid values
    ! are on stderr, i.e. in the test log), and the valid neighbours pass.
    ! Writes logs/clock/reject.nml (run from the repo root).

    use precision,  only : dp
    use fesm_clock

    implicit none

    character(len=*), parameter :: fnml = "logs/clock/reject.nml"
    type(clock_t)    :: clk
    type(cadence_t)  :: c1, c2
    type(schedule_t) :: s
    integer :: nfail, stat, u

    nfail = 0
    call execute_command_line("mkdir -p logs/clock")

    ! ----- dt = 0.3333: not a whole number of ticks --------------------------
    call clock_init(clk, 0.0_dp, 100.0_dp, 1.0_dp, "365_day", "tropical")
    call cadence_init(c1, clk, "geo", dt=0.3333_dp, stat=stat)
    call check("cadence dt = 0.3333 rejected", stat == 1)
    call cadence_init(c1, clk, "geo", per_year=3, stat=stat)
    call check("... per_year = 3 is not a multiple of dtt = 1: rejected", stat == 1)

    ! ----- dt_isos = 15 with dtt = 10 ----------------------------------------
    call clock_init(clk, 0.0_dp, 1000.0_dp, 10.0_dp, "365_day", "tropical")
    call cadence_init(c1, clk, "isos", dt=15.0_dp, stat=stat)
    call check("cadence dt = 15 with dtt = 10 rejected", stat == 1)
    call cadence_init(c1, clk, "isos", dt=20.0_dp, stat=stat)
    call check("... dt = 20 accepted", stat == 0)
    call cadence_init(c2, clk, "geo", dt=10.0_dp, phase=5.0_dp, stat=stat)
    call check("cadence phase 5 off the step grid (dtt = 10) rejected", stat == 1)

    ! ----- schedule time 105 off the step grid -------------------------------
    call schedule_init(s, clk, "out", times=[100.0_dp, 105.0_dp], stat=stat)
    call check("schedule time 105 with dtt = 10 rejected", stat == 1)
    call schedule_init(s, clk, "out", times=[100.0_dp, 110.0_dp, 5000.0_dp], stat=stat)
    call check("... 100, 110 (and out-of-range 5000) accepted", stat == 0)
    call schedule_init(s, clk, "out", dt=15.0_dp, stat=stat)
    call check("schedule const dt = 15 with dtt = 10 rejected", stat == 1)

    ! ----- time_end = 105 with geo = 10 --------------------------------------
    call clock_init(clk, 0.0_dp, 105.0_dp, 5.0_dp, "365_day", "tropical")
    call cadence_init(c1, clk, "geo", dt=10.0_dp)
    call schedule_init(s, clk, "restart", times=[50.0_dp], restart=.true.)
    call clock_start(clk, stat=stat)
    call check("time_end = 105 off the geo grid (10) rejected at start", stat == 1)
    call clock_init(clk, 0.0_dp, 100.0_dp, 5.0_dp, "365_day", "tropical")
    call cadence_init(c1, clk, "geo", dt=10.0_dp)
    call schedule_init(s, clk, "restart", times=[55.0_dp], restart=.true.)
    call clock_start(clk, stat=stat)
    call check("restart time 55 off the geo grid rejected at start", stat == 1)

    ! ----- disjoint phases ---------------------------------------------------
    call clock_init(clk, 0.0_dp, 100.0_dp, 5.0_dp, "365_day", "tropical")
    call cadence_init(c1, clk, "a", dt=10.0_dp, phase=0.0_dp)
    call cadence_init(c2, clk, "b", dt=10.0_dp, phase=5.0_dp)
    call schedule_init(s, clk, "restart", dt=50.0_dp, restart=.true.)
    call clock_start(clk, stat=stat)
    call check("cadences 10 yr with phases 0 and 5 (no common time) rejected", stat == 1)
    call clock_init(clk, 0.0_dp, 120.0_dp, 5.0_dp, "365_day", "tropical")
    call cadence_init(c1, clk, "a", dt=10.0_dp, phase=0.0_dp)
    call cadence_init(c2, clk, "b", dt=15.0_dp, phase=5.0_dp)
    call schedule_init(s, clk, "restart", dt=60.0_dp, restart=.true.)
    call clock_start(clk, stat=stat)
    call check("cadences 10/0 and 15/5 (common 30k + 20) with restarts at 60, 120 rejected", stat == 1)
    call clock_init(clk, 0.0_dp, 110.0_dp, 5.0_dp, "365_day", "tropical")
    call cadence_init(c1, clk, "a", dt=10.0_dp, phase=0.0_dp)
    call cadence_init(c2, clk, "b", dt=15.0_dp, phase=5.0_dp)
    call schedule_init(s, clk, "restart", times=[50.0_dp], restart=.true.)
    call clock_start(clk, stat=stat)
    call check("... restarts at 50 and time_end 110 (on 30k + 20) accepted", stat == 0)
    call cadence_init(c1, clk, "late", dt=10.0_dp, stat=stat)
    call check("cadence registered after clock_start rejected", stat == 1)

    ! ----- calendar, year, step grid, const ------------------------------------
    call clock_init(clk, 0.0_dp, 100.0_dp, 1.0_dp, "gregorian", "tropical", stat=stat)
    call check("calendar gregorian rejected", stat == 1)
    call clock_init(clk, 0.0_dp, 100.0_dp, 1.0_dp, "365_day", "sidereal", stat=stat)
    call check("year sidereal rejected", stat == 1)
    call clock_init(clk, 0.0_dp, 105.0_dp, 10.0_dp, "365_day", "tropical", stat=stat)
    call check("time_end = 105 off the step grid (dtt = 10) rejected", stat == 1)
    call clock_init(clk, 0.0_dp, 100.0_dp, 0.3333_dp, "365_day", "tropical", stat=stat)
    call check("dtt = 0.3333 rejected", stat == 1)
    call clock_init(clk, 0.0_dp, 100.0_dp, 1.0_dp, "365_day", "tropical", method="const", stat=stat)
    call check("method const without time_const rejected", stat == 1)
    call clock_init(clk, 0.0_dp, 100.0_dp, 1.0_dp, "365_day", "tropical", method="const", &
                    time_const=0.0_dp, time_const_cal=1950.0_dp, stat=stat)
    call check("method const with both time_const and time_const_cal rejected", stat == 1)
    call clock_init(clk, 0.0_dp, 100.0_dp, 1.0_dp, "365_day", "tropical", method="relative", stat=stat)
    call check("method relative rejected", stat == 1)
    call clock_init(clk, 0.0_dp, 100.0_dp, 1.0_dp, "365_day", "tropical", time_ref=1950.5_dp, stat=stat)
    call check("time_ref = 1950.5 (not whole years) rejected", stat == 1)
    call clock_init(clk, 100.0_dp, 0.0_dp, 1.0_dp, "365_day", "tropical", stat=stat)
    call check("time_end before time_init rejected", stat == 1)
    call clock_start(clk, stat=stat)
    call check("clock_start on a rejected clock rejected", stat == 1)

    ! ----- namelist: no defaults for calendar and year --------------------------
    open(newunit=u, file=fnml, status="replace", action="write")
    write(u,"(a)") "&ctrl"
    write(u,"(a)") "    time_init = 0"
    write(u,"(a)") "    time_end  = 100"
    write(u,"(a)") "    dtt       = 1"
    write(u,"(a)") '    year      = "tropical"'
    write(u,"(a)") "/"
    write(u,"(a)") "&tm_bad"
    write(u,"(a)") '    method = "every"'
    write(u,"(a)") "/"
    write(u,"(a)") "&tm_range"
    write(u,"(a)") '    method = "times"'
    write(u,"(a)") '    times  = "0:15:90"'
    write(u,"(a)") "/"
    close(u)
    call clock_init(clk, fnml, "ctrl", stat=stat)
    call check("namelist without calendar rejected", stat == 1)
    call clock_init(clk, 0.0_dp, 100.0_dp, 10.0_dp, "365_day", "tropical")
    call schedule_init(s, clk, "out", fnml, "tm_bad", stat=stat)
    call check("schedule method every rejected", stat == 1)
    call schedule_init(s, clk, "out", fnml, "tm_range", stat=stat)
    call check("schedule range 0:15:90 with dtt = 10 rejected", stat == 1)

    call summary("test_clock_reject")

contains

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

end program test_clock_reject
