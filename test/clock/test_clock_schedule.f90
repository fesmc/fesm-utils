program test_clock_schedule
    ! fesm_clock schedules:
    !   * file with values, t0:dt:t1, comments, duplicates, out-of-range times;
    !     time_end added; off-grid time rejected;
    !   * `times` string in a namelist (values and ranges), const (output
    !     includes time_init, restart excludes it), none;
    !   * due/cursor in a loop, t_next; kill switch (clock_finish): output due,
    !     restart due only where every cadence fired;
    !   * the example of docs/dev/clock-design.md section 6 from a namelist with
    !     a defaults file (calendar and year from the defaults), describe table.
    ! Writes logs/clock/sched* (run from the repo root).

    use precision,  only : dp
    use fesm_clock

    implicit none

    character(len=*), parameter :: ftimes = "logs/clock/sched_times.txt"
    character(len=*), parameter :: fnml   = "logs/clock/sched.nml"
    character(len=*), parameter :: fdef   = "logs/clock/sched_defaults.nml"
    type(clock_t)    :: clk
    type(schedule_t) :: s, s_out, s_rst
    type(cadence_t)  :: c_geo, c_ice
    integer :: nfail, stat, ndue, u, i
    real(dp), allocatable :: t(:)

    nfail = 0
    call execute_command_line("mkdir -p logs/clock")

    open(newunit=u, file=ftimes, status="replace", action="write")
    write(u,"(a)") "# output times"
    write(u,"(a)") "0"
    write(u,"(a)") "50     # inline comment"
    write(u,"(a)") "50"
    write(u,"(a)") ""
    write(u,"(a)") "100:100:500"
    write(u,"(a)") "-1000"
    write(u,"(a)") "5000"
    close(u)

    open(newunit=u, file=fnml, status="replace", action="write")
    write(u,"(a)") "&tm_times"
    write(u,"(a)") '    method = "times"'
    write(u,"(a)") "    times  = 0 500 600:100:800"
    write(u,"(a)") "/"
    write(u,"(a)") "&tm_file"
    write(u,"(a)") '    method = "file"'
    write(u,"(a)") '    file   = "'//ftimes//'"'
    write(u,"(a)") "/"
    write(u,"(a)") "&tm_const"
    write(u,"(a)") '    method = "const"'
    write(u,"(a)") "    dt     = 250"
    write(u,"(a)") "/"
    write(u,"(a)") "&tm_rst"
    write(u,"(a)") '    method = "const"'
    write(u,"(a)") "    dt     = 500"
    write(u,"(a)") "/"
    write(u,"(a)") "&tm_none"
    write(u,"(a)") '    method = "none"'
    write(u,"(a)") "/"
    write(u,"(a)") "&tm_bad"
    write(u,"(a)") '    method = "times"'
    write(u,"(a)") "    times  = 0, 105"
    write(u,"(a)") "/"
    ! Section 6 of the design note
    write(u,"(a)") "&ctrl"
    write(u,"(a)") '    method    = "rel"'
    write(u,"(a)") "    time_init = -21000"
    write(u,"(a)") "    time_end  = 0"
    write(u,"(a)") "    dtt       = 5"
    write(u,"(a)") "/"
    write(u,"(a)") "&tm_2D"
    write(u,"(a)") '    method = "const"'
    write(u,"(a)") "    dt     = 1000"
    write(u,"(a)") "/"
    write(u,"(a)") "&tm_restart"
    write(u,"(a)") '    method = "const"'
    write(u,"(a)") "    dt     = 5000"
    write(u,"(a)") "/"
    close(u)

    open(newunit=u, file=fdef, status="replace", action="write")
    write(u,"(a)") "&ctrl"
    write(u,"(a)") '    calendar = "365_day"'
    write(u,"(a)") '    year     = "tropical"'
    write(u,"(a)") "    time_ref = 1950"
    write(u,"(a)") "/"
    close(u)

    call clock_init(clk, 0.0_dp, 1000.0_dp, 10.0_dp, "365_day", "tropical")

    ! ----- file -----------------------------------------------------------------------
    call schedule_init(s, clk, "out_file", fnml, "tm_file")
    t = schedule_times(s, clk)
    call check("file: values, range, comments, duplicates, out-of-range, time_end", &
        same_list(t, [0.0_dp, 50.0_dp, 100.0_dp, 200.0_dp, 300.0_dp, 400.0_dp, 500.0_dp, 1000.0_dp]))

    ! ----- times string ------------------------------------------------------------------
    call schedule_init(s, clk, "out_times", fnml, "tm_times")
    t = schedule_times(s, clk)
    call check("times string: 0 500 600:100:800 + time_end", &
        same_list(t, [0.0_dp, 500.0_dp, 600.0_dp, 700.0_dp, 800.0_dp, 1000.0_dp]))

    ! ----- const output / restart, none ----------------------------------------------------
    call schedule_init(s, clk, "out_const", fnml, "tm_const")
    t = schedule_times(s, clk)
    call check("const output: multiples of 250 from 0", &
        same_list(t, [0.0_dp, 250.0_dp, 500.0_dp, 750.0_dp, 1000.0_dp]))
    call schedule_init(s, clk, "rst_const", fnml, "tm_rst", restart=.true.)
    t = schedule_times(s, clk)
    call check("const restart: excludes time_init", same_list(t, [500.0_dp, 1000.0_dp]))
    call schedule_init(s, clk, "none", fnml, "tm_none")
    call check("none: empty", schedule_size(s) == 0)

    call clock_init(clk, 100.0_dp, 1000.0_dp, 10.0_dp, "365_day", "tropical")
    call schedule_init(s, clk, "out_const", fnml, "tm_const")
    t = schedule_times(s, clk)
    call check("const output from time_init = 100: absolute multiples, time_init included", &
        same_list(t, [100.0_dp, 250.0_dp, 500.0_dp, 750.0_dp, 1000.0_dp]))
    call schedule_init(s, clk, "out_bad", fnml, "tm_bad", stat=stat)
    call check("times 0, 105 (off the step grid) rejected", stat == 1)

    ! ----- due loop, cursor, t_next ----------------------------------------------------------
    call clock_init(clk, 0.0_dp, 1000.0_dp, 10.0_dp, "365_day", "tropical")
    call schedule_init(s, clk, "out_file", fnml, "tm_file")
    call clock_start(clk)
    call check("t_next at 0 is 50", schedule_t_next(s, clk) == 50.0_dp)
    ndue = 0
    do
        if (schedule_due(s, clk)) ndue = ndue + 1
        if (schedule_due(s, clk) .and. clock_time(clk) == 300.0_dp) &
            call check("t_next at 300 is 400", schedule_t_next(s, clk) == 400.0_dp)
        if (clock_finished(clk)) exit
        call clock_advance(clk)
    end do
    call check("due exactly at the 8 listed times", ndue == 8)
    call check("t_next at the end: none left (huge)", schedule_t_next(s, clk) == huge(1.0_dp))

    ! ----- kill switch ------------------------------------------------------------------------
    call clock_init(clk, 0.0_dp, 1000.0_dp, 10.0_dp, "365_day", "tropical")
    call cadence_init(c_geo, clk, "geo", dt=100.0_dp)
    call schedule_init(s_out, clk, "out", times=[500.0_dp])
    call schedule_init(s_rst, clk, "restart", dt=500.0_dp, restart=.true.)
    call clock_start(clk)
    do i = 1, 15
        call clock_advance(clk)
    end do
    call clock_finish(clk)
    call check("kill at 150: finished", clock_finished(clk))
    call check("kill at 150: output due", schedule_due(s_out, clk))
    call check("kill at 150: restart not allowed (geo grid)", .not. clock_restart_allowed(clk))
    call check("kill at 150: restart not due", .not. schedule_due(s_rst, clk))

    call clock_init(clk, 0.0_dp, 1000.0_dp, 10.0_dp, "365_day", "tropical")
    call cadence_init(c_geo, clk, "geo", dt=100.0_dp)
    call schedule_init(s_rst, clk, "restart", dt=500.0_dp, restart=.true.)
    call clock_start(clk)
    do i = 1, 20
        call clock_advance(clk)
    end do
    call clock_finish(clk)
    call check("kill at 200: restart due (on every cadence grid)", schedule_due(s_rst, clk))

    ! ----- section 6 of the design note (describe in the log) ------------------------------
    call clock_init(clk, fnml, "ctrl", defaults_file=fdef)
    call cadence_init(c_ice, clk, "ice", every_step=.true.)
    call cadence_init(c_geo, clk, "geo", dt=10.0_dp)
    call schedule_init(s_out, clk, "out_2D", fnml, "tm_2D")
    call schedule_init(s_rst, clk, "restart", fnml, "tm_restart", restart=.true.)
    call clock_start(clk, stat=stat)
    call check("section 6 example: starts", stat == 0)
    call check("section 6 example: calendar and year from the defaults file", &
        trim(clk%calendar) == "365_day" .and. trim(clk%year) == "tropical")
    call check("section 6 example: 22 outputs, 5 restarts", schedule_size(s_out) == 22 .and. schedule_size(s_rst) == 5)

    call summary("test_clock_schedule")

contains

    logical function same_list(a, b)
        real(dp), intent(in) :: a(:), b(:)
        same_list = size(a) == size(b)
        if (same_list) same_list = all(a == b)
        if (.not. same_list) write(*,"(a,*(1x,g0))") "    got:", a
    end function same_list

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

end program test_clock_schedule
