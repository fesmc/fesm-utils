program test_clock_exact
    ! fesm_clock: model time is n*tick, never a sum.
    !   * 1e6 steps of dtt = 0.1 from -130000: time() equals to_years(n_init + k*dt)
    !     bit for bit at every step and ends exactly on time_end;
    !   * dtt = 1/12 (from 1.0/12.0) with a per_year = 12 cadence: due every step,
    !     intervals exact, ends exactly on time_end.
    ! Self-contained: no files.

    use, intrinsic :: iso_fortran_env, only : int64
    use precision,  only : dp
    use fesm_clock

    implicit none

    type(clock_t)   :: clk
    type(cadence_t) :: c_mon
    integer(int64)  :: n0, dt, k
    integer :: nfail, nfire
    logical :: exact, exact_int
    real(dp) :: t0, t1, tprev

    nfail = 0

    ! ----- dtt = 0.1, 1e6 steps -----------------------------------------------
    call clock_init(clk, -130000.0_dp, -30000.0_dp, 0.1_dp, "365_day", "tropical")
    call clock_start(clk)
    n0 = clock_to_ticks(clk, -130000.0_dp, "time_init")
    dt = clock_to_ticks(clk, 0.1_dp, "dtt")
    call check("dtt = 0.1 is 3600 ticks", dt == 3600_int64)
    call check("clock_dtt returns 0.1 exactly", clock_dtt(clk) == 0.1_dp)

    exact = .true.
    k = 0
    do while (.not. clock_finished(clk))
        call clock_advance(clk)
        k = clock_step_index(clk)
        if (clock_time(clk) /= clock_to_years(clk, n0 + k*dt)) exact = .false.
    end do
    call check("1e6 steps taken", k == 1000000_int64)
    call check("time() == to_years(n_init + k*dt) at every step", exact)
    call check("ends exactly on time_end = -30000", clock_time(clk) == -30000.0_dp)
    call check("time_elapsed = 1e5 exactly", clock_time_elapsed(clk) == 100000.0_dp)

    ! ----- dtt = 1/12 via a per_year cadence ----------------------------------
    call clock_init(clk, 0.0_dp, 100.0_dp, 1.0_dp/12.0_dp, "360_day", "360d")
    call cadence_init(c_mon, clk, "month", per_year=12)
    call clock_start(clk)
    call check("dtt = 1/12 is 3000 ticks", clock_to_ticks(clk, clock_dtt(clk), "dtt") == 3000_int64)

    nfire = 0
    exact = .true.
    exact_int = .true.
    tprev = 0.0_dp
    do while (.not. clock_finished(clk))
        call clock_advance(clk)
        k = clock_step_index(clk)
        if (clock_time(clk) /= clock_to_years(clk, k*3000_int64)) exact = .false.
        if (mod(k, 12_int64) == 0) then
            if (clock_time(clk) /= real(k/12, dp)) exact_int = .false.
        end if
        if (cadence_due(c_mon, clk)) then
            call cadence_interval(c_mon, clk, t0, t1)
            if (t0 /= tprev) exact = .false.
            tprev = t1
            call cadence_mark(c_mon, clk)
            nfire = nfire + 1
        end if
    end do
    call check("1200 monthly steps", clock_step_index(clk) == 1200_int64)
    call check("per_year = 12 cadence due every step", nfire == 1200)
    call check("time() exact and t1 of one interval == t0 of the next", exact)
    call check("whole years exact every 12 steps", exact_int)
    call check("ends exactly on time_end = 100", clock_time(clk) == 100.0_dp)
    call check("cadence_dt_years = 1/12", cadence_dt_years(c_mon, clk) == clock_to_years(clk, 3000_int64))

    call summary("test_clock_exact")

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

end program test_clock_exact
