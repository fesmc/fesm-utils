program test_clock_negative
    ! fesm_clock at large negative times: -1e8 -> -1e8 + 1000 yr, dtt = 0.5.
    !   * cadences 10 yr (phase 0) and 20 yr (phase 5) fire on the absolute
    !     grids (modulo with floor semantics), with exact intervals;
    !   * to_years(-n) == -to_years(n);
    !   * the grids 0 mod 10 and 5 mod 20 are disjoint: valid without a restart
    !     schedule (describe reports the empty restart grid), is_multiple_of_all.
    ! Self-contained: no files.

    use, intrinsic :: iso_fortran_env, only : int64
    use precision,  only : dp
    use fesm_clock

    implicit none

    type(clock_t)   :: clk
    type(cadence_t) :: c10, c20
    integer :: nfail, n10, n20
    integer(int64) :: n
    logical :: ok10, ok20, oksym
    real(dp) :: t0, t1, first20

    nfail = 0

    call clock_init(clk, -1.0e8_dp, -1.0e8_dp + 1000.0_dp, 0.5_dp, "365_day", "tropical")
    call cadence_init(c10, clk, "c10", dt=10.0_dp)
    call cadence_init(c20, clk, "c20", dt=20.0_dp, phase=5.0_dp)
    call clock_start(clk)
    call check("next c10 time after -1e8 is -1e8 + 10", cadence_t_next(c10, clk) == -99999990.0_dp)
    call check("next c20 time after -1e8 is -1e8 + 5", cadence_t_next(c20, clk) == -99999995.0_dp)

    n10 = 0
    n20 = 0
    ok10 = .true.
    ok20 = .true.
    first20 = 0.0_dp
    do while (.not. clock_finished(clk))
        call clock_advance(clk)
        if (cadence_due(c10, clk)) then
            call cadence_interval(c10, clk, t0, t1)
            if (modulo(t1, 10.0_dp) /= 0.0_dp) ok10 = .false.
            if (n10 > 0 .and. t1 - t0 /= 10.0_dp) ok10 = .false.
            call cadence_mark(c10, clk)
            n10 = n10 + 1
        end if
        if (cadence_due(c20, clk)) then
            call cadence_interval(c20, clk, t0, t1)
            if (modulo(t1 - 5.0_dp, 20.0_dp) /= 0.0_dp) ok20 = .false.
            if (n20 == 0) first20 = t1
            if (n20 > 0 .and. t1 - t0 /= 20.0_dp) ok20 = .false.
            call cadence_mark(c20, clk)
            n20 = n20 + 1
        end if
    end do
    call check("c10 fires 100 times, on multiples of 10, intervals 10", n10 == 100 .and. ok10)
    call check("c20 (phase 5) fires 50 times on 5 + 20k, intervals 20", n20 == 50 .and. ok20)
    call check("c20 first interval [-1e8, -1e8 + 5] (shorter than dt)", first20 == -99999995.0_dp)
    call check("ends exactly on -1e8 + 1000", clock_time(clk) == -99999000.0_dp)
    call check("is_multiple_of_all false at the end (c20 grid)", .not. is_multiple_of_all(clk%n, [c10, c20]))
    call check("is_multiple_of_all true for c10 alone", is_multiple_of_all(clk%n, [c10]))

    oksym = .true.
    do n = 1_int64, 3600000000000_int64, 999999937_int64
        if (clock_to_years(clk, -n) /= -clock_to_years(clk, n)) oksym = .false.
    end do
    call check("to_years(-n) == -to_years(n)", oksym)
    call check("to_ticks(-1e8) = -3.6e12", clock_to_ticks(clk, -1.0e8_dp, "t") == -3600000000000_int64)

    call summary("test_clock_negative")

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

end program test_clock_negative
