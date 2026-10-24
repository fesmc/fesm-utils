program test_clock_convert
    ! fesm_clock conversions (N = 36000): to_ticks accepts every decimal
    ! literal that denotes a whole number of ticks (0.1, 0.01, -21e3, 1/12 to
    ! 16 digits) and rejects the rest (1/12 to 15 digits, 0.3333, beyond 2^53
    ! ticks); to_years is correctly rounded and odd.
    ! Self-contained: no files.

    use, intrinsic :: iso_fortran_env, only : int64
    use precision,  only : dp
    use fesm_clock

    implicit none

    type(clock_t)  :: clk
    integer        :: nfail, stat
    integer(int64) :: n

    nfail = 0
    call clock_init(clk, 0.0_dp, 1.0_dp, 1.0_dp, "365_day", "tropical")
    call check("default ticks_per_year = 36000", clk%tpy == CLOCK_TICKS_PER_YEAR .and. clk%tpy == 36000_int64)

    call check("to_ticks(0.1) = 3600",      clock_to_ticks(clk, 0.1_dp, "x") == 3600_int64)
    call check("to_ticks(0.01) = 360",      clock_to_ticks(clk, 0.01_dp, "x") == 360_int64)
    call check("to_ticks(0.001) = 36",      clock_to_ticks(clk, 0.001_dp, "x") == 36_int64)
    call check("to_ticks(1/32) = 1125",     clock_to_ticks(clk, 0.03125_dp, "x") == 1125_int64)
    call check("to_ticks(-21e3) = -7.56e8", clock_to_ticks(clk, -21.0e3_dp, "x") == -756000000_int64)
    call check("to_ticks(0.0833333333333333) = 3000 (16 digits)", &
        clock_to_ticks(clk, 0.0833333333333333_dp, "x") == 3000_int64)
    call check("to_ticks(1.0/12.0) = 3000", clock_to_ticks(clk, 1.0_dp/12.0_dp, "x") == 3000_int64)
    call check("to_ticks(1.0/3.0) = 12000", clock_to_ticks(clk, 1.0_dp/3.0_dp, "x") == 12000_int64)

    n = clock_to_ticks(clk, 0.083333333333333_dp, "x", stat=stat)
    call check("to_ticks(0.083333333333333) rejected (15 digits)", stat == 1)
    n = clock_to_ticks(clk, 0.3333_dp, "dt", stat=stat)
    call check("to_ticks(0.3333) rejected", stat == 1)
    n = clock_to_ticks(clk, 1.0e12_dp, "t", stat=stat)
    call check("to_ticks(1e12) rejected (beyond 2^53 ticks)", stat == 1)
    n = clock_to_ticks(clk, 2.5e11_dp, "t", stat=stat)
    call check("to_ticks(2.5e11) accepted (below 2^53 ticks)", stat == 0 .and. n == 9000000000000000_int64)

    call check("is_representable(0.1)", clock_is_representable(clk, 0.1_dp))
    call check("not is_representable(0.3333)", .not. clock_is_representable(clk, 0.3333_dp))
    call check("not is_representable(1/365)", .not. clock_is_representable(clk, 1.0_dp/365.0_dp))
    call check("to_years(3600) == 0.1 (correctly rounded)", clock_to_years(clk, 3600_int64) == 0.1_dp)
    call check("to_years(3000) == 1.0/12.0", clock_to_years(clk, 3000_int64) == 1.0_dp/12.0_dp)
    call check("to_years(-3600) == -0.1", clock_to_years(clk, -3600_int64) == -0.1_dp)

    call clock_init(clk, 0.0_dp, 1.0_dp, 1.0_dp, "365_day", "tropical", ticks_per_year=1800)
    call check("N = 1800: 0.001 not representable", .not. clock_is_representable(clk, 0.001_dp))
    call check("N = 1800: 1/12 representable", clock_is_representable(clk, 1.0_dp/12.0_dp))

    call summary("test_clock_convert")

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

end program test_clock_convert
