program test_clock_calendar
    ! fesm_clock calendar service: day_of_year, month, day_of_month at known
    ! year fractions for 360_day and 365_day (also at negative times),
    ! sec_day; calendar none stops on query (here via stat).
    ! Self-contained: no files.

    use precision,      only : dp
    use phys_constants, only : sec_day_360d_tropical
    use fesm_clock

    implicit none

    type(clock_t) :: clk
    integer :: nfail, stat, doy, i
    logical :: ok

    nfail = 0

    ! ----- 360_day: one step per day ---------------------------------------------------
    call clock_init(clk, 0.0_dp, 2.0_dp, 1.0_dp/360.0_dp, "360_day", "tropical")
    call clock_start(clk)
    call check("360_day t = 0: day 1, Jan 1", clock_day_of_year(clk) == 1 .and. clock_month(clk) == 1 &
        .and. clock_day_of_month(clk) == 1)
    ok = .true.
    do i = 1, 720
        call clock_advance(clk)
        doy = modulo(i, 360) + 1
        if (clock_day_of_year(clk) /= doy) ok = .false.
        if (clock_month(clk) /= (doy - 1)/30 + 1) ok = .false.
        if (clock_day_of_month(clk) /= doy - 30*((doy - 1)/30)) ok = .false.
        if (i == 31) then
            call check("360_day day 32 = Feb 2", clock_month(clk) == 2 .and. clock_day_of_month(clk) == 2)
        end if
        if (i == 359) then
            call check("360_day day 360 = Dec 30", clock_day_of_year(clk) == 360 .and. clock_month(clk) == 12 &
                .and. clock_day_of_month(clk) == 30)
        end if
    end do
    call check("360_day: day_of_year, month, day_of_month at every day of two years", ok)
    call check("360_day on a tropical year: sec_day = sec_day_360d_tropical", &
        clock_sec_day(clk) == sec_day_360d_tropical)

    ! ----- 365_day: quarter years -------------------------------------------------------
    call clock_init(clk, -1.0_dp, 1.0_dp, 0.25_dp, "noleap", "365d")
    call clock_start(clk)
    call check("noleap is 365_day", trim(clk%calendar) == "365_day")
    call check("365_day t = -1: Jan 1", clock_day_of_year(clk) == 1 .and. clock_month(clk) == 1)
    call clock_advance(clk)
    call check("365_day t = -0.75 (fraction 0.25): day 92 = Apr 2", clock_year_fraction(clk) == 0.25_dp &
        .and. clock_day_of_year(clk) == 92 .and. clock_month(clk) == 4 .and. clock_day_of_month(clk) == 2)
    call clock_advance(clk)
    call check("365_day t = -0.5: day 183 = Jul 2", clock_day_of_year(clk) == 183 .and. &
        clock_month(clk) == 7 .and. clock_day_of_month(clk) == 2)
    call clock_advance(clk)
    call check("365_day t = -0.25: day 274 = Oct 1", clock_day_of_year(clk) == 274 .and. &
        clock_month(clk) == 10 .and. clock_day_of_month(clk) == 1)
    call clock_advance(clk)
    call check("365_day t = 0: Jan 1", clock_day_of_year(clk) == 1 .and. clock_year_fraction(clk) == 0.0_dp)
    call check("365_day on a 365d year: sec_day = 86400", clock_sec_day(clk) == 86400.0_dp)

    ! ----- none -------------------------------------------------------------------------
    call clock_init(clk, 0.0_dp, 1.0_dp, 0.5_dp, "none", "tropical")
    call clock_start(clk)
    doy = clock_day_of_year(clk, stat=stat)
    call check("calendar none: day_of_year stops (stat = 1)", stat == 1)
    doy = clock_month(clk, stat=stat)
    call check("calendar none: month stops (stat = 1)", stat == 1)
    call check("calendar none: sec_day stops (stat = 1)", clock_sec_day(clk, stat=stat) == 0.0_dp .and. stat == 1)
    call clock_advance(clk)
    call check("calendar none: year_fraction still defined (0.5)", clock_year_fraction(clk) == 0.5_dp)

    call summary("test_clock_calendar")

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

end program test_clock_calendar
