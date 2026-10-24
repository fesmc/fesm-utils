program test_clock_forcing
    ! fesm_clock views for method rel, cal, const (time_const and
    ! time_const_cal): time_rel, time_cal, and time_cal = time_rel + time_ref
    ! (exact in ticks, to an ulp in years) at every step.
    ! Self-contained: no files.

    use, intrinsic :: iso_fortran_env, only : int64
    use precision,  only : dp
    use fesm_clock

    implicit none

    type(clock_t) :: clk
    integer :: nfail
    logical :: ok_rel, ok_sum

    nfail = 0

    ! ----- rel ----------------------------------------------------------------------
    call clock_init(clk, -21000.0_dp, -20990.0_dp, 0.1_dp, "365_day", "tropical", method="rel")
    call clock_start(clk)
    call run_views(ok_rel, ok_sum, "rel")
    call check("rel: time_rel = time at every step", ok_rel)
    call check("rel: time_cal = time_rel + 1950 at every step", ok_sum)
    call check("rel: time_cal at the end = -19040", clock_time_cal(clk) == -19040.0_dp)

    ! ----- cal ----------------------------------------------------------------------
    call clock_init(clk, 1850.0_dp, 1860.0_dp, 0.1_dp, "365_day", "tropical", method="cal", time_ref=1950.0_dp)
    call clock_start(clk)
    call run_views(ok_rel, ok_sum, "cal")
    call check("cal: time_cal = time at every step", ok_rel)
    call check("cal: time_cal = time_rel + 1950 at every step", ok_sum)
    call check("cal: time_rel at the end = -90", clock_time_rel(clk) == -90.0_dp)

    ! ----- const, time_const ---------------------------------------------------------
    call clock_init(clk, 0.0_dp, 10.0_dp, 0.1_dp, "365_day", "tropical", method="const", time_const=-21000.0_dp)
    call clock_start(clk)
    call run_views(ok_rel, ok_sum, "const")
    call check("const (time_const): time_rel frozen at -21000", ok_rel)
    call check("const (time_const): time_cal = time_rel + 1950", ok_sum)
    call check("const (time_const): model time advances to 10", clock_time(clk) == 10.0_dp)
    call check("const (time_const): time_cal = -19050", clock_time_cal(clk) == -19050.0_dp)

    ! ----- const, time_const_cal ------------------------------------------------------
    call clock_init(clk, 0.0_dp, 10.0_dp, 0.1_dp, "365_day", "tropical", method="const", &
                    time_const_cal=1850.0_dp)
    call clock_start(clk)
    call check("const (time_const_cal): time_rel = -100", clock_time_rel(clk) == -100.0_dp)
    call check("const (time_const_cal): time_cal = 1850", clock_time_cal(clk) == 1850.0_dp)
    call clock_advance(clk)
    call check("const (time_const_cal): frozen after advance", clock_time_cal(clk) == 1850.0_dp &
        .and. clock_time(clk) == 0.1_dp)

    ! ----- time_ref ---------------------------------------------------------------------
    call clock_init(clk, 0.0_dp, 10.0_dp, 1.0_dp, "365_day", "tropical", time_ref=2000.0_dp)
    call check("rel, time_ref = 2000: time_cal = 2000 at 0", clock_time_cal(clk) == 2000.0_dp)

    call summary("test_clock_forcing")

contains

    subroutine run_views(ok_rel, ok_sum, method)
        ! Step to the end, checking the views at every step.
        logical,          intent(out) :: ok_rel, ok_sum
        character(len=*), intent(in)  :: method
        integer(int64) :: n_ref, n_rel, n_cal
        real(dp) :: frozen

        ok_rel = .true.
        ok_sum = .true.
        n_ref  = clock_to_ticks(clk, 1950.0_dp, "time_ref")
        frozen = clock_time_rel(clk)
        do
            select case (method)
            case ("rel")
                if (clock_time_rel(clk) /= clock_time(clk)) ok_rel = .false.
            case ("cal")
                if (clock_time_cal(clk) /= clock_time(clk)) ok_rel = .false.
            case ("const")
                if (clock_time_rel(clk) /= frozen .or. frozen /= -21000.0_dp) ok_rel = .false.
            end select
            n_rel = clock_to_ticks(clk, clock_time_rel(clk), "time_rel")
            n_cal = clock_to_ticks(clk, clock_time_cal(clk), "time_cal")
            if (n_cal /= n_rel + n_ref) ok_sum = .false.
            if (abs(clock_time_cal(clk) - (clock_time_rel(clk) + 1950.0_dp)) > spacing(clock_time_cal(clk))) &
                ok_sum = .false.
            if (clock_finished(clk)) exit
            call clock_advance(clk)
        end do
    end subroutine run_views

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

end program test_clock_forcing
