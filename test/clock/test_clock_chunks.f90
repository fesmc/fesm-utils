program test_clock_chunks
    ! fesm_clock: chunks with different dtt fire every cadence at the times of
    ! one long run.
    !   run A:  0 -> 200, dtt = 1;
    !   run B:  0 -> 100, dtt = 1, restart bundle at 100, then
    !           100 -> 200, dtt = 0.5, continued from the bundle.
    ! Cadences geo (10 yr) and clim (5 yr): identical firing times and
    ! identical (t0, t1) intervals, bit for bit; time_elapsed continues.
    ! Writes logs/clock/chunks_restart.nc (run from the repo root).

    use precision,  only : dp
    use fesm_clock

    implicit none

    character(len=*), parameter :: frst = "logs/clock/chunks_restart.nc"
    integer, parameter :: nmax = 100
    real(dp) :: geo_a(2,nmax), clim_a(2,nmax), geo_b(2,nmax), clim_b(2,nmax)
    integer  :: ngeo_a, nclim_a, ngeo_b, nclim_b, nfail
    real(dp) :: elapsed_a, elapsed_b
    logical  :: continued

    nfail = 0
    call execute_command_line("mkdir -p logs/clock && rm -f "//frst)

    ngeo_a = 0; nclim_a = 0; ngeo_b = 0; nclim_b = 0

    call run_chunk(0.0_dp, 200.0_dp, 1.0_dp, "", "", geo_a, ngeo_a, clim_a, nclim_a, elapsed_a, continued)
    call run_chunk(0.0_dp, 100.0_dp, 1.0_dp, "", frst, geo_b, ngeo_b, clim_b, nclim_b, elapsed_b, continued)
    call check("chunk 1 is not a continuation", .not. continued)
    call run_chunk(100.0_dp, 200.0_dp, 0.5_dp, frst, "", geo_b, ngeo_b, clim_b, nclim_b, elapsed_b, continued)
    call check("chunk 2 is a continuation", continued)

    call check("geo fires 20 times in both", ngeo_a == 20 .and. ngeo_b == 20)
    call check("clim fires 40 times in both", nclim_a == 40 .and. nclim_b == 40)
    if (ngeo_a == ngeo_b) call check("geo intervals identical (bitwise)", &
        all(geo_a(:,1:ngeo_a) == geo_b(:,1:ngeo_b)))
    if (nclim_a == nclim_b) call check("clim intervals identical (bitwise)", &
        all(clim_a(:,1:nclim_a) == clim_b(:,1:nclim_b)))
    call check("geo interval [100, 110] after the restart", any(geo_b(1,1:ngeo_b) == 100.0_dp .and. &
        geo_b(2,1:ngeo_b) == 110.0_dp))
    call check("time_elapsed = 200 at the end of both", elapsed_a == 200.0_dp .and. elapsed_b == 200.0_dp)

    call summary("test_clock_chunks")

contains

    subroutine run_chunk(ti, te, dtt, rst_in, rst_out, geo, ngeo, clim, nclim, elapsed, continued)
        real(dp),         intent(in)    :: ti, te, dtt
        character(len=*), intent(in)    :: rst_in, rst_out
        real(dp),         intent(inout) :: geo(:,:), clim(:,:)
        integer,          intent(inout) :: ngeo, nclim
        real(dp),         intent(out)   :: elapsed
        logical,          intent(out)   :: continued

        type(clock_t)    :: clk
        type(cadence_t)  :: c_geo, c_clim
        type(schedule_t) :: s_rst
        real(dp) :: t0, t1

        call clock_init(clk, ti, te, dtt, "365_day", "tropical")
        call cadence_init(c_geo,  clk, "geo",  dt=10.0_dp)
        call cadence_init(c_clim, clk, "clim", dt=5.0_dp)
        if (len(rst_out) > 0) then
            call schedule_init(s_rst, clk, "restart", times=[te], restart=.true.)
        else
            call schedule_init(s_rst, clk, "restart", restart=.true.)
        end if
        if (len(rst_in) > 0) then
            call clock_restart_read(clk, rst_in)
            call cadence_restart_read(c_geo,  clk, rst_in)
            call cadence_restart_read(c_clim, clk, rst_in)
        end if
        call clock_start(clk)

        do
            if (schedule_due(s_rst, clk)) then
                call clock_restart_write(clk, rst_out)
                call cadence_restart_write(c_geo, rst_out)
                call cadence_restart_write(c_clim, rst_out)
            end if
            if (clock_finished(clk)) exit

            call clock_advance(clk)

            if (cadence_due(c_geo, clk)) then
                call cadence_interval(c_geo, clk, t0, t1)
                ngeo = ngeo + 1
                geo(:,ngeo) = [t0, t1]
                call cadence_mark(c_geo, clk)
            end if
            if (cadence_due(c_clim, clk)) then
                call cadence_interval(c_clim, clk, t0, t1)
                nclim = nclim + 1
                clim(:,nclim) = [t0, t1]
                call cadence_mark(c_clim, clk)
            end if
        end do
        elapsed   = clock_time_elapsed(clk)
        continued = clk%continued

    end subroutine run_chunk

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

end program test_clock_chunks
