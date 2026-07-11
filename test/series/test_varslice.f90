program test_varslice
    ! Exercise varslice time-slicing against synthetic NetCDF data with known
    ! values, confirming the fixes to:
    !   #1  varslice_init_arg now allocates par%filenames / sets with_time_sub
    !   #2  slice_method="range_sum" returns data (previously all missing)
    !   #3  time axis longer than the old 10000-element work arrays
    !   #4  multi-file globbing via get_matching_files (+ nc_read_multifile)
    !
    ! Synthetic field temp(x,time) has the analytic value
    !     temp(i,t) = time(t) + (i-1)*100
    ! so every slice/aggregate has a hand-checkable expected value.

    use precision
    use ncio,     only : nc_create, nc_write_dim, nc_write
    use varslice, only : varslice_class, varslice_init_arg, varslice_update, varslice_end

    implicit none

    integer, parameter :: nx = 3
    real(wp), parameter :: rtol = 1.0e-4_wp

    type(varslice_class) :: vs
    integer :: nfail
    real(wp) :: exp1(nx)

    nfail = 0

    ! --- build synthetic input files ----------------------------------------
    call make_file("test_vs_single.nc",   nx, 1, 10)   ! one file, years 1..10
    call make_file("test_vs_multi_01.nc", nx, 1,  5)   ! first half
    call make_file("test_vs_multi_02.nc", nx, 6, 10)   ! second half
    call make_file("test_vs_big.nc",       1, 1, 12000) ! long time axis (>10000)

    ! ========================================================================
    ! Single-file cases (also exercises fix #1: init via varslice_init_arg)
    ! ========================================================================
    call varslice_init_arg(vs, "test_vs_single.nc", "temp", "u", "u", &
                           with_time=.TRUE., time_par=[1.0_wp,10.0_wp,1.0_wp,0.0_wp])

    ! exact @ time=5  -> [5,105,205]
    call varslice_update(vs, time=[5.0_wp], method="exact")
    exp1 = val_col(5.0_wp, nx)
    call check_shape("exact shape", vs, nx, 1, nfail)
    call check_col("exact @5", vs%var(:,1,1,1), exp1, nfail)

    ! range @ [3,5] -> 3 time slices (t=3,4,5)
    call varslice_update(vs, time=[3.0_wp,5.0_wp], method="range")
    call check_shape("range shape", vs, nx, 3, nfail)
    call check_col("range col t=3", vs%var(:,1,1,1), val_col(3.0_wp,nx), nfail)
    call check_col("range col t=4", vs%var(:,2,1,1), val_col(4.0_wp,nx), nfail)
    call check_col("range col t=5", vs%var(:,3,1,1), val_col(5.0_wp,nx), nfail)

    ! range_mean @ [1,10] -> mean over t=1..10 (per point): time-mean = 5.5
    call varslice_update(vs, time=[1.0_wp,10.0_wp], method="range_mean", rep=1)
    call check_shape("range_mean shape", vs, nx, 1, nfail)
    call check_col("range_mean @[1,10]", vs%var(:,1,1,1), val_col(5.5_wp,nx), nfail)

    ! range_sum @ [1,10] -> sum over t=1..10 (per point).  FIX #2:
    ! previously get_indices had no "range_sum" case, so this returned missing.
    ! point i: sum_{t=1..10}(t) + 10*(i-1)*100 = 55 + 1000*(i-1)
    call varslice_update(vs, time=[1.0_wp,10.0_wp], method="range_sum", rep=1)
    call check_shape("range_sum shape", vs, nx, 1, nfail)
    call check_col("range_sum @[1,10]", vs%var(:,1,1,1), &
                   [55.0_wp, 1055.0_wp, 2055.0_wp], nfail)

    ! interp @ 4.5 -> halfway between t=4 and t=5 -> [4.5,104.5,204.5]
    call varslice_update(vs, time=[4.5_wp], method="interp")
    call check_shape("interp shape", vs, nx, 1, nfail)
    call check_col("interp @4.5", vs%var(:,1,1,1), val_col(4.5_wp,nx), nfail)

    call varslice_end(vs)

    ! ========================================================================
    ! Multi-file case (fix #4: glob + get_matching_files + nc_read_multifile)
    ! ========================================================================
    call varslice_init_arg(vs, "test_vs_multi_*.nc", "temp", "u", "u", &
                           with_time=.TRUE., time_par=[1.0_wp,10.0_wp,1.0_wp,0.0_wp])

    ! exact @ 7 -> value from the SECOND file, proving the files were stitched
    call varslice_update(vs, time=[7.0_wp], method="exact")
    call check_col("multifile exact @7", vs%var(:,1,1,1), val_col(7.0_wp,nx), nfail)

    ! range @ [4,7] -> spans both files: t=4,5 (file1) and 6,7 (file2)
    call varslice_update(vs, time=[4.0_wp,7.0_wp], method="range")
    call check_shape("multifile range shape", vs, nx, 4, nfail)
    call check_col("multifile range t=4", vs%var(:,1,1,1), val_col(4.0_wp,nx), nfail)
    call check_col("multifile range t=6", vs%var(:,3,1,1), val_col(6.0_wp,nx), nfail)
    call check_col("multifile range t=7", vs%var(:,4,1,1), val_col(7.0_wp,nx), nfail)

    call varslice_end(vs)

    ! ========================================================================
    ! Long time axis (fix #3: work arrays sized to n, not fixed at 10000)
    ! ========================================================================
    call varslice_init_arg(vs, "test_vs_big.nc", "temp", "u", "u", &
                           with_time=.TRUE., time_par=[1.0_wp,12000.0_wp,1.0_wp,0.0_wp])

    ! exact @ 9000 (index 9000 > 10000-th only in count, but n=12000 array) -> 9000
    call varslice_update(vs, time=[9000.0_wp], method="exact")
    call check_col("big exact @9000", vs%var(:,1,1,1), [9000.0_wp], nfail)

    ! range_mean over the whole 12000-long axis -> time-mean = 6000.5
    call varslice_update(vs, time=[1.0_wp,12000.0_wp], method="range_mean", rep=1)
    call check_col("big range_mean", vs%var(:,1,1,1), [6000.5_wp], nfail)

    call varslice_end(vs)

    ! --- cleanup synthetic files --------------------------------------------
    call execute_command_line("rm -f test_vs_single.nc test_vs_multi_01.nc "// &
                              "test_vs_multi_02.nc test_vs_big.nc")

    ! --- report -------------------------------------------------------------
    write(*,*)
    if (nfail == 0) then
        write(*,*) "test_varslice: PASS"
    else
        write(*,*) "test_varslice: FAIL", nfail, "checks"
        stop 1
    end if

contains

    function val_col(t, n) result(col)
        ! Analytic column of temp at time t: temp(i,t) = t + (i-1)*100
        real(wp), intent(in) :: t
        integer,  intent(in) :: n
        real(wp) :: col(n)
        integer  :: i
        do i = 1, n
            col(i) = t + real(i-1,wp)*100.0_wp
        end do
    end function val_col

    subroutine make_file(filename, n, y0, y1)
        ! Write temp(x=n, time=y0..y1) with temp(i,t) = time(t) + (i-1)*100
        character(len=*), intent(in) :: filename
        integer,          intent(in) :: n, y0, y1
        integer  :: i, t, nt
        real(dp), allocatable :: x(:), time(:)
        real(wp), allocatable :: var(:,:)

        nt = y1 - y0 + 1
        allocate(x(n), time(nt), var(n,nt))
        do i = 1, n
            x(i) = real(i-1,dp)
        end do
        do t = 1, nt
            time(t) = real(y0 + (t-1),dp)
            do i = 1, n
                var(i,t) = real(time(t),wp) + real(i-1,wp)*100.0_wp
            end do
        end do

        call nc_create(filename)
        call nc_write_dim(filename, "x",    x=x)
        call nc_write_dim(filename, "time", x=time)
        call nc_write(filename, "temp", var, dim1="x", dim2="time")

        deallocate(x, time, var)
    end subroutine make_file

    subroutine check_shape(label, vs, nsp, ntime, nfail)
        character(len=*),     intent(in)    :: label
        type(varslice_class), intent(in)    :: vs
        integer,              intent(in)    :: nsp, ntime
        integer,              intent(inout) :: nfail
        if (size(vs%var,1) == nsp .and. size(vs%var,2) == ntime) then
            write(*,"(a,a)") "  PASS ", label
        else
            write(*,"(a,a,4i6)") "  FAIL ", label//"  shape=", shape(vs%var)
            nfail = nfail + 1
        end if
    end subroutine check_shape

    subroutine check_col(label, got, expected, nfail)
        character(len=*), intent(in)    :: label
        real(wp),         intent(in)    :: got(:), expected(:)
        integer,          intent(inout) :: nfail
        real(wp) :: err
        if (size(got) /= size(expected)) then
            write(*,"(a,a,2i6)") "  FAIL ", label//"  size mismatch", &
                                 size(got), size(expected)
            nfail = nfail + 1
            return
        end if
        err = maxval(abs(got - expected))
        if (err <= rtol*maxval(abs(expected)) + rtol) then
            write(*,"(a,a,es12.4)") "  PASS ", label//"  err=", err
        else
            write(*,"(a,a,es12.4)") "  FAIL ", label//"  err=", err
            write(*,"(a,10f12.3)")  "        got     =", got
            write(*,"(a,10f12.3)")  "        expected=", expected
            nfail = nfail + 1
        end if
    end subroutine check_col

end program test_varslice
