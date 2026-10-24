program test_varslice_interval
    ! slice_method = "interval_mean" (varslice V5): the overlap-weighted mean
    ! over [t0,t1] of the linear interpolant (irregular axis, partial
    ! segments, hold outside), of cell means over CF time bounds, a single
    ! time (the interpolated value), and rep = nsub on a monthly axis (each
    ! year standing for [year, year+1)).

    use precision
    use ncio
    use varslice, only : varslice_class, varslice_init_nml, varslice_update, varslice_end

    implicit none

    character(len=*), parameter :: fnml = "test_ds_interval.nml"
    real(dp), parameter :: a = 3.0_dp, b = 0.25_dp
    real(wp), parameter :: tol = 2e-5_wp
    type(varslice_class) :: vs
    real(dp) :: t4(4), tb(2,3), tm(36)
    real(wp) :: v(2,2,36)
    integer  :: k, nfail

    nfail = 0

    ! --- linear data on an irregular axis (yr CE) -----------------------------
    t4 = [0.0_dp, 10.0_dp, 30.0_dp, 60.0_dp]
    do k = 1, 4
        v(:,:,k) = real(a + b*t4(k), wp)
    end do
    call ncfile("test_ds_interval_lin.nc", v(:,:,1:4), t4, "yr CE")

    ! --- cell means over time bounds: cells [0,10], [10,30], [30,60] ----------
    tb = reshape([0.0_dp,10.0_dp, 10.0_dp,30.0_dp, 30.0_dp,60.0_dp], [2,3])
    v(:,:,1) = 1.0_wp; v(:,:,2) = 2.0_wp; v(:,:,3) = 4.0_wp
    call ncfile("test_ds_interval_bnds.nc", v(:,:,1:3), 0.5_dp*(tb(1,:)+tb(2,:)), "yr CE")
    call nc_write_dim("test_ds_interval_bnds.nc", "nv", x=1.0_dp, dx=1.0_dp, nx=2)
    call nc_write("test_ds_interval_bnds.nc", "time_bnds", tb, dim1="nv", dim2="time")
    call nc_write_attr("test_ds_interval_bnds.nc", "time", "bounds", "time_bnds")

    ! --- monthly, 3 years (noleap, mid-month); value = 100*year_index + month --
    do k = 1, 36
        tm(k) = 365.0_dp*((k-1)/12) + 15.0_dp + 30.0_dp*mod(k-1,12)
        v(:,:,k) = real(100*((k-1)/12 + 1) + mod(k-1,12) + 1, wp)
    end do
    call ncfile("test_ds_interval_mon.nc", v, tm, "days since 2000-01-01", cal="noleap")

    call write_nml()

    ! Linear interpolant
    call varslice_init_nml(vs, fnml, "lin", verbose=.false.)
    call varslice_update(vs, [5.0_dp, 45.0_dp], method="interval_mean")
    call check("linear data: mean = value at the midpoint", near(vs%var(1,1,1,1), a + b*25.0_dp))
    call varslice_update(vs, [12.0_dp, 18.0_dp], method="interval_mean")
    call check("inside one segment", near(vs%var(1,1,1,1), a + b*15.0_dp))
    call varslice_update(vs, [-10.0_dp, 5.0_dp], method="interval_mean")
    call check("hold before the axis + partial segment", &
               near(vs%var(1,1,1,1), (10.0_dp*a + 5.0_dp*(a + b*2.5_dp))/15.0_dp))
    call varslice_update(vs, [50.0_dp, 80.0_dp], method="interval_mean")
    call check("partial segment + hold after the axis", &
               near(vs%var(2,2,1,1), (10.0_dp*(a + b*55.0_dp) + 20.0_dp*(a + b*60.0_dp))/30.0_dp))
    call varslice_update(vs, [20.0_dp], method="interval_mean")
    call check("single time: interpolated value", near(vs%var(1,2,1,1), a + b*20.0_dp))
    call varslice_update(vs, [10.0_dp, 30.0_dp], method="range_mean")
    call check("range_mean unchanged (samples in the closed range)", &
               near(vs%var(1,1,1,1), a + b*20.0_dp))
    call varslice_end(vs)

    ! Cell means over the time bounds
    call varslice_init_nml(vs, fnml, "bnds", verbose=.false.)
    call check("bounds read", allocated(vs%time_bnds))
    if (allocated(vs%time_bnds)) call check("bounds values", all(vs%time_bnds .eq. tb))
    call varslice_update(vs, [5.0_dp, 40.0_dp], method="interval_mean")
    call check("overlap-weighted cell means", near(vs%var(1,1,1,1), (5.0_dp*1 + 20.0_dp*2 + 10.0_dp*4)/35.0_dp))
    call varslice_update(vs, [-5.0_dp, 5.0_dp], method="interval_mean")
    call check("cells: hold before", near(vs%var(1,1,1,1), 1.0_dp))
    call varslice_update(vs, [15.0_dp], method="interval_mean")
    call check("cells: single time in a cell", near(vs%var(1,1,1,1), 2.0_dp))
    call varslice_end(vs)

    ! Monthly with rep = 12
    call varslice_init_nml(vs, fnml, "mon", verbose=.false.)
    call varslice_update(vs, [2001.0_dp, 2002.0_dp], method="interval_mean", rep=12)
    call check("rep=12 over a calendar year: its cycle", size(vs%var,3) .eq. 12 .and. &
               near(vs%var(1,1,1,1), 201.0_dp) .and. near(vs%var(1,1,12,1), 212.0_dp))
    call varslice_update(vs, [2000.5_dp, 2002.5_dp], method="interval_mean", rep=12)
    call check("rep=12 over two years: overlap-weighted years", &
               near(vs%var(1,1,3,1), (0.5_dp*103 + 1.0_dp*203 + 0.5_dp*303)/2.0_dp))
    call varslice_update(vs, [2001.3_dp], method="interval_mean", rep=12)
    call check("rep=12 at a time: the cycle of its year", near(vs%var(1,1,7,1), 207.0_dp))
    call varslice_update(vs, [2001.0_dp, 2002.0_dp], method="interval_mean", rep=1)
    call check("rep=1 on a monthly axis: linear interpolant mean", abs(vs%var(1,1,1,1) - 206.5_wp) .lt. 1.0_wp)
    call varslice_end(vs)

    if (nfail .gt. 0) then
        write(*,*) "test_varslice_interval: ", nfail, " FAILED"
        error stop 1
    end if
    write(*,*) "test_varslice_interval: all passed"

contains

    logical function near(x, ref)
        real(wp), intent(in) :: x
        real(dp), intent(in) :: ref
        near = abs(real(x,dp) - ref) .le. tol*max(1.0_dp, abs(ref))
        if (.not. near) write(*,*) "      got ", x, " expected ", ref
    end function near

    subroutine ncfile(fname, v, t, units, cal)
        character(len=*), intent(in) :: fname, units
        real(wp),         intent(in) :: v(:,:,:)
        real(dp),         intent(in) :: t(:)
        character(len=*), intent(in), optional :: cal
        call nc_create(fname)
        call nc_write_dim(fname, "xc", x=1.0_dp, dx=1.0_dp, nx=size(v,1))
        call nc_write_dim(fname, "yc", x=1.0_dp, dx=1.0_dp, nx=size(v,2))
        call nc_write_dim(fname, "time", x=t, units=units, calendar=cal, unlimited=.false.)
        call nc_write(fname, "v", v, dim1="xc", dim2="yc", dim3="time")
    end subroutine ncfile

    subroutine write_nml()
        integer :: u
        open(newunit=u, file=fnml, status="replace", action="write")
        call grp(u, "lin",  "test_ds_interval_lin.nc")
        call grp(u, "bnds", "test_ds_interval_bnds.nc")
        call grp(u, "mon",  "test_ds_interval_mon.nc")
        close(u)
    end subroutine write_nml

    subroutine grp(u, name, file)
        integer,          intent(in) :: u
        character(len=*), intent(in) :: name, file
        write(u,"(a)") "&"//name
        write(u,"(a)") '    filename = "'//file//'"'
        write(u,"(a)") '    name     = "v"'
        write(u,"(a)") '    units    = "1" "1"'
        write(u,"(a)") '    scaling  = 1.0 0.0'
        write(u,"(a)") '    time     = -1 0 0 0 0'
        write(u,"(a)") "/"
    end subroutine grp

    subroutine check(label, ok)
        character(len=*), intent(in) :: label
        logical,          intent(in) :: ok
        write(*,"(a,a)") merge("  ok   ", "  FAIL ", ok), label
        if (.not. ok) nfail = nfail + 1
    end subroutine check

end program test_varslice_interval
