program test_anomaly
    ! fesm_datasrc anomalies: anomaly_ref = "slice" is exactly 0 at
    ! ref_time and the change since then elsewhere; anomaly_ref = "group"
    ! against a static group of the same quantity in other units (the
    ! conversion is applied before differencing).

    use precision
    use ncio
    use coordinates,  only : grid_class, grid_init
    use fesm_fields
    use fesm_datasrc

    implicit none

    character(len=*), parameter :: fpar = "test_ds_anom.nml", fdef = "par/datasrc_defaults.nml"
    integer, parameter :: nx = 3, ny = 2
    type(grid_class)    :: grid
    type(field_state_t) :: st
    type(datasrc_t)     :: ds
    real(wp) :: v(nx,ny,3), r(nx,ny)
    real(dp) :: a(nx,ny)
    integer  :: h, i, k, stat, nfail

    nfail = 0

    call grid_init(grid, name="TEST-DS-ANOM", mtype="cartesian", units="kilometers", &
                   x0=0.0_dp, dx=10.0_dp, nx=nx, y0=0.0_dp, dy=10.0_dp, ny=ny)

    ! T(t) = 5 + 0.1*t + i [degC] at t = 0, 100, 200 (yr CE); reference 278.15 + i [K]
    do k = 1, 3
        do i = 1, nx
            v(i,:,k) = real(5.0_dp + 0.1_dp*100.0_dp*(k-1) + i, wp)
        end do
    end do
    do i = 1, nx
        r(i,:) = real(278.15_dp + i, wp)
    end do
    call nc_create("test_ds_anom.nc")
    call nc_write_dim("test_ds_anom.nc", "xc", x=0.0_dp, dx=10.0_dp, nx=nx)
    call nc_write_dim("test_ds_anom.nc", "yc", x=0.0_dp, dx=10.0_dp, nx=ny)
    call nc_write_dim("test_ds_anom.nc", "time", x=[0.0_dp, 100.0_dp, 200.0_dp], units="yr CE", unlimited=.false.)
    call nc_write("test_ds_anom.nc", "T", v, dim1="xc", dim2="yc", dim3="time", units="degC")
    call nc_write("test_ds_anom.nc", "T_ref", r, dim1="xc", dim2="yc", units="K")

    call write_nml()

    call fields_init(st, "anom_test:exports")
    call fields_add(st, "dT",  "K", FIELD_INSTANT, FIELD_INTENSIVE, h, grid="TEST-DS-ANOM", dims=[nx,ny])
    call fields_add(st, "dTg", "K", FIELD_INSTANT, FIELD_INTENSIVE, h, grid="TEST-DS-ANOM", dims=[nx,ny])
    call fields_lock(st)

    call datasrc_init(ds, fpar, "anom", st, [character(len=4) :: "dT", "dTg"], grid, fdef, stat=stat)
    call check("init", stat .eq. 0)
    call datasrc_describe(ds)

    call datasrc_update(ds, 100.0_dp)
    call datasrc_get(ds, "dT", a)
    call check("slice: exactly 0 at ref_time", all(a .eq. 0.0_dp))
    call datasrc_get(ds, "dTg", a)
    call check("group: converted (degC -> K) before differencing", maxval(abs(a - 10.0_dp)) .lt. 1e-4_dp)

    call datasrc_update(ds, 150.0_dp)
    call datasrc_get(ds, "dT", a)
    call check("slice: the change since ref_time", maxval(abs(a - 5.0_dp)) .lt. 1e-4_dp)
    call datasrc_get(ds, "dTg", a)
    call check("group at another time", maxval(abs(a - 15.0_dp)) .lt. 1e-4_dp)
    call datasrc_end(ds)

    call datasrc_init(ds, fpar, "anom_bad", st, ["dT"], grid, fdef, stat=stat)
    call check("ref_time outside the axis: stat /= 0", stat .ne. 0)

    if (nfail .gt. 0) then
        write(*,*) "test_anomaly: ", nfail, " FAILED"
        error stop 1
    end if
    write(*,*) "test_anomaly: all passed"

contains

    subroutine write_nml()
        integer :: u
        open(newunit=u, file=fpar, status="replace", action="write")
        write(u,"(a)") "&anom"
        write(u,"(a)") '    dT  = "t_slice"'
        write(u,"(a)") '    dTg = "t_group"'
        write(u,"(a)") "/"
        write(u,"(a)") "&anom_bad"
        write(u,"(a)") '    dT  = "t_slice_out"'
        write(u,"(a)") "/"
        write(u,"(a)") "&t_slice"
        write(u,"(a)") '    filename    = "test_ds_anom.nc"'
        write(u,"(a)") '    name        = "T"'
        write(u,"(a)") '    anomaly_ref = "slice"'
        write(u,"(a)") '    ref_time    = 100.0'
        write(u,"(a)") "/"
        write(u,"(a)") "&t_slice_out"
        write(u,"(a)") '    filename    = "test_ds_anom.nc"'
        write(u,"(a)") '    name        = "T"'
        write(u,"(a)") '    anomaly_ref = "slice"'
        write(u,"(a)") '    ref_time    = 500.0'
        write(u,"(a)") "/"
        write(u,"(a)") "&t_group"
        write(u,"(a)") '    filename    = "test_ds_anom.nc"'
        write(u,"(a)") '    name        = "T"'
        write(u,"(a)") '    anomaly_ref = "group"'
        write(u,"(a)") '    ref_group   = "t_ref"'
        write(u,"(a)") "/"
        write(u,"(a)") "&t_ref"
        write(u,"(a)") '    filename = "test_ds_anom.nc"'
        write(u,"(a)") '    name     = "T_ref"'
        write(u,"(a)") "/"
        close(u)
    end subroutine write_nml

    subroutine check(label, ok)
        character(len=*), intent(in) :: label
        logical,          intent(in) :: ok
        write(*,"(a,a)") merge("  ok   ", "  FAIL ", ok), label
        if (.not. ok) nfail = nfail + 1
    end subroutine check

end program test_anomaly
