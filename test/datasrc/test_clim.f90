program test_clim
    ! fesm_datasrc, climatologies: a 12-month 2-D climatology served as
    ! (nx,ny,12), its annual mean (time_method = "mean") as a rank-2 field,
    ! a rank-4 climatology (x,y,depth,month) only as its mean (rank 3), a
    ! clim_db-style group (time = 0 0 0 0 12, degC -> K by scaling) accepted
    ! unchanged, and the cycle of a transient monthly axis for the year of t.

    use precision
    use ncio
    use coordinates,  only : grid_class, grid_init
    use fesm_fields
    use fesm_datasrc

    implicit none

    character(len=*), parameter :: fpar = "test_ds_clim.nml", fdef = "par/datasrc_defaults.nml"
    integer, parameter :: nx = 3, ny = 2, nm = 12, nz = 2
    type(grid_class)    :: grid
    type(field_state_t) :: st
    type(datasrc_t)     :: ds
    real(wp) :: c2(nx,ny,nm), c3(nx,ny,nz,nm), mon(nx,ny,24)
    real(dp) :: v2(nx,ny), v3(nx,ny,nm), w3(nx,ny,nz), tm(24)
    integer  :: h, m, k, stat, nfail

    nfail = 0

    call grid_init(grid, name="TEST-DS-CLIM", mtype="cartesian", units="kilometers", &
                   x0=0.0_dp, dx=10.0_dp, nx=nx, y0=0.0_dp, dy=10.0_dp, ny=ny)

    do m = 1, nm
        c2(:,:,m)   = real(m, wp)                 ! mean 6.5
        c3(:,:,1,m) = real(m, wp)
        c3(:,:,2,m) = real(m, wp) + 10.0_wp       ! mean 16.5
    end do
    call nc_create("test_ds_clim.nc")
    call nc_write_dim("test_ds_clim.nc", "xc", x=0.0_dp, dx=10.0_dp, nx=nx)
    call nc_write_dim("test_ds_clim.nc", "yc", x=0.0_dp, dx=10.0_dp, nx=ny)
    call nc_write_dim("test_ds_clim.nc", "depth", x=[5.0_dp, 50.0_dp])
    call nc_write_dim("test_ds_clim.nc", "month", x=1.0_dp, dx=1.0_dp, nx=nm)
    call nc_write("test_ds_clim.nc", "t2m", c2, dim1="xc", dim2="yc", dim3="month", units="degC")
    call nc_write("test_ds_clim.nc", "to",  c3, dim1="xc", dim2="yc", dim3="depth", dim4="month", units="degC")

    ! Transient monthly, 2 years (2000, 2001), mid-month, noleap; value = year index*100 + month
    do k = 1, 24
        tm(k) = 365.0_dp*((k-1)/12) + 15.0_dp + 30.0_dp*mod(k-1,12)
        mon(:,:,k) = real(100*((k-1)/12 + 1) + mod(k-1,12) + 1, wp)
    end do
    call nc_create("test_ds_clim_mon.nc")
    call nc_write_dim("test_ds_clim_mon.nc", "xc", x=0.0_dp, dx=10.0_dp, nx=nx)
    call nc_write_dim("test_ds_clim_mon.nc", "yc", x=0.0_dp, dx=10.0_dp, nx=ny)
    call nc_write_dim("test_ds_clim_mon.nc", "time", x=tm, units="days since 2000-01-01", calendar="noleap", &
                      unlimited=.false.)
    call nc_write("test_ds_clim_mon.nc", "tas", mon, dim1="xc", dim2="yc", dim3="time", units="K")

    call write_nml()

    call fields_init(st, "clim_test:exports")
    call fields_add(st, "tas_m",  "K", FIELD_INSTANT, FIELD_INTENSIVE, h, grid="TEST-DS-CLIM", dims=[nx,ny,nm], &
                    zdim="month")
    call fields_add(st, "tas_a",  "K", FIELD_INSTANT, FIELD_INTENSIVE, h, grid="TEST-DS-CLIM", dims=[nx,ny])
    call fields_add(st, "to_a",   "K", FIELD_INSTANT, FIELD_INTENSIVE, h, grid="TEST-DS-CLIM", dims=[nx,ny,nz], &
                    zdim="depth")
    call fields_add(st, "tas_db", "K", FIELD_INSTANT, FIELD_INTENSIVE, h, grid="TEST-DS-CLIM", dims=[nx,ny,nm], &
                    zdim="month")
    call fields_add(st, "tas_tr", "K", FIELD_INSTANT, FIELD_INTENSIVE, h, grid="TEST-DS-CLIM", dims=[nx,ny,nm], &
                    zdim="month")
    call fields_lock(st)

    call datasrc_init(ds, fpar, "clim", st, [character(len=8) :: "tas_m", "tas_a", "to_a", "tas_db", "tas_tr"], &
                      grid, fdef, stat=stat)
    call check("init", stat .eq. 0)
    call datasrc_describe(ds)
    call datasrc_update(ds, 2001.4_dp)
    call datasrc_fill(ds, st)

    call fields_get(st, fields_handle(st, "tas_m"), v3)
    call check("12-month climatology served (nx,ny,12), degC -> K", &
               abs(v3(1,1,3) - (3.0_dp + 273.15_dp)) .lt. 1e-4_dp .and. abs(v3(3,2,12) - 285.15_dp) .lt. 1e-4_dp)
    call fields_get(st, fields_handle(st, "tas_a"), v2)
    call check("annual mean (time_method = mean)", maxval(abs(v2 - (6.5_dp + 273.15_dp))) .lt. 1e-4_dp)
    call fields_get(st, fields_handle(st, "to_a"), w3)
    call check("rank-4 climatology as its mean", abs(w3(2,1,1) - 279.65_dp) .lt. 1e-4_dp .and. &
               abs(w3(2,1,2) - 289.65_dp) .lt. 1e-4_dp)
    call fields_get(st, fields_handle(st, "tas_db"), v3)
    call check("clim_db group (time = 0 0 0 0 12, scaling 273.15) unchanged", &
               abs(v3(2,2,7) - (7.0_dp + 273.15_dp)) .lt. 1e-4_dp)
    call fields_get(st, fields_handle(st, "tas_tr"), v3)
    call check("transient monthly: the cycle of the year of t", abs(v3(1,1,5) - 205.0_dp) .lt. 1e-4_dp)
    call datasrc_end(ds)

    call datasrc_init(ds, fpar, "clim_bad", st, ["to_a"], grid, fdef, stat=stat)
    call check("rank-4 climatology without mean: stat /= 0", stat .ne. 0)
    call datasrc_init(ds, fpar, "clim_bad2", st, ["tas_a"], grid, fdef, stat=stat)
    call check("2-D climatology into a rank-2 field without mean: stat /= 0", stat .ne. 0)

    if (nfail .gt. 0) then
        write(*,*) "test_clim: ", nfail, " FAILED"
        error stop 1
    end if
    write(*,*) "test_clim: all passed"

contains

    subroutine write_nml()
        integer :: u
        open(newunit=u, file=fpar, status="replace", action="write")
        write(u,"(a)") "&clim"
        write(u,"(a)") '    tas_m  = "c_t2m"'
        write(u,"(a)") '    tas_a  = "c_t2m_mean"'
        write(u,"(a)") '    to_a   = "c_to_mean"'
        write(u,"(a)") '    tas_db = "c_t2m_db"'
        write(u,"(a)") '    tas_tr = "c_tas_monthly"'
        write(u,"(a)") "/"
        write(u,"(a)") "&clim_bad"
        write(u,"(a)") '    to_a = "c_to"'
        write(u,"(a)") "/"
        write(u,"(a)") "&clim_bad2"
        write(u,"(a)") '    tas_a = "c_t2m"'
        write(u,"(a)") "/"
        write(u,"(a)") "&c_t2m"
        write(u,"(a)") '    filename = "test_ds_clim.nc"'
        write(u,"(a)") '    name     = "t2m"'
        write(u,"(a)") '    mode     = "clim"'
        write(u,"(a)") "/"
        write(u,"(a)") "&c_t2m_mean"
        write(u,"(a)") '    filename    = "test_ds_clim.nc"'
        write(u,"(a)") '    name        = "t2m"'
        write(u,"(a)") '    mode        = "clim"'
        write(u,"(a)") '    time_method = "mean"'
        write(u,"(a)") "/"
        write(u,"(a)") "&c_to"
        write(u,"(a)") '    filename = "test_ds_clim.nc"'
        write(u,"(a)") '    name     = "to"'
        write(u,"(a)") '    mode     = "clim"'
        write(u,"(a)") "/"
        write(u,"(a)") "&c_to_mean"
        write(u,"(a)") '    filename    = "test_ds_clim.nc"'
        write(u,"(a)") '    name        = "to"'
        write(u,"(a)") '    mode        = "clim"'
        write(u,"(a)") '    time_method = "mean"'
        write(u,"(a)") "/"
        write(u,"(a)") "! a clim_db.nml group, unchanged"
        write(u,"(a)") "&c_t2m_db"
        write(u,"(a)") '    filename = "test_ds_clim.nc"'
        write(u,"(a)") '    name     = "t2m"'
        write(u,"(a)") '    units    = "degC" "K"'
        write(u,"(a)") '    scaling  = 1.0 273.15'
        write(u,"(a)") '    time     = 0 0 0 0 12'
        write(u,"(a)") "/"
        write(u,"(a)") "&c_tas_monthly"
        write(u,"(a)") '    filename = "test_ds_clim_mon.nc"'
        write(u,"(a)") '    name     = "tas"'
        write(u,"(a)") "/"
        close(u)
    end subroutine write_nml

    subroutine check(label, ok)
        character(len=*), intent(in) :: label
        logical,          intent(in) :: ok
        write(*,"(a,a)") merge("  ok   ", "  FAIL ", ok), label
        if (.not. ok) nfail = nfail + 1
    end subroutine check

end program test_clim
