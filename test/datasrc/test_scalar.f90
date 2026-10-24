program test_scalar
    ! fesm_datasrc, 0-D fields through series: an ascii series (units and
    ! time_units required) and a netCDF series (ka BP, decoded with
    ! time_ref), INSTANT interp, MEAN over an interval, fill into a
    ! field_state_t with validity [t,t] / [t0,t1].

    use precision
    use ncio
    use coordinates,  only : grid_class, grid_init
    use fesm_fields
    use fesm_datasrc

    implicit none

    character(len=*), parameter :: fpar = "test_ds_scalar.nml", fdef = "par/datasrc_defaults.nml"
    type(grid_class)    :: grid
    type(field_state_t) :: st
    type(datasrc_t)     :: ds
    real(dp) :: x, t0, t1
    integer  :: h, u, stat, nfail

    nfail = 0

    call grid_init(grid, name="TEST-DS-0D", mtype="cartesian", units="kilometers", &
                   x0=0.0_dp, dx=1.0_dp, nx=2, y0=0.0_dp, dy=1.0_dp, ny=2)

    ! ascii: time (clock years) and sea level [m]
    open(newunit=u, file="test_ds_bsl.dat", status="replace", action="write")
    write(u,"(a)") "# time  bsl"
    write(u,"(a)") "-20000  -120.0"
    write(u,"(a)") "-10000   -60.0"
    write(u,"(a)") "     0     0.0"
    close(u)

    ! netCDF: CO2 [ppm] over ka BP
    call nc_create("test_ds_co2.nc")
    call nc_write_dim("test_ds_co2.nc", "time", x=[30.0_dp, 20.0_dp, 10.0_dp, 0.0_dp], units="ka BP", &
                      unlimited=.false.)
    call nc_write("test_ds_co2.nc", "co2", [190.0_dp, 200.0_dp, 260.0_dp, 280.0_dp], dim1="time", units="1")

    call write_nml()

    call fields_init(st, "forcing_test:exports")
    call fields_add(st, "bsl",  "m", FIELD_INSTANT, FIELD_INTENSIVE, h)
    call fields_add(st, "co2",  "1", FIELD_INSTANT, FIELD_INTENSIVE, h)
    call fields_add(st, "bslm", "m", FIELD_MEAN,    FIELD_INTENSIVE, h)
    call fields_add(st, "bsl_cm", "1", FIELD_INSTANT, FIELD_INTENSIVE, h)
    call fields_lock(st)

    ! --- INSTANT ------------------------------------------------------------------
    call datasrc_init(ds, fpar, "forcing", st, [character(len=4) :: "bsl", "co2"], grid, fdef, &
                      time_ref=1950.0_dp, time_span=[-20000.0_dp, 0.0_dp], stat=stat)
    call check("init ascii + netCDF series", stat .eq. 0)
    call datasrc_describe(ds)
    call datasrc_update(ds, -15000.0_dp)
    call datasrc_get(ds, "bsl", x)
    call check("ascii: interp", abs(x - (-90.0_dp)) .lt. 1e-9_dp)
    call datasrc_get(ds, "co2", x)
    call check("netCDF ka BP with time_ref: interp", abs(x - 230.0_dp) .lt. 1e-9_dp)
    call datasrc_fill(ds, st)
    call fields_get(st, fields_handle(st, "bsl"), x, t0, t1)
    call check("fill: value and validity [t,t]", abs(x + 90.0_dp) .lt. 1e-9_dp .and. &
               t0 .eq. -15000.0_dp .and. t1 .eq. -15000.0_dp)
    call datasrc_end(ds)

    ! --- MEAN over an interval: the mean of the interpolant ----------------------
    call datasrc_init(ds, fpar, "forcing_mean", st, ["bslm"], grid, fdef, stat=stat)
    call check("init MEAN scalar", stat .eq. 0)
    call datasrc_update(ds, -15000.0_dp, -5000.0_dp)
    call datasrc_get(ds, "bslm", x)
    call check("MEAN over [t0,t1]", abs(x - (-60.0_dp)) .lt. 1e-9_dp)
    call datasrc_fill(ds, st)
    call fields_finalize_interval(st, -5000.0_dp)
    call fields_get(st, fields_handle(st, "bslm"), x, t0, t1)
    call check("MEAN fill: validity [t0,t1]", abs(x + 60.0_dp) .lt. 1e-9_dp .and. &
               t0 .eq. -15000.0_dp .and. t1 .eq. -5000.0_dp)
    call datasrc_end(ds)

    ! --- requirements -------------------------------------------------------------
    call datasrc_init(ds, fpar, "no_units", st, ["bsl"], grid, fdef, stat=stat)
    call check("ascii without units(1): stat /= 0", stat .ne. 0)
    call datasrc_init(ds, fpar, "no_tunits", st, ["bsl"], grid, fdef, stat=stat)
    call check("ascii without time_units: stat /= 0", stat .ne. 0)
    call datasrc_init(ds, fpar, "conv", st, ["bsl_cm"], grid, fdef, stat=stat)
    call check("no conversion m -> 1: stat /= 0", stat .ne. 0)

    if (nfail .gt. 0) then
        write(*,*) "test_scalar: ", nfail, " FAILED"
        error stop 1
    end if
    write(*,*) "test_scalar: all passed"

contains

    subroutine write_nml()
        integer :: u
        open(newunit=u, file=fpar, status="replace", action="write")
        write(u,"(a)") "&forcing"
        write(u,"(a)") '    bsl = "bsl_test"'
        write(u,"(a)") '    co2 = "co2_test"'
        write(u,"(a)") "/"
        write(u,"(a)") "&forcing_mean"
        write(u,"(a)") '    bslm = "bsl_test"'
        write(u,"(a)") "/"
        write(u,"(a)") "&no_units"
        write(u,"(a)") '    bsl = "bsl_no_units"'
        write(u,"(a)") "/"
        write(u,"(a)") "&no_tunits"
        write(u,"(a)") '    bsl = "bsl_no_tunits"'
        write(u,"(a)") "/"
        write(u,"(a)") "&conv"
        write(u,"(a)") '    bsl_cm = "bsl_test"'
        write(u,"(a)") "/"
        write(u,"(a)") "&bsl_test"
        write(u,"(a)") '    filename   = "test_ds_bsl.dat"'
        write(u,"(a)") '    units      = "m" "m"'
        write(u,"(a)") '    time_units = "clock"'
        write(u,"(a)") "/"
        write(u,"(a)") "&bsl_no_units"
        write(u,"(a)") '    filename   = "test_ds_bsl.dat"'
        write(u,"(a)") '    time_units = "clock"'
        write(u,"(a)") "/"
        write(u,"(a)") "&bsl_no_tunits"
        write(u,"(a)") '    filename   = "test_ds_bsl.dat"'
        write(u,"(a)") '    units      = "m" "m"'
        write(u,"(a)") "/"
        write(u,"(a)") "&co2_test"
        write(u,"(a)") '    filename = "test_ds_co2.nc"'
        write(u,"(a)") '    name     = "co2"'
        write(u,"(a)") "/"
        close(u)
    end subroutine write_nml

    subroutine check(label, ok)
        character(len=*), intent(in) :: label
        logical,          intent(in) :: ok
        write(*,"(a,a)") merge("  ok   ", "  FAIL ", ok), label
        if (.not. ok) nfail = nfail + 1
    end subroutine check

end program test_scalar
