program test_examples
    ! The worked examples of the fesm_datasrc design note (§8), on synthetic
    ! files, with the field groups in a library file (var_defs):
    !   (a) background ice on the globe from an ICE-6G-like history: an
    !       irregular ka BP axis read from the file, time_ref, time_span,
    !       f_grnd filled with a constant
    !   (b) static geothermal: mW m-2 -> W m-2, fill nearest, {domain}; and
    !       the constant alternative
    !   (c) a scalar bsl curve from an ascii file in clock years
    !   (d) transient geo as anomalies of a relative sea-level history:
    !       dz_bed = -(rsl(t) - rsl(0)); dz_ss and bsl from one const group
    !       (2-D and 0-D)

    use precision
    use ncio
    use constants,    only : mv
    use coordinates,  only : grid_class, grid_init
    use fesm_fields
    use fesm_datasrc

    implicit none

    character(len=*), parameter :: fpar = "test_ds_ex_par.nml", flib = "test_ds_ex_lib.nml", &
                                   fdef = "par/datasrc_defaults.nml"
    integer, parameter :: gx = 6, gy = 4, rx = 4, ry = 3
    real(dp), parameter :: t_ref = 1950.0_dp
    type(grid_class)    :: globe, grl
    type(field_state_t) :: ice, geo, geo_c, forcing, grl_geo
    type(datasrc_t)     :: ds
    character(len=64)   :: var_defs(1)
    real(dp) :: tk(5), h2(gx,gy), r2(rx,ry), x, t0, t1
    real(wp) :: hv(gx,gy,5), fg(gx,gy,5), rsl(rx,ry,5), ghf(rx,ry)
    integer  :: h, k, stat, nfail

    nfail = 0
    var_defs(1) = flib

    call grid_init(globe, name="GLOBAL-TEST", mtype="latlon", units="degrees", &
                   x0=-150.0_dp, dx=60.0_dp, nx=gx, y0=-67.5_dp, dy=45.0_dp, ny=gy)
    call grid_init(grl, name="GRL-TEST", mtype="cartesian", units="kilometers", &
                   x0=0.0_dp, dx=100.0_dp, nx=rx, y0=0.0_dp, dy=100.0_dp, ny=ry)

    call write_files()
    call write_nml()

    ! === (a) background ice on the globe =========================================
    call fields_init(ice, "ice_bg:exports")
    call fields_add(ice, "H_ice",  "m", FIELD_INSTANT, FIELD_EXTENSIVE, h, grid="GLOBAL-TEST", dims=[gx,gy])
    call fields_add(ice, "z_srf",  "m", FIELD_INSTANT, FIELD_INTENSIVE, h, grid="GLOBAL-TEST", dims=[gx,gy])
    call fields_add(ice, "f_grnd", "1", FIELD_INSTANT, FIELD_INTENSIVE, h, grid="GLOBAL-TEST", dims=[gx,gy])
    call fields_lock(ice)

    call datasrc_init(ds, fpar, "ice_bg", ice, [character(len=8) :: "H_ice", "z_srf", "f_grnd"], globe, fdef, &
                      required=[character(len=8) :: "H_ice", "z_srf", "f_grnd"], var_defs=var_defs, &
                      time_ref=t_ref, time_span=[-3000.0_dp, 0.0_dp], stat=stat)
    call check("(a) init", stat .eq. 0)
    call datasrc_describe(ds)
    call datasrc_update(ds, -1500.0_dp, -750.0_dp)     ! geometry is INSTANT: the values at t1
    call datasrc_fill(ds, ice)
    call fields_get(ice, fields_handle(ice, "H_ice"), h2, t0, t1)
    call check("(a) H_ice at t1 on the irregular axis (500 yr step)", &
               abs(h2(2,3) - hval(2, 3, -750.0_dp)) .lt. 1e-3_dp .and. t0 .eq. -750.0_dp .and. t1 .eq. -750.0_dp)
    call fields_get(ice, fields_handle(ice, "f_grnd"), h2)
    call check("(a) f_grnd: missing cells filled with 0", h2(1,1) .eq. 0.0_dp .and. abs(h2(3,2) - 1.0_dp) .lt. 1e-6_dp)
    call datasrc_end(ds)

    ! === (b) static geothermal ======================================================
    call fields_init(geo, "geothermal_grl:exports")
    call fields_add(geo, "Q_geo", "W m-2", FIELD_INSTANT, FIELD_INTENSIVE, h, grid="GRL-TEST", dims=[rx,ry])
    call fields_lock(geo)
    call datasrc_init(ds, fpar, "geothermal_grl", geo, ["Q_geo"], grl, fdef, domain="Greenland", &
                      var_defs=var_defs, stat=stat)
    call check("(b) init", stat .eq. 0)
    call datasrc_describe(ds)
    call datasrc_update(ds, 0.0_dp)            ! nothing read; validity [t,t]
    call datasrc_fill(ds, geo)
    call fields_get(geo, fields_handle(geo, "Q_geo"), r2, t0, t1)
    call check("(b) mW m-2 -> W m-2, the missing cell filled", abs(r2(2,2) - 0.060_dp) .lt. 1e-7_dp .and. &
               r2(1,1) .gt. 0.0_dp .and. t0 .eq. 0.0_dp .and. t1 .eq. 0.0_dp)
    call datasrc_end(ds)

    call fields_init(geo_c, "geothermal_c:exports")
    call fields_add(geo_c, "Q_geo", "W m-2", FIELD_INSTANT, FIELD_INTENSIVE, h, grid="GRL-TEST", dims=[rx,ry])
    call fields_lock(geo_c)
    call datasrc_init(ds, fpar, "geothermal_const", geo_c, ["Q_geo"], grl, fdef, var_defs=var_defs, stat=stat)
    call datasrc_update(ds, 0.0_dp)
    call datasrc_get(ds, "Q_geo", r2)
    call check("(b) mode = const, value = 0.042", stat .eq. 0 .and. all(r2 .eq. 0.042_dp))
    call datasrc_end(ds)

    ! === (c) scalar bsl curve ========================================================
    call fields_init(forcing, "forcing:exports")
    call fields_add(forcing, "bsl", "m", FIELD_INSTANT, FIELD_INTENSIVE, h)
    call fields_lock(forcing)
    call datasrc_init(ds, fpar, "forcing", forcing, ["bsl"], globe, fdef, var_defs=var_defs, &
                      time_ref=t_ref, time_span=[-3000.0_dp, 0.0_dp], stat=stat)
    call check("(c) init", stat .eq. 0)
    call datasrc_update(ds, -2500.0_dp)        ! INSTANT: series value at t
    call datasrc_fill(ds, forcing)             ! put(h_bsl, value, t1=t)
    call fields_get(forcing, fields_handle(forcing, "bsl"), x, t0, t1)
    call check("(c) bsl at t", abs(x - (-75.0_dp)) .lt. 1e-9_dp .and. t1 .eq. -2500.0_dp)
    call datasrc_end(ds)

    ! === (d) geo as anomalies of a relative sea-level history =========================
    call fields_init(grl_geo, "geo_grl:exports")
    call fields_add(grl_geo, "dz_bed", "m", FIELD_INSTANT, FIELD_INTENSIVE, h, grid="GRL-TEST", dims=[rx,ry])
    call fields_add(grl_geo, "dz_ss",  "m", FIELD_INSTANT, FIELD_INTENSIVE, h, grid="GRL-TEST", dims=[rx,ry])
    call fields_add(grl_geo, "bsl",    "m", FIELD_INSTANT, FIELD_INTENSIVE, h)
    call fields_lock(grl_geo)
    call datasrc_init(ds, fpar, "geo_grl", grl_geo, [character(len=8) :: "dz_bed", "dz_ss", "bsl"], grl, fdef, &
                      domain="Greenland", var_defs=var_defs, time_ref=t_ref, time_span=[-3000.0_dp, 0.0_dp], &
                      stat=stat)
    call check("(d) init", stat .eq. 0)
    call datasrc_describe(ds)
    call datasrc_update(ds, 0.0_dp)
    call datasrc_get(ds, "dz_bed", r2)
    call check("(d) dz_bed = 0 at the reference", all(r2 .eq. 0.0_dp))
    call datasrc_update(ds, -2000.0_dp)
    call datasrc_fill(ds, grl_geo)
    call fields_get(grl_geo, fields_handle(grl_geo, "dz_bed"), r2)
    call check("(d) dz_bed = -(rsl(t) - rsl(0))", abs(r2(3,2) - (-(rval(3, 2, -2000.0_dp) - rval(3, 2, 0.0_dp)))) &
               .lt. 1e-3_dp)
    call fields_get(grl_geo, fields_handle(grl_geo, "dz_ss"), r2)
    call fields_get(grl_geo, fields_handle(grl_geo, "bsl"), x)
    call check("(d) one const group: 2-D dz_ss and 0-D bsl = 0", all(r2 .eq. 0.0_dp) .and. x .eq. 0.0_dp)
    call datasrc_end(ds)

    if (nfail .gt. 0) then
        write(*,*) "test_examples: ", nfail, " FAILED"
        error stop 1
    end if
    write(*,*) "test_examples: all passed"

contains

    real(dp) function hval(i, j, t)
        ! Ice thickness at clock time t (linear in t)
        integer,  intent(in) :: i, j
        real(dp), intent(in) :: t
        hval = 100.0_dp*i + 10.0_dp*j - 0.5_dp*t
    end function hval

    real(dp) function rval(i, j, t)
        ! Relative sea level at clock time t: falls as the bed rebounds
        integer,  intent(in) :: i, j
        real(dp), intent(in) :: t
        rval = -0.02_dp*t + i + 0.1_dp*j
    end function rval

    subroutine write_files()
        integer :: i, j, u

        tk = [3.0_dp, 2.0_dp, 1.0_dp, 0.5_dp, 0.0_dp]        ! ka BP: 1 kyr, then 500 yr
        do k = 1, 5
            do j = 1, gy; do i = 1, gx
                hv(i,j,k) = real(hval(i, j, -1000.0_dp*tk(k)), wp)
                fg(i,j,k) = 1.0_wp
            end do; end do
            fg(1,1,k) = mv
            do j = 1, ry; do i = 1, rx
                rsl(i,j,k) = real(rval(i, j, -1000.0_dp*tk(k)), wp)
            end do; end do
        end do

        call nc_create("test_ds_ex_GLOBAL-TEST_IceHistory.nc")
        call nc_write_dim("test_ds_ex_GLOBAL-TEST_IceHistory.nc", "lon", x=-150.0_dp, dx=60.0_dp, nx=gx)
        call nc_write_dim("test_ds_ex_GLOBAL-TEST_IceHistory.nc", "lat", x=-67.5_dp, dx=45.0_dp, nx=gy)
        call nc_write_dim("test_ds_ex_GLOBAL-TEST_IceHistory.nc", "time", x=tk, units="ka BP", unlimited=.false.)
        call nc_write("test_ds_ex_GLOBAL-TEST_IceHistory.nc", "H_ice", hv, dim1="lon", dim2="lat", dim3="time", &
                      units="m")
        call nc_write("test_ds_ex_GLOBAL-TEST_IceHistory.nc", "z_srf", hv, dim1="lon", dim2="lat", dim3="time", &
                      units="m")
        call nc_write("test_ds_ex_GLOBAL-TEST_IceHistory.nc", "f_grnd", fg, dim1="lon", dim2="lat", dim3="time", &
                      units="1", missing_value=mv)
        call nc_write_attr("test_ds_ex_GLOBAL-TEST_IceHistory.nc", "grid_name", "GLOBAL-TEST")

        ghf = 60.0_wp
        ghf(1,1) = mv
        call nc_create("test_ds_ex_Greenland_GRL-TEST_GHF.nc")
        call nc_write_dim("test_ds_ex_Greenland_GRL-TEST_GHF.nc", "xc", x=0.0_dp, dx=100.0_dp, nx=rx)
        call nc_write_dim("test_ds_ex_Greenland_GRL-TEST_GHF.nc", "yc", x=0.0_dp, dx=100.0_dp, nx=ry)
        call nc_write("test_ds_ex_Greenland_GRL-TEST_GHF.nc", "ghf", ghf, dim1="xc", dim2="yc", units="mW m-2", &
                      missing_value=mv)
        call nc_write_attr("test_ds_ex_Greenland_GRL-TEST_GHF.nc", "grid_name", "GRL-TEST")

        call nc_create("test_ds_ex_Greenland_GRL-TEST_RSL.nc")
        call nc_write_dim("test_ds_ex_Greenland_GRL-TEST_RSL.nc", "xc", x=0.0_dp, dx=100.0_dp, nx=rx)
        call nc_write_dim("test_ds_ex_Greenland_GRL-TEST_RSL.nc", "yc", x=0.0_dp, dx=100.0_dp, nx=ry)
        call nc_write_dim("test_ds_ex_Greenland_GRL-TEST_RSL.nc", "time", x=tk, units="ka BP", unlimited=.false.)
        call nc_write("test_ds_ex_Greenland_GRL-TEST_RSL.nc", "rsl", rsl, dim1="xc", dim2="yc", dim3="time", units="m")

        open(newunit=u, file="test_ds_ex_sealevel.dat", status="replace", action="write")
        write(u,"(a)") "# time [forcing yr]   bsl [m]"
        write(u,"(a)") "-3000   -90.0"
        write(u,"(a)") "-2000   -60.0"
        write(u,"(a)") "    0     0.0"
        close(u)
    end subroutine write_files

    subroutine write_nml()
        integer :: u
        open(newunit=u, file=fpar, status="replace", action="write")
        write(u,"(a)") "&ice_bg                     ! data ice on the globe"
        write(u,"(a)") '    backend = "data"'
        write(u,"(a)") '    grid    = "GLOBAL-TEST"'
        write(u,"(a)") '    dt      = 100.0'
        write(u,"(a)") '    clock   = "forcing"'
        write(u,"(a)") '    H_ice   = "ice6g_c_H_ice"'
        write(u,"(a)") '    z_srf   = "ice6g_c_z_srf"'
        write(u,"(a)") '    f_grnd  = "ice6g_c_f_grnd"'
        write(u,"(a)") "/"
        write(u,"(a)") "&geothermal_grl"
        write(u,"(a)") '    Q_geo = "q_geo_shapiro"'
        write(u,"(a)") "/"
        write(u,"(a)") "&geothermal_const"
        write(u,"(a)") '    Q_geo = "q_geo_const"'
        write(u,"(a)") "/"
        write(u,"(a)") "&forcing"
        write(u,"(a)") '    bsl = "bsl_spratt2016"'
        write(u,"(a)") "/"
        write(u,"(a)") "&geo_grl"
        write(u,"(a)") '    backend = "data"'
        write(u,"(a)") '    clock   = "forcing"'
        write(u,"(a)") '    dz_bed  = "ice6g_c_rsl_dz_bed"'
        write(u,"(a)") '    dz_ss   = "zero"'
        write(u,"(a)") '    bsl     = "zero"'
        write(u,"(a)") "/"
        close(u)

        open(newunit=u, file=flib, status="replace", action="write")
        write(u,"(a)") "&ice6g_c_H_ice"
        write(u,"(a)") '    filename = "test_ds_ex_{grid_name}_IceHistory.nc"'
        write(u,"(a)") '    name     = "H_ice"'
        write(u,"(a)") "/"
        write(u,"(a)") "&ice6g_c_z_srf"
        write(u,"(a)") '    filename = "test_ds_ex_{grid_name}_IceHistory.nc"'
        write(u,"(a)") '    name     = "z_srf"'
        write(u,"(a)") "/"
        write(u,"(a)") "&ice6g_c_f_grnd"
        write(u,"(a)") '    filename   = "test_ds_ex_{grid_name}_IceHistory.nc"'
        write(u,"(a)") '    name       = "f_grnd"'
        write(u,"(a)") '    fill       = "const"'
        write(u,"(a)") '    fill_value = 0.0'
        write(u,"(a)") "/"
        write(u,"(a)") "&q_geo_shapiro"
        write(u,"(a)") '    filename = "test_ds_ex_{domain}_{grid_name}_GHF.nc"'
        write(u,"(a)") '    name     = "ghf"'
        write(u,"(a)") '    fill     = "nearest"    ! file units "mW m-2" -> W m-2 (x 1e-3)'
        write(u,"(a)") "/"
        write(u,"(a)") "&q_geo_const"
        write(u,"(a)") '    mode  = "const"'
        write(u,"(a)") '    value = 0.042'
        write(u,"(a)") "/"
        write(u,"(a)") "&bsl_spratt2016"
        write(u,"(a)") '    filename   = "test_ds_ex_sealevel.dat"   ! ascii: time value'
        write(u,"(a)") '    units      = "m" "m"                      ! ascii has no attribute'
        write(u,"(a)") '    time_units = "clock"                      ! first column is forcing time'
        write(u,"(a)") "/"
        write(u,"(a)") "&ice6g_c_rsl_dz_bed"
        write(u,"(a)") '    filename    = "test_ds_ex_{domain}_{grid_name}_RSL.nc"'
        write(u,"(a)") '    name        = "rsl"'
        write(u,"(a)") '    scaling     = -1.0 0.0'
        write(u,"(a)") '    anomaly_ref = "slice"'
        write(u,"(a)") '    ref_time    = 0.0'
        write(u,"(a)") "/"
        write(u,"(a)") "&zero"
        write(u,"(a)") '    mode  = "const"'
        write(u,"(a)") '    value = 0.0'
        write(u,"(a)") "/"
        close(u)
    end subroutine write_nml

    subroutine check(label, ok)
        character(len=*), intent(in) :: label
        logical,          intent(in) :: ok
        write(*,"(a,a)") merge("  ok   ", "  FAIL ", ok), label
        if (.not. ok) nfail = nfail + 1
    end subroutine check

end program test_examples
