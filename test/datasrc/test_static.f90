program test_static
    ! fesm_datasrc, static fields: const, static remapped (remap = auto from
    ! the file's grid_name, mW m-2 -> W m-2), static on the instance grid
    ! with fill, a masked field with missing values, validity [t,t]; and
    ! the init errors (stat): shape, grid, units attribute, conversion.

    use precision
    use ncio
    use constants,    only : mv
    use coordinates,  only : grid_class, grid_init
    use fesm_fields
    use fesm_datasrc

    implicit none

    character(len=*), parameter :: fpar = "test_ds_static.nml", fdef = "par/datasrc_defaults.nml"
    integer, parameter :: nx = 4, ny = 3, ns = 6, ms = 5
    real(dp), parameter :: a = 40.0_dp, b = 0.5_dp, c = -0.2_dp
    type(grid_class)    :: grid
    type(field_state_t) :: st
    type(datasrc_t)     :: ds
    real(dp) :: v2(nx,ny), q, t0, t1
    logical  :: valid(nx,ny)
    real(wp) :: src(ns,ms), loc(nx,ny)
    integer  :: h, i, j, stat, nfail

    nfail = 0

    call grid_init(grid, name="TEST-DS-TGT", mtype="latlon", units="degrees", &
                   x0=-7.5_dp, dx=5.0_dp, nx=nx, y0=62.5_dp, dy=5.0_dp, ny=ny)

    ! Source grid: 6 x 5 cells of 5 degrees, a linear field (bilinear exact)
    do j = 1, ms; do i = 1, ns
        src(i,j) = real(a + b*(-15.0_dp + 5.0_dp*(i-1)) + c*(55.0_dp + 5.0_dp*(j-1)), wp)
    end do; end do
    call nc_create("test_ds_static_src.nc")
    call nc_write_dim("test_ds_static_src.nc", "lon", x=-15.0_dp, dx=5.0_dp, nx=ns)
    call nc_write_dim("test_ds_static_src.nc", "lat", x=55.0_dp, dx=5.0_dp, nx=ms)
    call nc_write("test_ds_static_src.nc", "ghf", src, dim1="lon", dim2="lat", units="mW m-2")
    call nc_write_attr("test_ds_static_src.nc", "grid_name", "TEST-DS-SRC")

    ! On the instance grid, two missing cells
    do j = 1, ny; do i = 1, nx
        loc(i,j) = real(10*i + j, wp)
    end do; end do
    loc(2,2) = mv
    loc(4,3) = mv
    call nc_create("test_ds_static_loc.nc")
    call nc_write_dim("test_ds_static_loc.nc", "lon", x=-7.5_dp, dx=5.0_dp, nx=nx)
    call nc_write_dim("test_ds_static_loc.nc", "lat", x=62.5_dp, dx=5.0_dp, nx=ny)
    call nc_write("test_ds_static_loc.nc", "h_sed", loc, dim1="lon", dim2="lat", units="m", missing_value=mv)
    call nc_write("test_ds_static_loc.nc", "t_bed", loc, dim1="lon", dim2="lat", units="K", missing_value=mv)
    call nc_write_attr("test_ds_static_loc.nc", "grid_name", "TEST-DS-TGT")

    call write_nml()

    ! --- the state ----------------------------------------------------------------
    call fields_init(st, "geo_test:exports")
    call fields_add(st, "Q_geo", "W m-2", FIELD_INSTANT, FIELD_INTENSIVE, h, grid="TEST-DS-TGT", dims=[nx,ny])
    call fields_add(st, "Q_c",   "W m-2", FIELD_INSTANT, FIELD_INTENSIVE, h, grid="TEST-DS-TGT", dims=[nx,ny])
    call fields_add(st, "H_sed", "m",     FIELD_INSTANT, FIELD_EXTENSIVE, h, grid="TEST-DS-TGT", dims=[nx,ny])
    call fields_add(st, "T_bed", "K",     FIELD_INSTANT, FIELD_INTENSIVE, h, grid="TEST-DS-TGT", dims=[nx,ny], &
                    masked=.true.)
    call fields_add(st, "q0",    "W m-2", FIELD_INSTANT, FIELD_INTENSIVE, h)
    call fields_add(st, "V3",    "m",     FIELD_INSTANT, FIELD_INTENSIVE, h, grid="TEST-DS-TGT", dims=[nx,ny,2])
    call fields_add(st, "W_oth", "m",     FIELD_INSTANT, FIELD_INTENSIVE, h, grid="OTHER-GRID", dims=[nx,ny])
    call fields_lock(st)

    call datasrc_init(ds, fpar, "geo", st, [character(len=8) :: "Q_geo", "Q_c", "H_sed", "T_bed", "q0"], &
                      grid, fdef, stat=stat)
    call check("init", stat .eq. 0)
    call datasrc_describe(ds)
    call datasrc_update(ds, 100.0_dp)
    call datasrc_fill(ds, st)

    ! const
    call fields_get(st, fields_handle(st, "Q_c"), v2, t0, t1)
    call check("const: value", all(v2 .eq. 0.042_dp))
    call check("const: validity [t,t]", t0 .eq. 100.0_dp .and. t1 .eq. 100.0_dp)
    call fields_get(st, fields_handle(st, "q0"), q)
    call check("const scalar from the same group", q .eq. 0.042_dp)

    ! static remapped from grid_name, converted
    call fields_get(st, fields_handle(st, "Q_geo"), v2)
    call check("static remapped (bilin) and converted to W m-2", &
               abs(v2(2,2) - 1e-3_dp*(a + b*(-2.5_dp) + c*67.5_dp)) .lt. 1e-7_dp .and. &
               abs(v2(4,1) - 1e-3_dp*(a + b*7.5_dp + c*62.5_dp)) .lt. 1e-7_dp)

    ! static on the grid, fill = nearest
    call fields_get(st, fields_handle(st, "H_sed"), v2)
    call check("static on the grid, filled", all(v2 .gt. 0.0_dp) .and. v2(1,1) .eq. 11.0_dp)

    ! masked field: missing cells invalid
    call fields_get(st, fields_handle(st, "T_bed"), v2, valid=valid)
    call check("masked: missing cells invalid", .not. valid(2,2) .and. .not. valid(4,3) .and. &
               count(valid) .eq. nx*ny-2 .and. v2(3,2) .eq. 32.0_dp)

    ! update(t) again: static fields keep their value, validity moves
    call datasrc_update(ds, 200.0_dp)
    call datasrc_fill(ds, st)
    call fields_get(st, fields_handle(st, "H_sed"), v2, t0, t1)
    call check("static: validity follows the update", t0 .eq. 200.0_dp .and. t1 .eq. 200.0_dp .and. &
               v2(3,1) .eq. 31.0_dp)
    call check("datasrc_has", datasrc_has(ds, "Q_geo") .and. .not. datasrc_has(ds, "V3"))
    call datasrc_end(ds)

    ! --- init errors with stat ----------------------------------------------------
    call datasrc_init(ds, fpar, "bad_shape", st, ["V3"], grid, fdef, stat=stat)
    call check("shape mismatch: stat /= 0", stat .ne. 0)
    call datasrc_init(ds, fpar, "bad_grid", st, ["W_oth"], grid, fdef, stat=stat)
    call check("field on another grid: stat /= 0", stat .ne. 0)
    call datasrc_init(ds, fpar, "bad_units", st, ["Q_geo"], grid, fdef, stat=stat)
    call check("units(1) differs from the attribute: stat /= 0", stat .ne. 0)
    call datasrc_init(ds, fpar, "bad_conv", st, ["Q_geo"], grid, fdef, stat=stat)
    call check("no conversion K -> W m-2: stat /= 0", stat .ne. 0)

    if (nfail .gt. 0) then
        write(*,*) "test_static: ", nfail, " FAILED"
        error stop 1
    end if
    write(*,*) "test_static: all passed"

contains

    subroutine write_nml()
        integer :: u
        open(newunit=u, file=fpar, status="replace", action="write")
        write(u,"(a)") "&geo"
        write(u,"(a)") '    Q_geo = "ghf_test"'
        write(u,"(a)") '    Q_c   = "ghf_const"'
        write(u,"(a)") '    H_sed = "sed_test"'
        write(u,"(a)") '    T_bed = "tbed_test"'
        write(u,"(a)") '    q0    = "ghf_const"'
        write(u,"(a)") "/"
        write(u,"(a)") "&bad_shape"
        write(u,"(a)") '    V3 = "sed_test"'
        write(u,"(a)") "/"
        write(u,"(a)") "&bad_grid"
        write(u,"(a)") '    W_oth = "sed_test"'
        write(u,"(a)") "/"
        write(u,"(a)") "&bad_units"
        write(u,"(a)") '    Q_geo = "ghf_wrong_units"'
        write(u,"(a)") "/"
        write(u,"(a)") "&bad_conv"
        write(u,"(a)") '    Q_geo = "tbed_test"'
        write(u,"(a)") "/"
        write(u,"(a)") "&ghf_test"
        write(u,"(a)") '    filename = "test_ds_static_src.nc"'
        write(u,"(a)") '    name     = "ghf"'
        write(u,"(a)") "/"
        write(u,"(a)") "&ghf_wrong_units"
        write(u,"(a)") '    filename = "test_ds_static_src.nc"'
        write(u,"(a)") '    name     = "ghf"'
        write(u,"(a)") '    units    = "W m-2" "auto"'
        write(u,"(a)") "/"
        write(u,"(a)") "&ghf_const"
        write(u,"(a)") '    mode  = "const"'
        write(u,"(a)") '    value = 0.042'
        write(u,"(a)") "/"
        write(u,"(a)") "&sed_test"
        write(u,"(a)") '    filename = "test_ds_static_loc.nc"'
        write(u,"(a)") '    name     = "h_sed"'
        write(u,"(a)") '    fill     = "nearest"'
        write(u,"(a)") "/"
        write(u,"(a)") "&tbed_test"
        write(u,"(a)") '    filename = "test_ds_static_loc.nc"'
        write(u,"(a)") '    name     = "t_bed"'
        write(u,"(a)") "/"
        close(u)
    end subroutine write_nml

    subroutine check(label, ok)
        character(len=*), intent(in) :: label
        logical,          intent(in) :: ok
        write(*,"(a,a)") merge("  ok   ", "  FAIL ", ok), label
        if (.not. ok) nfail = nfail + 1
    end subroutine check

end program test_static
