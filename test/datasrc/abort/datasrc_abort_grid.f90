program datasrc_abort_grid
    ! Expected abort (make check-datasrc): a field group whose file is not
    ! on the instance grid and has no grid to remap from (no grid_name, no
    ! grid_src): varslice stops on the size check.

    use, intrinsic :: iso_fortran_env, only : output_unit
    use precision,   only : wp, dp
    use ncio
    use coordinates, only : grid_class, grid_init
    use fesm_fields
    use fesm_datasrc

    implicit none

    type(grid_class)    :: grid
    type(field_state_t) :: st
    type(datasrc_t)     :: ds
    real(wp) :: v(5,4)
    integer  :: h, u

    write(output_unit, "(a)") "EXPECT: the variable is not on the target grid (set remap and grid_src to remap it)."
    flush(output_unit)

    call grid_init(grid, name="TEST-DS-ABORT", mtype="cartesian", units="kilometers", &
                   x0=0.0_dp, dx=10.0_dp, nx=3, y0=0.0_dp, dy=10.0_dp, ny=2)
    v = 1.0_wp
    call nc_create("test_ds_abort_grid.nc")
    call nc_write_dim("test_ds_abort_grid.nc", "xc", x=0.0_dp, dx=10.0_dp, nx=5)
    call nc_write_dim("test_ds_abort_grid.nc", "yc", x=0.0_dp, dx=10.0_dp, nx=4)
    call nc_write("test_ds_abort_grid.nc", "probe", v, dim1="xc", dim2="yc", units="m")

    open(newunit=u, file="test_ds_abort_grid.nml", status="replace", action="write")
    write(u,"(a)") "&inst"
    write(u,"(a)") '    probe = "probe_group"'
    write(u,"(a)") "/"
    write(u,"(a)") "&probe_group"
    write(u,"(a)") '    filename = "test_ds_abort_grid.nc"'
    write(u,"(a)") '    name     = "probe"'
    write(u,"(a)") "/"
    close(u)

    call fields_init(st, "abort_test")
    call fields_add(st, "probe", "m", FIELD_INSTANT, FIELD_INTENSIVE, h, grid="TEST-DS-ABORT", dims=[3,2])
    call fields_lock(st)

    call datasrc_init(ds, "test_ds_abort_grid.nml", "inst", st, ["probe"], grid, "par/datasrc_defaults.nml")

    write(output_unit, "(a)") "FAIL: datasrc_abort_grid did not abort"

end program datasrc_abort_grid
