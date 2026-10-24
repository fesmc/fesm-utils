program datasrc_abort_mode
    ! Expected abort (make check-datasrc): a varslice group whose mode
    ! contradicts its legacy time key (static data declared transient).

    use, intrinsic :: iso_fortran_env, only : output_unit
    use precision, only : wp, dp
    use ncio
    use varslice, only : varslice_class, varslice_init_nml

    implicit none

    type(varslice_class) :: vs
    integer :: u

    write(output_unit, "(a)") "EXPECT: mode conflicts with the time key."
    flush(output_unit)

    call nc_create("test_ds_abort_mode.nc")
    call nc_write_dim("test_ds_abort_mode.nc", "xc", x=1.0_dp, dx=1.0_dp, nx=3)
    call nc_write_dim("test_ds_abort_mode.nc", "yc", x=1.0_dp, dx=1.0_dp, nx=2)
    call nc_write("test_ds_abort_mode.nc", "probe", reshape([(real(u,wp), u=1,6)], [3,2]), dim1="xc", dim2="yc")

    open(newunit=u, file="test_ds_abort_mode.nml", status="replace", action="write")
    write(u,"(a)") "&probe"
    write(u,"(a)") '    filename = "test_ds_abort_mode.nc"'
    write(u,"(a)") '    name     = "probe"'
    write(u,"(a)") '    units    = "m" "m"'
    write(u,"(a)") '    scaling  = 1.0 0.0'
    write(u,"(a)") '    time     = 0 0 0 0 1'
    write(u,"(a)") '    mode     = "transient"'
    write(u,"(a)") "/"
    close(u)

    call varslice_init_nml(vs, "test_ds_abort_mode.nml", "probe", verbose=.false.)

    write(output_unit, "(a)") "FAIL: datasrc_abort_mode did not abort"

end program datasrc_abort_mode
