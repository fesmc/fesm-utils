program datasrc_abort_time_order
    ! Expected abort (make check-datasrc): a time axis from the file that is
    ! not strictly increasing (two files in the wrong order).

    use, intrinsic :: iso_fortran_env, only : output_unit
    use precision, only : wp, dp
    use ncio
    use varslice, only : varslice_class, varslice_init_nml

    implicit none

    type(varslice_class) :: vs
    integer :: u

    write(output_unit, "(a)") "EXPECT: the time axis is not strictly increasing."
    flush(output_unit)

    call ncfile("test_ds_abort_order_1.nc", [10.0_dp, 11.0_dp])
    call ncfile("test_ds_abort_order_2.nc", [ 5.0_dp,  6.0_dp])

    open(newunit=u, file="test_ds_abort_order.nml", status="replace", action="write")
    write(u,"(a)") "&probe"
    write(u,"(a)") '    filename = "test_ds_abort_order_*.nc"'
    write(u,"(a)") '    name     = "probe"'
    write(u,"(a)") '    units    = "m" "m"'
    write(u,"(a)") '    scaling  = 1.0 0.0'
    write(u,"(a)") '    time     = -1 0 0 0 0'
    write(u,"(a)") "/"
    close(u)

    call varslice_init_nml(vs, "test_ds_abort_order.nml", "probe", verbose=.false.)

    write(output_unit, "(a)") "FAIL: datasrc_abort_time_order did not abort"

contains

    subroutine ncfile(fname, t)
        character(len=*), intent(in) :: fname
        real(dp),         intent(in) :: t(:)
        real(wp) :: v(2,2,size(t))
        v = 1.0_wp
        call nc_create(fname)
        call nc_write_dim(fname, "xc", x=1.0_dp, dx=1.0_dp, nx=2)
        call nc_write_dim(fname, "yc", x=1.0_dp, dx=1.0_dp, nx=2)
        call nc_write_dim(fname, "time", x=t, units="yr CE", unlimited=.false.)
        call nc_write(fname, "probe", v, dim1="xc", dim2="yc", dim3="time")
    end subroutine ncfile

end program datasrc_abort_time_order
