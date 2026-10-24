program test_varslice_units_attr
    ! par%units_file (varslice V3): the units attribute of the variable, ""
    ! without one, through both init paths.

    use precision
    use ncio
    use varslice, only : varslice_class, varslice_init_nml, varslice_init_arg, varslice_end

    implicit none

    type(varslice_class) :: vs
    integer :: nfail, u

    nfail = 0

    call nc_create("test_ds_uattr.nc")
    call nc_write_dim("test_ds_uattr.nc", "xc", x=1.0_dp, dx=1.0_dp, nx=3)
    call nc_write_dim("test_ds_uattr.nc", "yc", x=1.0_dp, dx=1.0_dp, nx=2)
    call nc_write("test_ds_uattr.nc", "ghf",  reshape([(real(u,wp), u=1,6)], [3,2]), dim1="xc", dim2="yc", &
                  units="mW m-2")
    call nc_write("test_ds_uattr.nc", "bare", reshape([(real(u,wp), u=1,6)], [3,2]), dim1="xc", dim2="yc")

    open(newunit=u, file="test_ds_uattr.nml", status="replace", action="write")
    write(u,"(a)") "&ghf"
    write(u,"(a)") '    filename = "test_ds_uattr.nc"'
    write(u,"(a)") '    name     = "ghf"'
    write(u,"(a)") '    units    = "W m-2" "W m-2"'
    write(u,"(a)") '    scaling  = 1e-3 0.0'
    write(u,"(a)") '    time     = 0 0 0 0 1'
    write(u,"(a)") "/"
    close(u)

    call varslice_init_nml(vs, "test_ds_uattr.nml", "ghf", verbose=.false.)
    call check("units_file from the attribute (nml)", trim(vs%par%units_file) .eq. "mW m-2")
    call check("group units kept as metadata", trim(vs%par%units_in) .eq. "W m-2")
    call varslice_end(vs)

    call varslice_init_arg(vs, "test_ds_uattr.nc", "bare", "1", "1", with_time=.false.)
    call check("units_file empty without an attribute (arg)", len_trim(vs%par%units_file) .eq. 0)
    call varslice_end(vs)

    if (nfail .gt. 0) then
        write(*,*) "test_varslice_units_attr: ", nfail, " FAILED"
        error stop 1
    end if
    write(*,*) "test_varslice_units_attr: all passed"

contains

    subroutine check(label, ok)
        character(len=*), intent(in) :: label
        logical,          intent(in) :: ok
        write(*,"(a,a)") merge("  ok   ", "  FAIL ", ok), label
        if (.not. ok) nfail = nfail + 1
    end subroutine check

end program test_varslice_units_attr
