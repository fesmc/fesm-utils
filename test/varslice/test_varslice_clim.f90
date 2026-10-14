program test_varslice_clim
    ! varslice with climatologies (a sub-annual cycle without years, time =
    ! 0 ... 12): to(x,y,depth,month) and t2m(x,y,month). A plain update keeps
    ! the cycle (the last dimension), a range mean with rep=1 gives its mean,
    ! and varslice_sub_mean averages the cycle.

    use precision
    use ncio
    use constants, only : mv
    use varslice

    implicit none

    character(len=*), parameter :: fname = "test_varslice_clim.nc"
    character(len=*), parameter :: fnml  = "test/varslice/test_varslice_clim.nml"
    integer, parameter :: nx = 4, ny = 3, nz = 2, nm = 12
    real(wp), parameter :: tol = 1e-4_wp

    type(varslice_class) :: vo, va
    real(wp) :: to(nx,ny,nz,nm), t2m(nx,ny,nm)
    real(wp), allocatable :: v(:,:,:)
    integer :: m, nfail

    nfail = 0

    ! Month m: value m (plus 10*level for the ocean), so the cycle mean is 6.5
    do m = 1, nm
        to(:,:,1,m) = real(m,wp)
        to(:,:,2,m) = real(m,wp) + 10.0_wp
        t2m(:,:,m)  = 250.0_wp + real(m,wp)
    end do
    to(1,1,2,:) = mv        ! a cell without data at depth

    call nc_create(fname)
    call nc_write_dim(fname, "xc",    x=1.0_wp, dx=1.0_wp, nx=nx)
    call nc_write_dim(fname, "yc",    x=1.0_wp, dx=1.0_wp, nx=ny)
    call nc_write_dim(fname, "depth", x=[5.0_wp, 50.0_wp])
    call nc_write_dim(fname, "month", x=1.0_wp, dx=1.0_wp, nx=nm)
    call nc_write(fname, "to",  to,  dim1="xc", dim2="yc", dim3="depth", dim4="month", missing_value=mv)
    call nc_write(fname, "t2m", t2m, dim1="xc", dim2="yc", dim3="month")

    ! --- ocean (x,y,depth,month) ---
    call varslice_init_nml(vo, fnml, "ocn_clim")

    call varslice_update(vo)
    call check("ocean cycle kept: shape", all(shape(vo%var) .eq. [nx,ny,nz,nm]), nfail)
    call check("ocean cycle kept: July", abs(vo%var(2,2,1,7) - (7.0_wp+273.15_wp)) .lt. tol, nfail)
    call check("ocean depth axis", abs(vo%z(2) - 50.0_wp) .lt. tol, nfail)

    v = varslice_sub_mean(vo)
    call check("ocean sub_mean", abs(v(2,2,2) - (16.5_wp+273.15_wp)) .lt. tol, nfail)
    call check("ocean sub_mean missing", v(1,1,2) .eq. mv, nfail)

    call varslice_update(vo, [0.0_dp,0.0_dp], method="range_mean", rep=1)
    call check("ocean range_mean rep=1: shape", all(shape(vo%var) .eq. [nx,ny,nz,1]), nfail)
    call check("ocean range_mean rep=1: mean", abs(vo%var(3,2,1,1) - (6.5_wp+273.15_wp)) .lt. tol, nfail)

    call varslice_update(vo, [0.0_dp,0.0_dp], method="range", rep=12)
    call check("ocean range rep=12: cycle again", size(vo%var,4) .eq. nm, nfail)

    ! --- atmosphere (x,y,month): unchanged layout, month in the third slot ---
    call varslice_init_nml(va, fnml, "atm_clim")
    call varslice_update(va, [0.0_dp,0.0_dp], method="range_mean", rep=12)
    call check("atm cycle: shape", all(shape(va%var) .eq. [nx,ny,nm,1]), nfail)
    call check("atm cycle: March", abs(va%var(1,1,3,1) - 253.0_wp) .lt. tol, nfail)
    call varslice_update(va, [0.0_dp,0.0_dp], method="range_mean", rep=1)
    call check("atm range_mean rep=1", abs(va%var(1,1,1,1) - 256.5_wp) .lt. tol, nfail)

    if (nfail .gt. 0) then
        write(*,*) "test_varslice_clim: ", nfail, " FAILED"
        error stop 1
    end if
    write(*,*) "test_varslice_clim: all passed"

contains

    subroutine check(label, ok, nfail)
        character(len=*), intent(IN)    :: label
        logical,          intent(IN)    :: ok
        integer,          intent(INOUT) :: nfail
        write(*,"(a,a)") merge("  ok   ", "  FAIL ", ok), label
        if (.not. ok) nfail = nfail + 1
    end subroutine check

end program test_varslice_clim
