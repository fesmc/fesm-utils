program test_zremap
    ! Exercise the zremap module on a synthetic ice dome with a field linear in
    ! elevation, f = f0 + g*z, and end gradients dfdz_lo = dfdz_hi = g, so the
    ! remapped field must equal f0 + g*z_now wherever there are data:
    !   * reconstruction at z_ref and remapping to z_ref-200 m (interior);
    !   * extrapolation above the summit and below the dome margin;
    !   * the fallback f_ref + dfdz_lo*(z_now-z_ref) far from the dome;
    !   * the gradient diagnostic dfdz = g on the dome;
    !   * the same with stride 1 and stride 3 (irregular last table point).
    ! The namelist is written by the test itself, so no external files are needed.

    use precision, only: wp
    use constants, only: mv
    use zremap

    implicit none

    character(len=*), parameter :: nml_file = "test_zremap_tmp.nml"
    character(len=*), parameter :: out_file = "test_zremap_out.nc"
    integer,  parameter :: nx = 101, ny = 81
    real(wp), parameter :: dx = 10e3_wp             ! [m]
    real(wp), parameter :: f0 = 1.0_wp, g = 2e-3_wp ! f = f0 + g*z
    real(wp), parameter :: f_out = 5.0_wp           ! f_ref outside the dome
    real(wp), parameter :: R = 35.0_wp              ! [points] dome radius
    real(wp), parameter :: H = 3000.0_wp            ! [m] dome height
    real(wp), parameter :: tol = 1e-3_wp

    real(wp) :: z_ref(nx,ny), f_ref(nx,ny), rr(nx,ny)
    logical  :: mask(nx,ny)
    integer  :: nfail, i, j

    nfail = 0

    write(*,*) "==================================================="
    write(*,*) " zremap test suite"
    write(*,*) "==================================================="
    write(*,"(a)") ""

    do j = 1, ny
    do i = 1, nx
        rr(i,j) = sqrt(real(i-51,wp)**2 + real(j-41,wp)**2) / R
    end do
    end do
    mask  = rr < 1.0_wp
    z_ref = 0.0_wp
    where (mask) z_ref = H*sqrt(1.0_wp - rr**2)
    f_ref = f_out
    where (mask) f_ref = f0 + g*z_ref

    call run_case(1)
    call run_case(3)

    write(*,"(a)") ""
    if (nfail == 0) then
        write(*,"(a)") " All zremap tests PASSED."
    else
        write(*,"(a,i0,a)") " zremap tests FAILED: ", nfail, " case(s)."
    end if

    call remove_file(nml_file)

    if (nfail > 0) stop 1

contains

    subroutine run_case(stride)
        integer, intent(in) :: stride

        type(zremap_class) :: zr
        real(wp) :: z_now(nx,ny), f_now(nx,ny), err
        logical  :: inner(nx,ny)
        character(len=16) :: lbl

        write(lbl,"(a,i0)") "stride=", stride
        call write_namelist(stride)

        call zremap_init(zr, nml_file, "zremap", nx, ny, dx)
        call zremap_update(zr, f_ref, z_ref, mask)

        ! Points well inside the dome (all table points around them have data)
        inner = rr < 0.8_wp

        ! Reconstruction at the reference elevation
        call zremap_map(zr, f_now, z_ref)
        err = maxval(abs(f_now - (f0 + g*z_ref)), mask=inner)
        call check(trim(lbl)//" reconstruction", err)

        ! Lowering by 200 m
        z_now = z_ref - 200.0_wp
        call zremap_map(zr, f_now, z_now)
        err = maxval(abs(f_now - (f0 + g*z_now)), mask=inner)
        call check(trim(lbl)//" lowering 200 m", err)

        ! Extrapolation above the summit and below the dome margin
        z_now = z_ref
        z_now(51,41) = H + 500.0_wp
        z_now(51,41-30) = -300.0_wp     ! r = 30/35 of the radius
        call zremap_map(zr, f_now, z_now)
        call check(trim(lbl)//" above summit", abs(f_now(51,41) - (f0 + g*z_now(51,41))))
        call check(trim(lbl)//" below margin", abs(f_now(51,11) - (f0 + g*z_now(51,11))))

        ! Fallback far from the dome: f_ref + dfdz_lo*(z_now-z_ref)
        z_now(1,1) = 1000.0_wp
        call zremap_map(zr, f_now, z_now)
        call check(trim(lbl)//" fallback", abs(f_now(1,1) - (f_out + g*1000.0_wp)))

        ! Gradient diagnostic
        err = maxval(abs(zr%dfdz - g), mask=inner)
        call check(trim(lbl)//" dfdz (relative)", err/g)

        if (stride == 3) then
            call zremap_write(zr, out_file)
            write(*,"(a)") "  written: "//out_file
        end if

        call zremap_end(zr)

    end subroutine run_case

    subroutine check(name, err)
        character(len=*), intent(in) :: name
        real(wp),         intent(in) :: err

        if (err <= tol) then
            write(*,"(a,a,es10.2)") "  PASS  ", name//"  err=", err
        else
            write(*,"(a,a,es10.2)") "  FAIL  ", name//"  err=", err
            nfail = nfail + 1
        end if

    end subroutine check

    subroutine write_namelist(stride)
        integer, intent(in) :: stride
        integer :: u

        open(newunit=u, file=nml_file, status="replace", action="write")
        write(u,"(a)")        "&zremap"
        write(u,"(a)")        "    sigma   = 50e3"
        write(u,"(a)")        "    z_min   = 0.0"
        write(u,"(a)")        "    z_max   = 3200.0"
        write(u,"(a)")        "    dz      = 100.0"
        write(u,"(a)")        "    w_min   = 1e-3"
        write(u,"(a,i0)")     "    stride  = ", stride
        write(u,"(a)")        "    dfdz_lo = 2e-3"
        write(u,"(a)")        "    dfdz_hi = 2e-3"
        write(u,"(a)")        "/"
        close(u)

    end subroutine write_namelist

    subroutine remove_file(filename)
        character(len=*), intent(in) :: filename
        integer :: u, ios

        open(newunit=u, file=filename, status="old", iostat=ios)
        if (ios == 0) close(u, status="delete")

    end subroutine remove_file

end program test_zremap
