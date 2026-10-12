program test_zremap
    ! Exercise the zremap module on a synthetic ice dome with a field linear in
    ! elevation, f = f0 + g*z, and end gradients dfdz_lo = dfdz_hi = g, so the
    ! remapped field must equal f0 + g*z_now wherever there are data:
    !   * reconstruction at z_ref and remapping to z_ref-200 m (interior);
    !   * extrapolation above the summit and below the dome margin;
    !   * the fallback f_ref + dfdz_lo*(z_now-z_ref) far from the dome;
    !   * the gradient diagnostic dfdz = g on the dome;
    !   * the support w (~1 in the dome, 0 far away) and dz_ext (0 inside the
    !     data range, > 0 above the summit);
    !   * continuity from the dome to the fallback (off the dome f_ref rises
    !     with distance, off the linear relation): no jumps along a row at
    !     fixed z_now beyond the slope of the fallback itself;
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
    real(wp), parameter :: df_out = 0.1_wp          ! d(f_ref)/d(points) outside the dome
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
    f_ref = f0 + df_out*(rr - 1.0_wp)*R
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
        real(wp) :: z_now(nx,ny), f_now(nx,ny), dz_ext(nx,ny), err
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
        call zremap_map(zr, f_now, z_now, dz_ext)
        call check(trim(lbl)//" above summit", abs(f_now(51,41) - (f0 + g*z_now(51,41))))
        call check(trim(lbl)//" below margin", abs(f_now(51,11) - (f0 + g*z_now(51,11))))

        ! Extrapolation distance: within the data range on the dome (away from
        ! the summit), beyond it above the summit and below the margin
        call check(trim(lbl)//" dz_ext inside", maxval(abs(dz_ext), mask=inner .and. rr > 0.2_wp))
        call check(trim(lbl)//" dz_ext summit >= 500 m", max(0.0_wp, 500.0_wp - dz_ext(51,41)))
        call check(trim(lbl)//" dz_ext margin < 0", max(0.0_wp, dz_ext(51,11)))
        call check(trim(lbl)//" dz_ext no data = mv", abs(dz_ext(1,1) - mv))

        ! Support: ~1 at the centre, 0 far away
        call check(trim(lbl)//" w centre", abs(zr%w(51,41) - 1.0_wp)*1e-1_wp)
        call check(trim(lbl)//" w corner", zr%w(1,1))

        ! Fallback far from the dome: f_ref + dfdz_lo*(z_now-z_ref)
        z_now(1,1) = 1000.0_wp
        call zremap_map(zr, f_now, z_now)
        call check(trim(lbl)//" fallback", abs(f_now(1,1) - (f_ref(1,1) + g*1000.0_wp)))

        ! Continuity from the dome to the fallback along the centre row at
        ! z_now = 1000 m: profile f0+g*1000 inside, fallback f_ref+g*1000 far
        ! away; the largest step between neighbours stays within twice the
        ! fallback slope df_out (a hard switch would jump by ~1.5)
        z_now = 1000.0_wp
        call zremap_map(zr, f_now, z_now)
        err = maxval(abs(f_now(2:nx,41) - f_now(1:nx-1,41)))
        write(*,"(a,es10.2)") "        max step along the row: ", err
        call check(trim(lbl)//" continuity", max(0.0_wp, err - 2.0_wp*df_out))
        call check(trim(lbl)//" continuity centre", abs(f_now(51,41) - (f0 + g*1000.0_wp)))
        call check(trim(lbl)//" continuity edge", abs(f_now(1,41) - (f_ref(1,41) + g*1000.0_wp)))

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
        write(u,"(a)")        "    w_full  = 0.25"
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
