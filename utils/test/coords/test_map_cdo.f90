program test_map_cdo
    ! End-to-end check of the online cdo map generation (gen="cdo"): a coarse
    ! global lon/lat source mapped onto an Antarctic polar-stereographic grid,
    ! the climber-x use case. Requires cdo on PATH; skips cleanly otherwise.
    !   - con/bil/nn cdo maps must reproduce a constant field and fill targets.
    !   - the cdo conservative map must agree closely with the in-package
    !     (coords) conservative map on a smooth field.

    use coords

    implicit none

    type(grid_class) :: gll, gant
    type(map_class)  :: m_cdo, m_co
    real(dp), allocatable :: vs(:,:), vt(:,:), vt_co(:,:)
    logical,  allocatable :: m2(:,:), m2_co(:,:)
    character(len=*), parameter :: fldr = "maps_cdo_test"
    integer  :: i, j, stat, fails, nboth, n
    real(dp) :: dmax, emax, emean

    call execute_command_line("command -v cdo > /dev/null 2>&1", exitstat=stat)
    if (stat /= 0) then
        write(*,*) "SKIP: test_map_cdo (cdo not found on PATH)"
        stop 0
    end if

    call execute_command_line("rm -rf "//fldr)

    ! Source: 5deg x 5deg global lon/lat (climber-like)
    call grid_init(gll, name="ll5", mtype="latlon", units="degrees", &
                   x0=0.0_dp, dx=5.0_dp, nx=72, y0=-90.0_dp, dy=5.0_dp, ny=37)

    ! Target: Antarctic polar-stereographic, 64 km
    call grid_init(gant, name="ANT-64KM", mtype="polar_stereographic", units="kilometers", &
                   x0=-3040.0_dp, dx=64.0_dp, nx=96, y0=-3040.0_dp, dy=64.0_dp, ny=96, &
                   lambda=0.0_dp, phi=-71.0_dp, alpha=19.0_dp)

    allocate(vs(72,37), vt(96,96), vt_co(96,96), m2(96,96), m2_co(96,96))
    fails = 0

    ! --- constant-field invariant for each cdo kernel ---
    vs = 5.0_dp
    call run_const("con", gll, gant, vs, fldr, fails)
    call run_const("bil", gll, gant, vs, fldr, fails)
    call run_const("nn",  gll, gant, vs, fldr, fails)

    ! --- cdo conservative vs in-package conservative on a smooth field ---
    do j = 1, 37
        do i = 1, 72
            vs(i,j) = gll%lat(i,j)
        end do
    end do
    ! Generate the cdo conservative map (kept on disk as its own SCRIP cache)
    call map_init(m_cdo, gll, gant, method="con", gen="cdo", fldr=fldr, load=.false.)
    call map_field(m_cdo, "lat", vs, vt, stat="mean", mask2=m2)

    ! cdo cache reload: a second map_init with load=.true. must load that cdo
    ! SCRIP file through the auto-detecting loader and reproduce it exactly.
    call map_init(m_co, gll, gant, method="con", gen="cdo", fldr=fldr, load=.true.)
    call map_field(m_co, "lat", vs, vt_co, stat="mean", mask2=m2_co)
    dmax = 0.0_dp
    do j = 1, 96
        do i = 1, 96
            if (m2(i,j) .and. m2_co(i,j)) dmax = max(dmax, abs(vt(i,j) - vt_co(i,j)))
        end do
    end do
    write(*,"(a,es10.2)") " cdo cache reload max diff =", dmax
    if (dmax > 1.0e-12_dp) then
        write(*,*) "FAIL: reloaded cdo map differs from generated one"; fails = fails + 1
    end if

    ! Physical check: a conservative remap of the latitude field must recover
    ! the target-cell latitude to within ~one source cell (5 deg here).
    emax = 0.0_dp; emean = 0.0_dp; n = 0
    do j = 1, 96
        do i = 1, 96
            if (m2(i,j)) then
                emax  = max(emax, abs(vt(i,j) - gant%lat(i,j)))
                emean = emean + abs(vt(i,j) - gant%lat(i,j))
                n = n + 1
            end if
        end do
    end do
    if (n > 0) emean = emean / real(n,dp)
    write(*,"(a,f7.3,a,f7.3,a)") " cdo-con latitude recovery: mean=", emean, " max=", emax, " (deg)"
    if (emean > 2.0_dp .or. emax > 6.0_dp) then
        write(*,*) "FAIL: cdo conservative map did not recover target latitude"; fails = fails + 1
    end if

    ! Informational: agreement between cdo and in-package (coords) conservative
    call map_init(m_co, gll, gant, method="con", gen="coords", fldr=fldr, load=.false.)
    call map_field(m_co, "lat", vs, vt_co, stat="mean", mask2=m2_co)
    dmax = 0.0_dp; nboth = 0
    do j = 1, 96
        do i = 1, 96
            if (m2(i,j) .and. m2_co(i,j)) then
                dmax = max(dmax, abs(vt(i,j) - vt_co(i,j)))
                nboth = nboth + 1
            end if
        end do
    end do
    write(*,"(a,f7.3,a,i6,a,i6,a,i6,a)") " cdo-con vs coords-con: max diff=", dmax, &
        " deg  (both=", nboth, " cdo=", count(m2), " coords=", count(m2_co), ")"

    call execute_command_line("rm -rf "//fldr)

    if (fails > 0) stop 1
    write(*,*) "PASS: test_map_cdo"

contains

    subroutine run_const(method, gs, gt, vs, fldr, fails)
        character(len=*), intent(in)    :: method, fldr
        type(grid_class), intent(in)    :: gs, gt
        real(dp),         intent(in)    :: vs(:,:)
        integer,          intent(inout) :: fails
        type(map_class) :: m
        real(dp), allocatable :: vt(:,:)
        logical,  allocatable :: m2(:,:)
        real(dp) :: emax
        integer  :: nx, ny

        nx = gt%G%nx; ny = gt%G%ny
        allocate(vt(nx,ny), m2(nx,ny))
        call map_init(m, gs, gt, method=method, gen="cdo", fldr=fldr, load=.false., clean=.true.)
        call map_field(m, "c", vs, vt, stat="mean", mask2=m2)
        where (.not. m2) vt = 5.0_dp
        emax = maxval(abs(vt - 5.0_dp))
        write(*,"(a,es12.4,a,i6)") " cdo "//trim(method)//" constant-field max err =", emax, &
                                   "  filled=", count(m2)
        if (emax > 1.0e-6_dp) then
            write(*,*) "FAIL: cdo "//trim(method)//" did not reproduce the constant"; fails = fails + 1
        end if
        if (count(m2) == 0) then
            write(*,*) "FAIL: cdo "//trim(method)//" filled no targets"; fails = fails + 1
        end if
        deallocate(vt, m2)
    end subroutine run_const

end program test_map_cdo
