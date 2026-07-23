program test_conservative_rotated
    ! Conservative remapping FROM a rotated-pole source grid (RACMO-like) onto a
    ! polar-stereographic target, gen="coords". Guards conservative_planar's
    ! cross-system (do_project) path: a rotated source's axis values are ROTATED
    ! lon/lat (rlon/rlat), not geographic, so the cell corners must be inverse-
    ! rotated to geographic before projecting into the target plane. Without that
    ! the source polygons land at the wrong place (equatorial, for the ANT pole)
    ! and NO target cell is covered.
    !
    ! Self-contained (no cdo). Checks: (a) the source covers the whole target,
    ! (b) a constant field is preserved exactly, (c) a field equal to geographic
    ! latitude remaps back to the target cell-center latitude (correct placement).

    use coords
    implicit none

    type(grid_class) :: gs, gt
    type(map_class)  :: m
    integer, parameter :: nxs=101, nys=108, nxt=61, nyt=61
    real(dp), parameter :: lon_p=20.0_dp, lat_p=5.0_dp   ! RACMO ANT rotated N-pole
    real(dp), parameter :: rlon0=-20.0_dp, rlat0=-28.0_dp, drl=0.4_dp
    character(len=*), parameter :: fldr = "maps_con_rot_test"
    real(dp) :: vs(nxs,nys), vc(nxt,nyt)
    logical  :: msk(nxt,nyt)
    integer  :: i, j, ncov, fails
    real(dp) :: emax_const, emax_lat, d

    call execute_command_line("rm -rf "//fldr)
    fails = 0

    ! Rotated-pole source: axis = rotated lon/lat (degrees), geographic lon/lat
    ! recovered by inverse rotation. This footprint covers the south polar cap.
    call grid_init(gs, name="racmo-rot", mtype="rotated_latitude_longitude", units="degrees", &
                   x0=rlon0, dx=drl, nx=nxs, y0=rlat0, dy=drl, ny=nys, &
                   lambda=lon_p, phi=lat_p, lon180=.true.)

    ! Target: polar-stereographic over Antarctica, inside the source footprint.
    call grid_init(gt, name="ANT-stereo", mtype="polar_stereographic", units="kilometers", &
                   x0=-1500.0_dp, dx=50.0_dp, nx=nxt, y0=-1500.0_dp, dy=50.0_dp, ny=nyt, &
                   lambda=0.0_dp, phi=-71.0_dp, lon180=.true.)

    call map_init(m, gs, gt, method="con", gen="coords", fldr=fldr, load=.false.)

    ! (a) coverage: every target cell lies inside the source footprint
    call map_field(m, "c", vs*0.0_dp + 5.0_dp, vc, stat="mean", mask2=msk)
    ncov = count(msk)
    write(*,"(a,i6,a,i6)") " covered target cells:", ncov, " / ", nxt*nyt
    if (ncov /= nxt*nyt) then
        write(*,*) "FAIL: rotated source did not fully cover target"; fails = fails + 1
    end if

    ! (b) constant field preserved exactly on covered cells
    emax_const = maxval(abs(vc - 5.0_dp), mask=msk)
    write(*,"(a,es12.3)") " max|const-5| on covered:", emax_const
    if (emax_const > 1.0e-9_dp) then
        write(*,*) "FAIL: constant field not preserved"; fails = fails + 1
    end if

    ! (c) latitude field remaps to target cell-center latitude (placement check)
    do j = 1, nys
        do i = 1, nxs
            vs(i,j) = gs%lat(i,j)
        end do
    end do
    call map_field(m, "v", vs, vc, stat="mean", mask2=msk)
    emax_lat = 0.0_dp
    do j = 1, nyt
        do i = 1, nxt
            if (.not. msk(i,j)) cycle
            d = abs(vc(i,j) - gt%lat(i,j))
            emax_lat = max(emax_lat, d)
        end do
    end do
    write(*,"(a,f10.4)") " max|remapped_lat - target_lat| (deg):", emax_lat
    if (emax_lat > 1.0_dp) then   ! target cell ~0.45 deg; slack for cell averaging
        write(*,*) "FAIL: geographic placement wrong (lat mismatch)"; fails = fails + 1
    end if

    call execute_command_line("rm -rf "//fldr)

    if (fails > 0) stop 1
    write(*,*) "PASS: test_conservative_rotated"

end program test_conservative_rotated
