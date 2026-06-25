program test_conservative_xsys
    ! Cross-system conservative: global lat-lon source -> regional stereographic
    ! target. Source cell corners are projected into the target plane and clipped.
    ! Checks: (1) a constant field maps to a constant; (2) for interior target
    ! cells fully covered by the source, the projected overlap areas partition
    ! the target cell exactly (sum of link weights == planar cell area).

    use coords

    implicit none

    type(grid_class) :: gs, gt
    type(map_class)  :: map
    integer, parameter :: nxs=73, nys=40, nxt=20, nyt=30
    real(dp) :: vs(nxs,nys), vt(nxt,nyt)
    logical  :: m2(nxt,nyt)
    integer  :: i, j, c, j1, j2, fails
    real(dp) :: wsum, cell_area, emax_sum, emax_const

    ! Global-ish lat-lon source covering the target generously
    call grid_init(gs, name="ll-src", mtype="latlon", units="degrees", &
                   x0=-120.0_dp, dx=2.5_dp, nx=nxs, y0=50.0_dp, dy=1.0_dp, ny=nys)

    ! Regional oblique-stereographic target over central Greenland
    call grid_init(gt, name="grl-tgt", mtype="stereographic", units="kilometers", &
                   dx=40.0_dp, nx=nxt, dy=40.0_dp, ny=nyt, &
                   lambda=-40.0_dp, phi=72.0_dp, alpha=7.5_dp)

    call map_init_conservative(map, gs, gt)

    fails = 0

    ! (1) constant field
    vs = 7.0_dp
    call map_field(map, "const", vs, vt, method="mean", mask2=m2)
    emax_const = maxval(abs(vt - 7.0_dp), mask=m2)
    write(*,*) "covered cells:", count(m2), " of", size(m2)
    write(*,*) "constant-field max error =", emax_const
    if (count(m2) < size(m2)) then
        write(*,*) "FAIL: not all target cells covered"; fails = fails + 1
    end if
    if (emax_const > 1.0e-9_dp) then
        write(*,*) "FAIL: constant field not preserved"; fails = fails + 1
    end if

    ! (2) interior cells: sum of overlap weights == planar cell area (dx*dy)
    cell_area = gt%G%dx * gt%G%dy
    emax_sum = 0.0_dp
    do j = 2, nyt-1
        do i = 2, nxt-1
            c = (j-1)*nxt + i
            j1 = map%wm%dst_off(c); j2 = map%wm%dst_off(c+1) - 1
            wsum = 0.0_dp
            if (j2 >= j1) wsum = sum(map%wm%w(j1:j2))
            emax_sum = max(emax_sum, abs(wsum - cell_area))
        end do
    end do
    write(*,*) "interior cell area =", cell_area, " max |sum w - area| =", emax_sum
    if (emax_sum > 1.0e-6_dp * cell_area) then
        write(*,*) "FAIL: projected overlaps do not partition interior cells"; fails = fails + 1
    end if

    if (fails > 0) stop 1
    write(*,*) "PASS: test_conservative_xsys"

end program test_conservative_xsys
