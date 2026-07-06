program test_conservative_desclat
    ! Regression: analytic separable lat-lon conservative (conservative_latlon)
    ! with a target whose latitude axis is stored DESCENDING (north -> south),
    ! the standard ordering for Gaussian grids (e.g. VILMA's Gauss-Legendre grid).
    !
    ! The separable path builds sin(lat) bands from the axis edges; before the
    ! orientation fix, interval_overlap assumed ascending edges, so a descending
    ! target produced zero band overlaps -> zero links -> an all-missing map (the
    ! "n_links = 0" symptom on the geo(0.125 deg) -> vilma con map). Here the
    ! target is an exact 2x refinement of the source over the same domain but with
    ! its latitude axis reversed, so a correct map must cover every target cell and
    ! conserve the area integral.

    use coords

    implicit none

    type(grid_class) :: gs, gt
    type(map_class)  :: map
    real(dp) :: vs(4,3), vt(8,6)
    real(dp) :: ylat_desc(6)
    logical  :: m2(8,6)
    integer  :: i, j, fails
    real(dp) :: emax_const

    ! Source 4x3 cells: [0,40] lon x [30,60] lat, ascending latitude
    call grid_init(gs, name="ll-src", mtype="latlon", units="degrees", &
                   x0=5.0_dp, dx=10.0_dp, nx=4, y0=35.0_dp, dy=10.0_dp, ny=3)

    ! Target 8x6 cells: exact 2x refinement of the same domain, but latitude axis
    ! stored DESCENDING (57.5, 52.5, ... , 32.5) to mimic a Gaussian grid.
    do j = 1, 6
        ylat_desc(j) = 57.5_dp - real(j-1,dp)*5.0_dp
    end do
    call grid_init(gt, name="ll-tgt-desc", mtype="latlon", units="degrees", &
                   x=[(2.5_dp + real(i-1,dp)*5.0_dp, i=1,8)], y=ylat_desc)

    do j = 1, 3
        do i = 1, 4
            vs(i,j) = real(i,dp) + 10.0_dp*real(j,dp)
        end do
    end do

    call map_init_conservative(map, gs, gt)

    fails = 0

    ! (0) links were actually generated (the direct n_links = 0 regression)
    write(*,*) "n_links =", map%wm%n_links
    if (map%wm%n_links == 0) then
        write(*,*) "FAIL: conservative_latlon produced zero links for descending-lat target"
        fails = fails + 1
    end if

    ! (1) every target cell is covered (before the fix: none were)
    call map_field(map, "c", vs*0.0_dp + 3.0_dp, vt, stat="mean", mask2=m2)
    write(*,*) "covered:", count(m2), " of", size(m2)
    if (count(m2) < size(m2)) then
        write(*,*) "FAIL: not all target cells covered"; fails = fails + 1
    end if

    ! (2) constant field preserved exactly (partition normalizes to 1)
    emax_const = maxval(abs(vt - 3.0_dp), mask=m2)
    write(*,*) "const err =", emax_const
    if (emax_const > 1.0e-12_dp) then
        write(*,*) "FAIL: constant field not preserved"; fails = fails + 1
    end if

    ! (3) global area-integral conservation (physical grid areas). Tolerance is
    ! loose (matches test_conservative_sphere): the map weights use exact sin(lat)
    ! bands while grid%area uses the coordinates cell-area formula, so the two
    ! area conventions differ at ~1e-3. The exact conservation invariant is the
    ! constant-field check above (err = 0); this is a sanity bound.
    call map_field(map, "v", vs, vt, stat="mean", mask2=m2)
    block
        real(dp) :: int_s, int_t
        int_s = sum(vs * gs%area)
        int_t = sum(vt * gt%area)
        write(*,*) "source integral =", int_s, " target integral =", int_t, &
                   " rel diff =", abs(int_s-int_t)/abs(int_s)
        if (abs(int_s - int_t) > 1.0e-2_dp * abs(int_s)) then
            write(*,*) "FAIL: area integral not conserved"; fails = fails + 1
        end if
    end block

    if (fails > 0) stop 1
    write(*,*) "PASS: test_conservative_desclat"

end program test_conservative_desclat
