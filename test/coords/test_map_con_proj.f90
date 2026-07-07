program test_map_con_proj
    ! Regression test for conservative remapping onto a PROJECTED (polar-
    ! stereographic) target -- the climber geo/ice coupling case (lat-lon geo
    ! grid -> stereographic ice grid). The in-package (coords) conservative map
    ! must agree closely with cdo's spherical conservative map (gencon/YAC).
    !
    ! Guards the projected-source edge subsampling in conservative_planar: a
    ! rough field over a 1deg source shows the chord (straight-projected-edge)
    ! bias that subsampling removes -- without it the bulk max diff here is ~0.8,
    ! with it ~0.09. Requires cdo on PATH; skips cleanly otherwise.

    use coords

    implicit none

    type(grid_class) :: gs, gt
    type(map_class)  :: m_co, m_cdo
    integer, parameter :: nxs=360, nys=180, nxt=173, nyt=173
    real(dp), parameter :: d2r = 3.14159265358979_dp/180.0_dp
    real(dp), parameter :: tol_bulk = 0.30_dp     ! new code ~0.09; old (chord) ~0.79
    character(len=*), parameter :: fldr = "maps_con_proj_test"
    real(dp) :: vs(nxs,nys), vco(nxt,nyt), vcdo(nxt,nyt), vconst(nxt,nyt)
    logical  :: mco(nxt,nyt), mcdo(nxt,nyt), mc(nxt,nyt)
    integer  :: i, j, stat, fails, only_co, only_cdo
    real(dp) :: d, lat, dmax_bulk, dmax_pole, emax_const

    call execute_command_line("command -v cdo > /dev/null 2>&1", exitstat=stat)
    if (stat /= 0) then
        write(*,*) "SKIP: test_map_con_proj (cdo not found on PATH)"
        stop 0
    end if
    call execute_command_line("rm -rf "//fldr)

    ! Source: 1deg global lon/lat.  Target: NH polar-stereographic (climber-like).
    call grid_init(gs, name="ll1deg", mtype="latlon", units="degrees", &
                   x0=-179.5_dp, dx=1.0_dp, nx=nxs, y0=-89.5_dp, dy=1.0_dp, ny=nys)
    call grid_init(gt, name="NHr-64KM", mtype="polar_stereographic", units="kilometers", &
                   x0=-5504.0_dp, dx=64.0_dp, nx=nxt, y0=-5504.0_dp, dy=64.0_dp, ny=nyt, &
                   lambda=-45.0_dp, phi=70.0_dp, alpha=1.0e-3_dp, lon180=.true.)

    fails = 0

    ! rough source field (high-wavenumber, so cell-area weighting errors show)
    do j = 1, nys
        do i = 1, nxs
            vs(i,j) = 1000.0_dp * cos(9.0_dp*gs%lat(i,j)*d2r) &
                                * cos(13.0_dp*gs%lon(i,j)*d2r)
        end do
    end do

    call map_init(m_co, gs, gt, method="con", gen="coords", fldr=fldr, load=.false.)

    ! --- constant field must be preserved exactly (partition normalizes) ---
    call map_field(m_co, "c", vs*0.0_dp + 5.0_dp, vconst, stat="mean", mask2=mc)
    emax_const = maxval(abs(vconst - 5.0_dp), mask=mc)
    if (emax_const > 1.0e-9_dp) then
        write(*,*) "FAIL: constant field not preserved, err =", emax_const; fails = fails + 1
    end if

    ! --- rough field: coords-con vs cdo-con ---
    call map_field(m_co, "v", vs, vco, stat="mean", mask2=mco)

    call map_init(m_cdo, gs, gt, method="con", gen="cdo", fldr=fldr, load=.false.)
    call map_field(m_cdo, "v", vs, vcdo, stat="mean", mask2=mcdo)

    only_co = 0; only_cdo = 0; dmax_bulk = 0.0_dp; dmax_pole = 0.0_dp
    do j = 1, nyt
        do i = 1, nxt
            if (mco(i,j) .and. .not. mcdo(i,j)) only_co = only_co + 1
            if (mcdo(i,j) .and. .not. mco(i,j)) only_cdo = only_cdo + 1
            if (.not. (mco(i,j) .and. mcdo(i,j))) cycle
            d = abs(vco(i,j) - vcdo(i,j)); lat = gt%lat(i,j)
            if (lat < 85.0_dp) then
                dmax_bulk = max(dmax_bulk, d)
            else
                dmax_pole = max(dmax_pole, d)
            end if
        end do
    end do

    write(*,"(a,i5,a,i5)") " conservative coverage mismatch: only-coords=", only_co, &
                           "  only-cdo=", only_cdo
    write(*,"(a,f8.4,a,f8.4)") " max|coords-cdo| con: bulk(lat<85)=", dmax_bulk, &
                               "  cap(lat>85)=", dmax_pole

    ! conservative maps every target cell that its source covers -> coverage must match
    if (only_co /= 0 .or. only_cdo /= 0) then
        write(*,*) "FAIL: coords/cdo conservative coverage differs"; fails = fails + 1
    end if
    if (dmax_bulk > tol_bulk) then
        write(*,*) "FAIL: coords-con vs cdo-con bulk diff too large =", dmax_bulk; fails = fails + 1
    end if

    call execute_command_line("rm -rf "//fldr)

    if (fails > 0) stop 1
    write(*,*) "PASS: test_map_con_proj"

end program test_map_con_proj
