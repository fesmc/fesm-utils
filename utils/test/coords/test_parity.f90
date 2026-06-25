program test_parity
    ! Exercise the non-grid->grid map_field combos and their sp/int wrappers
    ! (the interpolation-parity gaps closed in coords-followups.md §2):
    !   points->grid (2D target: reset / mask_pack / fill / smoothing)
    !   grid->points, points->points (1D targets: reset / mask_pack)
    !
    ! Two kinds of check:
    !   - A constant source transfers exactly to every covered target, and
    !     reset / mask_pack behave as documented (robust, geometry-free).
    !   - fill / smoothing on points->grid must match the trusted grid->grid
    !     path bit-for-bit (same neighbour structure), proving the optional
    !     post-apply controls are forwarded identically.

    use coords

    implicit none

    type(grid_class)   :: gs, gt
    type(points_class) :: ps, pt
    type(map_class)    :: m_pg, m_gp, m_pp
    integer, parameter :: nxs=5, nys=5, ns=nxs*nys
    integer, parameter :: nxt=4, nyt=4
    integer, parameter :: ntp=6
    real(dp), parameter :: C = 7.0_dp, miss = -9999.0_dp, rtol = 1.0e-9_dp
    real(dp) :: slon(ns), slat(ns), tlon(ntp), tlat(ntp)
    integer  :: ix, iy, c0, fails

    real(dp) :: vs2(nxs,nys), vs1(ns)
    real(dp) :: vg(nxt,nyt), vp(ntp)
    logical  :: mg(nxt,nyt), mp(ntp), packg(nxt,nyt), packp(ntp)
    real(sp) :: vs1_sp(ns), vs2_sp(nxs,nys), vg_sp(nxt,nyt), vp_sp(ntp)
    integer  :: vs1_i(ns),  vs2_i(nxs,nys),  vg_i(nxt,nyt),  vp_i(ntp)

    ! larger grids for the fill/smooth forwarding comparison
    type(grid_class)   :: gsB, gtB
    type(points_class) :: psB
    type(map_class)    :: mGB, mPB
    integer, parameter :: nxsb=8, nysb=8, nsb=nxsb*nysb
    integer, parameter :: nxtb=12, nytb=12
    real(dp) :: slonB(nsb), slatB(nsb)
    real(dp) :: fB2d(nxsb,nysb), fB1d(nsb), fB2d_h(nxsb,nysb), fB1d_h(nsb)
    real(dp) :: vA(nxtb,nytb), vB(nxtb,nytb), vBase(nxtb,nytb)
    integer  :: nmiss0

    fails = 0

    ! --- source: 5x5 latlon grid + matching point set (grid flatten order) ----
    c0 = 0
    do iy = 1, nys
        do ix = 1, nxs
            c0 = c0 + 1
            slon(c0) = -10.0_dp + real(ix-1,dp)*5.0_dp
            slat(c0) =  40.0_dp + real(iy-1,dp)*5.0_dp
        end do
    end do
    call grid_init(gs, name="src", mtype="latlon", units="degrees", &
                   x0=-10.0_dp, dx=5.0_dp, nx=nxs, y0=40.0_dp, dy=5.0_dp, ny=nys)
    call points_init(ps, name="psrc", mtype="latlon", units="degrees", x=slon, y=slat)

    call grid_init(gt, name="tgt", mtype="latlon", units="degrees", &
                   x0=-8.0_dp, dx=5.0_dp, nx=nxt, y0=42.0_dp, dy=5.0_dp, ny=nyt)
    tlon = [-8.0_dp, -3.0_dp,  2.0_dp,  7.0_dp, -5.0_dp,  4.0_dp]
    tlat = [42.0_dp, 47.0_dp, 52.0_dp, 57.0_dp, 50.0_dp, 45.0_dp]
    call points_init(pt, name="ptgt", mtype="latlon", units="degrees", x=tlon, y=tlat)

    vs2 = C;            vs1 = C
    vs2_sp = real(C,sp); vs1_sp = real(C,sp)
    vs2_i  = nint(C);    vs1_i  = nint(C)

    call map_init(m_pg, ps, gt, max_neighbors=8)   ! points -> grid
    call map_init(m_gp, gs, pt, max_neighbors=8)   ! grid   -> points
    call map_init(m_pp, ps, pt, max_neighbors=8)   ! points -> points

    ! ===== points -> grid (2D target) =========================================
    call map_field(m_pg, "v", vs1, vg, method="shepard", missing_value=miss, mask2=mg)
    call chk("pg: constant transfer", all(mg) .and. maxval(abs(vg-C)) < rtol, fails)

    packg = .true.; packg(1,1) = .false.; packg(3,2) = .false.
    call map_field(m_pg, "v", vs1, vg, method="shepard", missing_value=miss, &
                   mask2=mg, mask_pack=packg)
    call chk("pg: mask_pack excludes", (.not. mg(1,1)) .and. (.not. mg(3,2)) .and. &
             vg(1,1) == miss .and. vg(3,2) == miss, fails)
    call chk("pg: mask_pack keeps rest", mg(2,2) .and. abs(vg(2,2)-C) < rtol, fails)

    vg = -1.0_dp
    packg = .true.; packg(4,4) = .false.
    call map_field(m_pg, "v", vs1, vg, method="shepard", missing_value=miss, &
                   mask_pack=packg, reset=.false.)
    call chk("pg: reset=F preserves", abs(vg(4,4)+1.0_dp) < rtol, fails)
    call chk("pg: reset=F maps rest",  abs(vg(1,1)-C)     < rtol, fails)

    ! ===== grid -> points (1D target) =========================================
    call map_field(m_gp, "v", vs2, vp, method="shepard", missing_value=miss, mask2=mp)
    call chk("gp: constant transfer", all(mp) .and. maxval(abs(vp-C)) < rtol, fails)

    packp = .true.; packp(2) = .false.
    call map_field(m_gp, "v", vs2, vp, method="shepard", missing_value=miss, &
                   mask2=mp, mask_pack=packp)
    call chk("gp: mask_pack", (.not. mp(2)) .and. vp(2) == miss .and. &
             abs(vp(1)-C) < rtol, fails)

    vp = -1.0_dp; packp = .true.; packp(3) = .false.
    call map_field(m_gp, "v", vs2, vp, method="shepard", missing_value=miss, &
                   mask_pack=packp, reset=.false.)
    call chk("gp: reset=F preserves", abs(vp(3)+1.0_dp) < rtol .and. abs(vp(1)-C) < rtol, fails)

    ! ===== points -> points (1D target) =======================================
    call map_field(m_pp, "v", vs1, vp, method="shepard", missing_value=miss, mask2=mp)
    call chk("pp: constant transfer", all(mp) .and. maxval(abs(vp-C)) < rtol, fails)

    ! ===== sp / integer wrappers (constant transfer) ==========================
    call map_field(m_pg, "v", vs1_sp, vg_sp, method="shepard", &
                   missing_value=real(miss,sp), mask2=mg)
    call chk("pg sp", all(mg) .and. maxval(abs(vg_sp-real(C,sp))) < 1.0e-5_sp, fails)
    call map_field(m_pg, "v", vs1_i, vg_i, method="shepard", missing_value=nint(miss), mask2=mg)
    call chk("pg int", all(mg) .and. all(vg_i == nint(C)), fails)

    call map_field(m_gp, "v", vs2_sp, vp_sp, method="shepard", &
                   missing_value=real(miss,sp), mask2=mp)
    call chk("gp sp", all(mp) .and. maxval(abs(vp_sp-real(C,sp))) < 1.0e-5_sp, fails)
    call map_field(m_gp, "v", vs2_i, vp_i, method="shepard", missing_value=nint(miss), mask2=mp)
    call chk("gp int", all(mp) .and. all(vp_i == nint(C)), fails)

    call map_field(m_pp, "v", vs1_sp, vp_sp, method="shepard", &
                   missing_value=real(miss,sp), mask2=mp)
    call chk("pp sp", all(mp) .and. maxval(abs(vp_sp-real(C,sp))) < 1.0e-5_sp, fails)
    call map_field(m_pp, "v", vs1_i, vp_i, method="shepard", missing_value=nint(miss), mask2=mp)
    call chk("pp int", all(mp) .and. all(vp_i == nint(C)), fails)

    ! ===== fill / smoothing forwarding: points->grid == grid->grid ============
    ! Build an 8x8 latlon source (and its point-set twin, same order) with a
    ! varying field, and a 12x12 target inside coverage (large enough for the
    ! nr=4 fill stencil and the gaussian kernel).
    c0 = 0
    do iy = 1, nysb
        do ix = 1, nxsb
            c0 = c0 + 1
            slonB(c0) = -14.0_dp + real(ix-1,dp)*4.0_dp
            slatB(c0) =  38.0_dp + real(iy-1,dp)*4.0_dp
            fB2d(ix,iy) = real(ix,dp) + 10.0_dp*real(iy,dp)
        end do
    end do
    fB1d = reshape(fB2d, [nsb])
    call grid_init(gsB, name="srcB", mtype="latlon", units="degrees", &
                   x0=-14.0_dp, dx=4.0_dp, nx=nxsb, y0=38.0_dp, dy=4.0_dp, ny=nysb)
    call points_init(psB, name="psrcB", mtype="latlon", units="degrees", x=slonB, y=slatB)
    call grid_init(gtB, name="tgtB", mtype="latlon", units="degrees", &
                   x0=-12.0_dp, dx=2.0_dp, nx=nxtb, y0=40.0_dp, dy=2.0_dp, ny=nytb)
    call map_init(mGB, gsB, gtB, max_neighbors=8)
    call map_init(mPB, psB, gtB, max_neighbors=8)

    ! Smoothing: gaussian-smoothed grid->grid vs points->grid must be identical.
    call map_field(mGB, "f", fB2d, vA, method="shepard", missing_value=miss, &
                   filt_method="gaussian", filt_par=[1.0_dp, 1.0_dp])
    call map_field(mPB, "f", fB1d, vB, method="shepard", missing_value=miss, &
                   filt_method="gaussian", filt_par=[1.0_dp, 1.0_dp])
    call chk("smooth: points->grid == grid->grid", maxval(abs(vA-vB)) < rtol, fails)

    ! Fill: blank an interior 3x3 source block so a cluster of interior target
    ! cells is missing, then nn-fill. grid->grid and points->grid must agree,
    ! and the fill must actually reduce the missing count.
    fB2d_h = fB2d
    fB2d_h(4:6,4:6) = miss
    fB1d_h = reshape(fB2d_h, [nsb])
    call map_field(mGB, "f", fB2d_h, vBase, method="shepard", missing_value=miss)
    nmiss0 = count(vBase == miss)
    call map_field(mGB, "f", fB2d_h, vA, method="shepard", missing_value=miss, fill_method="nn")
    call map_field(mPB, "f", fB1d_h, vB, method="shepard", missing_value=miss, fill_method="nn")
    call chk("fill: created interior holes", nmiss0 > 0, fails)
    call chk("fill: reduced missing count", count(vA == miss) < nmiss0, fails)
    call chk("fill: points->grid == grid->grid", maxval(abs(vA-vB)) < rtol, fails)

    if (fails > 0) then
        write(*,*) "test_parity FAILED with", fails, "errors"; stop 1
    end if
    write(*,*) "PASS: test_parity (points/grid combos + sp/int wrappers; reset/mask_pack/fill/smooth)"

contains

    subroutine chk(label, cond, nfail)
        character(len=*), intent(in)    :: label
        logical,          intent(in)    :: cond
        integer,          intent(inout) :: nfail
        if (cond) then
            write(*,"(a,a)") "  ok   ", label
        else
            write(*,"(a,a)") "  FAIL ", label
            nfail = nfail + 1
        end if
    end subroutine chk

end program test_parity
