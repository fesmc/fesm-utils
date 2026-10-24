program test_coupler_global
    ! Regional -> global conservative exchange (docs/dev/coupler-global-design.md,
    ! 1 and 8). A polar-stereographic disk domain (a load inside a disk of
    ! radius R, missing outside; centred 500 km from the pole, so it covers the
    ! pole and crosses the 0/360 and 180 meridians) is mapped onto a global
    ! 1 deg lon-lat grid (WGS84) and onto a Gaussian T64 grid (nlat = 96,
    ! spherical).
    !   Global target: sum_k phi_jk = 1 for every source cell (the whole
    !       domain is covered), to CLOSE_SPH on the sphere and CLOSE_ELL on
    !       WGS84 (see there).
    !   (b) frac2 in [0,1]; sum(wfull_dst*frac2) = sum_valid(wcov_src) on an
    !       unpruned map; frac2 = 1 for targets inside the disk, 0 outside.
    !   (f) a pruned and an unpruned map give the same norm="dst" integral
    !       (with areas) to 1e-13.
    !   frac2 present vs absent gives an identical var2 (norm="frac").
    !   norm="dst" without areas = frac2 * mean (map metric).
    ! Through the coupler (registered model areas, 2):
    !   (a) remap norm="dst": coupler_integral on the global grid equals the
    !       covered part of the regional integral to 1e-12, and the whole
    !       regional integral to 1e-12 (sphere) / 1e-8 (WGS84).
    !   (c) two overlapping regional grids (polar stereographic + regional
    !       lon-lat), owners from coupler_owner_reduce: merged integral =
    !       background part + owned regional integrals; merge(A,B) == merge(B,A)
    !       bit for bit; without owners max(frac_tot) > 1 + 1e-2.
    !   (d) an intensive merge of a constant returns the constant to 1e-14.
    !   (e) split with var_src == ref_src returns ref_dst exactly.
    ! Pole and seam geometry (check_pole_seam), "con" to both global grids:
    !   an ANT-like polar-stereographic grid over the pole and Greenland-like
    !   oblique stereographic grids across the 0/360 meridian and across the
    !   dateline: sum_k phi_jk = 1 for every source cell; a constant maps to 1
    !   (norm="frac", 1e-10) and, with norm="dst", to 1 on every target cell
    !   inside the domain (CLOSE_SPH / CLOSE_ELL): no cell is lost or
    !   double-counted along the seam or around the pole.

    use coords
    use coupler

    implicit none

    character(len=*), parameter :: fldr = "maps_coupler_global_test"
    real(dp), parameter :: R_DISK = 1500.0_dp      ! km
    real(dp), parameter :: X_DISK = 300.0_dp       ! km, disk centre (X_DISK, Y_DISK): 500 km from
    real(dp), parameter :: Y_DISK = 400.0_dp       ! the pole; the disk holds the pole and crosses
                                                   ! the 0/360 (+y axis) and 180 (-y) meridians
    ! closure of sum_k phi_jk = 1 on a global target. Sphere: ~1e-12, up to
    ! ~2e-10 for a cell holding the pole (round-off of ~1500 geodesic areas of
    ! polar strip pieces). WGS84: ~1e-8 for 100 km cells (the pieces have
    ! great-circle edges, their geodesic areas on the ellipsoid do not add up
    ! exactly).
    real(dp), parameter :: CLOSE_SPH = 1.0e-9_dp
    real(dp), parameter :: CLOSE_ELL = 1.0e-7_dp
    type(grid_class) :: gps, gll, gga
    integer :: fails

    fails = 0
    call execute_command_line("rm -rf "//fldr)

    call grid_init(gps, name="glob-ant", mtype="polar_stereographic", units="kilometers", planet="WGS84", &
                   lon180=.false., x0=-2000.0_dp, dx=50.0_dp, nx=81, y0=-2000.0_dp, dy=50.0_dp, ny=81, &
                   lambda=0.0_dp, phi=-71.0_dp)
    call grid_init(gll, name="glob-ll1", mtype="latlon", units="degrees", planet="WGS84", &
                   x0=0.5_dp, dx=1.0_dp, nx=360, y0=-89.5_dp, dy=1.0_dp, ny=180)
    ! Gauss bounds from the weights reach the poles (grid_init_gaussian)
    call grid_init_gaussian(gga, "glob-t64", nlon=192, nlat=96, planet="Spherical Earth")

    call check_target(gps, gll)
    call check_target(gps, gga)
    call check_coupler(gps, gll, gga)
    call check_pole_seam(gll, gga)

    call execute_command_line("rm -rf "//fldr)
    if (fails > 0) then
        write(*,*) "FAIL: test_coupler_global (", fails, " failure(s))"
        stop 1
    end if
    write(*,*) "PASS: test_coupler_global"

contains

    subroutine check_target(g1, g2)
        type(grid_class), intent(in) :: g1, g2
        type(map_class) :: mp, mu
        real(dp), allocatable :: v1(:,:), v2(:,:), v2b(:,:), f2(:,:), a1(:,:), a2(:,:), w2(:,:)
        logical,  allocatable :: m2(:,:), inside(:,:)
        real(dp) :: lhs, rhs, ip, iu, emax, rad, fmin, fmax
        real(dp), allocatable :: cov(:)
        integer  :: i, j, nx2, ny2

        nx2 = g2%G%nx; ny2 = g2%G%ny
        call disk_field(g1, v1)
        allocate(v2(nx2,ny2), v2b(nx2,ny2), f2(nx2,ny2), m2(nx2,ny2), w2(nx2,ny2))
        allocate(a1(g1%G%nx,g1%G%ny), a2(nx2,ny2))
        a1 = g1%G%dx*g1%G%dy*1.0e6_dp                 ! plane model area [m^2]
        a2 = g2%area*g2%cs%xy_conv**2                  ! [m^2]

        call map_init(mp, g1, g2, method="con", fldr=fldr)                       ! pruned (default)
        call map_init(mu, g1, g2, method="con", fldr=fldr//"/unpruned", weight_cutoff=0.0_dp)
        write(*,*) trim(g2%name), ": links pruned / unpruned =", mp%wm%n_links, mu%wm%n_links

        ! frac2 present vs absent: identical var2
        call map_field(mp, "v", v1, v2)
        call map_field(mp, "v", v1, v2b, frac2=f2, mask2=m2)
        if (any(v2 /= v2b)) then
            write(*,*) "FAIL: frac2 changes var2 (", trim(g2%name), ")"; fails = fails + 1
        end if

        ! (b) frac2 range, coverage identity (unpruned), inside / outside
        call map_field(mu, "v", v1, v2b, frac2=f2, mask2=m2)
        fmin = minval(f2); fmax = maxval(f2)
        lhs = sum(reshape(mu%wm%wfull_dst, [nx2,ny2])*f2)
        rhs = sum(pack(mu%wm%wcov_src, reshape(v1, [size(v1)]) /= mv_dp))
        write(*,*) "  frac2 range =", fmin, fmax, " coverage rel err =", abs(lhs-rhs)/rhs
        if (fmin < 0.0_dp .or. fmax > 1.0_dp) then
            write(*,*) "FAIL: frac2 outside [0,1]"; fails = fails + 1
        end if
        if (abs(lhs - rhs) > 1.0e-10_dp*rhs) then
            write(*,*) "FAIL: sum(wfull_dst*frac2) /= sum_valid(wcov_src)"; fails = fails + 1
        end if
        allocate(inside(nx2,ny2))
        emax = 0.0_dp
        do j = 1, ny2
            do i = 1, nx2
                rad = disk_radius(g1, g2%lon(i,j), g2%lat(i,j))
                if (rad < R_DISK - 300.0_dp) then
                    emax = max(emax, abs(f2(i,j) - 1.0_dp))
                else if (rad > R_DISK + 300.0_dp .and. f2(i,j) /= 0.0_dp) then
                    write(*,*) "FAIL: frac2 /= 0 outside the disk at", i, j, f2(i,j); fails = fails + 1
                    exit
                end if
            end do
        end do
        write(*,*) "  max |frac2 - 1| inside the disk =", emax
        if (emax > close_tol(g2)) then
            write(*,*) "FAIL: frac2 /= 1 inside the disk"; fails = fails + 1
        end if

        ! norm="dst" without areas = frac2 * (covered mean), where frac2 < 1;
        ! on a fully covered target the coverage may exceed 1 by the closure
        ! error of sum_k phi_jk, which frac2 clips and norm="dst" keeps
        call map_field(mu, "v", v1, v2, norm="dst")
        call map_field(mu, "v", v1, v2b)
        emax = maxval(abs(v2 - f2*v2b), mask=m2 .and. f2 < 1.0_dp)
        write(*,*) "  max |dst - frac2*mean| (frac2 < 1) =", emax
        if (emax > 1.0e-12_dp*maxval(abs(v1))) then
            write(*,*) "FAIL: norm=dst (no areas) /= frac2*mean"; fails = fails + 1
        end if
        emax = maxval(abs(v2/v2b - 1.0_dp), mask=m2 .and. f2 == 1.0_dp)
        write(*,*) "  max |dst/mean - 1| (frac2 = 1) =", emax
        if (emax > close_tol(g2)) then
            write(*,*) "FAIL: norm=dst (no areas) /= mean on covered targets"; fails = fails + 1
        end if

        ! (f) pruned vs unpruned: same model-metric integral (norm="dst" with areas)
        call map_field(mp, "v", v1, v2, norm="dst", area1=a1, area2=a2)
        ip = sum(a2*v2, mask=v2 /= mv_dp)
        call map_field(mu, "v", v1, v2, norm="dst", area1=a1, area2=a2)
        iu = sum(a2*v2, mask=v2 /= mv_dp)
        ! source amount on the covered part: sum_j a1_j v_j sum_k phi_jk
        allocate(cov(size(v1)))
        cov = 0.0_dp
        where (mu%wm%wfull_src > 0.0_dp) cov = mu%wm%wcov_src/mu%wm%wfull_src
        rhs = sum(reshape(a1*v1, [size(v1)])*cov, mask=reshape(v1, [size(v1)]) /= mv_dp)
        write(*,*) "  integral pruned / unpruned / covered source =", ip, iu, rhs
        write(*,*) "  covered fraction of the domain - 1 =", rhs/sum(a1*v1, mask=v1 /= mv_dp) - 1.0_dp
        if (abs(rhs/sum(a1*v1, mask=v1 /= mv_dp) - 1.0_dp) > close_tol(g2)) then
            write(*,*) "FAIL: covered fraction of a global target /= 1"; fails = fails + 1
        end if
        if (abs(ip - iu) > 1.0e-13_dp*abs(iu)) then
            write(*,*) "FAIL: pruned and unpruned integrals differ"; fails = fails + 1
        end if
        if (abs(iu - rhs) > 1.0e-12_dp*abs(rhs)) then
            write(*,*) "FAIL: norm=dst integral /= covered source integral"; fails = fails + 1
        end if

        ! norm="sum" with area1: the amounts per target add up to the covered
        ! source amount
        call map_field(mp, "v", v1, w2, norm="sum", area1=a1)
        ip = sum(w2, mask=w2 /= mv_dp)
        if (abs(ip - rhs) > 1.0e-12_dp*abs(rhs)) then
            write(*,*) "FAIL: norm=sum amounts /= covered source amount", ip, rhs; fails = fails + 1
        end if
    end subroutine check_target

    subroutine check_coupler(gps, gll, gga)
        type(grid_class), intent(in) :: gps, gll, gga
        type(coupler_class) :: cpl
        type(coupler_merge_t) :: acc
        type(grid_class) :: gb
        real(dp), allocatable :: va(:,:), vb(:,:), v2(:,:), oa(:,:), ob(:,:), bg(:,:), m_ab(:,:), m_ba(:,:)
        real(dp), allocatable :: ftot(:,:), c2(:,:), ref_s(:,:), ref_d(:,:), sp_d(:,:)
        real(dp) :: lonc, latc, i_glob, i_cov, i_bg, i_own, i_merged, i_reg
        integer  :: im, k

        call coupler_init(cpl, map_fldr=fldr//"/cpl")
        call coupler_add_grid(cpl, "glob-ant", gps)
        call coupler_add_grid(cpl, "glob-ll1", gll)
        call coupler_add_grid(cpl, "glob-t64", gga)
        call disk_field(gps, va)

        ! (a) regional -> global, norm="dst" with the registered areas
        do k = 1, 2
            if (k == 1) then
                call remap(cpl, va, "glob-ant", v2, "glob-ll1", method="con", norm="dst")
                i_glob = coupler_integral(cpl, "glob-ll1", v2, missing_value=mv_dp, describe=.true.)
                im = map_index(cpl, "glob-ant", "glob-ll1")
            else
                call remap(cpl, va, "glob-ant", v2, "glob-t64", method="con", norm="dst")
                i_glob = coupler_integral(cpl, "glob-t64", v2, missing_value=mv_dp, describe=.true.)
                im = map_index(cpl, "glob-ant", "glob-t64")
            end if
            i_cov = covered_integral(cpl, im, va, "glob-ant")
            i_reg = coupler_integral(cpl, "glob-ant", va, missing_value=mv_dp, describe=.true.)
            write(*,*) "(a) global / covered regional integral, rel diff =", i_glob, i_cov, &
                       abs(i_glob - i_cov)/abs(i_cov), "  regional total =", i_reg, &
                       " rel diff =", abs(i_glob - i_reg)/abs(i_reg)
            if (abs(i_glob - i_cov) > 1.0e-12_dp*abs(i_cov)) then
                write(*,*) "FAIL: (a) global integral /= covered regional integral"; fails = fails + 1
            end if
            ! the global target covers the whole domain: sum_k phi_jk = 1
            if (abs(i_glob - i_reg) > merge(1.0e-12_dp, 1.0e-8_dp, k == 2)*abs(i_reg)) then
                write(*,*) "FAIL: (a) global integral /= regional integral"; fails = fails + 1
            end if
        end do

        ! second regional domain: lon-lat box north of the disk centre, partly
        ! inside the disk; its 0.25 deg cells nest in the 1 deg global cells, so
        ! the owner split (per B cell) matches on the global grid
        call oblimap_projection_inverse(X_DISK*1000.0_dp, Y_DISK*1000.0_dp, lonc, latc, gps%cs%proj)
        lonc = real(floor(lonc), dp); latc = real(floor(latc), dp)
        call grid_init(gb, name="glob-box", mtype="latlon", units="degrees", planet="WGS84", &
                       x0=lonc - 10.0_dp + 0.125_dp, dx=0.25_dp, nx=80, &
                       y0=latc + 6.0_dp + 0.125_dp, dy=0.25_dp, ny=48)
        call coupler_add_grid(cpl, "glob-box", gb)
        allocate(vb(gb%G%nx, gb%G%ny))
        vb = 700.0_dp + 2.0_dp*(gb%x - lonc)

        ! owners: A owns its valid disk; B owns what A leaves
        allocate(oa(gps%G%nx, gps%G%ny), ob(gb%G%nx, gb%G%ny))
        oa = merge(1.0_dp, 0.0_dp, va /= mv_dp)
        ob = 1.0_dp
        call coupler_owner_reduce(cpl, "glob-ant", oa, "glob-box", ob)
        write(*,*) "(c) owner_reduce: B owner range", minval(ob), maxval(ob), " cells owned by A (ob < 0.01):", &
                   count(ob < 0.01_dp), " of", size(ob)
        if (minval(ob) > 0.01_dp .or. maxval(ob) < 0.99_dp) then
            write(*,*) "FAIL: (c) the domains do not partly overlap"; fails = fails + 1
        end if
        allocate(bg(gll%G%nx, gll%G%ny), ftot(gll%G%nx, gll%G%ny))
        bg = 10.0_dp

        ! (c) merge A then B, and B then A
        call coupler_merge_begin(acc, cpl, "glob-ll1")
        call coupler_merge_add(acc, cpl, "glob-ant", va, owner=oa)
        call coupler_merge_add(acc, cpl, "glob-box", vb, owner=ob)
        call coupler_merge_end(acc, cpl, bg, m_ab, frac_tot=ftot)
        call coupler_merge_begin(acc, cpl, "glob-ll1")
        call coupler_merge_add(acc, cpl, "glob-box", vb, owner=ob)
        call coupler_merge_add(acc, cpl, "glob-ant", va, owner=oa)
        call coupler_merge_end(acc, cpl, bg, m_ba)
        if (any(m_ab /= m_ba)) then
            write(*,*) "FAIL: (c) merge(A,B) /= merge(B,A)"; fails = fails + 1
        end if
        i_merged = coupler_integral(cpl, "glob-ll1", m_ab)
        i_bg     = coupler_integral(cpl, "glob-ll1", (1.0_dp - min(ftot, 1.0_dp))*bg)
        i_own    = covered_integral(cpl, map_index(cpl, "glob-ant", "glob-ll1"), merge(oa*va, mv_dp, va /= mv_dp), &
                                    "glob-ant") &
                 + covered_integral(cpl, map_index(cpl, "glob-box", "glob-ll1"), ob*vb, "glob-box")
        write(*,*) "(c) merged / background + owned, rel diff =", i_merged, i_bg + i_own, &
                   abs(i_merged - i_bg - i_own)/abs(i_merged), "  max frac_tot =", maxval(ftot)
        if (abs(i_merged - i_bg - i_own) > 1.0e-12_dp*abs(i_merged)) then
            write(*,*) "FAIL: (c) merged integral /= background + owned integrals"; fails = fails + 1
        end if
        ! without owners the domains double-count
        call coupler_merge_begin(acc, cpl, "glob-ll1")
        call coupler_merge_add(acc, cpl, "glob-ant", va)
        call coupler_merge_add(acc, cpl, "glob-box", vb)
        call coupler_merge_end(acc, cpl, bg, v2, frac_tot=ftot, overlap_tol=10.0_dp)
        write(*,*) "(c) without owners: max frac_tot =", maxval(ftot)
        if (maxval(ftot) <= 1.01_dp) then
            write(*,*) "FAIL: (c) overlapping domains without owners not detected"; fails = fails + 1
        end if

        ! (d) intensive merge of a constant
        allocate(c2(gll%G%nx, gll%G%ny))
        call coupler_merge_begin(acc, cpl, "glob-ll1", kind="intensive")
        call coupler_merge_add(acc, cpl, "glob-ant", merge(3.5_dp, mv_dp, va /= mv_dp), owner=oa)
        call coupler_merge_add(acc, cpl, "glob-box", vb*0.0_dp + 3.5_dp, owner=ob)
        call coupler_merge_end(acc, cpl, bg*0.0_dp + 3.5_dp, c2)
        write(*,*) "(d) intensive merge of 3.5: max error =", maxval(abs(c2 - 3.5_dp))
        if (maxval(abs(c2 - 3.5_dp)) > 1.0e-14_dp*3.5_dp) then
            write(*,*) "FAIL: (d) intensive merge does not preserve a constant"; fails = fails + 1
        end if

        ! (e) split with var_src == ref_src returns ref_dst exactly
        allocate(ref_s(gll%G%nx, gll%G%ny), ref_d(gps%G%nx, gps%G%ny))
        ref_s = 5.0_dp + gll%y*0.01_dp
        ref_d = 1.0_dp + gps%x*0.001_dp
        call coupler_split(cpl, "glob-ll1", ref_s, "glob-ant", sp_d, ref_src=ref_s, ref_dst=ref_d)
        if (any(sp_d /= ref_d)) then
            write(*,*) "FAIL: (e) split of the reference /= ref_dst"; fails = fails + 1
        else
            write(*,*) "(e) split(var_src == ref_src) == ref_dst"
        end if
    end subroutine check_coupler

    real(dp) function close_tol(g2) result(tol)
        ! Closure tolerance of sum_k phi_jk = 1 on the target planet of g2.
        type(grid_class), intent(in) :: g2
        tol = merge(CLOSE_SPH, CLOSE_ELL, g2%cs%planet%is_sphere)
    end function close_tol

    subroutine check_pole_seam(gll, gga)
        ! Projected cells over the pole and across the 0/360 meridian / the
        ! dateline onto global lon-lat (WGS84) and Gauss (sphere) targets.
        type(grid_class), intent(in) :: gll, gga
        type(grid_class) :: gant, ggrl
        real(dp) :: lam(2) = [0.0_dp, 180.0_dp]
        integer  :: k
        call grid_init(gant, name="pole-ant", mtype="polar_stereographic", units="kilometers", &
                       planet="WGS84", lon180=.false., x0=-3000.0_dp, dx=100.0_dp, nx=61, &
                       y0=-3000.0_dp, dy=100.0_dp, ny=61, lambda=0.0_dp, phi=-71.0_dp)
        call check_closure(gant, gll)
        call check_closure(gant, gga)
        do k = 1, 2
            call grid_init(ggrl, name="seam-grl"//merge("000", "180", k == 1), mtype="stereographic", &
                           units="kilometers", planet="WGS84", x0=-800.0_dp, dx=40.0_dp, nx=41, &
                           y0=-1500.0_dp, dy=40.0_dp, ny=76, lambda=lam(k), phi=72.0_dp, alpha=7.5_dp)
            call check_closure(ggrl, gll)
            call check_closure(ggrl, gga)
        end do
    end subroutine check_pole_seam

    subroutine check_closure(g1, g2)
        ! A global target covers every source cell once: sum_k phi_jk = 1;
        ! a constant maps to 1 (norm="frac") and its norm="dst" map is 1 on
        ! every target cell lying inside the source domain.
        type(grid_class), intent(in) :: g1, g2
        type(map_class) :: mp
        real(dp), allocatable :: one(:,:), v2(:,:)
        logical,  allocatable :: m2(:,:)
        real(dp) :: tol, emax_src, emax_frac, emax_dst
        integer  :: i, j, nin

        tol = close_tol(g2)
        call map_init(mp, g1, g2, method="con", fldr=fldr)
        emax_src = maxval(abs(mp%wm%wcov_src/mp%wm%wfull_src - 1.0_dp))

        allocate(one(g1%G%nx, g1%G%ny), v2(g2%G%nx, g2%G%ny), m2(g2%G%nx, g2%G%ny))
        one = 1.0_dp
        call map_field(mp, "one", one, v2, mask2=m2)
        emax_frac = maxval(abs(v2 - 1.0_dp), mask=m2)
        call map_field(mp, "one", one, v2, norm="dst", mask2=m2)
        emax_dst = 0.0_dp
        nin = 0
        do j = 1, g2%G%ny
            do i = 1, g2%G%nx
                if (.not. cell_in_domain(g1, g2, i, j)) cycle
                nin = nin + 1
                emax_dst = max(emax_dst, abs(v2(i,j) - 1.0_dp))
            end do
        end do
        write(*,"(1x,a,' -> ',a,': max|1-sum phi| =',es9.2,'  const frac err =',es9.2, &
             & '  const dst err on',i6,' inner targets =',es9.2)") trim(g1%name), trim(g2%name), &
             emax_src, emax_frac, nin, emax_dst
        if (any(mp%wm%wfull_src <= 0.0_dp)) then
            write(*,*) "FAIL: source cells without links on a global target"; fails = fails + 1
        end if
        if (emax_src > tol) then
            write(*,*) "FAIL: sum_k phi_jk /= 1 on a global target"; fails = fails + 1
        end if
        if (emax_frac > 1.0e-10_dp) then
            write(*,*) "FAIL: a constant does not map to 1 (norm=frac)"; fails = fails + 1
        end if
        if (nin == 0 .or. emax_dst > tol) then
            write(*,*) "FAIL: norm=dst constant /= 1 on targets inside the domain"; fails = fails + 1
        end if
    end subroutine check_closure

    logical function cell_in_domain(g1, g2, i, j) result(inside)
        ! Target cell (i,j) of lon-lat g2 lies inside the projected domain of
        ! g1: its boundary, sampled 4 points per edge, projects into g1's cells.
        type(grid_class), intent(in) :: g1, g2
        integer,          intent(in) :: i, j
        real(dp) :: lo, la, x, y, t, xb(4), yb(4)
        integer  :: e, s
        xb = [g2%G%x_bnds(i-1), g2%G%x_bnds(i), g2%G%x_bnds(i), g2%G%x_bnds(i-1)]
        yb = [g2%G%y_bnds(j-1), g2%G%y_bnds(j-1), g2%G%y_bnds(j), g2%G%y_bnds(j)]
        yb = max(-90.0_dp, min(90.0_dp, yb))
        inside = .false.
        do e = 1, 4
            do s = 0, 3
                t  = real(s, dp)/4.0_dp
                lo = xb(e) + t*(xb(mod(e,4)+1) - xb(e))
                la = yb(e) + t*(yb(mod(e,4)+1) - yb(e))
                if (la*g1%cs%proj%phi <= 0.0_dp) return        ! other hemisphere
                call oblimap_projection(lo, la, x, y, g1%cs%proj)
                x = x/g1%cs%xy_conv; y = y/g1%cs%xy_conv
                if (x < g1%G%x_bnds(0) .or. x > g1%G%x_bnds(g1%G%nx) .or. &
                    y < g1%G%y_bnds(0) .or. y > g1%G%y_bnds(g1%G%ny)) return
            end do
        end do
        inside = .true.
    end function cell_in_domain

    integer function map_index(cpl, src, dst) result(im)
        type(coupler_class), intent(in) :: cpl
        character(len=*),    intent(in) :: src, dst
        do im = 1, cpl%nmaps
            if (trim(cpl%maps(im)%src) == src .and. trim(cpl%maps(im)%dst) == dst .and. &
                trim(cpl%maps(im)%method) == "con") return
        end do
        im = -1
        write(*,*) "FAIL: no con map ", src, " -> ", dst
        stop 1
    end function map_index

    real(dp) function covered_integral(cpl, im, v, src) result(tot)
        ! sum_j area_j v_j sum_k phi_jk: the source integral over the part of
        ! each source cell covered by the target (registered areas).
        type(coupler_class), intent(in) :: cpl
        integer,             intent(in) :: im
        real(dp),            intent(in) :: v(:,:)
        character(len=*),    intent(in) :: src
        real(dp), allocatable :: cov(:,:)
        integer :: ig
        associate (wm => cpl%maps(im)%map%wm)
            allocate(cov(size(v,1), size(v,2)))
            cov = 0.0_dp
            where (reshape(wm%wfull_src, shape(v)) > 0.0_dp) &
                cov = reshape(wm%wcov_src, shape(v))/reshape(wm%wfull_src, shape(v))
        end associate
        do ig = 1, cpl%ngrids
            if (trim(cpl%grids(ig)%name) == src) exit
        end do
        tot = sum(cpl%grids(ig)%area*v*cov, mask=v /= mv_dp)
    end function covered_integral

    subroutine disk_field(g, v)
        ! Load inside the disk (smooth, positive), missing outside.
        type(grid_class), intent(in) :: g
        real(dp), allocatable, intent(out) :: v(:,:)
        integer :: i, j
        real(dp) :: r
        allocate(v(g%G%nx, g%G%ny))
        do j = 1, g%G%ny
            do i = 1, g%G%nx
                r = sqrt((g%x(i,j) - X_DISK)**2 + (g%y(i,j) - Y_DISK)**2)
                if (r <= R_DISK) then
                    v(i,j) = 1000.0_dp + 500.0_dp*cos(0.5_dp*3.14159_dp*r/R_DISK) + 0.1_dp*g%x(i,j)
                else
                    v(i,j) = mv_dp
                end if
            end do
        end do
    end subroutine disk_field

    real(dp) function disk_radius(g, lon, lat) result(r)
        ! Distance [km] of (lon,lat) from the disk centre, in the source plane.
        type(grid_class), intent(in) :: g
        real(dp), intent(in) :: lon, lat
        real(dp) :: x, y
        if (lat > -40.0_dp) then
            r = huge(1.0_dp)
            return
        end if
        call oblimap_projection(lon, lat, x, y, g%cs%proj)
        r = sqrt((x/1000.0_dp - X_DISK)**2 + (y/1000.0_dp - Y_DISK)**2)
    end function disk_radius

end program test_coupler_global
