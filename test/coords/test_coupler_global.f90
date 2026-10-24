program test_coupler_global
    ! Regional -> global conservative exchange (docs/dev/coupler-global-design.md,
    ! 1 and 8). A polar-stereographic disk domain (a load inside a disk of
    ! radius R, missing outside; centred 2500 km from the pole on 90E) is mapped
    ! onto a global 1 deg lon-lat grid and onto a Gaussian T64 grid (nlat = 96).
    ! The disk avoids the pole and the 0/360 meridian: projected source cells
    ! that straddle that meridian are mis-built by conservative_spherical
    ! (corner longitudes interpolated the long way round; a pre-existing
    ! issue, not changed here).
    !   (b) frac2 in [0,1]; sum(wfull_dst*frac2) = sum_valid(wcov_src) on an
    !       unpruned map; frac2 = 1 for targets inside the disk, 0 outside.
    !   (f) a pruned and an unpruned map give the same norm="dst" integral
    !       (with areas) to 1e-13.
    !   frac2 present vs absent gives an identical var2 (norm="frac").
    !   norm="dst" without areas = frac2 * mean (map metric).
    ! Through the coupler (registered model areas, 2):
    !   (a) remap norm="dst": coupler_integral on the global grid equals the
    !       covered part of the regional integral to 1e-12.
    !   (c) two overlapping regional grids (polar stereographic + regional
    !       lon-lat), owners from coupler_owner_reduce: merged integral =
    !       background part + owned regional integrals; merge(A,B) == merge(B,A)
    !       bit for bit; without owners max(frac_tot) > 1 + 1e-2.
    !   (d) an intensive merge of a constant returns the constant to 1e-14.
    !   (e) split with var_src == ref_src returns ref_dst exactly.

    use coords
    use coupler

    implicit none

    character(len=*), parameter :: fldr = "maps_coupler_global_test"
    real(dp), parameter :: R_DISK = 1500.0_dp      ! km
    real(dp), parameter :: X_DISK = 2500.0_dp      ! km, disk centre (X_DISK, 0): >= 1000 km from
                                                   ! the pole, away from the 0/360 meridian
    type(grid_class) :: gps, gll, gga
    real(dp), allocatable :: glat(:)
    integer :: fails

    fails = 0
    call execute_command_line("rm -rf "//fldr)

    call grid_init(gps, name="glob-ant", mtype="polar_stereographic", units="kilometers", planet="WGS84", &
                   lon180=.false., x0=500.0_dp, dx=50.0_dp, nx=81, y0=-2000.0_dp, dy=50.0_dp, ny=81, &
                   lambda=0.0_dp, phi=-71.0_dp)
    call grid_init(gll, name="glob-ll1", mtype="latlon", units="degrees", planet="WGS84", &
                   x0=0.5_dp, dx=1.0_dp, nx=360, y0=-89.5_dp, dy=1.0_dp, ny=180)
    allocate(glat(96))
    call gaussian_latitudes_calc(96, glat)
    call grid_init(gga, name="glob-t64", mtype="gaussian", units="degrees", planet="Spherical Earth", &
                   x0=0.0_dp, dx=1.875_dp, nx=192, y=glat)

    call check_target(gps, gll)
    call check_target(gps, gga)
    call check_coupler(gps, gll, gga)

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
        ! Tolerance: conservative_spherical clips against the target cell, whose
        ! poleward parallel edge (NSUB great-circle chords) is concave; the
        ! Sutherland-Hodgman clip then loses a sliver (~0.2% at 70S for 1 deg
        ! cells; a pre-existing limitation of the generator, not changed here).
        if (emax > 1.0e-2_dp) then
            write(*,*) "FAIL: frac2 /= 1 inside the disk"; fails = fails + 1
        end if

        ! norm="dst" without areas = frac2 * (covered mean)
        call map_field(mu, "v", v1, v2, norm="dst")
        call map_field(mu, "v", v1, v2b)
        emax = maxval(abs(v2 - f2*v2b), mask=m2)
        write(*,*) "  max |dst - frac2*mean| =", emax
        if (emax > 1.0e-12_dp*maxval(abs(v1))) then
            write(*,*) "FAIL: norm=dst (no areas) /= frac2*mean"; fails = fails + 1
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
        write(*,*) "  covered fraction of the domain =", rhs/sum(a1*v1, mask=v1 /= mv_dp)
        ! 1 for a global target up to the clip sliver loss (see above)
        if (rhs/sum(a1*v1, mask=v1 /= mv_dp) < 0.99_dp) then
            write(*,*) "FAIL: covered fraction of a global target < 0.99"; fails = fails + 1
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
                       abs(i_glob - i_cov)/abs(i_cov), "  regional total =", i_reg
            if (abs(i_glob - i_cov) > 1.0e-12_dp*abs(i_cov)) then
                write(*,*) "FAIL: (a) global integral /= covered regional integral"; fails = fails + 1
            end if
        end do

        ! second regional domain: lon-lat box north of the disk centre, partly
        ! inside the disk; its 0.25 deg cells nest in the 1 deg global cells, so
        ! the owner split (per B cell) matches on the global grid
        call oblimap_projection_inverse(X_DISK*1000.0_dp, 0.0_dp, lonc, latc, gps%cs%proj)
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
                r = sqrt((g%x(i,j) - X_DISK)**2 + g%y(i,j)**2)
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
        r = sqrt((x/1000.0_dp - X_DISK)**2 + (y/1000.0_dp)**2)
    end function disk_radius

end program test_coupler_global
