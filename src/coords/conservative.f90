module conservative
    ! Analytic conservative remapping via polygon clipping (improvement-doc §2,
    ! stage B: planar Sutherland-Hodgman). For each target cell, candidate source
    ! cells are found with a k-d-tree radius query, then each source cell polygon
    ! is clipped against the target cell and the overlap area becomes the link
    ! weight. The result is a MAP_WEIGHT applied as an area-weighted mean.
    !
    ! Two planar modes are supported:
    !   - same-system: source and target share a Cartesian/projected plane;
    !     cells are clipped directly as axis-aligned rectangles.
    !   - lat-lon source -> projected target: each source cell edge is subsampled
    !     into NSUB segments and projected into the target plane, so the clipped
    !     source polygon follows the curved parallels/meridians (chording straight
    !     across them biases the overlap area, growing with source cell size and
    !     projection distortion; subsampling brings a 1deg source to sub-metre
    !     agreement with cdo, matching a fine source).
    ! Lat-lon / Gaussian *targets* need great-circle clipping (conservative_spherical).

    use precision,       only: dp
    use constants, only: degrees_to_radians
    use coordinates,     only: grid_class, compare_coord, axis_midpoint_bounds
    use oblimap_projection_module, only: oblimap_projection, oblimap_projection_inverse
    use planet,          only: planet_area
    use kdtree
    use weight_map
    use polygons
    !$ use omp_lib

    implicit none
    private

    integer, parameter :: NSUB = 8   ! great-circle segments per parallel edge (stage C)

    ! Per-thread link accumulator: each thread appends its target cells' links
    ! here (no shared counter), then bufs_to_wm merges them into the CSR weight
    ! store in destination order.
    type :: link_buf
        integer,  allocatable :: src(:), dst(:)
        real(dp), allocatable :: w(:)
        integer :: n = 0
    end type link_buf

    public :: conservative_weights

contains

    subroutine buf_push(b, s, d, wv)
        ! Append one link to a per-thread buffer, doubling capacity as needed.
        type(link_buf), intent(inout) :: b
        integer,        intent(in)    :: s, d
        real(dp),       intent(in)    :: wv
        integer,  allocatable :: ti(:)
        real(dp), allocatable :: tr(:)
        integer :: cap
        if (.not. allocated(b%src)) allocate(b%src(1024), b%dst(1024), b%w(1024))
        cap = size(b%src)
        if (b%n >= cap) then
            allocate(ti(2*cap)); ti(1:b%n) = b%src(1:b%n); call move_alloc(ti, b%src)
            allocate(ti(2*cap)); ti(1:b%n) = b%dst(1:b%n); call move_alloc(ti, b%dst)
            allocate(tr(2*cap)); tr(1:b%n) = b%w(1:b%n);   call move_alloc(tr, b%w)
        end if
        b%n = b%n + 1
        b%src(b%n) = s; b%dst(b%n) = d; b%w(b%n) = wv
    end subroutine buf_push

    subroutine bufs_to_wm(wm, bufs, n1, n2)
        ! Merge per-thread link buffers into the weight store, ordered by
        ! destination (a counting sort over dst), as weight_map_index expects.
        type(weight_map_t), intent(inout) :: wm
        type(link_buf),     intent(in)    :: bufs(0:)
        integer,            intent(in)    :: n1, n2
        integer :: t, c, d, p, nl, nt
        integer, allocatable :: off(:), pos(:)
        nt = size(bufs)
        nl = 0
        do t = 0, nt-1
            nl = nl + bufs(t)%n
        end do
        call weight_map_alloc(wm, MAP_WEIGHT, n1, n2, nl)
        allocate(off(n2+1)); off = 0
        do t = 0, nt-1
            do c = 1, bufs(t)%n
                d = bufs(t)%dst(c); off(d+1) = off(d+1) + 1
            end do
        end do
        off(1) = 1
        do d = 1, n2
            off(d+1) = off(d+1) + off(d)
        end do
        allocate(pos(n2)); pos(1:n2) = off(1:n2)
        do t = 0, nt-1
            do c = 1, bufs(t)%n
                d = bufs(t)%dst(c)
                p = pos(d); pos(d) = p + 1
                wm%src(p) = bufs(t)%src(c); wm%dst(p) = d; wm%w(p) = bufs(t)%w(c)
            end do
        end do
        deallocate(off, pos)
        call weight_map_index(wm)
    end subroutine bufs_to_wm

    subroutine conservative_weights(wm, grid1, grid2)
        ! Build the conservative MAP_WEIGHT weight store for grid1 -> grid2 via
        ! analytic polygon clipping. Dispatch: planar clip for Cartesian/projected
        ! targets, spherical clip (great-circle arcs) for lat-lon / Gaussian
        ! targets. The map_class target metadata is assembled by the caller
        ! (mapping%map_init_conservative); this routine only fills the weight store,
        ! which keeps `conservative` independent of `mapping` (no circular use).
        type(weight_map_t), intent(inout) :: wm
        type(grid_class),   intent(in)    :: grid1, grid2
        if (grid2%cs%is_cartesian) then
            ! target is Cartesian/projected: planar clip. Use the analytic separable
            ! overlap when source and target share one Cartesian plane (both are
            ! axis-aligned rectangles); otherwise the general Sutherland-Hodgman
            ! clip (e.g. lat-lon source -> projected target).
            if (same_cartesian_system(grid1, grid2)) then
                call conservative_planar_reg(wm, grid1, grid2)
            else
                call conservative_planar(wm, grid1, grid2)
            end if
        else if (is_latlon(grid1) .and. is_latlon(grid2)) then
            ! both lat-lon/Gaussian: axis-aligned in (lon,lat), so the cell overlap
            ! is separable and analytic -- no great-circle polygon clipping needed.
            call conservative_latlon(wm, grid1, grid2)
        else
            call conservative_spherical(wm, grid1, grid2)
        end if
    end subroutine conservative_weights

    logical function is_latlon(grid)
        ! A lat-lon / Gaussian grid: axes are (lon,lat) in degrees, neither a
        ! Cartesian nor a projected system (both flags false; see coordinates).
        type(grid_class), intent(in) :: grid
        is_latlon = (.not. grid%cs%is_cartesian) .and. (.not. grid%cs%is_projection)
    end function is_latlon

    logical function same_cartesian_system(grid1, grid2)
        ! True when both grids live in the same Cartesian/projected plane, so their
        ! cells are axis-aligned rectangles in shared xy coordinates.
        type(grid_class), intent(in) :: grid1, grid2
        same_cartesian_system = compare_coord(grid1, grid2) .and. grid2%cs%is_cartesian &
                                .and. grid1%cs%is_cartesian
    end function same_cartesian_system

    subroutine conservative_planar(wm, grid1, grid2)
        type(weight_map_t), intent(inout) :: wm
        type(grid_class),   intent(in)    :: grid1, grid2

        type(kdtree_t) :: tree
        logical  :: same_sys, do_project
        integer  :: n1, n2, nx1, ny1, nx2, ny2
        integer  :: i, j, ic, jc, c, cs, cc, k, nf, no, nthreads, tid
        integer  :: e, en, s, nsv
        real(dp) :: rad, diag2, area, maxsrcdiag, d, xyc, px, py
        real(dp) :: cx, cy, cz, srcdeg, tgtdeg, tt, plon, plat, glon, glat, src_xyc
        real(dp) :: lx(4), ly(4), tcx(4), tcy(4), qpt(2)
        real(dp), allocatable :: emb(:,:), ox(:), oy(:)
        real(dp), allocatable :: scornx(:,:), scorny(:,:)   ! source cell boundary in target plane
        integer,  allocatable :: cidx(:)
        real(dp), allocatable :: cd2(:)
        logical,  allocatable :: need(:)                    ! source cell is a candidate of some target
        real(dp), allocatable :: wfull_src(:)               ! full source-cell area in the target plane
        type(link_buf), allocatable :: bufs(:)

        same_sys   = compare_coord(grid1, grid2) .and. grid2%cs%is_cartesian
        do_project = (.not. same_sys) .and. grid2%cs%is_projection .and. (.not. grid1%cs%is_cartesian)

        if (.not. (same_sys .or. do_project)) then
            write(*,*) "conservative:: map_init_conservative: stage-B planar clip supports"
            write(*,*) "  (a) same Cartesian/projected system, or"
            write(*,*) "  (b) lat-lon source -> projected target."
            write(*,*) "  Lat-lon/Gaussian targets need the spherical clip (stage C)."
            stop
        end if

        nx1 = grid1%G%nx; ny1 = grid1%G%ny; n1 = nx1*ny1
        nx2 = grid2%G%nx; ny2 = grid2%G%ny; n2 = nx2*ny2
        xyc = grid2%cs%xy_conv

        ! Express the source cell corners in the TARGET plane (the planar clip
        ! runs there in both modes). The candidate-search tree differs by mode:
        !   - same-system: planar tree over projected source centers.
        !   - cross-system (lat-lon -> projection): tree over source centers on the
        !     UNIT SPHERE. Projecting source centers would distort enormously far
        !     from the projection origin (stereographic diverges near the antipode),
        !     inflating the search radius until every target pulls nearly every
        !     source and the k-d tree degrades to O(n_src*n_tgt). The sphere is
        !     undistorted, so the candidate set per target stays O(1).
        !     (see docs/coords-performance.md)
        ! same-system cells are exact rectangles (4 corners); projected source
        ! cells trace their curved lon/lat edges with NSUB samples per edge.
        if (do_project) then
            nsv = 4*NSUB
        else
            nsv = 4
        end if
        allocate(scornx(nsv, n1), scorny(nsv, n1), wfull_src(n1))
        wfull_src = 0.0_dp
        if (same_sys) then
            allocate(emb(2, n1))
        else
            allocate(emb(3, n1))
        end if
        maxsrcdiag = 0.0_dp
        do jc = 1, ny1
            do ic = 1, nx1
                cs = (jc-1)*nx1 + ic
                if (same_sys) then
                    call cell_corners_grid(grid1, ic, jc, lx, ly)
                    scornx(:,cs) = lx; scorny(:,cs) = ly
                    wfull_src(cs) = polygon_area(lx, ly)
                    emb(1,cs) = grid1%x(ic,jc); emb(2,cs) = grid1%y(ic,jc)
                    d = sqrt((maxval(scornx(:,cs))-minval(scornx(:,cs)))**2 &
                           + (maxval(scorny(:,cs))-minval(scorny(:,cs)))**2)
                    maxsrcdiag = max(maxsrcdiag, d)
                else
                    ! candidate search uses the undistorted unit-sphere center
                    call lonlat_to_xyz(grid1%lon(ic,jc), grid1%lat(ic,jc), cx, cy, cz)
                    emb(1,cs) = cx; emb(2,cs) = cy; emb(3,cs) = cz
                end if
            end do
        end do
        call kdtree_build(tree, emb)

        ! candidate radius
        if (same_sys) then
            ! planar: half target diagonal + largest source-cell diagonal
            diag2 = sqrt(grid2%G%dx**2 + grid2%G%dy**2)
            rad   = 0.5_dp*diag2 + maxsrcdiag
        else
            ! chord on the unit sphere from source + target angular cell size.
            ! source is lon/lat (degrees); the projected target's cell size (grid
            ! units) is converted to degrees via the planet radius. Mirrors the
            ! radius in conservative_spherical: chord(src + tgt + 0.5*tgt).
            srcdeg = max(grid1%G%dx, grid1%G%dy)
            tgtdeg = (max(grid2%G%dx, grid2%G%dy)*xyc/grid2%cs%planet%a) &
                     / degrees_to_radians
            rad    = chord(srcdeg + 1.5_dp*tgtdeg)
        end if

        if (do_project) then
            ! Project only the source cells that are a candidate of some target
            ! cell. Far from the target the projection is ill-defined, not just
            ! distorted: stereographic sends the antipode of its center to
            ! infinity (a pole corner of a global source divides by zero). The
            ! candidate sets come from the sphere tree alone, so marking them
            ! first leaves the weights unchanged.
            allocate(need(n1), cidx(n1), cd2(n1))
            need = .false.
            do j = 1, ny2
                do i = 1, nx2
                    call lonlat_to_xyz(grid2%lon(i,j), grid2%lat(i,j), cx, cy, cz)
                    call kdtree_radius(tree, [cx,cy,cz], rad, cidx, cd2, nf)
                    need(cidx(1:nf)) = .true.
                end do
            end do
            deallocate(cidx, cd2)

            ! lx,ly are the source cell corners in the SOURCE axis system:
            ! geographic lon/lat for a plain lat-lon source, or ROTATED
            ! lon/lat (rlon/rlat) for a rotated-pole source. Subsample each
            ! cell edge into NSUB segments IN THAT AXIS SYSTEM (so the
            ! straight edge is traced along the true cell boundary), recover
            ! geographic lon/lat when the source is projected (rotated pole),
            ! then project every sample into the target plane so the source
            ! polygon follows the curved parallels/meridians instead of
            ! chording across them. The chord error grows with source cell
            ! size and projection distortion; subsampling cuts it from ~30 m
            ! to <1 m vs cdo for a 1deg source. (A separate ~O(100 m)
            ! planar-vs-spherical residual remains within a few cells of the
            ! geographic pole.)
            src_xyc = grid1%cs%xy_conv
            do jc = 1, ny1
                do ic = 1, nx1
                    cs = (jc-1)*nx1 + ic
                    if (.not. need(cs)) cycle
                    call cell_corners_grid(grid1, ic, jc, lx, ly)
                    k = 0
                    do e = 1, 4
                        en = mod(e,4) + 1
                        do s = 0, NSUB-1
                            tt   = real(s,dp)/real(NSUB,dp)
                            plon = lx(e) + tt*(lx(en) - lx(e))
                            plat = ly(e) + tt*(ly(en) - ly(e))
                            if (grid1%cs%is_projection) then
                                ! source axis is rotated/projected -> geographic
                                call oblimap_projection_inverse(plon*src_xyc, plat*src_xyc, &
                                                                glon, glat, grid1%cs%proj)
                                plon = glon; plat = glat
                            end if
                            call oblimap_projection(plon, plat, px, py, grid2%cs%proj)
                            k = k + 1
                            scornx(k,cs) = px/xyc; scorny(k,cs) = py/xyc
                        end do
                    end do
                    wfull_src(cs) = polygon_area(scornx(:,cs), scorny(:,cs))
                end do
            end do
            deallocate(need)
        end if

        nthreads = 1
        !$ nthreads = omp_get_max_threads()
        allocate(bufs(0:nthreads-1))

        !$omp parallel default(shared) &
        !$omp   private(i, j, c, cc, cs, nf, no, tcx, tcy, qpt, cx, cy, cz, area, cidx, cd2, ox, oy, tid)
        tid = 0
        !$ tid = omp_get_thread_num()
        allocate(cidx(n1), cd2(n1))
        !$omp do schedule(guided) collapse(2)
        do j = 1, ny2
            do i = 1, nx2
                c = (j-1)*nx2 + i                 ! target cell (column-major)
                call cell_corners_grid(grid2, i, j, tcx, tcy)
                if (same_sys) then
                    qpt = [grid2%x(i,j), grid2%y(i,j)]
                    call kdtree_radius(tree, qpt, rad, cidx, cd2, nf)
                else
                    call lonlat_to_xyz(grid2%lon(i,j), grid2%lat(i,j), cx, cy, cz)
                    call kdtree_radius(tree, [cx,cy,cz], rad, cidx, cd2, nf)
                end if
                do cc = 1, nf
                    cs = cidx(cc)
                    call polygon_clip(scornx(:,cs), scorny(:,cs), tcx, tcy, ox, oy, no)
                    if (no >= 3) then
                        area = polygon_area(ox(1:no), oy(1:no))
                        if (area > 0.0_dp) call buf_push(bufs(tid), cs, c, area)
                    end if
                end do
            end do
        end do
        !$omp end do
        deallocate(cidx, cd2)
        if (allocated(ox)) deallocate(ox, oy)
        !$omp end parallel

        call bufs_to_wm(wm, bufs, n1, n2)

        ! Coverage payload: full cell areas in the target plane (the map metric)
        allocate(wm%wfull_src(n1), wm%wfull_dst(n2))
        wm%wfull_src = wfull_src
        do j = 1, ny2
            do i = 1, nx2
                call cell_corners_grid(grid2, i, j, tcx, tcy)
                wm%wfull_dst((j-1)*nx2 + i) = polygon_area(tcx, tcy)
            end do
        end do
        call weight_map_coverage(wm, set_cov=.true.)

        call kdtree_free(tree)
        deallocate(emb, scornx, scorny, bufs, wfull_src)
    end subroutine conservative_planar

    subroutine conservative_spherical(wm, grid1, grid2)
        ! Stage C: clip on the unit sphere for lat-lon / Gaussian targets. Cell
        ! parallels (constant-latitude edges) are subsampled into NSUB great-
        ! circle segments; overlap area is the geodesic polygon area.
        type(weight_map_t), intent(inout) :: wm
        type(grid_class),   intent(in)    :: grid1, grid2

        type(kdtree_t) :: tree
        integer  :: n1, n2, nx1, ny1, nx2, ny2
        integer  :: i, j, ic, jc, c, cs, cc, nf, no, nv, nthreads, tid
        real(dp) :: rad, area, a, f, srcdeg, tgtdeg
        real(dp), allocatable :: emb(:,:)
        real(dp), allocatable :: slon(:,:), slat(:,:)     ! source cell corners (lon/lat)
        real(dp) :: tlon(4), tlat(4)
        integer,  allocatable :: cidx(:)
        real(dp), allocatable :: cd2(:)
        real(dp), allocatable :: sv(:,:), cv(:,:), ov(:,:)
        real(dp) :: olon(64), olat(64)
        type(link_buf), allocatable :: bufs(:)
        logical,  allocatable :: linked(:)
        real(dp), allocatable :: srad(:)                   ! source cap radii [rad]
        real(dp) :: cx, cy, cz, trad, ang
        real(dp), parameter :: CAP_MARGIN = 1.0e-9_dp      ! [rad], ~6 mm

        nx1 = grid1%G%nx; ny1 = grid1%G%ny; n1 = nx1*ny1
        nx2 = grid2%G%nx; ny2 = grid2%G%ny; n2 = nx2*ny2
        a = grid2%cs%planet%a; f = grid2%cs%planet%f

        ! Source cell corners as lon/lat, plus 3D centers for the tree
        allocate(slon(4,n1), slat(4,n1), emb(3,n1))
        do jc = 1, ny1
            do ic = 1, nx1
                cs = (jc-1)*nx1 + ic
                call cell_corners_lonlat(grid1, ic, jc, slon(:,cs), slat(:,cs))
                call lonlat_to_xyz(grid1%lon(ic,jc), grid1%lat(ic,jc), cx, cy, cz)
                emb(1,cs) = cx; emb(2,cs) = cy; emb(3,cs) = cz
            end do
        end do
        call kdtree_build(tree, emb)

        ! Spherical-cap pre-reject: each cell polygon lies in the cap around its
        ! centre whose radius reaches its farthest vertex (caps < 90 deg are
        ! convex, so they contain the great-circle edges too). A candidate pair
        ! whose caps do not meet cannot overlap and is not clipped; the weights
        ! do not change (docs/dev/coupler-global-design.md, 7 i).
        allocate(srad(n1))
        !$omp parallel do default(shared) schedule(guided) private(cs, ic, jc, sv)
        do cs = 1, n1
            ic = mod(cs-1, nx1) + 1; jc = (cs-1)/nx1 + 1
            call build_sphere_poly(slon(:,cs), slat(:,cs), grid1%lon(ic,jc), grid1%lat(ic,jc), sv)
            srad(cs) = cap_radius(emb(:,cs), sv)
        end do
        !$omp end parallel do

        ! candidate radius (chord on unit sphere): source + target angular cell
        ! size, both in DEGREES. grid2 is the lat-lon/Gaussian target, so its dx/dy
        ! are already degrees. grid1 (source) may be projected or cartesian, whose
        ! dx/dy are in axis units (km/m), NOT degrees -- convert via the planet
        ! radius (as the cross-system planar path does at conservative_planar).
        ! Using the raw axis value as degrees inflates the radius by ~100x, pulling
        ! nearly every source cell per target and degrading toward O(n_src*n_tgt).
        if (grid1%cs%is_cartesian .or. grid1%cs%is_projection) then
            srcdeg = (max(grid1%G%dx, grid1%G%dy)*grid1%cs%xy_conv/grid1%cs%planet%a) &
                     / degrees_to_radians
        else
            srcdeg = max(grid1%G%dx, grid1%G%dy)
        end if
        tgtdeg = max(grid2%G%dx, grid2%G%dy)
        rad = chord(srcdeg + 1.5_dp*tgtdeg)

        nthreads = 1
        !$ nthreads = omp_get_max_threads()
        allocate(bufs(0:nthreads-1))

        !$omp parallel default(shared) &
        !$omp   private(i, j, c, cc, cs, ic, jc, nf, no, nv, tlon, tlat, cx, cy, cz, area, &
        !$omp           cidx, cd2, sv, cv, ov, olon, olat, tid, trad, ang)
        tid = 0
        !$ tid = omp_get_thread_num()
        allocate(cidx(n1), cd2(n1))
        !$omp do schedule(guided) collapse(2)
        do j = 1, ny2
            do i = 1, nx2
                c = (j-1)*nx2 + i
                call cell_corners_lonlat(grid2, i, j, tlon, tlat)
                call build_sphere_poly(tlon, tlat, grid2%lon(i,j), grid2%lat(i,j), cv)
                call lonlat_to_xyz(grid2%lon(i,j), grid2%lat(i,j), cx, cy, cz)
                trad = cap_radius([cx,cy,cz], cv)
                call kdtree_radius(tree, [cx,cy,cz], rad, cidx, cd2, nf)
                do cc = 1, nf
                    cs = cidx(cc)
                    ang = acos(max(-1.0_dp, min(1.0_dp, emb(1,cs)*cx + emb(2,cs)*cy + emb(3,cs)*cz)))
                    if (ang > srad(cs) + trad + CAP_MARGIN) cycle
                    ic = mod(cs-1, nx1) + 1; jc = (cs-1)/nx1 + 1
                    call build_sphere_poly(slon(:,cs), slat(:,cs), &
                                           grid1%lon(ic,jc), grid1%lat(ic,jc), sv)
                    call polygon_clip_sphere(sv, cv, ov, no)
                    if (no >= 3 .and. no <= 64) then
                        do nv = 1, no
                            call xyz_to_lonlat(ov(1,nv), ov(2,nv), ov(3,nv), olon(nv), olat(nv))
                        end do
                        area = planet_area(a, f, olon(1:no), olat(1:no))
                        if (area > 0.0_dp) call buf_push(bufs(tid), cs, c, area)
                    end if
                end do
            end do
        end do
        !$omp end do
        deallocate(cidx, cd2)
        if (allocated(sv)) deallocate(sv)
        if (allocated(cv)) deallocate(cv)
        if (allocated(ov)) deallocate(ov)
        !$omp end parallel

        call bufs_to_wm(wm, bufs, n1, n2)

        ! Coverage payload: full cell areas (planet_area of the same great-circle
        ! polygons that are clipped; the map metric). Sources only where linked.
        allocate(wm%wfull_src(n1), wm%wfull_dst(n2), linked(n1))
        wm%wfull_src = 0.0_dp
        linked = .false.
        linked(wm%src(1:wm%n_links)) = .true.
        !$omp parallel do default(shared) schedule(guided) private(cs, ic, jc, sv)
        do cs = 1, n1
            if (.not. linked(cs)) cycle
            ic = mod(cs-1, nx1) + 1; jc = (cs-1)/nx1 + 1
            call build_sphere_poly(slon(:,cs), slat(:,cs), grid1%lon(ic,jc), grid1%lat(ic,jc), sv)
            wm%wfull_src(cs) = sphere_poly_area(sv, a, f)
        end do
        !$omp end parallel do
        !$omp parallel do default(shared) schedule(guided) collapse(2) private(i, j, tlon, tlat, cv)
        do j = 1, ny2
            do i = 1, nx2
                call cell_corners_lonlat(grid2, i, j, tlon, tlat)
                call build_sphere_poly(tlon, tlat, grid2%lon(i,j), grid2%lat(i,j), cv)
                wm%wfull_dst((j-1)*nx2 + i) = sphere_poly_area(cv, a, f)
            end do
        end do
        !$omp end parallel do
        call weight_map_coverage(wm, set_cov=.true.)

        call kdtree_free(tree)
        deallocate(emb, slon, slat, bufs, linked, srad)
    end subroutine conservative_spherical

    real(dp) function cap_radius(ctr, verts) result(r)
        ! Angular radius [rad] of the cap around unit vector ctr that holds
        ! all polygon vertices.
        real(dp), intent(in) :: ctr(3), verts(:,:)
        real(dp) :: dmin
        integer  :: k
        dmin = 1.0_dp
        do k = 1, size(verts, 2)
            dmin = min(dmin, ctr(1)*verts(1,k) + ctr(2)*verts(2,k) + ctr(3)*verts(3,k))
        end do
        r = acos(max(-1.0_dp, min(1.0_dp, dmin)))
    end function cap_radius

    real(dp) function sphere_poly_area(verts, a, f) result(area)
        ! Area on the planet (a, f) of a polygon given as unit vectors.
        real(dp), intent(in) :: verts(:,:), a, f
        real(dp) :: plon(size(verts,2)), plat(size(verts,2))
        integer  :: nv
        do nv = 1, size(verts,2)
            call xyz_to_lonlat(verts(1,nv), verts(2,nv), verts(3,nv), plon(nv), plat(nv))
        end do
        area = planet_area(a, f, plon, plat)
    end function sphere_poly_area

    subroutine lon_edges(grid, xe)
        ! Longitude cell edges (0:nx) of a lat-lon / Gaussian grid for the
        ! separable lat-lon path: the grid's x_bnds when the axis is ascending.
        ! A grid built from 0..360 centres with lon180=.true. stores them
        ! remapped in place (0..179, -179..-0.7): its explicit bounds (shifted
        ! with the centres) are unwrapped to a monotonic sequence; its default
        ! bounds hold a wrong seam edge (~0 instead of ~180), so the centres are
        ! unwrapped first and their midpoint edges taken. (For an ascending
        ! axis unwrap_lon is the identity, so all branches agree there.)
        type(grid_class), intent(in)  :: grid
        real(dp),         intent(out) :: xe(0:)
        real(dp), allocatable :: xc(:)
        integer :: i
        logical :: ascending
        ascending = .true.
        do i = 2, grid%G%nx
            if (grid%G%x(i) < grid%G%x(i-1)) ascending = .false.
        end do
        if (ascending) then
            xe = grid%G%x_bnds
        else if (grid%G%explicit_bnds) then
            call unwrap_lon(grid%G%x_bnds, xe)
        else
            allocate(xc(grid%G%nx))
            call unwrap_lon(grid%G%x, xc)
            call axis_midpoint_bounds(xc, xe, grid%G%dx)
        end if
    end subroutine lon_edges

    subroutine unwrap_lon(xin, xout)
        ! Copy longitude centers to a monotonic ascending sequence, undoing the
        ! in-place -180..180 remap (lon180=.true. on a 0..360 grid) that leaves the
        ! stored centers non-monotonic across the seam. Indices are preserved; only
        ! the values are shifted by whole turns so adjacent-midpoint edges are sane.
        real(dp), intent(in)  :: xin(:)
        real(dp), intent(out) :: xout(:)
        integer :: i
        xout(1) = xin(1)
        do i = 2, size(xin)
            xout(i) = xin(i)
            do while (xout(i) < xout(i-1)); xout(i) = xout(i) + 360.0_dp; end do
        end do
    end subroutine unwrap_lon

    real(dp) function interval_overlap(a, b, c, d) result(ov)
        ! Length of the overlap of intervals [a,b] and [c,d] (0 if disjoint).
        ! Orientation-independent: each interval is normalized to [min,max] first,
        ! so descending axes (e.g. a Gaussian grid stored north->south, or a
        ! sin(lat) band from such an axis) overlap correctly. Without this the
        ! separable lat/planar paths return zero links for descending targets.
        real(dp), intent(in) :: a, b, c, d
        ov = max(0.0_dp, min(max(a,b), max(c,d)) - max(min(a,b), min(c,d)))
    end function interval_overlap

    real(dp) function lon_overlap(a, b, c, d) result(ov)
        ! Longitude overlap of target [a,b] and source [c,d], allowing a +-360 deg
        ! wrap so grids with different longitude conventions (0..360 vs -180..180)
        ! and the periodic seam are handled. Both intervals are < 360 wide, so at
        ! most one shift yields a positive overlap; take the largest.
        ! Orientation-independent: each interval is normalized to [min,max] first
        ! so a descending longitude axis overlaps correctly (see interval_overlap).
        real(dp), intent(in) :: a, b, c, d
        real(dp) :: alo, ahi, clo, chi
        integer :: k
        alo = min(a,b); ahi = max(a,b)
        clo = min(c,d); chi = max(c,d)
        ov = 0.0_dp
        do k = -1, 1
            ov = max(ov, min(ahi, chi + 360.0_dp*k) - max(alo, clo + 360.0_dp*k))
        end do
        ov = max(0.0_dp, ov)
    end function lon_overlap

    subroutine conservative_latlon(wm, grid1, grid2)
        ! Analytic conservative weights for lat-lon/Gaussian source -> lat-lon/
        ! Gaussian target. Cells are axis-aligned in (lon,lat), so the spherical
        ! overlap area of a source and target cell is separable:
        !     area = R^2 * (lon overlap, rad) * (sin(lat_hi) - sin(lat_lo))
        ! The R^2 and the deg->rad factor are global constants that cancel in the
        ! per-target weighted mean, so the stored weight is
        !     w = lon_overlap(deg) * sinlat_overlap.
        ! The candidate search is separable too: overlapping source rows/columns
        ! per target row/column are found by 1-D interval tests (O(1) amortized),
        ! never a k-d tree. This replaces the great-circle polygon clip
        ! (conservative_spherical) for the both-lat-lon case.
        type(weight_map_t), intent(inout) :: wm
        type(grid_class),   intent(in)    :: grid1, grid2

        integer  :: nx1, ny1, nx2, ny2, n1, n2
        integer  :: i1, j1, i2, j2, nthreads, tid, p, q
        real(dp) :: wlo, wla
        real(dp), allocatable :: sxe(:), sye(:), txe(:), tye(:)   ! 0:n edges
        real(dp), allocatable :: ssin(:), tsin(:)                 ! 0:ny sin(lat edge)
        ! CSR of per-target-column lon overlaps and per-target-row lat overlaps
        integer,  allocatable :: lon_off(:), lat_off(:)
        integer,  allocatable :: lon_i1(:),  lat_j1(:)
        real(dp), allocatable :: lon_w(:),   lat_w(:)
        type(link_buf), allocatable :: bufs(:)

        nx1 = grid1%G%nx; ny1 = grid1%G%ny; n1 = nx1*ny1
        nx2 = grid2%G%nx; ny2 = grid2%G%ny; n2 = nx2*ny2

        ! Cell edges from the grids' bounds; longitude edges via lon_edges (which
        ! unwraps a lon180-shifted axis, see there). lon_overlap handles the
        ! +-360 convention offset between source and target, so only
        ! monotonicity matters.
        allocate(sxe(0:nx1), sye(0:ny1), txe(0:nx2), tye(0:ny2))
        call lon_edges(grid1, sxe); sye = grid1%G%y_bnds
        call lon_edges(grid2, txe); tye = grid2%G%y_bnds

        ! Clamp latitude edges to the poles (polar cells are half cells) and take
        ! the sine, in which latitude bands overlap linearly.
        allocate(ssin(0:ny1), tsin(0:ny2))
        do j1 = 0, ny1
            ssin(j1) = sin(max(-90.0_dp, min(90.0_dp, sye(j1)))*degrees_to_radians)
        end do
        do j2 = 0, ny2
            tsin(j2) = sin(max(-90.0_dp, min(90.0_dp, tye(j2)))*degrees_to_radians)
        end do

        ! Build the 1-D overlap tables (two-pass CSR: count, then fill).
        call build_axis_overlap_lon(txe, sxe, nx2, nx1, lon_off, lon_i1, lon_w)
        call build_axis_overlap_lat(tsin, ssin, ny2, ny1, lat_off, lat_j1, lat_w)

        nthreads = 1
        !$ nthreads = omp_get_max_threads()
        allocate(bufs(0:nthreads-1))

        ! Assemble links as the (row overlap) x (column overlap) cross product per
        ! target cell. Parallel over target rows; per-thread buffers avoid a shared
        ! append counter (bufs_to_wm merges them in destination order).
        !$omp parallel default(shared) &
        !$omp   private(i2, j2, i1, j1, p, q, wlo, wla, tid)
        tid = 0
        !$ tid = omp_get_thread_num()
        !$omp do schedule(guided)
        do j2 = 1, ny2
            do q = lat_off(j2), lat_off(j2+1)-1
                j1  = lat_j1(q); wla = lat_w(q)
                do i2 = 1, nx2
                    do p = lon_off(i2), lon_off(i2+1)-1
                        i1  = lon_i1(p); wlo = lon_w(p)
                        call buf_push(bufs(tid), (j1-1)*nx1 + i1, (j2-1)*nx2 + i2, wlo*wla)
                    end do
                end do
            end do
        end do
        !$omp end do
        !$omp end parallel

        call bufs_to_wm(wm, bufs, n1, n2)

        ! Coverage payload: full cells in the same metric (deg * sin-band)
        allocate(wm%wfull_src(n1), wm%wfull_dst(n2))
        do j1 = 1, ny1
            do i1 = 1, nx1
                wm%wfull_src((j1-1)*nx1 + i1) = abs(sxe(i1) - sxe(i1-1))*abs(ssin(j1) - ssin(j1-1))
            end do
        end do
        do j2 = 1, ny2
            do i2 = 1, nx2
                wm%wfull_dst((j2-1)*nx2 + i2) = abs(txe(i2) - txe(i2-1))*abs(tsin(j2) - tsin(j2-1))
            end do
        end do
        call weight_map_coverage(wm, set_cov=.true.)

        deallocate(sxe, sye, txe, tye, ssin, tsin, bufs)
        deallocate(lon_off, lon_i1, lon_w, lat_off, lat_j1, lat_w)
    end subroutine conservative_latlon

    subroutine build_axis_overlap_lon(txe, sxe, nt, ns, off, idx, w)
        ! Per-target-column CSR list of overlapping source columns and the
        ! (wrap-aware) longitude overlap width. Two-pass to size exactly.
        real(dp), intent(in)  :: txe(0:), sxe(0:)
        integer,  intent(in)  :: nt, ns
        integer,  allocatable, intent(out) :: off(:), idx(:)
        real(dp), allocatable, intent(out) :: w(:)
        integer  :: i2, i1, cnt, p
        real(dp) :: ov
        allocate(off(nt+1)); off = 0
        do i2 = 1, nt
            cnt = 0
            do i1 = 1, ns
                if (lon_overlap(txe(i2-1), txe(i2), sxe(i1-1), sxe(i1)) > 0.0_dp) cnt = cnt + 1
            end do
            off(i2+1) = cnt
        end do
        off(1) = 1
        do i2 = 1, nt
            off(i2+1) = off(i2+1) + off(i2)
        end do
        allocate(idx(off(nt+1)-1), w(off(nt+1)-1))
        do i2 = 1, nt
            p = off(i2)
            do i1 = 1, ns
                ov = lon_overlap(txe(i2-1), txe(i2), sxe(i1-1), sxe(i1))
                if (ov > 0.0_dp) then
                    idx(p) = i1; w(p) = ov; p = p + 1
                end if
            end do
        end do
    end subroutine build_axis_overlap_lon

    subroutine build_axis_overlap_lat(tsin, ssin, nt, ns, off, idx, w)
        ! Per-target-row CSR list of overlapping source rows and the sin(lat)-band
        ! overlap. Non-periodic; edges are monotonic so ranges are contiguous.
        real(dp), intent(in)  :: tsin(0:), ssin(0:)
        integer,  intent(in)  :: nt, ns
        integer,  allocatable, intent(out) :: off(:), idx(:)
        real(dp), allocatable, intent(out) :: w(:)
        integer  :: j2, j1, cnt, p
        real(dp) :: ov
        allocate(off(nt+1)); off = 0
        do j2 = 1, nt
            cnt = 0
            do j1 = 1, ns
                if (interval_overlap(tsin(j2-1), tsin(j2), ssin(j1-1), ssin(j1)) > 0.0_dp) cnt = cnt + 1
            end do
            off(j2+1) = cnt
        end do
        off(1) = 1
        do j2 = 1, nt
            off(j2+1) = off(j2+1) + off(j2)
        end do
        allocate(idx(off(nt+1)-1), w(off(nt+1)-1))
        do j2 = 1, nt
            p = off(j2)
            do j1 = 1, ns
                ov = interval_overlap(tsin(j2-1), tsin(j2), ssin(j1-1), ssin(j1))
                if (ov > 0.0_dp) then
                    idx(p) = j1; w(p) = ov; p = p + 1
                end if
            end do
        end do
    end subroutine build_axis_overlap_lat

    subroutine conservative_planar_reg(wm, grid1, grid2)
        ! Analytic conservative weights when source and target share one Cartesian
        ! plane: cells are axis-aligned rectangles, so the overlap area is the
        ! separable product of the x- and y-interval overlaps (in shared axis
        ! units). No polygon clip, no k-d tree. Exact to rounding.
        type(weight_map_t), intent(inout) :: wm
        type(grid_class),   intent(in)    :: grid1, grid2

        integer  :: nx1, ny1, nx2, ny2, n1, n2
        integer  :: i1, j1, i2, j2, nthreads, tid, p, q
        real(dp) :: wx, wy
        real(dp), allocatable :: sxe(:), sye(:), txe(:), tye(:)
        integer,  allocatable :: x_off(:), y_off(:), x_i1(:), y_j1(:)
        real(dp), allocatable :: x_w(:), y_w(:)
        type(link_buf), allocatable :: bufs(:)

        nx1 = grid1%G%nx; ny1 = grid1%G%ny; n1 = nx1*ny1
        nx2 = grid2%G%nx; ny2 = grid2%G%ny; n2 = nx2*ny2

        allocate(sxe(0:nx1), sye(0:ny1), txe(0:nx2), tye(0:ny2))
        sxe = grid1%G%x_bnds; sye = grid1%G%y_bnds
        txe = grid2%G%x_bnds; tye = grid2%G%y_bnds

        call build_axis_overlap_lat(txe, sxe, nx2, nx1, x_off, x_i1, x_w)
        call build_axis_overlap_lat(tye, sye, ny2, ny1, y_off, y_j1, y_w)

        nthreads = 1
        !$ nthreads = omp_get_max_threads()
        allocate(bufs(0:nthreads-1))

        !$omp parallel default(shared) private(i2, j2, i1, j1, p, q, wx, wy, tid)
        tid = 0
        !$ tid = omp_get_thread_num()
        !$omp do schedule(guided)
        do j2 = 1, ny2
            do q = y_off(j2), y_off(j2+1)-1
                j1 = y_j1(q); wy = y_w(q)
                do i2 = 1, nx2
                    do p = x_off(i2), x_off(i2+1)-1
                        i1 = x_i1(p); wx = x_w(p)
                        call buf_push(bufs(tid), (j1-1)*nx1 + i1, (j2-1)*nx2 + i2, wx*wy)
                    end do
                end do
            end do
        end do
        !$omp end do
        !$omp end parallel

        call bufs_to_wm(wm, bufs, n1, n2)

        ! Coverage payload: full cell areas in the shared plane
        allocate(wm%wfull_src(n1), wm%wfull_dst(n2))
        do j1 = 1, ny1
            do i1 = 1, nx1
                wm%wfull_src((j1-1)*nx1 + i1) = abs(sxe(i1) - sxe(i1-1))*abs(sye(j1) - sye(j1-1))
            end do
        end do
        do j2 = 1, ny2
            do i2 = 1, nx2
                wm%wfull_dst((j2-1)*nx2 + i2) = abs(txe(i2) - txe(i2-1))*abs(tye(j2) - tye(j2-1))
            end do
        end do
        call weight_map_coverage(wm, set_cov=.true.)

        deallocate(sxe, sye, txe, tye, bufs)
        deallocate(x_off, x_i1, x_w, y_off, y_j1, y_w)
    end subroutine conservative_planar_reg

    real(dp) function chord(deg) result(c)
        ! chord length on the unit sphere subtended by an angle of `deg` degrees
        real(dp), intent(in) :: deg
        c = 2.0_dp*sin(0.5_dp*deg*degrees_to_radians)
    end function chord

    subroutine cell_corners_lonlat(grid, i, j, clon, clat)
        ! Cell corners as lon/lat: from axis midpoints for lat-lon/Gaussian grids,
        ! or by inverse-projecting the projected-axis corners for projected grids.
        type(grid_class), intent(in)  :: grid
        integer,          intent(in)  :: i, j
        real(dp),         intent(out) :: clon(4), clat(4)
        real(dp) :: cx(4), cy(4), xyc, lo, la
        integer  :: k
        call cell_corners_grid(grid, i, j, cx, cy)
        if (grid%cs%is_projection) then
            xyc = grid%cs%xy_conv
            do k = 1, 4
                call oblimap_projection_inverse(cx(k)*xyc, cy(k)*xyc, lo, la, grid%cs%proj)
                clon(k) = lo; clat(k) = la
            end do
        else
            clon = cx; clat = cy     ! axis values are lon/lat
        end if
    end subroutine cell_corners_lonlat

    subroutine build_sphere_poly(clon, clat, ctrlon, ctrlat, verts)
        ! Build a spherical polygon (3D unit vectors) from 4 lon/lat corners,
        ! subsampling each edge into NSUB great-circle segments, oriented so the
        ! cell center is on the interior (+normal) side of every edge.
        real(dp), intent(in)  :: clon(4), clat(4), ctrlon, ctrlat
        real(dp), allocatable, intent(out) :: verts(:,:)
        real(dp) :: ll_lon(4*NSUB), ll_lat(4*NSUB)
        real(dp) :: t, cxv, cyv, czv, n1v(3), ev(3)
        integer  :: e, s, nv, k, en
        ! walk the 4 edges, subsampling
        nv = 0
        do e = 1, 4
            en = mod(e,4) + 1
            do s = 0, NSUB-1
                t = real(s,dp)/real(NSUB,dp)
                nv = nv + 1
                ll_lon(nv) = clon(e) + t*(clon(en) - clon(e))
                ll_lat(nv) = clat(e) + t*(clat(en) - clat(e))
            end do
        end do
        if (allocated(verts)) deallocate(verts)
        allocate(verts(3, nv))
        do k = 1, nv
            call lonlat_to_xyz(ll_lon(k), ll_lat(k), verts(1,k), verts(2,k), verts(3,k))
        end do
        ! Drop repeated vertices (a polar cell has two corners at the pole, so
        ! its polar edge collapses to one point): degenerate zero-length edges
        ! break the great-circle clip and the orientation test below.
        call drop_repeated(verts)
        nv = size(verts, 2)
        ! orient so center is interior: check edge 1 normal vs center
        call lonlat_to_xyz(ctrlon, ctrlat, cxv, cyv, czv)
        n1v(1) = verts(2,1)*verts(3,2) - verts(3,1)*verts(2,2)
        n1v(2) = verts(3,1)*verts(1,2) - verts(1,1)*verts(3,2)
        n1v(3) = verts(1,1)*verts(2,2) - verts(2,1)*verts(1,2)
        if (n1v(1)*cxv + n1v(2)*cyv + n1v(3)*czv < 0.0_dp) then
            ! reverse
            block
                real(dp) :: tmp(3)
                do k = 1, nv/2
                    tmp = verts(:,k); verts(:,k) = verts(:,nv-k+1); verts(:,nv-k+1) = tmp
                end do
            end block
        end if
    end subroutine build_sphere_poly

    subroutine drop_repeated(verts)
        ! Remove vertices that coincide with their predecessor (cyclically), to
        ! 1e-13 on the unit sphere (points at a pole differ by ~1e-17 with the
        ! longitude, as cos(90 deg) is not exactly 0 in floating point).
        real(dp), allocatable, intent(inout) :: verts(:,:)
        real(dp), parameter :: TOL = 1.0e-13_dp
        real(dp), allocatable :: tmp(:,:)
        logical,  allocatable :: keep(:)
        integer :: k, n, m, kp
        n = size(verts, 2)
        allocate(keep(n))
        keep = .true.
        kp = 1                                   ! last kept vertex
        do k = 2, n
            if (maxval(abs(verts(:,k) - verts(:,kp))) < TOL) then
                keep(k) = .false.
            else
                kp = k
            end if
        end do
        if (kp > 1) then
            if (maxval(abs(verts(:,kp) - verts(:,1))) < TOL) keep(kp) = .false.
        end if
        if (all(keep)) return
        m = count(keep)
        allocate(tmp(3, m))
        tmp = reshape(pack(verts, spread(keep, 1, 3)), [3, m])
        call move_alloc(tmp, verts)
    end subroutine drop_repeated

    subroutine lonlat_to_xyz(lon, lat, x, y, z)
        real(dp), intent(in)  :: lon, lat
        real(dp), intent(out) :: x, y, z
        real(dp) :: la, lo
        la = lat*degrees_to_radians; lo = lon*degrees_to_radians
        x = cos(la)*cos(lo); y = cos(la)*sin(lo); z = sin(la)
    end subroutine lonlat_to_xyz

    subroutine xyz_to_lonlat(x, y, z, lon, lat)
        real(dp), intent(in)  :: x, y, z
        real(dp), intent(out) :: lon, lat
        lon = atan2(y, x) / degrees_to_radians
        lat = asin(max(-1.0_dp, min(1.0_dp, z))) / degrees_to_radians
    end subroutine xyz_to_lonlat

    subroutine cell_corners_grid(grid, i, j, cx, cy)
        ! Four corners (CCW) of cell (i,j) in axis units, from the grid's cell
        ! edges (grid%G%x_bnds, grid%G%y_bnds).
        type(grid_class), intent(in)  :: grid
        integer,          intent(in)  :: i, j
        real(dp),         intent(out) :: cx(4), cy(4)
        real(dp) :: x1, x2, y1, y2

        x1 = grid%G%x_bnds(i-1); x2 = grid%G%x_bnds(i)
        y1 = grid%G%y_bnds(j-1); y2 = grid%G%y_bnds(j)

        cx = [x1, x2, x2, x1]
        cy = [y1, y1, y2, y2]
    end subroutine cell_corners_grid

end module conservative
