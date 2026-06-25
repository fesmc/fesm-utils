module conservative
    ! Analytic conservative remapping via polygon clipping (improvement-doc §2,
    ! stage B: planar Sutherland-Hodgman). For each target cell, candidate source
    ! cells are found with a k-d-tree radius query, then each source cell polygon
    ! is clipped against the target cell and the overlap area becomes the link
    ! weight. The result is a MAP_WEIGHT applied as an area-weighted mean.
    !
    ! Two planar modes are supported:
    !   - same-system: source and target share a Cartesian/projected plane;
    !     cells are clipped directly.
    !   - lat-lon source -> projected target: source cell corners (lon/lat) are
    !     projected into the target projection plane, then clipped there.
    ! Lat-lon / Gaussian *targets* need great-circle clipping (stage C, later).

    use precision,       only: dp
    use constants, only: degrees_to_radians
    use coordinates,     only: grid_class, compare_coord
    use oblimap_projection_module, only: oblimap_projection, oblimap_projection_inverse
    use planet,          only: planet_area
    use kdtree
    use weight_map
    use polygons
    use mapping,         only: map_class

    implicit none
    private

    integer, parameter :: NSUB = 8   ! great-circle segments per parallel edge (stage C)

    public :: map_init_conservative

contains

    subroutine map_init_conservative(map, grid1, grid2)
        ! Dispatch: planar clip for Cartesian/projected targets, spherical clip
        ! (great-circle arcs) for lat-lon / Gaussian targets.
        type(map_class),  intent(inout) :: map
        type(grid_class), intent(in)    :: grid1, grid2
        if (grid2%cs%is_cartesian) then
            call conservative_planar(map, grid1, grid2)
        else
            call conservative_spherical(map, grid1, grid2)
        end if
    end subroutine map_init_conservative

    subroutine conservative_planar(map, grid1, grid2)
        type(map_class),  intent(inout) :: map
        type(grid_class), intent(in)    :: grid1, grid2

        type(kdtree_t) :: tree
        logical  :: same_sys, do_project
        integer  :: n1, n2, nx1, ny1, nx2, ny2
        integer  :: i, j, ic, jc, c, cs, cc, k, nf, nl, nlmax, no
        real(dp) :: rad, diag2, area, maxsrcdiag, d, xyc, px, py
        real(dp) :: cx, cy, cz, srcdeg, tgtdeg
        real(dp) :: lx(4), ly(4), tcx(4), tcy(4), qpt(2)
        real(dp), allocatable :: emb(:,:), ox(:), oy(:)
        real(dp), allocatable :: scornx(:,:), scorny(:,:)   ! source corners in target plane
        integer,  allocatable :: cidx(:)
        real(dp), allocatable :: cd2(:)
        integer,  allocatable :: tsrc(:), tdst(:)
        real(dp), allocatable :: tw(:)

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

        map%name1 = grid1%name
        map%name2 = grid2%name
        map%cs    = grid2%cs
        map%is_grid = .true.
        map%G     = grid2%G
        map%npts  = n2
        map%nmax  = 0
        map%is_same_map = same_sys
        if (allocated(map%x)) deallocate(map%x, map%y, map%lon, map%lat)
        allocate(map%x(n2), map%y(n2), map%lon(n2), map%lat(n2))
        map%x   = reshape(grid2%x,   [n2])
        map%y   = reshape(grid2%y,   [n2])
        map%lon = reshape(grid2%lon, [n2])
        map%lat = reshape(grid2%lat, [n2])

        ! Express every source cell's corners in the TARGET plane (the planar clip
        ! runs there in both modes). The candidate-search tree differs by mode:
        !   - same-system: planar tree over projected source centers.
        !   - cross-system (lat-lon -> projection): tree over source centers on the
        !     UNIT SPHERE. Projecting source centers would distort enormously far
        !     from the projection origin (stereographic diverges near the antipode),
        !     inflating the search radius until every target pulls nearly every
        !     source and the k-d tree degrades to O(n_src*n_tgt). The sphere is
        !     undistorted, so the candidate set per target stays O(1).
        !     (see docs/coords-performance.md)
        allocate(scornx(4, n1), scorny(4, n1))
        if (same_sys) then
            allocate(emb(2, n1))
        else
            allocate(emb(3, n1))
        end if
        maxsrcdiag = 0.0_dp
        do jc = 1, ny1
            do ic = 1, nx1
                cs = (jc-1)*nx1 + ic
                call cell_corners_grid(grid1, ic, jc, lx, ly)
                if (same_sys) then
                    scornx(:,cs) = lx; scorny(:,cs) = ly
                    emb(1,cs) = grid1%x(ic,jc); emb(2,cs) = grid1%y(ic,jc)
                    d = sqrt((maxval(scornx(:,cs))-minval(scornx(:,cs)))**2 &
                           + (maxval(scorny(:,cs))-minval(scorny(:,cs)))**2)
                    maxsrcdiag = max(maxsrcdiag, d)
                else
                    ! lx,ly are lon,lat -> project corners into the target plane
                    do k = 1, 4
                        call oblimap_projection(lx(k), ly(k), px, py, grid2%cs%proj)
                        scornx(k,cs) = px/xyc; scorny(k,cs) = py/xyc
                    end do
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

        allocate(cidx(n1), cd2(n1))
        nlmax = max(n2*16, n1)
        allocate(tsrc(nlmax), tdst(nlmax), tw(nlmax))

        nl = 0
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
                        if (area > 0.0_dp) then
                            if (nl >= nlmax) call grow_links(tsrc, tdst, tw, nlmax)
                            nl = nl + 1
                            tsrc(nl) = cs
                            tdst(nl) = c
                            tw(nl)   = area
                        end if
                    end if
                end do
            end do
        end do

        call weight_map_alloc(map%wm, MAP_WEIGHT, n1, n2, nl)
        map%wm%src(1:nl) = tsrc(1:nl)
        map%wm%dst(1:nl) = tdst(1:nl)
        map%wm%w(1:nl)   = tw(1:nl)
        call weight_map_index(map%wm)

        call kdtree_free(tree)
        deallocate(emb, scornx, scorny, cidx, cd2, tsrc, tdst, tw)
        if (allocated(ox)) deallocate(ox, oy)
    end subroutine conservative_planar

    subroutine conservative_spherical(map, grid1, grid2)
        ! Stage C: clip on the unit sphere for lat-lon / Gaussian targets. Cell
        ! parallels (constant-latitude edges) are subsampled into NSUB great-
        ! circle segments; overlap area is the geodesic polygon area.
        type(map_class),  intent(inout) :: map
        type(grid_class), intent(in)    :: grid1, grid2

        type(kdtree_t) :: tree
        integer  :: n1, n2, nx1, ny1, nx2, ny2
        integer  :: i, j, ic, jc, c, cs, cc, nf, nl, nlmax, no, nv
        real(dp) :: rad, area, a, f
        real(dp), allocatable :: emb(:,:)
        real(dp), allocatable :: slon(:,:), slat(:,:)     ! source cell corners (lon/lat)
        real(dp) :: tlon(4), tlat(4)
        real(dp), allocatable :: cidx_r(:)
        integer,  allocatable :: cidx(:)
        real(dp), allocatable :: cd2(:)
        real(dp), allocatable :: sv(:,:), cv(:,:), ov(:,:)
        real(dp) :: olon(64), olat(64)
        integer,  allocatable :: tsrc(:), tdst(:)
        real(dp), allocatable :: tw(:)
        real(dp) :: cx, cy, cz

        nx1 = grid1%G%nx; ny1 = grid1%G%ny; n1 = nx1*ny1
        nx2 = grid2%G%nx; ny2 = grid2%G%ny; n2 = nx2*ny2
        a = grid2%cs%planet%a; f = grid2%cs%planet%f

        map%name1 = grid1%name; map%name2 = grid2%name
        map%cs = grid2%cs; map%is_grid = .true.; map%G = grid2%G
        map%npts = n2; map%nmax = 0; map%is_same_map = compare_coord(grid1, grid2)
        if (allocated(map%x)) deallocate(map%x, map%y, map%lon, map%lat)
        allocate(map%x(n2), map%y(n2), map%lon(n2), map%lat(n2))
        map%x = reshape(grid2%x,[n2]); map%y = reshape(grid2%y,[n2])
        map%lon = reshape(grid2%lon,[n2]); map%lat = reshape(grid2%lat,[n2])

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

        ! candidate radius (chord on unit sphere): target + source angular size
        rad = chord(max(grid1%G%dx, grid1%G%dy) + max(grid2%G%dx, grid2%G%dy) + &
                    0.5_dp*max(grid2%G%dx, grid2%G%dy))

        allocate(cidx(n1), cd2(n1), cidx_r(3))
        nlmax = max(n2*16, n1)
        allocate(tsrc(nlmax), tdst(nlmax), tw(nlmax))

        nl = 0
        do j = 1, ny2
            do i = 1, nx2
                c = (j-1)*nx2 + i
                call cell_corners_lonlat(grid2, i, j, tlon, tlat)
                call build_sphere_poly(tlon, tlat, grid2%lon(i,j), grid2%lat(i,j), cv)
                call lonlat_to_xyz(grid2%lon(i,j), grid2%lat(i,j), cx, cy, cz)
                call kdtree_radius(tree, [cx,cy,cz], rad, cidx, cd2, nf)
                do cc = 1, nf
                    cs = cidx(cc)
                    ic = mod(cs-1, nx1) + 1; jc = (cs-1)/nx1 + 1
                    call build_sphere_poly(slon(:,cs), slat(:,cs), &
                                           grid1%lon(ic,jc), grid1%lat(ic,jc), sv)
                    call polygon_clip_sphere(sv, cv, ov, no)
                    if (no >= 3) then
                        if (no <= 64) then
                            do nv = 1, no
                                call xyz_to_lonlat(ov(1,nv), ov(2,nv), ov(3,nv), olon(nv), olat(nv))
                            end do
                            area = planet_area(a, f, olon(1:no), olat(1:no))
                            if (area > 0.0_dp) then
                                if (nl >= nlmax) call grow_links(tsrc, tdst, tw, nlmax)
                                nl = nl + 1
                                tsrc(nl) = cs; tdst(nl) = c; tw(nl) = area
                            end if
                        end if
                    end if
                end do
            end do
        end do

        call weight_map_alloc(map%wm, MAP_WEIGHT, n1, n2, nl)
        map%wm%src(1:nl) = tsrc(1:nl); map%wm%dst(1:nl) = tdst(1:nl); map%wm%w(1:nl) = tw(1:nl)
        call weight_map_index(map%wm)

        call kdtree_free(tree)
        deallocate(emb, slon, slat, cidx, cd2, cidx_r, tsrc, tdst, tw)
        if (allocated(sv)) deallocate(sv)
        if (allocated(cv)) deallocate(cv)
        if (allocated(ov)) deallocate(ov)
    end subroutine conservative_spherical

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

    subroutine grow_links(tsrc, tdst, tw, cap)
        ! Double the temporary link buffers (safety for dense overlaps).
        integer,  allocatable, intent(inout) :: tsrc(:), tdst(:)
        real(dp), allocatable, intent(inout) :: tw(:)
        integer,               intent(inout) :: cap
        integer,  allocatable :: itmp(:)
        real(dp), allocatable :: rtmp(:)
        integer :: newcap
        newcap = cap*2
        allocate(itmp(newcap)); itmp(1:cap) = tsrc; call move_alloc(itmp, tsrc)
        allocate(itmp(newcap)); itmp(1:cap) = tdst; call move_alloc(itmp, tdst)
        allocate(rtmp(newcap)); rtmp(1:cap) = tw;   call move_alloc(rtmp, tw)
        cap = newcap
    end subroutine grow_links

    subroutine cell_corners_grid(grid, i, j, cx, cy)
        ! Four corners (CCW) of cell (i,j) in axis units, from adjacent-axis
        ! midpoints (so uniform and non-uniform/Gaussian grids share one path).
        type(grid_class), intent(in)  :: grid
        integer,          intent(in)  :: i, j
        real(dp),         intent(out) :: cx(4), cy(4)
        real(dp) :: x1, x2, y1, y2
        integer  :: nx, ny
        nx = grid%G%nx; ny = grid%G%ny

        if (i == 1) then
            x1 = grid%G%x(i) - 0.5_dp*(grid%G%x(i+1)-grid%G%x(i))
            x2 = 0.5_dp*(grid%G%x(i)+grid%G%x(i+1))
        else if (i == nx) then
            x1 = 0.5_dp*(grid%G%x(i)+grid%G%x(i-1))
            x2 = grid%G%x(i) + 0.5_dp*(grid%G%x(i)-grid%G%x(i-1))
        else
            x1 = 0.5_dp*(grid%G%x(i)+grid%G%x(i-1))
            x2 = 0.5_dp*(grid%G%x(i)+grid%G%x(i+1))
        end if

        if (j == 1) then
            y1 = grid%G%y(j) - 0.5_dp*(grid%G%y(j+1)-grid%G%y(j))
            y2 = 0.5_dp*(grid%G%y(j)+grid%G%y(j+1))
        else if (j == ny) then
            y1 = 0.5_dp*(grid%G%y(j)+grid%G%y(j-1))
            y2 = grid%G%y(j) + 0.5_dp*(grid%G%y(j)-grid%G%y(j-1))
        else
            y1 = 0.5_dp*(grid%G%y(j)+grid%G%y(j-1))
            y2 = 0.5_dp*(grid%G%y(j)+grid%G%y(j+1))
        end if

        cx = [x1, x2, x2, x1]
        cy = [y1, y1, y2, y2]
    end subroutine cell_corners_grid

end module conservative
