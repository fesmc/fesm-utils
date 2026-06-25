module mapping
    ! Build and apply mappings between two domains (the coords core).
    !
    ! map_init finds, for every target point, its nearest source points via a
    ! k-d tree and stores them as a MAP_DISTANCE weight_map (indices, distances,
    ! quadrants, neighbor positions). map_field applies that map to a field,
    ! choosing the per-point combination (nn / quadrant / shepard) at call time.
    !
    ! Metric selection mirrors the original library: if the source and target
    ! share a coordinate system and it is Cartesian/projected, neighbor search
    ! and distances are Euclidean in the projected plane; otherwise lon/lat are
    ! embedded on the unit sphere (chord distance is monotonic in geodesic
    ! distance) and exact distances are geodesic.

    use precision,       only: dp
    use coord_constants, only: degrees_to_radians, MISSING_VALUE_DEFAULT, ERR_DIST
    use coordinates,     only: points_class, grid_class, grid_axis_class, coord_system, &
                               grid_to_points, compare_coord
    use planet,          only: planet_distance, cartesian_distance, &
                               quadrant_latlon, quadrant_cartesian
    use kdtree
    use weight_map

    implicit none
    private

    type map_class
        character(len=256) :: name1, name2     ! source, target names
        type(coord_system) :: cs               ! target coordinate system
        logical :: is_grid = .false.
        type(grid_axis_class) :: G             ! target grid axes (when is_grid)
        integer :: npts = 0                    ! number of target points
        integer :: nmax = 0                    ! max neighbors stored per target
        logical :: is_same_map = .false.
        real(dp), allocatable :: x(:), y(:), lon(:), lat(:)   ! target coords (packed)
        type(weight_map_t) :: wm
    end type

    interface map_init
        module procedure map_init_grid_grid, map_init_grid_points
        module procedure map_init_points_grid, map_init_points_points
    end interface

    interface map_field
        module procedure map_field_grid_grid, map_field_points_points
        module procedure map_field_grid_points, map_field_points_grid
    end interface

    public :: map_class, map_init, map_field

contains

    ! ----- map_init overloads (convert to point sets, then delegate) ----------

    subroutine map_init_points_points(map, pts1, pts2, max_neighbors, dist_max)
        type(map_class),    intent(inout) :: map
        type(points_class), intent(in)    :: pts1, pts2
        integer,            intent(in)    :: max_neighbors
        real(dp), optional, intent(in)    :: dist_max
        call map_init_internal(map, pts1, pts2, max_neighbors, dist_max)
        map%is_grid = .false.
    end subroutine map_init_points_points

    subroutine map_init_grid_grid(map, grid1, grid2, max_neighbors, dist_max)
        type(map_class), intent(inout) :: map
        type(grid_class), intent(in)   :: grid1, grid2
        integer,         intent(in)    :: max_neighbors
        real(dp), optional, intent(in) :: dist_max
        type(points_class) :: pts1, pts2
        call grid_to_points(grid1, pts1, define=.true.)
        call grid_to_points(grid2, pts2, define=.true.)
        call map_init_internal(map, pts1, pts2, max_neighbors, dist_max)
        map%is_grid = .true.
        map%G = grid2%G
    end subroutine map_init_grid_grid

    subroutine map_init_grid_points(map, grid1, pts2, max_neighbors, dist_max)
        type(map_class), intent(inout)  :: map
        type(grid_class), intent(in)    :: grid1
        type(points_class), intent(in)  :: pts2
        integer,         intent(in)     :: max_neighbors
        real(dp), optional, intent(in)  :: dist_max
        type(points_class) :: pts1
        call grid_to_points(grid1, pts1, define=.true.)
        call map_init_internal(map, pts1, pts2, max_neighbors, dist_max)
        map%is_grid = .false.
    end subroutine map_init_grid_points

    subroutine map_init_points_grid(map, pts1, grid2, max_neighbors, dist_max)
        type(map_class), intent(inout)  :: map
        type(points_class), intent(in)  :: pts1
        type(grid_class), intent(in)    :: grid2
        integer,         intent(in)     :: max_neighbors
        real(dp), optional, intent(in)  :: dist_max
        type(points_class) :: pts2
        call grid_to_points(grid2, pts2, define=.true.)
        call map_init_internal(map, pts1, pts2, max_neighbors, dist_max)
        map%is_grid = .true.
        map%G = grid2%G
    end subroutine map_init_points_grid

    ! ----- core neighbor search ----------------------------------------------

    subroutine map_init_internal(map, pts1, pts2, max_neighbors, dist_max)
        type(map_class),    intent(inout) :: map
        type(points_class), intent(in)    :: pts1, pts2
        integer,            intent(in)    :: max_neighbors
        real(dp), optional, intent(in)    :: dist_max

        type(kdtree_t) :: tree
        logical  :: use_cart
        integer  :: n1, n2, ndim, i, c, j, nf, nl, nlmax
        integer  :: nq, nsel, s, best, tmpi
        real(dp) :: xyc, a, f, dmax, dist, q(3), tmpd
        real(dp), allocatable :: emb(:,:)
        integer,  allocatable :: knn_idx(:)
        real(dp), allocatable :: knn_d2(:)
        integer,  allocatable :: cand_idx(:)    ! candidate set, re-ranked by exact distance
        real(dp), allocatable :: cand_dist(:)
        ! temporary link store
        integer,  allocatable :: tsrc(:), tdst(:), tquad(:)
        real(dp), allocatable :: tdist(:), txn(:), tyn(:)

        n1 = pts1%npts
        n2 = pts2%npts
        map%name1 = pts1%name
        map%name2 = pts2%name
        map%cs    = pts2%cs
        map%npts  = n2
        map%nmax  = max_neighbors
        map%is_same_map = compare_coord(pts1, pts2)
        use_cart = map%is_same_map .and. pts2%cs%is_cartesian

        xyc = pts2%cs%xy_conv
        a   = pts2%cs%planet%a
        f   = pts2%cs%planet%f
        dmax = ERR_DIST
        if (present(dist_max)) dmax = dist_max

        ! Target coordinates (packed)
        if (allocated(map%x)) deallocate(map%x, map%y, map%lon, map%lat)
        allocate(map%x(n2), map%y(n2), map%lon(n2), map%lat(n2))
        map%x = pts2%x; map%y = pts2%y; map%lon = pts2%lon; map%lat = pts2%lat

        ! Build the source embedding for the k-d tree
        if (use_cart) then
            ndim = 2
            allocate(emb(2, n1))
            emb(1,:) = pts1%x * xyc
            emb(2,:) = pts1%y * xyc
        else
            ndim = 3
            allocate(emb(3, n1))
            call sphere_embed(pts1%lon, pts1%lat, emb)
        end if
        call kdtree_build(tree, emb)

        ! Over-fetch candidates by chord distance, then re-rank by the exact
        ! metric and keep the true max_neighbors nearest. This is needed because
        ! the spherical-chord embedding is not perfectly monotonic in the
        ! ellipsoidal geodesic distance, so the chord-kNN order/boundary can
        ! differ slightly from the exact result near ties.
        nq = min(n1, 2*max_neighbors + 8)
        allocate(knn_idx(nq), knn_d2(nq), cand_idx(nq), cand_dist(nq))
        nlmax = n2 * max_neighbors
        allocate(tsrc(nlmax), tdst(nlmax), tquad(nlmax), tdist(nlmax), txn(nlmax), tyn(nlmax))

        nl = 0
        do i = 1, n2
            if (use_cart) then
                q(1:2) = [pts2%x(i)*xyc, pts2%y(i)*xyc]
                call kdtree_knn(tree, q(1:2), nq, knn_idx, knn_d2, nf)
            else
                call sphere_embed_pt(pts2%lon(i), pts2%lat(i), q)
                call kdtree_knn(tree, q, nq, knn_idx, knn_d2, nf)
            end if

            ! exact distance for each candidate
            do c = 1, nf
                j = knn_idx(c)
                cand_idx(c) = j
                if (use_cart) then
                    cand_dist(c) = cartesian_distance(pts2%x(i)*xyc, pts2%y(i)*xyc, &
                                                      pts1%x(j)*xyc, pts1%y(j)*xyc)
                else
                    cand_dist(c) = planet_distance(a, f, pts2%lon(i), pts2%lat(i), &
                                                   pts1%lon(j), pts1%lat(j))
                end if
            end do

            ! selection-sort the first nsel candidates by exact distance
            nsel = min(nf, max_neighbors)
            do s = 1, nsel
                best = s
                do c = s+1, nf
                    if (cand_dist(c) < cand_dist(best)) best = c
                end do
                tmpd = cand_dist(s); cand_dist(s) = cand_dist(best); cand_dist(best) = tmpd
                tmpi = cand_idx(s);  cand_idx(s)  = cand_idx(best);  cand_idx(best)  = tmpi
            end do

            do s = 1, nsel
                j    = cand_idx(s)
                dist = cand_dist(s)
                if (dist > dmax) cycle
                nl = nl + 1
                tsrc(nl)  = j
                tdst(nl)  = i
                tdist(nl) = dist
                if (use_cart) then
                    tquad(nl) = quadrant_cartesian(pts2%x(i)*xyc, pts2%y(i)*xyc, &
                                                   pts1%x(j)*xyc, pts1%y(j)*xyc)
                else
                    tquad(nl) = quadrant_latlon(pts2%lon(i), pts2%lat(i), &
                                                pts1%lon(j), pts1%lat(j))
                end if
                ! Store neighbor position in the map's metric space (for bilinear):
                ! projected meters when Cartesian, otherwise lon/lat degrees.
                if (use_cart) then
                    txn(nl) = pts1%x(j)*xyc
                    tyn(nl) = pts1%y(j)*xyc
                else
                    txn(nl) = pts1%lon(j)
                    tyn(nl) = pts1%lat(j)
                end if
            end do
        end do

        ! Assemble the weight_map (links already grouped by target i)
        call weight_map_alloc(map%wm, MAP_DISTANCE, n1, n2, nl)
        map%wm%src(1:nl)      = tsrc(1:nl)
        map%wm%dst(1:nl)      = tdst(1:nl)
        map%wm%dist(1:nl)     = tdist(1:nl)
        map%wm%quadrant(1:nl) = tquad(1:nl)
        map%wm%xn(1:nl)       = txn(1:nl)
        map%wm%yn(1:nl)       = tyn(1:nl)
        call weight_map_index(map%wm)

        call kdtree_free(tree)
        deallocate(emb, knn_idx, knn_d2, cand_idx, cand_dist, tsrc, tdst, tquad, tdist, txn, tyn)
    end subroutine map_init_internal

    subroutine sphere_embed(lon, lat, emb)
        real(dp), intent(in)  :: lon(:), lat(:)
        real(dp), intent(out) :: emb(:,:)
        integer  :: i
        real(dp) :: la, lo
        do i = 1, size(lon)
            la = lat(i)*degrees_to_radians
            lo = lon(i)*degrees_to_radians
            emb(1,i) = cos(la)*cos(lo)
            emb(2,i) = cos(la)*sin(lo)
            emb(3,i) = sin(la)
        end do
    end subroutine sphere_embed

    subroutine sphere_embed_pt(lon, lat, v)
        real(dp), intent(in)  :: lon, lat
        real(dp), intent(out) :: v(3)
        real(dp) :: la, lo
        la = lat*degrees_to_radians
        lo = lon*degrees_to_radians
        v(1) = cos(la)*cos(lo)
        v(2) = cos(la)*sin(lo)
        v(3) = sin(la)
    end subroutine sphere_embed_pt

    ! ----- map_field overloads ------------------------------------------------

    subroutine map_field_points_points(map, name, var1, var2, method, missing_value, radius, mask2)
        type(map_class),  intent(in)  :: map
        character(len=*), intent(in)  :: name
        real(dp),         intent(in)  :: var1(:)
        real(dp),         intent(out) :: var2(:)
        character(len=*), intent(in)  :: method
        real(dp), optional, intent(in)  :: missing_value
        real(dp), optional, intent(in)  :: radius
        logical,  optional, intent(out) :: mask2(:)
        call map_apply_vec(map, var1, var2, method, missing_value, radius, mask2)
    end subroutine map_field_points_points

    subroutine map_field_grid_grid(map, name, var1, var2, method, missing_value, radius, mask2)
        type(map_class),  intent(in)  :: map
        character(len=*), intent(in)  :: name
        real(dp),         intent(in)  :: var1(:,:)
        real(dp),         intent(out) :: var2(:,:)
        character(len=*), intent(in)  :: method
        real(dp), optional, intent(in)  :: missing_value
        real(dp), optional, intent(in)  :: radius
        logical,  optional, intent(out) :: mask2(:,:)
        real(dp), allocatable :: v1(:), v2(:)
        logical,  allocatable :: m2(:)
        integer :: n1, n2
        n1 = size(var1); n2 = size(var2)
        allocate(v1(n1), v2(n2), m2(n2))
        v1 = reshape(var1, [n1])
        call map_apply_vec(map, v1, v2, method, missing_value, radius, m2)
        var2 = reshape(v2, shape(var2))
        if (present(mask2)) mask2 = reshape(m2, shape(mask2))
        deallocate(v1, v2, m2)
    end subroutine map_field_grid_grid

    subroutine map_field_grid_points(map, name, var1, var2, method, missing_value, radius, mask2)
        ! source on a grid (2D), target a point set (1D)
        type(map_class),  intent(in)  :: map
        character(len=*), intent(in)  :: name
        real(dp),         intent(in)  :: var1(:,:)
        real(dp),         intent(out) :: var2(:)
        character(len=*), intent(in)  :: method
        real(dp), optional, intent(in)  :: missing_value
        real(dp), optional, intent(in)  :: radius
        logical,  optional, intent(out) :: mask2(:)
        real(dp), allocatable :: v1(:)
        integer :: n1
        n1 = size(var1)
        allocate(v1(n1))
        v1 = reshape(var1, [n1])
        call map_apply_vec(map, v1, var2, method, missing_value, radius, mask2)
        deallocate(v1)
    end subroutine map_field_grid_points

    subroutine map_field_points_grid(map, name, var1, var2, method, missing_value, radius, mask2)
        ! source a point set (1D), target on a grid (2D)
        type(map_class),  intent(in)  :: map
        character(len=*), intent(in)  :: name
        real(dp),         intent(in)  :: var1(:)
        real(dp),         intent(out) :: var2(:,:)
        character(len=*), intent(in)  :: method
        real(dp), optional, intent(in)  :: missing_value
        real(dp), optional, intent(in)  :: radius
        logical,  optional, intent(out) :: mask2(:,:)
        real(dp), allocatable :: v2(:)
        logical,  allocatable :: m2(:)
        integer :: n2
        n2 = size(var2)
        allocate(v2(n2), m2(n2))
        call map_apply_vec(map, var1, v2, method, missing_value, radius, m2)
        var2 = reshape(v2, shape(var2))
        if (present(mask2)) mask2 = reshape(m2, shape(mask2))
        deallocate(v2, m2)
    end subroutine map_field_points_grid

    ! ----- apply dispatch (bilinear handled here; rest delegate to weight_map) -

    subroutine map_apply_vec(map, var1, var2, method, missing_value, radius, mask2)
        type(map_class),  intent(in)  :: map
        real(dp),         intent(in)  :: var1(:)
        real(dp),         intent(out) :: var2(:)
        character(len=*), intent(in)  :: method
        real(dp), optional, intent(in)  :: missing_value
        real(dp), optional, intent(in)  :: radius
        logical,  optional, intent(out) :: mask2(:)
        if (trim(method) == "bilinear" .or. trim(method) == "bilin") then
            call bilinear_apply(map, var1, var2, missing_value, mask2)
        else
            call weight_map_apply(map%wm, var1, var2, method, missing_value, radius, mask2)
        end if
    end subroutine map_apply_vec

    subroutine bilinear_apply(map, var1, var2, missing_value, mask2)
        ! Bilinear from the 4 quadrant corners. If all 4 are present and valid,
        ! true bilinear; otherwise graceful degradation to quadrant inverse-
        ! distance over the surviving corners (default per coords-design Q7).
        type(map_class),  intent(in)  :: map
        real(dp),         intent(in)  :: var1(:)
        real(dp),         intent(out) :: var2(:)
        real(dp), optional, intent(in)  :: missing_value
        logical,  optional, intent(out) :: mask2(:)

        real(dp) :: miss, a, f, xyc
        logical  :: use_cart
        integer  :: k, j, j1, j2, q
        integer  :: iq(4)                  ! link position of nearest valid neighbor per quadrant
        real(dp) :: dq(4)
        real(dp) :: z1, z2, z3, z4, tx, ty
        real(dp) :: alpha1, alpha2, alpha3, ymid1, ymid0
        real(dp) :: dx1, dx1t, dx2, dx2t, dy1, dy1t, p1, p0
        real(dp) :: wsum, vsum, w

        miss = MISSING_VALUE_DEFAULT
        if (present(missing_value)) miss = missing_value
        use_cart = map%is_same_map .and. map%cs%is_cartesian
        xyc = map%cs%xy_conv
        a   = map%cs%planet%a
        f   = map%cs%planet%f

        var2 = miss
        if (present(mask2)) mask2 = .false.

        do k = 1, map%wm%n_dst
            j1 = map%wm%dst_off(k); j2 = map%wm%dst_off(k+1) - 1
            if (j2 < j1) cycle

            ! nearest valid neighbor in each quadrant
            iq = 0; dq = huge(1.0_dp)
            do j = j1, j2
                q = map%wm%quadrant(j)
                if (q >= 1 .and. q <= 4) then
                    if (var1(map%wm%src(j)) /= miss .and. map%wm%dist(j) < dq(q)) then
                        dq(q) = map%wm%dist(j); iq(q) = j
                    end if
                end if
            end do

            if (all(iq > 0)) then
                ! all four corners -> true bilinear
                z1 = var1(map%wm%src(iq(1)))
                z2 = var1(map%wm%src(iq(2)))
                z3 = var1(map%wm%src(iq(3)))
                z4 = var1(map%wm%src(iq(4)))
                if (use_cart) then
                    tx = map%x(k)*xyc; ty = map%y(k)*xyc
                else
                    tx = map%lon(k);   ty = map%lat(k)
                end if
                ymid1 = 0.5_dp*(map%wm%yn(iq(2)) + map%wm%yn(iq(1)))
                ymid0 = 0.5_dp*(map%wm%yn(iq(3)) + map%wm%yn(iq(4)))
                if (use_cart) then
                    dx1  = cartesian_distance(tx, ymid1, map%wm%xn(iq(2)), map%wm%yn(iq(2)))
                    dx1t = cartesian_distance(map%wm%xn(iq(1)), map%wm%yn(iq(1)), map%wm%xn(iq(2)), map%wm%yn(iq(2)))
                    dx2  = cartesian_distance(tx, ymid0, map%wm%xn(iq(3)), map%wm%yn(iq(3)))
                    dx2t = cartesian_distance(map%wm%xn(iq(4)), map%wm%yn(iq(4)), map%wm%xn(iq(3)), map%wm%yn(iq(3)))
                    dy1  = cartesian_distance(tx, ty, tx, ymid0)
                    dy1t = cartesian_distance(tx, ymid1, tx, ymid0)
                else
                    dx1  = planet_distance(a, f, tx, ymid1, map%wm%xn(iq(2)), map%wm%yn(iq(2)))
                    dx1t = planet_distance(a, f, map%wm%xn(iq(1)), map%wm%yn(iq(1)), map%wm%xn(iq(2)), map%wm%yn(iq(2)))
                    dx2  = planet_distance(a, f, tx, ymid0, map%wm%xn(iq(3)), map%wm%yn(iq(3)))
                    dx2t = planet_distance(a, f, map%wm%xn(iq(4)), map%wm%yn(iq(4)), map%wm%xn(iq(3)), map%wm%yn(iq(3)))
                    dy1  = planet_distance(a, f, tx, ty, tx, ymid0)
                    dy1t = planet_distance(a, f, tx, ymid1, tx, ymid0)
                end if
                alpha1 = 0.0_dp; if (dx1t > 0.0_dp) alpha1 = dx1/dx1t
                alpha2 = 0.0_dp; if (dx2t > 0.0_dp) alpha2 = dx2/dx2t
                alpha3 = 0.0_dp; if (dy1t > 0.0_dp) alpha3 = dy1/dy1t
                p1 = z2 + alpha1*(z1 - z2)
                p0 = z3 + alpha2*(z4 - z3)
                var2(k) = p0 + alpha3*(p1 - p0)
                if (present(mask2)) mask2(k) = .true.
            else
                ! graceful degradation: quadrant inverse-distance over survivors
                wsum = 0.0_dp; vsum = 0.0_dp
                do q = 1, 4
                    if (iq(q) > 0) then
                        w = 1.0_dp / max(dq(q), 1.0e-12_dp)**2
                        wsum = wsum + w
                        vsum = vsum + w*var1(map%wm%src(iq(q)))
                    end if
                end do
                if (wsum > 0.0_dp) then
                    var2(k) = vsum/wsum
                    if (present(mask2)) mask2(k) = .true.
                end if
            end if
        end do
    end subroutine bilinear_apply

end module mapping
