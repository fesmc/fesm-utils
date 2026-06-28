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

    use precision,       only: dp, sp
    use constants, only: degrees_to_radians, mv_dp, ERR_DIST
    use coordinates,     only: points_class, grid_class, grid_axis_class, coord_system, &
                               grid_to_points, compare_coord, grid_write
    use planet,          only: planet_distance, cartesian_distance, &
                               quadrant_latlon, quadrant_cartesian
    use kdtree
    use weight_map
    use mapping_scrip,   only: map_scrip_class, map_scrip_load, map_scrip_to_weight_map, &
                               gen_map_filename
    use conservative,    only: conservative_weights
    use interp2D,        only: fill_weighted, fill_nearest, filter_poisson
    use gaussian_filter, only: filter_gaussian, filter_gaussian_fast
    use map_io,          only: weight_map_write, weight_map_read
    use ncio,            only: nc_exists_var
    use grid_to_cdo,     only: grid_cdo_write_desc_short, call_system_cdo

    implicit none
    private

    type map_class
        character(len=256) :: name1, name2     ! source, target names
        type(coord_system) :: cs               ! target coordinate system
        logical :: is_grid = .false.
        type(grid_axis_class) :: G             ! target grid axes (when is_grid)
        integer :: npts = 0                    ! number of target points
        integer :: nmax = 0                    ! max neighbors stored per target
        character(len=32) :: method = "nn"     ! interpolation kernel, fixed at init
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
        module procedure map_field_grid_grid_sp,    map_field_grid_grid_int
        module procedure map_field_points_grid_sp,  map_field_points_grid_int
        module procedure map_field_grid_points_sp,  map_field_grid_points_int
        module procedure map_field_points_points_sp, map_field_points_points_int
    end interface

    public :: map_class, map_init, map_field, map_read, map_init_conservative

contains

    ! ----- map_init overloads (convert to point sets, then delegate) ----------

    subroutine map_init_points_points(map, pts1, pts2, max_neighbors, dist_max, method)
        type(map_class),    intent(inout) :: map
        type(points_class), intent(in)    :: pts1, pts2
        integer,            intent(in)    :: max_neighbors
        real(dp), optional, intent(in)    :: dist_max
        character(len=*), optional, intent(in) :: method
        call map_init_internal(map, pts1, pts2, max_neighbors, dist_max, method)
        map%is_grid = .false.
    end subroutine map_init_points_points

    subroutine map_init_grid_grid(map, grid1, grid2, max_neighbors, dist_max, method, gen, fldr, load, clean)
        ! Build (or load) a grid -> grid map and cache it as a SCRIP(-superset)
        ! file under `fldr` (default "maps").
        !   method - interpolation kernel: "con" (conservative) or a distance
        !            kernel nn/shepard/quadrant/bilinear (default "nn").
        !   gen    - generator used when the weights must be made: "coords"
        !            (in-package, default) or "cdo" (external cdo system call).
        !   load   - if .true. (default) and a cache file already exists, it is
        !            loaded instead of regenerated; freshly generated maps are
        !            always saved to the cache. load=.false. forces regeneration.
        ! The loader auto-detects the file format, so a cdo-generated SCRIP file
        ! and an in-package (coords) cache interoperate at the same filename.
        type(map_class), intent(inout) :: map
        type(grid_class), intent(in)   :: grid1, grid2
        integer,  optional, intent(in) :: max_neighbors
        real(dp), optional, intent(in) :: dist_max
        character(len=*), optional, intent(in) :: method
        character(len=*), optional, intent(in) :: gen
        character(len=*), optional, intent(in) :: fldr
        logical,  optional, intent(in) :: load
        logical,  optional, intent(in) :: clean

        type(points_class) :: pts1, pts2
        character(len=32)  :: mtd, g
        character(len=256) :: mfldr, filename
        logical :: load_file, file_exists
        integer :: nmax

        mtd = "nn"
        if (present(method)) mtd = trim(method)
        g = "coords"
        if (present(gen)) g = trim(gen)
        mfldr = "maps"
        if (present(fldr)) mfldr = trim(fldr)
        load_file = .true.
        if (present(load)) load_file = load

        filename = gen_map_filename(grid1%name, grid2%name, mfldr, mtd)
        inquire(file=trim(filename), exist=file_exists)

        ! Load a cached map if one exists and loading is allowed
        if (load_file .and. file_exists) then
            call map_set_target_from_grid(map, grid1, grid2)
            map%method = trim(mtd)
            call map_load_wm(map, grid1, grid2, mtd, mfldr, filename)
            return
        end if

        ! Otherwise generate the weights
        select case (trim(g))
            case ("coords")
                if (trim(mtd) == "con") then
                    call map_init_conservative(map, grid1, grid2)
                else
                    nmax = 10
                    if (present(max_neighbors)) nmax = max_neighbors
                    call grid_to_points(grid1, pts1, define=.true.)
                    call grid_to_points(grid2, pts2, define=.true.)
                    call map_init_internal(map, pts1, pts2, nmax, dist_max, mtd)
                    map%is_grid = .true.
                    map%G = grid2%G
                end if
                ! Cache the freshly generated map
                call map_save_wm(map, mfldr, filename)
            case ("cdo")
                ! Generate the SCRIP map with an external cdo call (it writes
                ! `filename` itself, so no separate save step is needed).
                call map_init_cdo(map, grid1, grid2, mtd, mfldr, filename, clean)
            case default
                write(*,*) "mapping:: map_init: unknown gen '"//trim(g)//"'"
                write(*,*) "  expected 'coords' (in-package) or 'cdo' (external cdo call)."
                stop
        end select
    end subroutine map_init_grid_grid

    subroutine map_set_target_from_grid(map, grid1, grid2)
        ! Assemble the target (grid2) metadata on a grid -> grid map_class.
        ! Shared by the conservative generator and the cache loader.
        type(map_class),  intent(inout) :: map
        type(grid_class), intent(in)    :: grid1, grid2
        integer :: n2

        n2 = grid2%G%nx * grid2%G%ny
        map%name1 = grid1%name
        map%name2 = grid2%name
        map%cs    = grid2%cs
        map%is_grid = .true.
        map%G     = grid2%G
        map%npts  = n2
        map%nmax  = 0
        map%is_same_map = compare_coord(grid1, grid2)
        if (allocated(map%x)) deallocate(map%x, map%y, map%lon, map%lat)
        allocate(map%x(n2), map%y(n2), map%lon(n2), map%lat(n2))
        map%x   = reshape(grid2%x,   [n2])
        map%y   = reshape(grid2%y,   [n2])
        map%lon = reshape(grid2%lon, [n2])
        map%lat = reshape(grid2%lat, [n2])
    end subroutine map_set_target_from_grid

    subroutine map_init_conservative(map, grid1, grid2)
        ! Build a conservative (area-weighted) grid -> grid map as a MAP_WEIGHT
        ! using the in-package analytic polygon clip.
        type(map_class),  intent(inout) :: map
        type(grid_class), intent(in)    :: grid1, grid2

        call map_set_target_from_grid(map, grid1, grid2)
        map%method = "con"
        call conservative_weights(map%wm, grid1, grid2)
    end subroutine map_init_conservative

    subroutine map_load_wm(map, grid1, grid2, method, fldr, filename)
        ! Load a cached weight_map into map%wm, auto-detecting the file format:
        ! an in-package cache (weight_map_write, carries a "meta" variable) vs a
        ! standard/cdo SCRIP file (read via the mapping_scrip bridge).
        type(map_class),  intent(inout) :: map
        type(grid_class), intent(in)    :: grid1, grid2
        character(len=*), intent(in)    :: method, fldr, filename
        type(map_scrip_class) :: mps

        if (nc_exists_var(trim(filename), "meta")) then
            call weight_map_read(map%wm, trim(filename))
        else
            call map_scrip_load(mps, trim(grid1%name), trim(grid2%name), trim(fldr), trim(method))
            call map_scrip_to_weight_map(mps, map%wm)
        end if
    end subroutine map_load_wm

    subroutine map_save_wm(map, fldr, filename)
        ! Write map%wm to its cache file, creating the folder if needed.
        type(map_class),  intent(in) :: map
        character(len=*), intent(in) :: fldr, filename
        call execute_command_line("mkdir -p '"//trim(fldr)//"'")
        call weight_map_write(map%wm, trim(filename))
    end subroutine map_save_wm

    subroutine map_init_cdo(map, grid1, grid2, method, fldr, filename, clean)
        ! Generate a grid -> grid SCRIP map via an external cdo call and load it.
        ! This reproduces the climber-x flow: write short grid descriptions for
        ! both grids, write a source-grid NetCDF for cdo input, then
        !   cdo gen<method>,<dst.txt> -setgrid,<src.txt> <src.nc> <filename>
        ! and load the resulting SCRIP file. cdo writes `filename` directly, so
        ! the map is its own cache.
        type(map_class),  intent(inout) :: map
        type(grid_class), intent(in)    :: grid1, grid2
        character(len=*), intent(in)    :: method, fldr, filename
        logical, optional, intent(in)   :: clean

        type(map_scrip_class) :: mps
        character(len=512)     :: src_nc, desc1, desc2, cmd
        character(len=12)      :: xnm, ynm
        logical                :: do_clean

        do_clean = .false.
        if (present(clean)) do_clean = clean

        ! Target metadata (grid2); the kernel is baked into the cdo weights
        call map_set_target_from_grid(map, grid1, grid2)
        map%method = trim(method)

        call execute_command_line("mkdir -p '"//trim(fldr)//"'")

        ! Short grid descriptions for source and target
        call grid_cdo_write_desc_short(grid1, trim(fldr))
        call grid_cdo_write_desc_short(grid2, trim(fldr))

        ! Source-grid NetCDF for cdo input (xc/yc for projected/cartesian grids,
        ! lon/lat otherwise); -setgrid below reattaches the source description.
        src_nc = trim(fldr)//"/grid_"//trim(grid1%name)//".nc"
        if ((.not. grid1%cs%is_projection) .and. (.not. grid1%cs%is_cartesian)) then
            xnm = "lon"; ynm = "lat"
        else
            xnm = "xc";  ynm = "yc"
        end if
        call grid_write(grid1, trim(src_nc), trim(xnm), trim(ynm), .true.)

        ! Build and run the cdo command
        desc1 = trim(fldr)//"/grid_"//trim(grid1%name)//".txt"
        desc2 = trim(fldr)//"/grid_"//trim(grid2%name)//".txt"
        cmd = "cdo gen"//trim(method)//","//trim(desc2)//" -setgrid,"//trim(desc1)// &
              " "//trim(src_nc)//" "//trim(filename)
        call call_system_cdo(cmd)

        ! Load the generated SCRIP map into the weight store
        call map_scrip_load(mps, trim(grid1%name), trim(grid2%name), trim(fldr), trim(method))
        call map_scrip_to_weight_map(mps, map%wm)

        ! Optionally remove the intermediate grid files
        if (do_clean) then
            call execute_command_line("rm -f '"//trim(src_nc)//"' '"// &
                                      trim(desc1)//"' '"//trim(desc2)//"'")
        end if
    end subroutine map_init_cdo

    subroutine map_init_grid_points(map, grid1, pts2, max_neighbors, dist_max, method)
        type(map_class), intent(inout)  :: map
        type(grid_class), intent(in)    :: grid1
        type(points_class), intent(in)  :: pts2
        integer,         intent(in)     :: max_neighbors
        real(dp), optional, intent(in)  :: dist_max
        character(len=*), optional, intent(in) :: method
        type(points_class) :: pts1
        call grid_to_points(grid1, pts1, define=.true.)
        call map_init_internal(map, pts1, pts2, max_neighbors, dist_max, method)
        map%is_grid = .false.
    end subroutine map_init_grid_points

    subroutine map_init_points_grid(map, pts1, grid2, max_neighbors, dist_max, method)
        type(map_class), intent(inout)  :: map
        type(points_class), intent(in)  :: pts1
        type(grid_class), intent(in)    :: grid2
        integer,         intent(in)     :: max_neighbors
        real(dp), optional, intent(in)  :: dist_max
        character(len=*), optional, intent(in) :: method
        type(points_class) :: pts2
        call grid_to_points(grid2, pts2, define=.true.)
        call map_init_internal(map, pts1, pts2, max_neighbors, dist_max, method)
        map%is_grid = .true.
        map%G = grid2%G
    end subroutine map_init_points_grid

    ! ----- core neighbor search ----------------------------------------------

    subroutine map_init_internal(map, pts1, pts2, max_neighbors, dist_max, method)
        type(map_class),    intent(inout) :: map
        type(points_class), intent(in)    :: pts1, pts2
        integer,            intent(in)    :: max_neighbors
        real(dp), optional, intent(in)    :: dist_max
        character(len=*), optional, intent(in) :: method   ! distance kernel (default "shepard")

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
        map%method = "shepard"
        if (present(method)) map%method = trim(method)
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

    ! ----- load a map from a SCRIP file --------------------------------------

    subroutine map_read(map, src_name, dst_name, fldr, method)
        ! Stage A: load a pre-existing (CDO) SCRIP map file into a map_class as a
        ! MAP_WEIGHT, via the mapping_scrip reader + the weight_map bridge. The
        ! target grid dimensions are recovered from dst_grid_dims so the map can
        ! be applied with map_field and used to size target arrays (map%G%nx/ny).
        type(map_class),  intent(inout) :: map
        character(len=*), intent(in)    :: src_name, dst_name, fldr, method
        type(map_scrip_class) :: mps

        call map_scrip_load(mps, src_name, dst_name, fldr, method)
        call map_scrip_to_weight_map(mps, map%wm)

        map%name1 = trim(src_name)
        map%name2 = trim(dst_name)
        map%npts  = map%wm%n_dst
        map%nmax  = 0
        map%method = trim(method)

        if (allocated(mps%dst_grid_dims) .and. size(mps%dst_grid_dims) >= 2) then
            map%is_grid = .true.
            map%G%nx = mps%dst_grid_dims(1)
            map%G%ny = mps%dst_grid_dims(2)
        else
            map%is_grid = .false.
        end if
    end subroutine map_read

    ! ----- map_field overloads ------------------------------------------------

    subroutine map_field_to_grid(map, name, v1, var2, stat, method, missing_value, radius, mask2, &
                                 reset, mask_pack, fill_method, filt_method, filt_par, verbose)
        ! Core 2D-target worker: the source is already flattened to the vector v1,
        ! so this serves both grid (2D) and points (1D) sources. It carries the
        ! full post-apply surface (matching the legacy map_scrip_field):
        !   stat        - aggregation over each target's links ("mean" [default]/
        !                 "count"/"stdev").
        !   method      - optional interpolation-kernel override (distance maps
        !                 only: nn/shepard/quadrant/bilinear); defaults to the
        !                 kernel fixed at map_init (map%method).
        !   reset       - .true. (default) blanks var2 to missing first; .false.
        !                 leaves non-interpolated target points untouched (var2 is
        !                 read on entry, hence intent(inout)).
        !   mask_pack   - only interpolate where .true.; others keep the reset value.
        !   fill_method - fill remaining missing target points ("weighted"/"nn"/"none").
        !   filt_method - smooth the result ("gaussian"/"gaussian-fast"/"poisson"/"none").
        type(map_class),  intent(in)    :: map
        character(len=*), intent(in)    :: name
        real(dp),         intent(in)    :: v1(:)
        real(dp),         intent(inout) :: var2(:,:)
        character(len=*), optional, intent(in)  :: stat
        character(len=*), optional, intent(in)  :: method
        real(dp),         optional, intent(in)  :: missing_value
        real(dp),         optional, intent(in)  :: radius
        logical,          optional, intent(out) :: mask2(:,:)
        logical,          optional, intent(in)  :: reset
        logical,          optional, intent(in)  :: mask_pack(:,:)
        character(len=*), optional, intent(in)  :: fill_method
        character(len=*), optional, intent(in)  :: filt_method
        real(dp),         optional, intent(in)  :: filt_par(:)
        logical,          optional, intent(in)  :: verbose

        character(len=32) :: sta, knl
        real(dp) :: miss
        logical  :: reset_pts
        integer  :: n2, nx, ny
        real(dp), allocatable :: v2(:), var2_save(:,:)
        logical,  allocatable :: m2(:), m2d(:,:), maskp2d(:,:)

        sta = "mean"
        if (present(stat)) sta = trim(stat)
        knl = map_resolve_kernel(map, method)
        miss = mv_dp
        if (present(missing_value)) miss = missing_value
        reset_pts = .true.
        if (present(reset)) reset_pts = reset

        nx = size(var2,1); ny = size(var2,2)
        n2 = size(var2)

        ! Target-point packing mask (which points to interpolate)
        allocate(maskp2d(nx,ny))
        maskp2d = .true.
        if (present(mask_pack)) maskp2d = mask_pack

        ! Preserve existing target values when not resetting
        if (.not. reset_pts) then
            allocate(var2_save(nx,ny))
            var2_save = var2
        end if

        ! Apply the map over all targets, then restrict to the packing mask
        allocate(v2(n2), m2(n2), m2d(nx,ny))
        call map_apply_vec(map, v1, v2, knl, sta, miss, radius, m2)
        var2 = reshape(v2, [nx,ny])
        m2d  = reshape(m2, [nx,ny]) .and. maskp2d

        if (reset_pts) then
            where (.not. m2d) var2 = miss
        else
            where (.not. m2d) var2 = var2_save
        end if

        if (present(mask2)) mask2 = m2d

        ! Optional fill + filter post-processing
        if (present(fill_method) .or. present(filt_method)) then
            call map_postprocess(var2, miss, maskp2d, name, fill_method, filt_method, filt_par, verbose)
        end if

        deallocate(v2, m2, m2d, maskp2d)
        if (allocated(var2_save)) deallocate(var2_save)
    end subroutine map_field_to_grid

    subroutine map_field_to_points(map, name, v1, var2, stat, method, missing_value, radius, mask2, &
                                   reset, mask_pack)
        ! Core 1D-target worker: the source is already flattened to v1. 1D targets
        ! support reset/mask_pack only (fill/filter are inherently 2D grid ops):
        !   stat      - aggregation over each target's links ("mean"/"count"/"stdev").
        !   method    - optional interpolation-kernel override (distance maps only);
        !               defaults to the kernel fixed at map_init (map%method).
        !   reset     - .true. (default) blanks var2 to missing first; .false.
        !               leaves non-interpolated target points untouched.
        !   mask_pack - only interpolate where .true.; others keep the reset value.
        type(map_class),  intent(in)    :: map
        character(len=*), intent(in)    :: name
        real(dp),         intent(in)    :: v1(:)
        real(dp),         intent(inout) :: var2(:)
        character(len=*), optional, intent(in)  :: stat
        character(len=*), optional, intent(in)  :: method
        real(dp),         optional, intent(in)  :: missing_value
        real(dp),         optional, intent(in)  :: radius
        logical,          optional, intent(out) :: mask2(:)
        logical,          optional, intent(in)  :: reset
        logical,          optional, intent(in)  :: mask_pack(:)

        character(len=32) :: sta, knl
        real(dp) :: miss
        logical  :: reset_pts
        integer  :: n2
        real(dp), allocatable :: v2(:), var2_save(:)
        logical,  allocatable :: m2(:), maskp(:)

        sta = "mean"
        if (present(stat)) sta = trim(stat)
        knl = map_resolve_kernel(map, method)
        miss = mv_dp
        if (present(missing_value)) miss = missing_value
        reset_pts = .true.
        if (present(reset)) reset_pts = reset

        n2 = size(var2)

        allocate(maskp(n2))
        maskp = .true.
        if (present(mask_pack)) maskp = mask_pack

        if (.not. reset_pts) then
            allocate(var2_save(n2))
            var2_save = var2
        end if

        allocate(v2(n2), m2(n2))
        call map_apply_vec(map, v1, v2, knl, sta, miss, radius, m2)
        var2 = v2
        m2   = m2 .and. maskp

        if (reset_pts) then
            where (.not. m2) var2 = miss
        else
            where (.not. m2) var2 = var2_save
        end if

        if (present(mask2)) mask2 = m2

        deallocate(v2, m2, maskp)
        if (allocated(var2_save)) deallocate(var2_save)
    end subroutine map_field_to_points

    subroutine map_field_grid_grid(map, name, var1, var2, stat, method, missing_value, radius, mask2, &
                                   reset, mask_pack, fill_method, filt_method, filt_par, verbose)
        ! Map a 2D source field var1 onto the 2D target var2 (dp).
        type(map_class),  intent(in)    :: map
        character(len=*), intent(in)    :: name
        real(dp),         intent(in)    :: var1(:,:)
        real(dp),         intent(inout) :: var2(:,:)
        character(len=*), optional, intent(in)  :: stat
        character(len=*), optional, intent(in)  :: method
        real(dp),         optional, intent(in)  :: missing_value
        real(dp),         optional, intent(in)  :: radius
        logical,          optional, intent(out) :: mask2(:,:)
        logical,          optional, intent(in)  :: reset
        logical,          optional, intent(in)  :: mask_pack(:,:)
        character(len=*), optional, intent(in)  :: fill_method
        character(len=*), optional, intent(in)  :: filt_method
        real(dp),         optional, intent(in)  :: filt_par(:)
        logical,          optional, intent(in)  :: verbose

        real(dp), allocatable :: v1(:)
        allocate(v1(size(var1)))
        v1 = reshape(var1, [size(var1)])
        call map_field_to_grid(map, name, v1, var2, stat, method, missing_value, radius, mask2, &
                               reset, mask_pack, fill_method, filt_method, filt_par, verbose)
        deallocate(v1)
    end subroutine map_field_grid_grid

    subroutine map_field_points_points(map, name, var1, var2, stat, method, missing_value, radius, mask2, &
                                       reset, mask_pack)
        ! source a point set (1D), target a point set (1D)
        type(map_class),  intent(in)    :: map
        character(len=*), intent(in)    :: name
        real(dp),         intent(in)    :: var1(:)
        real(dp),         intent(inout) :: var2(:)
        character(len=*), optional, intent(in)  :: stat
        character(len=*), optional, intent(in)  :: method
        real(dp),         optional, intent(in)  :: missing_value
        real(dp),         optional, intent(in)  :: radius
        logical,          optional, intent(out) :: mask2(:)
        logical,          optional, intent(in)  :: reset
        logical,          optional, intent(in)  :: mask_pack(:)
        call map_field_to_points(map, name, var1, var2, stat, method, missing_value, radius, mask2, &
                                 reset, mask_pack)
    end subroutine map_field_points_points

    subroutine map_field_grid_grid_sp(map, name, var1, var2, stat, method, missing_value, radius, mask2, &
                                      reset, mask_pack, fill_method, filt_method, filt_par, verbose)
        ! Single-precision wrapper: accumulate in dp, store back in sp.
        type(map_class),  intent(in)    :: map
        character(len=*), intent(in)    :: name
        real(sp),         intent(in)    :: var1(:,:)
        real(sp),         intent(inout) :: var2(:,:)
        character(len=*), optional, intent(in)  :: stat
        character(len=*), optional, intent(in)  :: method
        real(sp),         optional, intent(in)  :: missing_value
        real(dp),         optional, intent(in)  :: radius
        logical,          optional, intent(out) :: mask2(:,:)
        logical,          optional, intent(in)  :: reset
        logical,          optional, intent(in)  :: mask_pack(:,:)
        character(len=*), optional, intent(in)  :: fill_method
        character(len=*), optional, intent(in)  :: filt_method
        real(sp),         optional, intent(in)  :: filt_par(:)
        logical,          optional, intent(in)  :: verbose

        real(dp), allocatable :: v1(:,:), v2(:,:), fpar(:)
        logical,  allocatable :: m2(:,:)
        real(dp) :: miss
        integer  :: nx, ny

        nx = size(var2,1); ny = size(var2,2)
        allocate(v1(size(var1,1),size(var1,2)), v2(nx,ny), m2(nx,ny))
        v1 = real(var1, dp)
        v2 = real(var2, dp)
        miss = mv_dp
        if (present(missing_value)) miss = real(missing_value, dp)
        if (present(filt_par)) then
            allocate(fpar(size(filt_par)))
            fpar = real(filt_par, dp)
        end if

        call map_field_grid_grid(map, name, v1, v2, stat=stat, method=method, missing_value=miss, radius=radius, &
                                 mask2=m2, reset=reset, mask_pack=mask_pack, fill_method=fill_method, &
                                 filt_method=filt_method, filt_par=fpar, verbose=verbose)

        var2 = real(v2, sp)
        if (present(mask2)) mask2 = m2
        deallocate(v1, v2, m2)
        if (allocated(fpar)) deallocate(fpar)
    end subroutine map_field_grid_grid_sp

    subroutine map_field_grid_grid_int(map, name, var1, var2, stat, method, missing_value, radius, mask2, &
                                       reset, mask_pack, fill_method, filt_method, filt_par, verbose)
        ! Integer wrapper: accumulate in dp, round back to integer.
        type(map_class),  intent(in)    :: map
        character(len=*), intent(in)    :: name
        integer,          intent(in)    :: var1(:,:)
        integer,          intent(inout) :: var2(:,:)
        character(len=*), optional, intent(in)  :: stat
        character(len=*), optional, intent(in)  :: method
        integer,          optional, intent(in)  :: missing_value
        real(dp),         optional, intent(in)  :: radius
        logical,          optional, intent(out) :: mask2(:,:)
        logical,          optional, intent(in)  :: reset
        logical,          optional, intent(in)  :: mask_pack(:,:)
        character(len=*), optional, intent(in)  :: fill_method
        character(len=*), optional, intent(in)  :: filt_method
        real(dp),         optional, intent(in)  :: filt_par(:)
        logical,          optional, intent(in)  :: verbose

        real(dp), allocatable :: v1(:,:), v2(:,:)
        logical,  allocatable :: m2(:,:)
        real(dp) :: miss
        integer  :: nx, ny

        nx = size(var2,1); ny = size(var2,2)
        allocate(v1(size(var1,1),size(var1,2)), v2(nx,ny), m2(nx,ny))
        v1 = real(var1, dp)
        v2 = real(var2, dp)
        miss = mv_dp
        if (present(missing_value)) miss = real(missing_value, dp)

        call map_field_grid_grid(map, name, v1, v2, stat=stat, method=method, missing_value=miss, radius=radius, &
                                 mask2=m2, reset=reset, mask_pack=mask_pack, fill_method=fill_method, &
                                 filt_method=filt_method, filt_par=filt_par, verbose=verbose)

        var2 = nint(v2)
        if (present(mask2)) mask2 = m2
        deallocate(v1, v2, m2)
    end subroutine map_field_grid_grid_int

    subroutine map_postprocess(var2, miss, maskp2d, name, fill_method, filt_method, filt_par, verbose)
        ! Shared fill + filter pass applied to a 2D target field after mapping
        ! (lifted from the legacy map_scrip_field so map_field matches its surface).
        real(dp),         intent(inout) :: var2(:,:)
        real(dp),         intent(in)    :: miss
        logical,          intent(in)    :: maskp2d(:,:)
        character(len=*), intent(in)    :: name
        character(len=*), optional, intent(in) :: fill_method
        character(len=*), optional, intent(in) :: filt_method
        real(dp),         optional, intent(in) :: filt_par(:)
        logical,          optional, intent(in) :: verbose

        logical,  allocatable :: mask2d(:,:)
        logical  :: verb
        integer  :: npts_apply
        real(dp) :: mean2, mean2b

        verb = .false.
        if (present(verbose)) verb = verbose

        allocate(mask2d(size(var2,1),size(var2,2)))
        mask2d = maskp2d

        ! === Filling ===
        if (present(fill_method)) then
            select case(trim(fill_method))
                case("weighted")
                    call fill_weighted(var2, miss, n=6, mask=mask2d)
                case("nn")
                    call fill_nearest(var2, miss, mask=mask2d)
                case("none")
                    ! Pass - no filling applied
                case default
                    write(*,*) "map_field:: Error: fill method not recognized: "//trim(fill_method)
                    write(*,*) "  fill_method = [weighted,nn,none]."
                    stop
            end select
        end if

        ! === Filtering ===
        if (present(filt_method)) then
            mask2d = (maskp2d .and. var2 /= miss)
            npts_apply = count(maskp2d)
            if (verb .and. npts_apply > 0) mean2 = sum(var2, mask=mask2d) / real(npts_apply,dp)

            select case(trim(filt_method))
                case("gaussian")
                    call filter_gaussian(var2, sigma=filt_par(1), dx=filt_par(2), mask=mask2d)
                case("gaussian-fast")
                    call filter_gaussian_fast(var2, sigma=filt_par(1), dx=filt_par(2), mask=mask2d)
                case("poisson")
                    call filter_poisson(var2, mask=mask2d, tol=filt_par(1), &
                                        missing_value=miss, wrapx=.false., verbose=.false.)
                case("none")
                    ! Pass - no filtering applied
                case default
                    write(*,*) "map_field:: Error: filtering method not recognized: "//trim(filt_method)
                    write(*,*) "  filt_method = [gaussian,gaussian-fast,poisson,none]."
                    stop
            end select

            if (verb .and. npts_apply > 0) then
                mean2b = sum(var2, mask=mask2d) / real(npts_apply,dp)
                write(*,"(4a,2g14.5)") trim(name), " - ", trim(filt_method), &
                                       ": mean[orig,filtered]: ", mean2, mean2b
            end if
        end if

        deallocate(mask2d)
    end subroutine map_postprocess

    subroutine map_field_grid_points(map, name, var1, var2, stat, method, missing_value, radius, mask2, &
                                     reset, mask_pack)
        ! source on a grid (2D), target a point set (1D)
        type(map_class),  intent(in)    :: map
        character(len=*), intent(in)    :: name
        real(dp),         intent(in)    :: var1(:,:)
        real(dp),         intent(inout) :: var2(:)
        character(len=*), optional, intent(in)  :: stat
        character(len=*), optional, intent(in)  :: method
        real(dp),         optional, intent(in)  :: missing_value
        real(dp),         optional, intent(in)  :: radius
        logical,          optional, intent(out) :: mask2(:)
        logical,          optional, intent(in)  :: reset
        logical,          optional, intent(in)  :: mask_pack(:)
        real(dp), allocatable :: v1(:)
        allocate(v1(size(var1)))
        v1 = reshape(var1, [size(var1)])
        call map_field_to_points(map, name, v1, var2, stat, method, missing_value, radius, mask2, &
                                 reset, mask_pack)
        deallocate(v1)
    end subroutine map_field_grid_points

    subroutine map_field_points_grid(map, name, var1, var2, stat, method, missing_value, radius, mask2, &
                                     reset, mask_pack, fill_method, filt_method, filt_par, verbose)
        ! source a point set (1D), target on a grid (2D)
        type(map_class),  intent(in)    :: map
        character(len=*), intent(in)    :: name
        real(dp),         intent(in)    :: var1(:)
        real(dp),         intent(inout) :: var2(:,:)
        character(len=*), optional, intent(in)  :: stat
        character(len=*), optional, intent(in)  :: method
        real(dp),         optional, intent(in)  :: missing_value
        real(dp),         optional, intent(in)  :: radius
        logical,          optional, intent(out) :: mask2(:,:)
        logical,          optional, intent(in)  :: reset
        logical,          optional, intent(in)  :: mask_pack(:,:)
        character(len=*), optional, intent(in)  :: fill_method
        character(len=*), optional, intent(in)  :: filt_method
        real(dp),         optional, intent(in)  :: filt_par(:)
        logical,          optional, intent(in)  :: verbose
        call map_field_to_grid(map, name, var1, var2, stat, method, missing_value, radius, mask2, &
                               reset, mask_pack, fill_method, filt_method, filt_par, verbose)
    end subroutine map_field_points_grid

    ! ----- sp / integer wrappers for the points/grid combos (accumulate in dp) -

    subroutine map_field_points_grid_sp(map, name, var1, var2, stat, method, missing_value, radius, mask2, &
                                        reset, mask_pack, fill_method, filt_method, filt_par, verbose)
        ! source a point set (1D), target on a grid (2D); sp store.
        type(map_class),  intent(in)    :: map
        character(len=*), intent(in)    :: name
        real(sp),         intent(in)    :: var1(:)
        real(sp),         intent(inout) :: var2(:,:)
        character(len=*), optional, intent(in)  :: stat
        character(len=*), optional, intent(in)  :: method
        real(sp),         optional, intent(in)  :: missing_value
        real(dp),         optional, intent(in)  :: radius
        logical,          optional, intent(out) :: mask2(:,:)
        logical,          optional, intent(in)  :: reset
        logical,          optional, intent(in)  :: mask_pack(:,:)
        character(len=*), optional, intent(in)  :: fill_method
        character(len=*), optional, intent(in)  :: filt_method
        real(sp),         optional, intent(in)  :: filt_par(:)
        logical,          optional, intent(in)  :: verbose

        real(dp), allocatable :: v1(:), v2(:,:), fpar(:)
        logical,  allocatable :: m2(:,:)
        real(dp) :: miss
        integer  :: nx, ny

        nx = size(var2,1); ny = size(var2,2)
        allocate(v1(size(var1)), v2(nx,ny), m2(nx,ny))
        v1 = real(var1, dp)
        v2 = real(var2, dp)
        miss = mv_dp
        if (present(missing_value)) miss = real(missing_value, dp)
        if (present(filt_par)) then
            allocate(fpar(size(filt_par)))
            fpar = real(filt_par, dp)
        end if

        call map_field_to_grid(map, name, v1, v2, stat=stat, method=method, missing_value=miss, radius=radius, &
                               mask2=m2, reset=reset, mask_pack=mask_pack, fill_method=fill_method, &
                               filt_method=filt_method, filt_par=fpar, verbose=verbose)

        var2 = real(v2, sp)
        if (present(mask2)) mask2 = m2
        deallocate(v1, v2, m2)
        if (allocated(fpar)) deallocate(fpar)
    end subroutine map_field_points_grid_sp

    subroutine map_field_points_grid_int(map, name, var1, var2, stat, method, missing_value, radius, mask2, &
                                         reset, mask_pack, fill_method, filt_method, filt_par, verbose)
        ! source a point set (1D), target on a grid (2D); integer store.
        type(map_class),  intent(in)    :: map
        character(len=*), intent(in)    :: name
        integer,          intent(in)    :: var1(:)
        integer,          intent(inout) :: var2(:,:)
        character(len=*), optional, intent(in)  :: stat
        character(len=*), optional, intent(in)  :: method
        integer,          optional, intent(in)  :: missing_value
        real(dp),         optional, intent(in)  :: radius
        logical,          optional, intent(out) :: mask2(:,:)
        logical,          optional, intent(in)  :: reset
        logical,          optional, intent(in)  :: mask_pack(:,:)
        character(len=*), optional, intent(in)  :: fill_method
        character(len=*), optional, intent(in)  :: filt_method
        real(dp),         optional, intent(in)  :: filt_par(:)
        logical,          optional, intent(in)  :: verbose

        real(dp), allocatable :: v1(:), v2(:,:)
        logical,  allocatable :: m2(:,:)
        real(dp) :: miss
        integer  :: nx, ny

        nx = size(var2,1); ny = size(var2,2)
        allocate(v1(size(var1)), v2(nx,ny), m2(nx,ny))
        v1 = real(var1, dp)
        v2 = real(var2, dp)
        miss = mv_dp
        if (present(missing_value)) miss = real(missing_value, dp)

        call map_field_to_grid(map, name, v1, v2, stat=stat, method=method, missing_value=miss, radius=radius, &
                               mask2=m2, reset=reset, mask_pack=mask_pack, fill_method=fill_method, &
                               filt_method=filt_method, filt_par=filt_par, verbose=verbose)

        var2 = nint(v2)
        if (present(mask2)) mask2 = m2
        deallocate(v1, v2, m2)
    end subroutine map_field_points_grid_int

    subroutine map_field_grid_points_sp(map, name, var1, var2, stat, method, missing_value, radius, mask2, &
                                        reset, mask_pack)
        ! source on a grid (2D), target a point set (1D); sp store.
        type(map_class),  intent(in)    :: map
        character(len=*), intent(in)    :: name
        real(sp),         intent(in)    :: var1(:,:)
        real(sp),         intent(inout) :: var2(:)
        character(len=*), optional, intent(in)  :: stat
        character(len=*), optional, intent(in)  :: method
        real(sp),         optional, intent(in)  :: missing_value
        real(dp),         optional, intent(in)  :: radius
        logical,          optional, intent(out) :: mask2(:)
        logical,          optional, intent(in)  :: reset
        logical,          optional, intent(in)  :: mask_pack(:)

        real(dp), allocatable :: v1(:), v2(:)
        logical,  allocatable :: m2(:)
        real(dp) :: miss

        allocate(v1(size(var1)), v2(size(var2)), m2(size(var2)))
        v1 = reshape(real(var1, dp), [size(var1)])
        v2 = real(var2, dp)
        miss = mv_dp
        if (present(missing_value)) miss = real(missing_value, dp)

        call map_field_to_points(map, name, v1, v2, stat=stat, method=method, missing_value=miss, radius=radius, &
                                 mask2=m2, reset=reset, mask_pack=mask_pack)

        var2 = real(v2, sp)
        if (present(mask2)) mask2 = m2
        deallocate(v1, v2, m2)
    end subroutine map_field_grid_points_sp

    subroutine map_field_grid_points_int(map, name, var1, var2, stat, method, missing_value, radius, mask2, &
                                         reset, mask_pack)
        ! source on a grid (2D), target a point set (1D); integer store.
        type(map_class),  intent(in)    :: map
        character(len=*), intent(in)    :: name
        integer,          intent(in)    :: var1(:,:)
        integer,          intent(inout) :: var2(:)
        character(len=*), optional, intent(in)  :: stat
        character(len=*), optional, intent(in)  :: method
        integer,          optional, intent(in)  :: missing_value
        real(dp),         optional, intent(in)  :: radius
        logical,          optional, intent(out) :: mask2(:)
        logical,          optional, intent(in)  :: reset
        logical,          optional, intent(in)  :: mask_pack(:)

        real(dp), allocatable :: v1(:), v2(:)
        logical,  allocatable :: m2(:)
        real(dp) :: miss

        allocate(v1(size(var1)), v2(size(var2)), m2(size(var2)))
        v1 = reshape(real(var1, dp), [size(var1)])
        v2 = real(var2, dp)
        miss = mv_dp
        if (present(missing_value)) miss = real(missing_value, dp)

        call map_field_to_points(map, name, v1, v2, stat=stat, method=method, missing_value=miss, radius=radius, &
                                 mask2=m2, reset=reset, mask_pack=mask_pack)

        var2 = nint(v2)
        if (present(mask2)) mask2 = m2
        deallocate(v1, v2, m2)
    end subroutine map_field_grid_points_int

    subroutine map_field_points_points_sp(map, name, var1, var2, stat, method, missing_value, radius, mask2, &
                                          reset, mask_pack)
        ! source a point set (1D), target a point set (1D); sp store.
        type(map_class),  intent(in)    :: map
        character(len=*), intent(in)    :: name
        real(sp),         intent(in)    :: var1(:)
        real(sp),         intent(inout) :: var2(:)
        character(len=*), optional, intent(in)  :: stat
        character(len=*), optional, intent(in)  :: method
        real(sp),         optional, intent(in)  :: missing_value
        real(dp),         optional, intent(in)  :: radius
        logical,          optional, intent(out) :: mask2(:)
        logical,          optional, intent(in)  :: reset
        logical,          optional, intent(in)  :: mask_pack(:)

        real(dp), allocatable :: v1(:), v2(:)
        logical,  allocatable :: m2(:)
        real(dp) :: miss

        allocate(v1(size(var1)), v2(size(var2)), m2(size(var2)))
        v1 = real(var1, dp)
        v2 = real(var2, dp)
        miss = mv_dp
        if (present(missing_value)) miss = real(missing_value, dp)

        call map_field_to_points(map, name, v1, v2, stat=stat, method=method, missing_value=miss, radius=radius, &
                                 mask2=m2, reset=reset, mask_pack=mask_pack)

        var2 = real(v2, sp)
        if (present(mask2)) mask2 = m2
        deallocate(v1, v2, m2)
    end subroutine map_field_points_points_sp

    subroutine map_field_points_points_int(map, name, var1, var2, stat, method, missing_value, radius, mask2, &
                                           reset, mask_pack)
        ! source a point set (1D), target a point set (1D); integer store.
        type(map_class),  intent(in)    :: map
        character(len=*), intent(in)    :: name
        integer,          intent(in)    :: var1(:)
        integer,          intent(inout) :: var2(:)
        character(len=*), optional, intent(in)  :: stat
        character(len=*), optional, intent(in)  :: method
        integer,          optional, intent(in)  :: missing_value
        real(dp),         optional, intent(in)  :: radius
        logical,          optional, intent(out) :: mask2(:)
        logical,          optional, intent(in)  :: reset
        logical,          optional, intent(in)  :: mask_pack(:)

        real(dp), allocatable :: v1(:), v2(:)
        logical,  allocatable :: m2(:)
        real(dp) :: miss

        allocate(v1(size(var1)), v2(size(var2)), m2(size(var2)))
        v1 = real(var1, dp)
        v2 = real(var2, dp)
        miss = mv_dp
        if (present(missing_value)) miss = real(missing_value, dp)

        call map_field_to_points(map, name, v1, v2, stat=stat, method=method, missing_value=miss, radius=radius, &
                                 mask2=m2, reset=reset, mask_pack=mask_pack)

        var2 = nint(v2)
        if (present(mask2)) mask2 = m2
        deallocate(v1, v2, m2)
    end subroutine map_field_points_points_int

    ! ----- apply dispatch (bilinear handled here; rest delegate to weight_map) -

    subroutine map_apply_vec(map, var1, var2, kernel, stat, missing_value, radius, mask2)
        ! Apply the map. `kernel` is the (already resolved) interpolation kernel
        ! and `stat` the aggregation. A baked MAP_WEIGHT map applies its stored
        ! weights (kernel only validated, not re-applied); a MAP_DISTANCE map
        ! builds weights from `kernel` over its neighbor store.
        type(map_class),  intent(in)  :: map
        real(dp),         intent(in)  :: var1(:)
        real(dp),         intent(out) :: var2(:)
        character(len=*), intent(in)  :: kernel
        character(len=*), intent(in)  :: stat
        real(dp), optional, intent(in)  :: missing_value
        real(dp), optional, intent(in)  :: radius
        logical,  optional, intent(out) :: mask2(:)

        if (map%wm%kind == MAP_WEIGHT) then
            ! Weights are baked (conservative / cdo): the kernel cannot change.
            call weight_map_apply(map%wm, var1, var2, stat=stat, &
                                  missing_value=missing_value, radius=radius, mask2=mask2)
        else
            ! Distance store: kernel selects how per-link weights are formed.
            if (trim(kernel) == "bilinear" .or. trim(kernel) == "bilin") then
                call bilinear_apply(map, var1, var2, missing_value, mask2)
            else
                call weight_map_apply(map%wm, var1, var2, kernel=kernel, stat=stat, &
                                      missing_value=missing_value, radius=radius, mask2=mask2)
            end if
        end if
    end subroutine map_apply_vec

    function map_resolve_kernel(map, method) result(kernel)
        ! Resolve the effective interpolation kernel for an apply call: the
        ! optional `method` override if present, otherwise the map's init-time
        ! kernel (map%method). Invalid combinations stop with a clear message
        ! rather than silently doing the wrong thing.
        type(map_class),  intent(in) :: map
        character(len=*), optional, intent(in) :: method
        character(len=32) :: kernel

        kernel = trim(map%method)
        if (.not. present(method)) return
        if (len_trim(method) == 0)  return

        if (map%wm%kind == MAP_WEIGHT) then
            ! Baked weights: an override that disagrees with how they were
            ! generated is a usage error.
            if (trim(method) /= trim(map%method)) then
                write(*,*) "mapping:: map_field: cannot override the kernel of a baked weight map."
                write(*,*) "  this map was generated with method='"//trim(map%method)//"'"
                write(*,*) "  (conservative/cdo weights are fixed); requested method='"//trim(method)//"'."
                write(*,*) "  Pass stat=[mean,count,stdev] to change the aggregation, or re-init the map."
                stop
            end if
            kernel = trim(map%method)
        else
            ! Distance store: only the distance kernels are available.
            select case (trim(method))
                case ("nn", "nearest", "shepard", "radius", "quadrant", "bilinear", "bilin")
                    kernel = trim(method)
                case default
                    write(*,*) "mapping:: map_field: kernel '"//trim(method)//"' is not available on a distance map."
                    write(*,*) "  valid kernels: nn, shepard, quadrant, bilinear."
                    stop
            end select
        end if
    end function map_resolve_kernel

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

        miss = mv_dp
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
