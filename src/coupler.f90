module coupler
    ! A thin multigrid coupling layer over the coords mapping core.
    !
    ! A coupler_class holds two things:
    !   * a registry of named grids (grid slots, keyed by name), and
    !   * a lazily-populated cache of grid->grid maps (map slots).
    !
    ! Grids are nodes; maps are directed, method-typed edges. remap(cpl, ...)
    ! looks up (or builds+caches on first use) the map for a given
    ! (source, target, method) and applies it, sizing the destination array
    ! itself. Each map is built once, stored once, and shared across all callers.
    !
    ! Storage: both registries are allocatable slot arrays without a fixed
    ! capacity. A slot holds its grid / map as an allocatable component; when
    ! the slot array grows (doubling), each component is moved with move_alloc,
    ! so the large grid and weight arrays are never copied and slot indices stay
    ! valid. Maps are built directly from the registered grids (no copies). A
    ! grid that is not registered is read once from grid_<name>.txt in the map
    ! folder into a "from_disk" slot and reused.
    !
    ! method selects the interpolation kernel and is baked into the map at build
    ! time, so it is part of the cache key: "con" (conservative/area-weighted,
    ! the default), a distance kernel (nn / shepard / quadrant / bilin), or
    ! "auto" (resolved per call, see resolve_method). stat (mean / count /
    ! stdev) is applied when the field is mapped and is NOT part of the key --
    ! it is forwarded to map_field per call, as are missing_value, mask2,
    ! reset, mask_pack, fill_method, filt_method and filt_par.
    !
    ! remap with src == dst copies (no map, no registration needed);
    ! mask_pack, reset, fill_method and filt_* still apply, and mask2 is
    ! var /= missing_value.
    !
    ! gen selects the map generator used when a map's weights must actually be
    ! built: "coords" (in-package, the default) or "cdo" (external cdo call). It
    ! is a build-time-only concern (irrelevant once the map is cached), so it is
    ! NOT part of the cache key. A coupler carries a default gen (set in
    ! coupler_init); a specific map can override it via coupler_prime(gen=),
    ! which is the point at which that map's weights are generated.

    use coords, only : dp, sp, mv_dp, grid_class, map_class, map_init, map_field, &
                       map_postprocess, grid_cdo_read_desc, grid_signature_str, grid_signature_match
    use fesm_log, only : fesm_error, to_str

    implicit none
    private

    integer, parameter :: SLOTS_INIT = 8   ! initial slot-array size (doubles on growth)

    type grid_slot
        character(len=256) :: name = ""
        type(grid_class), allocatable :: grid
        logical  :: from_disk = .false.   ! resolved from grid_<name>.txt, not registered
        real(dp) :: area_med  = -1.0_dp   ! median model cell area [m^2], once computed
    end type

    type map_slot
        character(len=256) :: src = "", dst = ""
        character(len=32)  :: method = ""
        type(map_class), allocatable :: map
    end type

    type coupler_class
        type(grid_slot), allocatable :: grids(:)
        integer            :: ngrids = 0
        type(map_slot),  allocatable :: maps(:)
        integer            :: nmaps  = 0
        character(len=512) :: map_fldr = "maps"   ! disk cache dir for map_init
        character(len=32)  :: gen = "coords"       ! default map generator; override
                                                   ! per map via coupler_prime(gen=)
        character(len=32)  :: on_mismatch = "error" ! cache file of other grids: error | rebuild
        logical            :: hash_names = .false.  ! grid hashes in cache file names
    end type

    interface remap
        ! One wrapper per (kind x rank); each sizes var_dst and calls the generic
        ! map_field, which dispatches on kind (dp / sp / int).
        !
        ! Note: the concrete `dp` and `sp` overloads are intentional and are NOT
        ! tied to the compile-time working precision `wp`. The coupler accepts
        ! both real kinds from callers regardless of how `wp` is set, so these
        ! must stay `real(dp)` / `real(sp)` — do not "genericize" the `_sp`
        ! wrappers to `real(wp)`: if `wp == dp` that collides with the `_dp`
        ! overload (ambiguous generic) and drops support for `sp` input.
        module procedure remap_2D_dp, remap_2D_sp, remap_2D_int
        module procedure remap_3D_dp, remap_3D_sp, remap_3D_int
    end interface

    public :: coupler_class
    public :: coupler_init, coupler_add_grid, coupler_take_grid, coupler_prime
    public :: remap

contains

    subroutine coupler_init(cpl, map_fldr, gen, on_mismatch, hash_names)
        ! Reset a coupler to an empty registry + cache. `gen` sets the coupler-wide
        ! default map generator ("coords" / "cdo"); it can be overridden per map at
        ! prime time via coupler_prime(gen=). `on_mismatch` ("error" default |
        ! "rebuild") and `hash_names` (default .false.) are passed to map_init
        ! for every map (cache file identity, see map_init).
        type(coupler_class), intent(out) :: cpl
        character(len=*), intent(in), optional :: map_fldr
        character(len=*), intent(in), optional :: gen
        character(len=*), intent(in), optional :: on_mismatch
        logical,          intent(in), optional :: hash_names

        allocate(cpl%grids(SLOTS_INIT), cpl%maps(SLOTS_INIT))
        cpl%ngrids = 0
        cpl%nmaps  = 0
        cpl%map_fldr = "maps"
        if (present(map_fldr)) cpl%map_fldr = trim(map_fldr)
        cpl%gen = "coords"
        if (present(gen)) cpl%gen = trim(gen)
        cpl%on_mismatch = "error"
        if (present(on_mismatch)) cpl%on_mismatch = trim(on_mismatch)
        cpl%hash_names = .false.
        if (present(hash_names)) cpl%hash_names = hash_names
    end subroutine coupler_init

    subroutine coupler_add_grid(cpl, name, grid)
        ! Register a grid under `name` (one copy, at init). Grid identity (and
        ! the map cache/disk filenames) flow through grid%name, so the stored
        ! copy is renamed to `name`. A name already resolved from disk is
        ! accepted if the grids agree (no-op); otherwise a repeated name is an
        ! error.
        type(coupler_class), intent(inout) :: cpl
        character(len=*),    intent(in)    :: name
        type(grid_class),    intent(in)    :: grid
        integer :: ig

        ig = find_grid(cpl, name)
        if (ig > 0) then
            call check_disk_slot(cpl, ig, grid, "coupler_add_grid")
            return
        end if
        ig = new_grid_slot(cpl, name)
        allocate(cpl%grids(ig)%grid, source=grid)
        cpl%grids(ig)%grid%name = trim(name)
    end subroutine coupler_add_grid

    subroutine coupler_take_grid(cpl, name, grid)
        ! Register a grid under `name` by moving it into the coupler (no copy;
        ! `grid` is deallocated on return). For grids where a second copy costs
        ! GBs (e.g. 1 km Antarctica). Same name rules as coupler_add_grid.
        type(coupler_class),           intent(inout) :: cpl
        character(len=*),              intent(in)    :: name
        type(grid_class), allocatable, intent(inout) :: grid
        integer :: ig

        if (.not. allocated(grid)) &
            call fesm_error("coupler_take_grid", "grid '"//trim(name)//"' is not allocated.")
        ig = find_grid(cpl, name)
        if (ig > 0) then
            call check_disk_slot(cpl, ig, grid, "coupler_take_grid")
            deallocate(grid)
            return
        end if
        ig = new_grid_slot(cpl, name)
        call move_alloc(grid, cpl%grids(ig)%grid)
        cpl%grids(ig)%grid%name = trim(name)
    end subroutine coupler_take_grid

    subroutine coupler_prime(cpl, src, dst, method, gen)
        ! Eagerly build a map up front (fail fast, pay the build/disk-load cost
        ! at init rather than mid-timestep). Same find-or-build path as remap.
        ! `gen` overrides the coupler-wide default generator for this one map --
        ! prime is where a specific map's weights get generated, so it is the
        ! natural override point. A no-op for src == dst (remap copies).
        ! method="auto" primes the map remap would choose for norm="frac".
        type(coupler_class), intent(inout) :: cpl
        character(len=*),    intent(in)    :: src, dst
        character(len=*), intent(in), optional :: method, gen

        integer :: im
        if (trim(src) == trim(dst)) return
        im = get_map(cpl, src, dst, resolve_method(cpl, src, dst, method), gen)
    end subroutine coupler_prime

    ! ----- registry internals --------------------------------------------------

    function find_grid(cpl, name) result(ig)
        ! Registry index of `name`, or -1 if not registered.
        type(coupler_class), intent(in) :: cpl
        character(len=*),    intent(in) :: name
        integer :: ig, i

        ig = -1
        do i = 1, cpl%ngrids
            if (trim(cpl%grids(i)%name) == trim(name)) then
                ig = i
                return
            end if
        end do
    end function find_grid

    function new_grid_slot(cpl, name) result(ig)
        ! Append an empty grid slot named `name`, doubling the slot array when
        ! full (components moved, never copied).
        type(coupler_class), intent(inout) :: cpl
        character(len=*),    intent(in)    :: name
        integer :: ig
        type(grid_slot), allocatable :: tmp(:)
        integer :: i

        if (.not. allocated(cpl%grids)) allocate(cpl%grids(SLOTS_INIT))
        if (cpl%ngrids >= size(cpl%grids)) then
            allocate(tmp(2*size(cpl%grids)))
            do i = 1, cpl%ngrids
                tmp(i)%name      = cpl%grids(i)%name
                tmp(i)%from_disk = cpl%grids(i)%from_disk
                tmp(i)%area_med  = cpl%grids(i)%area_med
                call move_alloc(cpl%grids(i)%grid, tmp(i)%grid)
            end do
            call move_alloc(tmp, cpl%grids)
        end if
        cpl%ngrids = cpl%ngrids + 1
        ig = cpl%ngrids
        cpl%grids(ig)%name = trim(name)
    end function new_grid_slot

    function new_map_slot(cpl, src, dst, method) result(im)
        ! Append an empty map slot for (src, dst, method), doubling when full.
        type(coupler_class), intent(inout) :: cpl
        character(len=*),    intent(in)    :: src, dst, method
        integer :: im
        type(map_slot), allocatable :: tmp(:)
        integer :: i

        if (.not. allocated(cpl%maps)) allocate(cpl%maps(SLOTS_INIT))
        if (cpl%nmaps >= size(cpl%maps)) then
            allocate(tmp(2*size(cpl%maps)))
            do i = 1, cpl%nmaps
                tmp(i)%src    = cpl%maps(i)%src
                tmp(i)%dst    = cpl%maps(i)%dst
                tmp(i)%method = cpl%maps(i)%method
                call move_alloc(cpl%maps(i)%map, tmp(i)%map)
            end do
            call move_alloc(tmp, cpl%maps)
        end if
        cpl%nmaps = cpl%nmaps + 1
        im = cpl%nmaps
        cpl%maps(im)%src    = trim(src)
        cpl%maps(im)%dst    = trim(dst)
        cpl%maps(im)%method = trim(method)
    end function new_map_slot

    subroutine check_disk_slot(cpl, ig, grid, caller)
        ! A second registration of an existing name: allowed only over a slot
        ! resolved from disk, and only if the grids agree.
        type(coupler_class), intent(in) :: cpl
        integer,             intent(in) :: ig
        type(grid_class),    intent(in) :: grid
        character(len=*),    intent(in) :: caller
        character(len=128) :: why

        if (.not. cpl%grids(ig)%from_disk) &
            call fesm_error(caller, "grid '"//trim(cpl%grids(ig)%name)//"' already registered.")
        if (.not. grid_signature_match(grid_signature_str(cpl%grids(ig)%grid), &
                                       grid_signature_str(grid), why)) &
            call fesm_error(caller, "grid '"//trim(cpl%grids(ig)%name)//"' differs from the grid "// &
                "already read from "//trim(cpl%map_fldr)//"/grid_"//trim(cpl%grids(ig)%name)//".txt ("// &
                trim(why)//").")
    end subroutine check_disk_slot

    function resolve_grid(cpl, name) result(ig)
        ! Registry index of `name`; an unregistered name is read once from
        ! grid_<name>.txt in the map folder into a from_disk slot.
        type(coupler_class), intent(inout) :: cpl
        character(len=*),    intent(in)    :: name
        integer :: ig

        ig = find_grid(cpl, name)
        if (ig > 0) return
        ig = new_grid_slot(cpl, name)
        allocate(cpl%grids(ig)%grid)
        call grid_cdo_read_desc(cpl%grids(ig)%grid, trim(name), trim(cpl%map_fldr))
        cpl%grids(ig)%from_disk = .true.
    end function resolve_grid

    function get_map(cpl, src, dst, method, gen) result(im)
        ! Find the cached (src->dst, method) map, or build+cache it on miss,
        ! directly from the registry's grids. `gen` selects the generator on a
        ! build miss; when absent the coupler-wide default (cpl%gen) is used.
        type(coupler_class), intent(inout) :: cpl
        character(len=*),    intent(in)    :: src, dst, method
        character(len=*), intent(in), optional :: gen
        integer :: im

        integer :: i, is, id
        character(len=32) :: g

        do i = 1, cpl%nmaps
            if (trim(cpl%maps(i)%src)    == trim(src)    .and. &
                trim(cpl%maps(i)%dst)    == trim(dst)    .and. &
                trim(cpl%maps(i)%method) == trim(method)) then
                im = i
                return
            end if
        end do

        is = resolve_grid(cpl, src)
        id = resolve_grid(cpl, dst)

        g = cpl%gen
        if (present(gen)) g = trim(gen)

        im = new_map_slot(cpl, src, dst, method)
        allocate(cpl%maps(im)%map)
        call map_init(cpl%maps(im)%map, cpl%grids(is)%grid, cpl%grids(id)%grid, &
                      method=trim(method), gen=trim(g), fldr=trim(cpl%map_fldr), load=.true., &
                      on_mismatch=trim(cpl%on_mismatch), hash_names=cpl%hash_names)
    end function get_map

    function resolve_method(cpl, src, dst, method, norm) result(mtd)
        ! The map method of a remap call: `method` as given (default "con").
        ! "auto" resolves to "con" for norm /= "frac" (amounts are remapped
        ! conservatively both ways); otherwise to "bilin" when the target is
        ! finer than the source (smaller median model cell area: refining) and
        ! to "con" when it is coarser or equal (coarsening).
        type(coupler_class), intent(inout) :: cpl
        character(len=*),    intent(in)    :: src, dst
        character(len=*), intent(in), optional :: method, norm
        character(len=32) :: mtd
        real(dp) :: a_src, a_dst

        mtd = "con"
        if (present(method)) mtd = trim(method)
        if (trim(mtd) /= "auto") return

        mtd = "con"
        if (present(norm)) then
            if (trim(norm) /= "frac") return
        end if
        a_src = median_area(cpl, resolve_grid(cpl, src))
        a_dst = median_area(cpl, resolve_grid(cpl, dst))
        if (a_dst < a_src) mtd = "bilin"
    end function resolve_method

    function median_area(cpl, ig) result(amed)
        ! Median model cell area [m^2] of grid slot ig, computed once. The
        ! default model area is the plane dx*dy for a projected grid and
        ! grid%area*xy_conv^2 otherwise (docs/dev/coupler-global-design.md 2.1).
        type(coupler_class), intent(inout) :: cpl
        integer,             intent(in)    :: ig
        real(dp) :: amed
        real(dp), allocatable :: a(:)

        if (cpl%grids(ig)%area_med >= 0.0_dp) then
            amed = cpl%grids(ig)%area_med
            return
        end if
        associate (g => cpl%grids(ig)%grid)
            if (g%cs%is_projection) then
                amed = abs(g%G%dx*g%G%dy)*g%cs%xy_conv**2
            else
                allocate(a(size(g%area)))
                a = reshape(abs(g%area), [size(g%area)])*g%cs%xy_conv**2
                amed = median(a)
            end if
        end associate
        cpl%grids(ig)%area_med = amed
    end function median_area

    function median(a) result(m)
        ! Median (upper middle element for even n) of a; a is reordered.
        ! Hoare quickselect, O(n) expected.
        real(dp), intent(inout) :: a(:)
        real(dp) :: m
        integer  :: k, lo, hi, i, j
        real(dp) :: pivot, t

        k = (size(a) + 1)/2
        lo = 1; hi = size(a)
        do while (lo < hi)
            pivot = a((lo + hi)/2)
            i = lo; j = hi
            do while (i <= j)
                do while (a(i) < pivot); i = i + 1; end do
                do while (a(j) > pivot); j = j - 1; end do
                if (i <= j) then
                    t = a(i); a(i) = a(j); a(j) = t
                    i = i + 1; j = j - 1
                end if
            end do
            if (k <= j) then
                hi = j
            else if (k >= i) then
                lo = i
            else
                exit
            end if
        end do
        m = a(k)
    end function median

    ! ----- remap: 2D ---------------------------------------------------------

    subroutine remap_2D_dp(cpl, var_src, src, var_dst, dst, method, stat, missing_value, mask2, &
                           reset, mask_pack, fill_method, filt_method, filt_par)
        type(coupler_class),   intent(inout) :: cpl
        real(dp),              intent(in)    :: var_src(:,:)
        character(len=*),      intent(in)    :: src, dst
        real(dp), allocatable, intent(inout) :: var_dst(:,:)
        character(len=*), intent(in),  optional :: method, stat
        real(dp),         intent(in),  optional :: missing_value
        logical,          intent(out), optional :: mask2(:,:)
        logical,          intent(in),  optional :: reset
        logical,          intent(in),  optional :: mask_pack(:,:)
        character(len=*), intent(in),  optional :: fill_method, filt_method
        real(dp),         intent(in),  optional :: filt_par(:)
        integer :: im

        if (trim(src) == trim(dst)) then
            call size_dst(size(var_src,1), size(var_src,2), 1, reset, var_dst_2D=var_dst)
            call identity_2D(var_src, var_dst, stat, missing_value, mask2, reset, mask_pack, &
                             fill_method, filt_method, filt_par)
            return
        end if
        im = get_map(cpl, src, dst, resolve_method(cpl, src, dst, method))
        call size_dst(cpl%maps(im)%map%G%nx, cpl%maps(im)%map%G%ny, 1, reset, var_dst_2D=var_dst)
        call map_field(cpl%maps(im)%map, "coupler", var_src, var_dst, stat=stat, &
                       missing_value=missing_value, mask2=mask2, reset=reset, mask_pack=mask_pack, &
                       fill_method=fill_method, filt_method=filt_method, filt_par=filt_par)
    end subroutine remap_2D_dp

    subroutine remap_2D_sp(cpl, var_src, src, var_dst, dst, method, stat, missing_value, mask2, &
                           reset, mask_pack, fill_method, filt_method, filt_par)
        type(coupler_class),   intent(inout) :: cpl
        real(sp),              intent(in)    :: var_src(:,:)
        character(len=*),      intent(in)    :: src, dst
        real(sp), allocatable, intent(inout) :: var_dst(:,:)
        character(len=*), intent(in),  optional :: method, stat
        real(sp),         intent(in),  optional :: missing_value
        logical,          intent(out), optional :: mask2(:,:)
        logical,          intent(in),  optional :: reset
        logical,          intent(in),  optional :: mask_pack(:,:)
        character(len=*), intent(in),  optional :: fill_method, filt_method
        real(sp),         intent(in),  optional :: filt_par(:)
        integer :: im

        if (trim(src) == trim(dst)) then
            call size_dst(size(var_src,1), size(var_src,2), 1, reset, var_dst_2D_sp=var_dst)
            call identity_2D_sp(var_src, var_dst, stat, missing_value, mask2, reset, mask_pack, &
                                fill_method, filt_method, filt_par)
            return
        end if
        im = get_map(cpl, src, dst, resolve_method(cpl, src, dst, method))
        call size_dst(cpl%maps(im)%map%G%nx, cpl%maps(im)%map%G%ny, 1, reset, var_dst_2D_sp=var_dst)
        call map_field(cpl%maps(im)%map, "coupler", var_src, var_dst, stat=stat, &
                       missing_value=missing_value, mask2=mask2, reset=reset, mask_pack=mask_pack, &
                       fill_method=fill_method, filt_method=filt_method, filt_par=filt_par)
    end subroutine remap_2D_sp

    subroutine remap_2D_int(cpl, var_src, src, var_dst, dst, method, stat, missing_value, mask2, &
                            reset, mask_pack)
        type(coupler_class),  intent(inout) :: cpl
        integer,              intent(in)    :: var_src(:,:)
        character(len=*),     intent(in)    :: src, dst
        integer, allocatable, intent(inout) :: var_dst(:,:)
        character(len=*), intent(in),  optional :: method, stat
        integer,          intent(in),  optional :: missing_value
        logical,          intent(out), optional :: mask2(:,:)
        logical,          intent(in),  optional :: reset
        logical,          intent(in),  optional :: mask_pack(:,:)
        integer :: im

        if (trim(src) == trim(dst)) then
            call size_dst(size(var_src,1), size(var_src,2), 1, reset, var_dst_2D_int=var_dst)
            call identity_2D_int(var_src, var_dst, missing_value, mask2, reset, mask_pack)
            return
        end if
        im = get_map(cpl, src, dst, resolve_method(cpl, src, dst, method))
        call size_dst(cpl%maps(im)%map%G%nx, cpl%maps(im)%map%G%ny, 1, reset, var_dst_2D_int=var_dst)
        call map_field(cpl%maps(im)%map, "coupler", var_src, var_dst, stat=stat, &
                       missing_value=missing_value, mask2=mask2, reset=reset, mask_pack=mask_pack)
    end subroutine remap_2D_int

    ! ----- remap: 3D (the 2D map applied per layer) ------------------------------

    subroutine remap_3D_dp(cpl, var_src, src, var_dst, dst, method, stat, missing_value, mask2, &
                           reset, mask_pack, fill_method, filt_method, filt_par)
        type(coupler_class),   intent(inout) :: cpl
        real(dp),              intent(in)    :: var_src(:,:,:)
        character(len=*),      intent(in)    :: src, dst
        real(dp), allocatable, intent(inout) :: var_dst(:,:,:)
        character(len=*), intent(in),  optional :: method, stat
        real(dp),         intent(in),  optional :: missing_value
        logical,          intent(out), optional :: mask2(:,:,:)
        logical,          intent(in),  optional :: reset
        logical,          intent(in),  optional :: mask_pack(:,:)
        character(len=*), intent(in),  optional :: fill_method, filt_method
        real(dp),         intent(in),  optional :: filt_par(:)
        integer :: im, nz, k
        logical :: identity

        nz = size(var_src, 3)
        identity = trim(src) == trim(dst)
        if (identity) then
            call size_dst(size(var_src,1), size(var_src,2), nz, reset, var_dst_3D=var_dst)
        else
            im = get_map(cpl, src, dst, resolve_method(cpl, src, dst, method))
            call size_dst(cpl%maps(im)%map%G%nx, cpl%maps(im)%map%G%ny, nz, reset, var_dst_3D=var_dst)
        end if
        do k = 1, nz
            if (identity .and. present(mask2)) then
                call identity_2D(var_src(:,:,k), var_dst(:,:,k), stat, missing_value, mask2(:,:,k), &
                                 reset, mask_pack, fill_method, filt_method, filt_par)
            else if (identity) then
                call identity_2D(var_src(:,:,k), var_dst(:,:,k), stat, missing_value, &
                                 reset=reset, mask_pack=mask_pack, fill_method=fill_method, &
                                 filt_method=filt_method, filt_par=filt_par)
            else if (present(mask2)) then
                call map_field(cpl%maps(im)%map, "coupler", var_src(:,:,k), var_dst(:,:,k), stat=stat, &
                               missing_value=missing_value, mask2=mask2(:,:,k), reset=reset, &
                               mask_pack=mask_pack, fill_method=fill_method, filt_method=filt_method, &
                               filt_par=filt_par)
            else
                call map_field(cpl%maps(im)%map, "coupler", var_src(:,:,k), var_dst(:,:,k), stat=stat, &
                               missing_value=missing_value, reset=reset, mask_pack=mask_pack, &
                               fill_method=fill_method, filt_method=filt_method, filt_par=filt_par)
            end if
        end do
    end subroutine remap_3D_dp

    subroutine remap_3D_sp(cpl, var_src, src, var_dst, dst, method, stat, missing_value, mask2, &
                           reset, mask_pack, fill_method, filt_method, filt_par)
        type(coupler_class),   intent(inout) :: cpl
        real(sp),              intent(in)    :: var_src(:,:,:)
        character(len=*),      intent(in)    :: src, dst
        real(sp), allocatable, intent(inout) :: var_dst(:,:,:)
        character(len=*), intent(in),  optional :: method, stat
        real(sp),         intent(in),  optional :: missing_value
        logical,          intent(out), optional :: mask2(:,:,:)
        logical,          intent(in),  optional :: reset
        logical,          intent(in),  optional :: mask_pack(:,:)
        character(len=*), intent(in),  optional :: fill_method, filt_method
        real(sp),         intent(in),  optional :: filt_par(:)
        integer :: im, nz, k
        logical :: identity

        nz = size(var_src, 3)
        identity = trim(src) == trim(dst)
        if (identity) then
            call size_dst(size(var_src,1), size(var_src,2), nz, reset, var_dst_3D_sp=var_dst)
        else
            im = get_map(cpl, src, dst, resolve_method(cpl, src, dst, method))
            call size_dst(cpl%maps(im)%map%G%nx, cpl%maps(im)%map%G%ny, nz, reset, var_dst_3D_sp=var_dst)
        end if
        do k = 1, nz
            if (identity .and. present(mask2)) then
                call identity_2D_sp(var_src(:,:,k), var_dst(:,:,k), stat, missing_value, mask2(:,:,k), &
                                    reset, mask_pack, fill_method, filt_method, filt_par)
            else if (identity) then
                call identity_2D_sp(var_src(:,:,k), var_dst(:,:,k), stat, missing_value, &
                                    reset=reset, mask_pack=mask_pack, fill_method=fill_method, &
                                    filt_method=filt_method, filt_par=filt_par)
            else if (present(mask2)) then
                call map_field(cpl%maps(im)%map, "coupler", var_src(:,:,k), var_dst(:,:,k), stat=stat, &
                               missing_value=missing_value, mask2=mask2(:,:,k), reset=reset, &
                               mask_pack=mask_pack, fill_method=fill_method, filt_method=filt_method, &
                               filt_par=filt_par)
            else
                call map_field(cpl%maps(im)%map, "coupler", var_src(:,:,k), var_dst(:,:,k), stat=stat, &
                               missing_value=missing_value, reset=reset, mask_pack=mask_pack, &
                               fill_method=fill_method, filt_method=filt_method, filt_par=filt_par)
            end if
        end do
    end subroutine remap_3D_sp

    subroutine remap_3D_int(cpl, var_src, src, var_dst, dst, method, stat, missing_value, mask2, &
                            reset, mask_pack)
        type(coupler_class),  intent(inout) :: cpl
        integer,              intent(in)    :: var_src(:,:,:)
        character(len=*),     intent(in)    :: src, dst
        integer, allocatable, intent(inout) :: var_dst(:,:,:)
        character(len=*), intent(in),  optional :: method, stat
        integer,          intent(in),  optional :: missing_value
        logical,          intent(out), optional :: mask2(:,:,:)
        logical,          intent(in),  optional :: reset
        logical,          intent(in),  optional :: mask_pack(:,:)
        integer :: im, nz, k
        logical :: identity

        nz = size(var_src, 3)
        identity = trim(src) == trim(dst)
        if (identity) then
            call size_dst(size(var_src,1), size(var_src,2), nz, reset, var_dst_3D_int=var_dst)
        else
            im = get_map(cpl, src, dst, resolve_method(cpl, src, dst, method))
            call size_dst(cpl%maps(im)%map%G%nx, cpl%maps(im)%map%G%ny, nz, reset, var_dst_3D_int=var_dst)
        end if
        do k = 1, nz
            if (identity .and. present(mask2)) then
                call identity_2D_int(var_src(:,:,k), var_dst(:,:,k), missing_value, mask2(:,:,k), &
                                     reset, mask_pack)
            else if (identity) then
                call identity_2D_int(var_src(:,:,k), var_dst(:,:,k), missing_value, &
                                     reset=reset, mask_pack=mask_pack)
            else if (present(mask2)) then
                call map_field(cpl%maps(im)%map, "coupler", var_src(:,:,k), var_dst(:,:,k), stat=stat, &
                               missing_value=missing_value, mask2=mask2(:,:,k), reset=reset, &
                               mask_pack=mask_pack)
            else
                call map_field(cpl%maps(im)%map, "coupler", var_src(:,:,k), var_dst(:,:,k), stat=stat, &
                               missing_value=missing_value, reset=reset, mask_pack=mask_pack)
            end if
        end do
    end subroutine remap_3D_int

    ! ----- identity shortcut (src == dst) --------------------------------------

    subroutine identity_2D(var_src, var_dst, stat, missing_value, mask2, reset, mask_pack, &
                           fill_method, filt_method, filt_par)
        ! remap of a grid onto itself: copy the valid values into var_dst
        ! (already sized). A self-map has one link per target, so mean and
        ! count are the value and stdev is 0. mask_pack, reset and the fill /
        ! filter post-processing apply as in map_field.
        real(dp),         intent(in)    :: var_src(:,:)
        real(dp),         intent(inout) :: var_dst(:,:)
        character(len=*), intent(in),  optional :: stat
        real(dp),         intent(in),  optional :: missing_value
        logical,          intent(out), optional :: mask2(:,:)
        logical,          intent(in),  optional :: reset
        logical,          intent(in),  optional :: mask_pack(:,:)
        character(len=*), intent(in),  optional :: fill_method, filt_method
        real(dp),         intent(in),  optional :: filt_par(:)

        real(dp) :: miss
        logical  :: reset_pts
        logical, allocatable :: m(:,:), maskp(:,:)

        miss = mv_dp
        if (present(missing_value)) miss = missing_value
        reset_pts = .true.
        if (present(reset)) reset_pts = reset

        allocate(maskp(size(var_src,1), size(var_src,2)), m(size(var_src,1), size(var_src,2)))
        maskp = .true.
        if (present(mask_pack)) maskp = mask_pack
        m = (var_src /= miss) .and. maskp

        if (reset_pts) then
            where (.not. m) var_dst = miss
        end if
        where (m) var_dst = var_src
        if (present(stat)) then
            if (trim(stat) == "stdev") then
                where (m) var_dst = 0.0_dp
            end if
        end if
        if (present(mask2)) mask2 = m

        if (present(fill_method) .or. present(filt_method)) &
            call map_postprocess(var_dst, miss, maskp, "coupler", fill_method=fill_method, &
                                 filt_method=filt_method, filt_par=filt_par)
    end subroutine identity_2D

    subroutine identity_2D_sp(var_src, var_dst, stat, missing_value, mask2, reset, mask_pack, &
                              fill_method, filt_method, filt_par)
        ! sp identity: in dp, stored back in sp (as map_field's sp wrappers do).
        real(sp),         intent(in)    :: var_src(:,:)
        real(sp),         intent(inout) :: var_dst(:,:)
        character(len=*), intent(in),  optional :: stat
        real(sp),         intent(in),  optional :: missing_value
        logical,          intent(out), optional :: mask2(:,:)
        logical,          intent(in),  optional :: reset
        logical,          intent(in),  optional :: mask_pack(:,:)
        character(len=*), intent(in),  optional :: fill_method, filt_method
        real(sp),         intent(in),  optional :: filt_par(:)

        real(dp), allocatable :: v1(:,:), v2(:,:), miss, fpar(:)

        allocate(v1(size(var_src,1), size(var_src,2)), v2(size(var_src,1), size(var_src,2)))
        v1 = real(var_src, dp)
        v2 = 0.0_dp
        if (present(reset)) then
            if (.not. reset) v2 = real(var_dst, dp)
        end if
        if (present(missing_value)) miss = real(missing_value, dp)
        if (present(filt_par)) fpar = real(filt_par, dp)
        ! unallocated miss / fpar are passed on as absent
        call identity_2D(v1, v2, stat, miss, mask2, reset, mask_pack, fill_method, filt_method, fpar)
        var_dst = real(v2, sp)
    end subroutine identity_2D_sp

    subroutine identity_2D_int(var_src, var_dst, missing_value, mask2, reset, mask_pack)
        ! Integer identity: copy where valid and packed.
        integer, intent(in)    :: var_src(:,:)
        integer, intent(inout) :: var_dst(:,:)
        integer, intent(in),  optional :: missing_value
        logical, intent(out), optional :: mask2(:,:)
        logical, intent(in),  optional :: reset
        logical, intent(in),  optional :: mask_pack(:,:)

        integer :: miss
        logical :: reset_pts
        logical, allocatable :: m(:,:)

        miss = nint(mv_dp)
        if (present(missing_value)) miss = missing_value
        reset_pts = .true.
        if (present(reset)) reset_pts = reset

        allocate(m(size(var_src,1), size(var_src,2)))
        m = (var_src /= miss)
        if (present(mask_pack)) m = m .and. mask_pack
        if (reset_pts) then
            where (.not. m) var_dst = miss
        end if
        where (m) var_dst = var_src
        if (present(mask2)) mask2 = m
    end subroutine identity_2D_int

    ! ----- destination sizing ----------------------------------------------------

    subroutine size_dst(nx, ny, nz, reset, var_dst_2D, var_dst_2D_sp, var_dst_2D_int, &
                        var_dst_3D, var_dst_3D_sp, var_dst_3D_int)
        ! Allocate the (one present) destination at (nx,ny[,nz]), reallocating
        ! on a shape change. With reset=.false. its values are kept where
        ! nothing is mapped, so it must already be allocated at that shape.
        integer, intent(in) :: nx, ny, nz
        logical, intent(in), optional :: reset
        real(dp), allocatable, intent(inout), optional :: var_dst_2D(:,:)
        real(sp), allocatable, intent(inout), optional :: var_dst_2D_sp(:,:)
        integer,  allocatable, intent(inout), optional :: var_dst_2D_int(:,:)
        real(dp), allocatable, intent(inout), optional :: var_dst_3D(:,:,:)
        real(sp), allocatable, intent(inout), optional :: var_dst_3D_sp(:,:,:)
        integer,  allocatable, intent(inout), optional :: var_dst_3D_int(:,:,:)
        logical :: keep, ok
        integer :: want(3), shp(3)

        keep = .false.
        if (present(reset)) keep = .not. reset
        want = [nx, ny, nz]
        shp  = 0
        ok   = .true.

        if (present(var_dst_2D)) then
            if (allocated(var_dst_2D)) shp = [shape(var_dst_2D), 1]
            if (shp(1) /= nx .or. shp(2) /= ny) then
                ok = .false.
                if (allocated(var_dst_2D)) deallocate(var_dst_2D)
                if (.not. keep) allocate(var_dst_2D(nx,ny))
            end if
        else if (present(var_dst_2D_sp)) then
            if (allocated(var_dst_2D_sp)) shp = [shape(var_dst_2D_sp), 1]
            if (shp(1) /= nx .or. shp(2) /= ny) then
                ok = .false.
                if (allocated(var_dst_2D_sp)) deallocate(var_dst_2D_sp)
                if (.not. keep) allocate(var_dst_2D_sp(nx,ny))
            end if
        else if (present(var_dst_2D_int)) then
            if (allocated(var_dst_2D_int)) shp = [shape(var_dst_2D_int), 1]
            if (shp(1) /= nx .or. shp(2) /= ny) then
                ok = .false.
                if (allocated(var_dst_2D_int)) deallocate(var_dst_2D_int)
                if (.not. keep) allocate(var_dst_2D_int(nx,ny))
            end if
        else if (present(var_dst_3D)) then
            if (allocated(var_dst_3D)) shp = shape(var_dst_3D)
            if (any(shp /= want)) then
                ok = .false.
                if (allocated(var_dst_3D)) deallocate(var_dst_3D)
                if (.not. keep) allocate(var_dst_3D(nx,ny,nz))
            end if
        else if (present(var_dst_3D_sp)) then
            if (allocated(var_dst_3D_sp)) shp = shape(var_dst_3D_sp)
            if (any(shp /= want)) then
                ok = .false.
                if (allocated(var_dst_3D_sp)) deallocate(var_dst_3D_sp)
                if (.not. keep) allocate(var_dst_3D_sp(nx,ny,nz))
            end if
        else if (present(var_dst_3D_int)) then
            if (allocated(var_dst_3D_int)) shp = shape(var_dst_3D_int)
            if (any(shp /= want)) then
                ok = .false.
                if (allocated(var_dst_3D_int)) deallocate(var_dst_3D_int)
                if (.not. keep) allocate(var_dst_3D_int(nx,ny,nz))
            end if
        end if

        if (keep .and. .not. ok) &
            call fesm_error("remap", "reset=.false. needs var_dst allocated at the target shape.", &
                            "var_dst shape = "//to_str(shp)//new_line("a")//"target shape  = "//to_str(want))
    end subroutine size_dst

end module coupler
