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
    use fesm_log, only : fesm_error, log_msg, LOG_INFO, to_str

    implicit none
    private

    real(dp), parameter :: MERGE_LOG_TOL = 1.0e-9_dp  ! coverage above 1 + this is logged by coupler_merge_end (round-off is not)
    integer, parameter :: SLOTS_INIT = 8   ! initial slot-array size (doubles on growth)

    type grid_slot
        character(len=256) :: name = ""
        type(grid_class), allocatable :: grid
        real(dp), allocatable :: area(:,:)  ! model cell area [m^2] (coupler_set_area)
        logical  :: from_disk = .false.   ! resolved from grid_<name>.txt, not registered
        real(dp) :: area_med  = -1.0_dp   ! median model cell area [m^2], once computed
    end type

    type coupler_merge_t
        ! Accumulator of coupler_merge_begin / _add / _end: per added source
        ! domain d, its coverage F_d and its amount (or frac*mean) A_d on the
        ! target grid. Opaque to callers.
        character(len=256) :: dst  = ""
        character(len=16)  :: kind = "extensive"
        integer :: nd = 0
        character(len=256), allocatable :: names(:)
        real(dp), allocatable :: cov(:,:,:), amt(:,:,:)   ! (nx, ny, domain)
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

    public :: coupler_class, coupler_merge_t
    public :: coupler_init, coupler_add_grid, coupler_take_grid, coupler_prime, coupler_method
    public :: coupler_set_area, coupler_integral
    public :: coupler_merge_begin, coupler_merge_add, coupler_merge_end
    public :: coupler_split, coupler_owner_reduce
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
        call default_area(cpl%grids(ig))
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
        call default_area(cpl%grids(ig))
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

    function coupler_method(cpl, src, dst, method, norm) result(resolved)
        ! The map method remap(cpl, ..., src, ..., dst, method=method, norm=norm)
        ! would use, without building or applying a map: `method` as given
        ! (default "con"), "auto" resolved as in remap (see resolve_method).
        ! Lets a host log the resolved method. An unregistered grid is read
        ! from grid_<name>.txt as in remap. remap with src == dst copies and
        ! uses no map whatever this returns.
        type(coupler_class), intent(inout) :: cpl
        character(len=*),    intent(in)    :: src, dst
        character(len=*), intent(in), optional :: method, norm
        character(len=32) :: resolved

        resolved = resolve_method(cpl, src, dst, method, norm)
    end function coupler_method

    subroutine coupler_set_area(cpl, name, area)
        ! Set the model cell area [m^2] of grid `name` (registered, or resolved
        ! from grid_<name>.txt): what the model on that grid integrates, e.g.
        ! the Gaussian quadrature R^2*dlon*w_j for a spectral model. remap with
        ! norm="dst"/"sum", coupler_integral and the merge use it.
        type(coupler_class), intent(inout) :: cpl
        character(len=*),    intent(in)    :: name
        real(dp),            intent(in)    :: area(:,:)
        integer :: ig

        ig = resolve_grid(cpl, name)
        if (size(area,1) /= cpl%grids(ig)%grid%G%nx .or. size(area,2) /= cpl%grids(ig)%grid%G%ny) &
            call fesm_error("coupler_set_area", "area of grid '"//trim(name)//"' has the wrong shape.", &
                "shape(area) = "//to_str(shape(area))//new_line("a")// &
                "grid shape  = "//to_str([cpl%grids(ig)%grid%G%nx, cpl%grids(ig)%grid%G%ny]))
        cpl%grids(ig)%area = area
        cpl%grids(ig)%area_med = -1.0_dp
    end subroutine coupler_set_area

    function coupler_integral(cpl, name, var, mask, missing_value, describe) result(total)
        ! Integral sum(area*var) over grid `name` with its model cell areas
        ! (Neumaier-compensated). Cells where mask is .false. or var equals
        ! missing_value are skipped. With describe=.true. the integral is
        ! logged together with the one over the ellipsoidal cell areas
        ! (grid%area*xy_conv^2), as a diagnostic of the metric difference.
        type(coupler_class), intent(in) :: cpl
        character(len=*),    intent(in) :: name
        real(dp),            intent(in) :: var(:,:)
        logical,  optional,  intent(in) :: mask(:,:)
        real(dp), optional,  intent(in) :: missing_value
        logical,  optional,  intent(in) :: describe
        real(dp) :: total
        logical, allocatable :: m(:,:)
        real(dp) :: tot_ell
        integer  :: ig

        ig = find_grid(cpl, name)
        if (ig < 0) call fesm_error("coupler_integral", "grid '"//trim(name)//"' is not registered.")
        if (size(var,1) /= size(cpl%grids(ig)%area,1) .or. size(var,2) /= size(cpl%grids(ig)%area,2)) &
            call fesm_error("coupler_integral", "var has the wrong shape for grid '"//trim(name)//"'.", &
                "shape(var) = "//to_str(shape(var))//new_line("a")// &
                "grid shape = "//to_str(shape(cpl%grids(ig)%area)))
        allocate(m(size(var,1), size(var,2)))
        m = .true.
        if (present(mask)) m = mask
        if (present(missing_value)) m = m .and. var /= missing_value
        total = neumaier_sum(cpl%grids(ig)%area, var, m)

        if (present(describe)) then
            if (describe) then
                associate (g => cpl%grids(ig)%grid)
                    tot_ell = neumaier_sum(g%area*g%cs%xy_conv**2, var, m)
                end associate
                call log_msg(LOG_INFO, "coupler_integral", trim(name)//": model-area integral = "// &
                    to_str(total)//", ellipsoidal-area integral = "//to_str(tot_ell)// &
                    " (relative difference "//to_str((tot_ell - total)/max(abs(total), tiny(1.0_dp)))//")")
            end if
        end if
    end function coupler_integral

    function neumaier_sum(a, v, m) result(s)
        ! sum(a*v, mask=m) with Neumaier compensation.
        real(dp), intent(in) :: a(:,:), v(:,:)
        logical,  intent(in) :: m(:,:)
        real(dp) :: s, c, t, x
        integer  :: i, j
        s = 0.0_dp; c = 0.0_dp
        do j = 1, size(v,2)
            do i = 1, size(v,1)
                if (.not. m(i,j)) cycle
                x = a(i,j)*v(i,j)
                t = s + x
                if (abs(s) >= abs(x)) then
                    c = c + ((s - t) + x)
                else
                    c = c + ((x - t) + s)
                end if
                s = t
            end do
        end do
        s = s + c
    end function neumaier_sum

    ! ----- merge: N regional fields + background -> global -----------------------

    subroutine coupler_merge_begin(acc, cpl, dst, kind)
        ! Start a merge onto grid `dst`. kind = "extensive" (default; amounts per
        ! area, conserved) or "intensive" (coverage-weighted means).
        type(coupler_merge_t), intent(out)   :: acc
        type(coupler_class),   intent(inout) :: cpl
        character(len=*),      intent(in)    :: dst
        character(len=*), optional, intent(in) :: kind
        integer :: ig

        ig = resolve_grid(cpl, dst)
        acc%dst = trim(dst)
        acc%kind = "extensive"
        if (present(kind)) acc%kind = trim(kind)
        if (trim(acc%kind) /= "extensive" .and. trim(acc%kind) /= "intensive") &
            call fesm_error("coupler_merge_begin", "unknown kind '"//trim(acc%kind)//"'.", &
                            "expected 'extensive' or 'intensive'.")
        acc%nd = 0
        allocate(acc%names(0))
        allocate(acc%cov(cpl%grids(ig)%grid%G%nx, cpl%grids(ig)%grid%G%ny, 0))
        allocate(acc%amt(cpl%grids(ig)%grid%G%nx, cpl%grids(ig)%grid%G%ny, 0))
    end subroutine coupler_merge_begin

    subroutine coupler_merge_add(acc, cpl, src, var_src, owner, missing_value)
        ! Add source domain `src` (its field var_src, owner weights in [0,1],
        ! default 1) to the merge. Per target cell k, with the src -> dst "con"
        ! map and valid source cells j (var_src /= missing_value):
        !   F_k = r_k sum_j w_jk o_j / wfull_k                       coverage
        !   A_k = sum_j phi_jk o_j area_j v_j / area_k     (extensive)
        !   A_k = r_k sum_j w_jk o_j v_j / wfull_k         (intensive)
        type(coupler_merge_t), intent(inout) :: acc
        type(coupler_class),   intent(inout) :: cpl
        character(len=*),      intent(in)    :: src
        real(dp),              intent(in)    :: var_src(:,:)
        real(dp), optional,    intent(in)    :: owner(:,:)
        real(dp), optional,    intent(in)    :: missing_value
        integer  :: im, is, id, nx, ny
        real(dp) :: miss
        real(dp), allocatable :: o(:,:), f(:,:), cov(:,:), amt(:,:), tmp3(:,:,:)
        character(len=256), allocatable :: tnames(:)

        if (acc%dst == "") call fesm_error("coupler_merge_add", "call coupler_merge_begin first.")
        if (any(acc%names == src)) &
            call fesm_error("coupler_merge_add", "domain '"//trim(src)//"' was already added.")
        miss = mv_dp
        if (present(missing_value)) miss = missing_value

        im = get_map(cpl, src, acc%dst, "con")
        is = find_grid(cpl, src); id = find_grid(cpl, acc%dst)
        nx = size(acc%cov,1); ny = size(acc%cov,2)

        allocate(o(size(var_src,1), size(var_src,2)), f(size(var_src,1), size(var_src,2)))
        o = 1.0_dp
        if (present(owner)) o = owner
        ! coverage: the owner weight where the value is valid
        f = merge(o, miss, var_src /= miss)
        allocate(cov(nx,ny), amt(nx,ny))
        call map_field(cpl%maps(im)%map, "merge", f, cov, missing_value=miss, norm="dst")
        where (cov == miss) cov = 0.0_dp
        ! amount (extensive, model metric) or frac*mean (intensive, map metric)
        f = merge(o*var_src, miss, var_src /= miss)
        if (trim(acc%kind) == "extensive") then
            call map_field(cpl%maps(im)%map, "merge", f, amt, missing_value=miss, norm="dst", &
                           area1=cpl%grids(is)%area, area2=cpl%grids(id)%area)
        else
            call map_field(cpl%maps(im)%map, "merge", f, amt, missing_value=miss, norm="dst")
        end if
        where (amt == miss) amt = 0.0_dp

        ! append the domain
        acc%nd = acc%nd + 1
        allocate(tnames(acc%nd)); tnames(1:acc%nd-1) = acc%names; tnames(acc%nd) = trim(src)
        call move_alloc(tnames, acc%names)
        allocate(tmp3(nx, ny, acc%nd)); tmp3(:,:,1:acc%nd-1) = acc%cov; tmp3(:,:,acc%nd) = cov
        call move_alloc(tmp3, acc%cov)
        allocate(tmp3(nx, ny, acc%nd)); tmp3(:,:,1:acc%nd-1) = acc%amt; tmp3(:,:,acc%nd) = amt
        call move_alloc(tmp3, acc%amt)
    end subroutine coupler_merge_add

    subroutine coupler_merge_end(acc, cpl, background, var_dst, frac_tot, overlap_tol)
        ! Combine the added domains with the background on the target grid:
        !   extensive: var = (1 - min(F,1))*bg + sum_d A_d
        !   intensive: var = (1 - min(F,1))*bg + sum_d A_d / max(F,1)
        ! with F = sum_d F_d (returned unclipped in frac_tot). Domains are summed
        ! in the order of their names, so the result does not depend on the
        ! order of coupler_merge_add calls. max(F) > 1 + overlap_tol (default
        ! 1e-2) stops, naming the overlapping domains (they need owner weights);
        ! a smaller excess is accepted, and logged with the double-counted amount
        ! when it exceeds 1e-9 (round-off in the sum of the coverages is not logged).
        type(coupler_merge_t),  intent(inout) :: acc
        type(coupler_class),    intent(inout) :: cpl
        real(dp),               intent(in)    :: background(:,:)
        real(dp), allocatable,  intent(inout) :: var_dst(:,:)
        real(dp), optional,     intent(out)   :: frac_tot(:,:)
        real(dp), optional,     intent(in)    :: overlap_tol
        real(dp), allocatable :: ftot(:,:), atot(:,:)
        integer,  allocatable :: ord(:)
        real(dp) :: tol, dbl
        integer  :: nx, ny, d, id, i, j
        character(len=:), allocatable :: who

        if (acc%dst == "") call fesm_error("coupler_merge_end", "call coupler_merge_begin first.")
        nx = size(acc%cov,1); ny = size(acc%cov,2)
        if (size(background,1) /= nx .or. size(background,2) /= ny) &
            call fesm_error("coupler_merge_end", "background has the wrong shape for grid '"// &
                trim(acc%dst)//"'.", "shape(background) = "//to_str(shape(background))// &
                new_line("a")//"grid shape        = "//to_str([nx, ny]))
        tol = 1.0e-2_dp
        if (present(overlap_tol)) tol = overlap_tol

        ord = name_order(acc%names)
        allocate(ftot(nx,ny), atot(nx,ny))
        ftot = 0.0_dp; atot = 0.0_dp
        do d = 1, acc%nd
            ftot = ftot + acc%cov(:,:,ord(d))
            atot = atot + acc%amt(:,:,ord(d))
        end do

        if (maxval(ftot) > 1.0_dp + tol) then
            who = ""
            do d = 1, acc%nd
                if (any(acc%cov(:,:,ord(d)) > 0.0_dp .and. ftot > 1.0_dp + tol)) &
                    who = who//" "//trim(acc%names(ord(d)))
            end do
            call fesm_error("coupler_merge_end", "domains overlap on grid '"//trim(acc%dst)// &
                "' (max coverage "//to_str(maxval(ftot))//" > 1 + overlap_tol): give them owner weights.", &
                "overlapping domains:"//who)
        end if
        if (maxval(ftot) > 1.0_dp + MERGE_LOG_TOL) then
            id = find_grid(cpl, acc%dst)
            dbl = 0.0_dp
            do j = 1, ny
                do i = 1, nx
                    if (ftot(i,j) > 1.0_dp) dbl = dbl + cpl%grids(id)%area(i,j)*atot(i,j)*(1.0_dp - 1.0_dp/ftot(i,j))
                end do
            end do
            call log_msg(LOG_INFO, "coupler_merge_end", trim(acc%dst)//": coverage up to "// &
                to_str(maxval(ftot))//" (within overlap_tol); double-counted amount ~ "//to_str(dbl))
        end if

        if (allocated(var_dst)) then
            if (size(var_dst,1) /= nx .or. size(var_dst,2) /= ny) deallocate(var_dst)
        end if
        if (.not. allocated(var_dst)) allocate(var_dst(nx,ny))
        if (trim(acc%kind) == "extensive") then
            var_dst = (1.0_dp - min(ftot, 1.0_dp))*background + atot
        else
            var_dst = (1.0_dp - min(ftot, 1.0_dp))*background + atot/max(ftot, 1.0_dp)
        end if
        if (present(frac_tot)) frac_tot = ftot
    end subroutine coupler_merge_end

    function name_order(names) result(ord)
        ! Indices that sort `names` (insertion sort; N is small).
        character(len=*), intent(in) :: names(:)
        integer, allocatable :: ord(:)
        integer :: i, j, t
        ord = [(i, i = 1, size(names))]
        do i = 2, size(names)
            t = ord(i); j = i - 1
            do while (j >= 1)
                if (llt(names(t), names(ord(j)))) then
                    ord(j+1) = ord(j); j = j - 1
                else
                    exit
                end if
            end do
            ord(j+1) = t
        end do
    end function name_order

    ! ----- split: global -> one regional grid -------------------------------------

    subroutine coupler_split(cpl, src, var_src, dst, var_dst, method, ref_src, ref_dst, &
                             missing_value, fill_method, mask2)
        ! Global field var_src on `src` onto the regional grid `dst`. Without
        ! references: remap (default method "bilin"). With ref_src and ref_dst:
        ! var_dst = ref_dst + remap(var_src - ref_src), where a missing anomaly
        ! counts as 0 (the reference is kept).
        type(coupler_class),   intent(inout) :: cpl
        character(len=*),      intent(in)    :: src, dst
        real(dp),              intent(in)    :: var_src(:,:)
        real(dp), allocatable, intent(inout) :: var_dst(:,:)
        character(len=*), optional, intent(in)  :: method
        real(dp),         optional, intent(in)  :: ref_src(:,:), ref_dst(:,:)
        real(dp),         optional, intent(in)  :: missing_value
        character(len=*), optional, intent(in)  :: fill_method
        logical,          optional, intent(out) :: mask2(:,:)
        character(len=32) :: mtd
        real(dp) :: miss
        real(dp), allocatable :: anom(:,:)

        mtd = "bilin"
        if (present(method)) mtd = trim(method)
        miss = mv_dp
        if (present(missing_value)) miss = missing_value

        if (present(ref_src) .neqv. present(ref_dst)) &
            call fesm_error("coupler_split", "give both ref_src and ref_dst, or neither.")

        if (.not. present(ref_src)) then
            call remap(cpl, var_src, src, var_dst, dst, method=mtd, missing_value=miss, &
                       fill_method=fill_method, mask2=mask2)
            return
        end if

        anom = merge(var_src - ref_src, miss, var_src /= miss .and. ref_src /= miss)
        call remap(cpl, anom, src, var_dst, dst, method=mtd, missing_value=miss, &
                   fill_method=fill_method, mask2=mask2)
        if (size(ref_dst,1) /= size(var_dst,1) .or. size(ref_dst,2) /= size(var_dst,2)) &
            call fesm_error("coupler_split", "ref_dst has the wrong shape for grid '"//trim(dst)//"'.")
        where (var_dst == miss)
            var_dst = ref_dst
        elsewhere
            var_dst = ref_dst + var_dst
        end where
    end subroutine coupler_split

    subroutine coupler_owner_reduce(cpl, src_hi, owner_hi, dst, owner)
        ! Priority-ordered ownership: remove from `owner` (on grid dst) the
        ! part of each cell already owned by the higher-priority domain src_hi:
        !   owner := max(0, owner - (covered fraction of owner_hi on dst))
        ! via the cached src_hi -> dst "con" map (norm="dst" in the map metric).
        ! Call once at init, in priority order.
        type(coupler_class), intent(inout) :: cpl
        character(len=*),    intent(in)    :: src_hi, dst
        real(dp),            intent(in)    :: owner_hi(:,:)
        real(dp),            intent(inout) :: owner(:,:)
        integer :: im
        real(dp), allocatable :: f(:,:)

        im = get_map(cpl, src_hi, dst, "con")
        if (size(owner,1) /= cpl%maps(im)%map%G%nx .or. size(owner,2) /= cpl%maps(im)%map%G%ny) &
            call fesm_error("coupler_owner_reduce", "owner has the wrong shape for grid '"//trim(dst)//"'.")
        allocate(f(size(owner,1), size(owner,2)))
        call map_field(cpl%maps(im)%map, "owner", owner_hi, f, norm="dst")
        where (f /= mv_dp) owner = max(0.0_dp, owner - f)
    end subroutine coupler_owner_reduce

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
                call move_alloc(cpl%grids(i)%area, tmp(i)%area)
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
        call default_area(cpl%grids(ig))
    end function resolve_grid

    subroutine default_area(slot)
        ! Default model cell area [m^2]: the plane dx*dy for a projected grid
        ! (what an ice-sheet model integrates), grid%area*xy_conv^2 otherwise
        ! (docs/dev/coupler-global-design.md, 2.1). coupler_set_area overrides it.
        type(grid_slot), intent(inout) :: slot
        associate (g => slot%grid)
            if (allocated(slot%area)) deallocate(slot%area)
            allocate(slot%area(g%G%nx, g%G%ny))
            if (g%cs%is_projection) then
                slot%area = abs(g%G%dx*g%G%dy)*g%cs%xy_conv**2
            else
                slot%area = g%area*g%cs%xy_conv**2
            end if
        end associate
        slot%area_med = -1.0_dp
    end subroutine default_area

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
        ! Median model cell area [m^2] of grid slot ig (computed once, reset by
        ! coupler_set_area).
        type(coupler_class), intent(inout) :: cpl
        integer,             intent(in)    :: ig
        real(dp) :: amed
        real(dp), allocatable :: a(:)

        if (cpl%grids(ig)%area_med < 0.0_dp) then
            a = reshape(abs(cpl%grids(ig)%area), [size(cpl%grids(ig)%area)])
            cpl%grids(ig)%area_med = median(a)
        end if
        amed = cpl%grids(ig)%area_med
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
                           reset, mask_pack, fill_method, filt_method, filt_par, norm, frac2)
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
        character(len=*), intent(in),  optional :: norm
        real(dp),         intent(out), optional :: frac2(:,:)
        integer :: im
        real(dp), allocatable :: scl(:,:)

        if (trim(src) == trim(dst)) then
            call identity_scale(cpl, src, norm, scl)
            call size_dst(size(var_src,1), size(var_src,2), 1, reset, var_dst_2D=var_dst)
            call identity_2D(var_src, var_dst, stat, missing_value, mask2, reset, mask_pack, &
                             fill_method, filt_method, filt_par, frac2, scl)
            return
        end if
        im = get_map(cpl, src, dst, resolve_method(cpl, src, dst, method, norm))
        call size_dst(cpl%maps(im)%map%G%nx, cpl%maps(im)%map%G%ny, 1, reset, var_dst_2D=var_dst)
        if (amount_norm(norm)) then
            ! norm="dst"/"sum": the registered model areas of both grids
            call map_field(cpl%maps(im)%map, "coupler", var_src, var_dst, stat=stat, &
                           missing_value=missing_value, mask2=mask2, reset=reset, mask_pack=mask_pack, &
                           fill_method=fill_method, filt_method=filt_method, filt_par=filt_par, &
                           norm=norm, frac2=frac2, area1=cpl%grids(find_grid(cpl, src))%area, &
                           area2=cpl%grids(find_grid(cpl, dst))%area)
        else
            call map_field(cpl%maps(im)%map, "coupler", var_src, var_dst, stat=stat, &
                           missing_value=missing_value, mask2=mask2, reset=reset, mask_pack=mask_pack, &
                           fill_method=fill_method, filt_method=filt_method, filt_par=filt_par, &
                           norm=norm, frac2=frac2)
        end if
    end subroutine remap_2D_dp

    subroutine remap_2D_sp(cpl, var_src, src, var_dst, dst, method, stat, missing_value, mask2, &
                           reset, mask_pack, fill_method, filt_method, filt_par, norm, frac2)
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
        character(len=*), intent(in),  optional :: norm
        real(dp),         intent(out), optional :: frac2(:,:)
        integer :: im
        real(dp), allocatable :: scl(:,:)

        if (trim(src) == trim(dst)) then
            call identity_scale(cpl, src, norm, scl)
            call size_dst(size(var_src,1), size(var_src,2), 1, reset, var_dst_2D_sp=var_dst)
            call identity_2D_sp(var_src, var_dst, stat, missing_value, mask2, reset, mask_pack, &
                                fill_method, filt_method, filt_par, frac2, scl)
            return
        end if
        im = get_map(cpl, src, dst, resolve_method(cpl, src, dst, method, norm))
        call size_dst(cpl%maps(im)%map%G%nx, cpl%maps(im)%map%G%ny, 1, reset, var_dst_2D_sp=var_dst)
        if (amount_norm(norm)) then
            ! norm="dst"/"sum": the registered model areas of both grids
            call map_field(cpl%maps(im)%map, "coupler", var_src, var_dst, stat=stat, &
                           missing_value=missing_value, mask2=mask2, reset=reset, mask_pack=mask_pack, &
                           fill_method=fill_method, filt_method=filt_method, filt_par=filt_par, &
                           norm=norm, frac2=frac2, area1=cpl%grids(find_grid(cpl, src))%area, &
                           area2=cpl%grids(find_grid(cpl, dst))%area)
        else
            call map_field(cpl%maps(im)%map, "coupler", var_src, var_dst, stat=stat, &
                           missing_value=missing_value, mask2=mask2, reset=reset, mask_pack=mask_pack, &
                           fill_method=fill_method, filt_method=filt_method, filt_par=filt_par, &
                           norm=norm, frac2=frac2)
        end if
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
    !
    ! Each layer goes through the 2D specific (one map lookup per layer, after
    ! the first a cache hit); mask2 / frac2 are gathered per layer through 2D
    ! work arrays (unallocated when absent, so they are passed on as absent).

    subroutine remap_3D_dp(cpl, var_src, src, var_dst, dst, method, stat, missing_value, mask2, &
                           reset, mask_pack, fill_method, filt_method, filt_par, norm, frac2)
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
        character(len=*), intent(in),  optional :: norm
        real(dp),         intent(out), optional :: frac2(:,:,:)
        integer :: nx, ny, nz, k
        real(dp), allocatable :: v2(:,:), f2(:,:)
        logical,  allocatable :: m2(:,:)

        call target_shape(cpl, src, dst, size(var_src,1), size(var_src,2), method, norm, nx, ny)
        nz = size(var_src, 3)
        call size_dst(nx, ny, nz, reset, var_dst_3D=var_dst)
        allocate(v2(nx,ny))
        if (present(mask2)) allocate(m2(nx,ny))
        if (present(frac2)) allocate(f2(nx,ny))
        do k = 1, nz
            v2 = var_dst(:,:,k)
            call remap_2D_dp(cpl, var_src(:,:,k), src, v2, dst, method, stat, missing_value, m2, &
                             reset, mask_pack, fill_method, filt_method, filt_par, norm, f2)
            var_dst(:,:,k) = v2
            if (present(mask2)) mask2(:,:,k) = m2
            if (present(frac2)) frac2(:,:,k) = f2
        end do
    end subroutine remap_3D_dp

    subroutine remap_3D_sp(cpl, var_src, src, var_dst, dst, method, stat, missing_value, mask2, &
                           reset, mask_pack, fill_method, filt_method, filt_par, norm, frac2)
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
        character(len=*), intent(in),  optional :: norm
        real(dp),         intent(out), optional :: frac2(:,:,:)
        integer :: nx, ny, nz, k
        real(sp), allocatable :: v2(:,:)
        real(dp), allocatable :: f2(:,:)
        logical,  allocatable :: m2(:,:)

        call target_shape(cpl, src, dst, size(var_src,1), size(var_src,2), method, norm, nx, ny)
        nz = size(var_src, 3)
        call size_dst(nx, ny, nz, reset, var_dst_3D_sp=var_dst)
        allocate(v2(nx,ny))
        if (present(mask2)) allocate(m2(nx,ny))
        if (present(frac2)) allocate(f2(nx,ny))
        do k = 1, nz
            v2 = var_dst(:,:,k)
            call remap_2D_sp(cpl, var_src(:,:,k), src, v2, dst, method, stat, missing_value, m2, &
                             reset, mask_pack, fill_method, filt_method, filt_par, norm, f2)
            var_dst(:,:,k) = v2
            if (present(mask2)) mask2(:,:,k) = m2
            if (present(frac2)) frac2(:,:,k) = f2
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
        integer :: nx, ny, nz, k
        integer, allocatable :: v2(:,:)
        logical, allocatable :: m2(:,:)

        call target_shape(cpl, src, dst, size(var_src,1), size(var_src,2), method, nx=nx, ny=ny)
        nz = size(var_src, 3)
        call size_dst(nx, ny, nz, reset, var_dst_3D_int=var_dst)
        allocate(v2(nx,ny))
        if (present(mask2)) allocate(m2(nx,ny))
        do k = 1, nz
            v2 = var_dst(:,:,k)
            call remap_2D_int(cpl, var_src(:,:,k), src, v2, dst, method, stat, missing_value, m2, &
                              reset, mask_pack)
            var_dst(:,:,k) = v2
            if (present(mask2)) mask2(:,:,k) = m2
        end do
    end subroutine remap_3D_int

    subroutine target_shape(cpl, src, dst, nxs, nys, method, norm, nx, ny)
        ! Target shape of a remap: the source shape (nxs, nys) for src == dst,
        ! else the target grid of the (src, dst, method) map (built if needed).
        type(coupler_class), intent(inout) :: cpl
        character(len=*),    intent(in)    :: src, dst
        integer,             intent(in)    :: nxs, nys
        character(len=*), intent(in), optional :: method, norm
        integer,             intent(out)   :: nx, ny
        integer :: im
        if (trim(src) == trim(dst)) then
            nx = nxs; ny = nys
        else
            im = get_map(cpl, src, dst, resolve_method(cpl, src, dst, method, norm))
            nx = cpl%maps(im)%map%G%nx; ny = cpl%maps(im)%map%G%ny
        end if
    end subroutine target_shape

    ! ----- identity shortcut (src == dst) --------------------------------------

    subroutine identity_2D(var_src, var_dst, stat, missing_value, mask2, reset, mask_pack, &
                           fill_method, filt_method, filt_par, frac2, scale)
        ! remap of a grid onto itself: copy the valid values into var_dst
        ! (already sized). A self-map has one link per target, so mean and
        ! count are the value and stdev is 0. mask_pack, reset and the fill /
        ! filter post-processing apply as in map_field; frac2 = merge(1,0,mask2).
        ! `scale` multiplies the copied values (the cell area for norm="sum").
        real(dp),         intent(in)    :: var_src(:,:)
        real(dp),         intent(inout) :: var_dst(:,:)
        character(len=*), intent(in),  optional :: stat
        real(dp),         intent(in),  optional :: missing_value
        logical,          intent(out), optional :: mask2(:,:)
        logical,          intent(in),  optional :: reset
        logical,          intent(in),  optional :: mask_pack(:,:)
        character(len=*), intent(in),  optional :: fill_method, filt_method
        real(dp),         intent(in),  optional :: filt_par(:)
        real(dp),         intent(out), optional :: frac2(:,:)
        real(dp),         intent(in),  optional :: scale(:,:)

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
        if (present(scale)) then
            where (m) var_dst = var_dst*scale
        end if
        if (present(stat)) then
            if (trim(stat) == "stdev") then
                where (m) var_dst = 0.0_dp
            end if
        end if
        if (present(mask2)) mask2 = m
        if (present(frac2)) frac2 = merge(1.0_dp, 0.0_dp, m)

        if (present(fill_method) .or. present(filt_method)) &
            call map_postprocess(var_dst, miss, maskp, "coupler", fill_method=fill_method, &
                                 filt_method=filt_method, filt_par=filt_par)
    end subroutine identity_2D

    subroutine identity_2D_sp(var_src, var_dst, stat, missing_value, mask2, reset, mask_pack, &
                              fill_method, filt_method, filt_par, frac2, scale)
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
        real(dp),         intent(out), optional :: frac2(:,:)
        real(dp),         intent(in),  optional :: scale(:,:)

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
        call identity_2D(v1, v2, stat, miss, mask2, reset, mask_pack, fill_method, filt_method, fpar, &
                         frac2, scale)
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

    logical function amount_norm(norm)
        ! norm present and "dst" / "sum" (the coupler then supplies model areas)
        character(len=*), intent(in), optional :: norm
        amount_norm = .false.
        if (present(norm)) amount_norm = trim(norm) == "dst" .or. trim(norm) == "sum"
    end function amount_norm

    subroutine identity_scale(cpl, name, norm, scale)
        ! The identity remap's value scale: the model cell area for norm="sum"
        ! (amount per cell); left unallocated (passed on as absent) otherwise.
        type(coupler_class), intent(in) :: cpl
        character(len=*),    intent(in) :: name
        character(len=*), intent(in), optional :: norm
        real(dp), allocatable, intent(out) :: scale(:,:)
        integer :: ig
        if (.not. present(norm)) return
        if (trim(norm) /= "sum") return
        ig = find_grid(cpl, name)
        if (ig < 0) call fesm_error("remap", "norm='sum' needs grid '"//trim(name)//"' registered (its area).")
        scale = cpl%grids(ig)%area
    end subroutine identity_scale

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
