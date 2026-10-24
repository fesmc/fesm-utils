module weight_map
    ! The single sparse source->target weight store for the coords library.
    !
    ! A weight_map_t is a link-centric (SCRIP-style) sparse matrix: parallel
    ! arrays src(:)/dst(:) addressing source and target points, plus a payload
    ! that depends on `kind`:
    !
    !   MAP_DISTANCE : links carry the neighbor distance (+ position, quadrant);
    !                  the weight is computed at apply time from the chosen
    !                  method, so one cached map serves nn/quadrant/shepard/
    !                  bilinear and adapts to a time-varying source mask.
    !   MAP_WEIGHT   : links carry a final weight w(:); apply just reduces it
    !                  (conservative mean/count/stdev, or an imported SCRIP map).
    !
    ! Links are stored grouped by destination (ascending). dst_off(:) is a CSR
    ! offset array giving, for each target k, the link range dst_off(k):dst_off(k+1)-1.
    !
    ! Coverage (conservative maps from the coords generator only; see
    ! docs/dev/coupler-global-design.md, 1): wfull_src/wfull_dst hold the full
    ! source/target cell areas and wcov_src/wcov_dst the link-weight sums before
    ! pruning, all in the map metric (the unit of w). q_src = wcov_src/s_src and
    ! r_dst = wcov_dst/s_dst (s: post-prune link sums) restore the coverage that
    ! weight_map_prune removed; weight_map_coverage derives them. With them,
    ! weight_map_apply offers norm="dst"/"sum" and the covered fraction frac2.

    use precision,       only: dp
    use constants, only: mv_dp
    use fesm_log,  only: fesm_error

    implicit none
    private

    integer, parameter, public :: MAP_DISTANCE = 1
    integer, parameter, public :: MAP_WEIGHT   = 2

    type weight_map_t
        integer :: kind    = MAP_DISTANCE
        integer :: n_src   = 0
        integer :: n_dst   = 0
        integer :: n_links = 0

        integer,  allocatable :: src(:), dst(:)     ! (n_links) link addresses
        real(dp), allocatable :: w(:)               ! (n_links) final weight (MAP_WEIGHT)

        ! MAP_DISTANCE payload
        real(dp), allocatable :: dist(:)            ! (n_links) neighbor distance
        real(dp), allocatable :: xn(:), yn(:)       ! (n_links) neighbor position (for bilinear)
        integer,  allocatable :: quadrant(:)        ! (n_links) 1..4

        ! CSR index by destination
        integer,  allocatable :: dst_off(:)         ! (n_dst+1)

        ! Coverage payload (conservative maps; map metric)
        real(dp), allocatable :: wfull_src(:)       ! (n_src) full source-cell area
        real(dp), allocatable :: wfull_dst(:)       ! (n_dst) full target-cell area
        real(dp), allocatable :: wcov_src(:)        ! (n_src) sum_k w_jk before pruning
        real(dp), allocatable :: wcov_dst(:)        ! (n_dst) sum_j w_jk before pruning
        real(dp), allocatable :: q_src(:)           ! (n_src) wcov_src / post-prune sum
        real(dp), allocatable :: r_dst(:)           ! (n_dst) wcov_dst / post-prune sum
    end type

    public :: weight_map_t
    public :: weight_map_alloc, weight_map_free, weight_map_index, weight_map_apply
    public :: weight_map_prune, weight_map_coverage, weight_map_has_coverage

contains

    subroutine weight_map_alloc(wm, kind, n_src, n_dst, n_links)
        type(weight_map_t), intent(inout) :: wm
        integer,            intent(in)    :: kind, n_src, n_dst, n_links
        call weight_map_free(wm)
        wm%kind    = kind
        wm%n_src   = n_src
        wm%n_dst   = n_dst
        wm%n_links = n_links
        allocate(wm%src(n_links), wm%dst(n_links))
        if (kind == MAP_WEIGHT) then
            allocate(wm%w(n_links))
        else
            allocate(wm%dist(n_links), wm%xn(n_links), wm%yn(n_links), wm%quadrant(n_links))
        end if
        allocate(wm%dst_off(n_dst+1))
        wm%dst_off = 1
    end subroutine weight_map_alloc

    subroutine weight_map_free(wm)
        type(weight_map_t), intent(inout) :: wm
        if (allocated(wm%src))      deallocate(wm%src)
        if (allocated(wm%dst))      deallocate(wm%dst)
        if (allocated(wm%w))        deallocate(wm%w)
        if (allocated(wm%dist))     deallocate(wm%dist)
        if (allocated(wm%xn))       deallocate(wm%xn)
        if (allocated(wm%yn))       deallocate(wm%yn)
        if (allocated(wm%quadrant)) deallocate(wm%quadrant)
        if (allocated(wm%dst_off))  deallocate(wm%dst_off)
        if (allocated(wm%wfull_src)) deallocate(wm%wfull_src)
        if (allocated(wm%wfull_dst)) deallocate(wm%wfull_dst)
        if (allocated(wm%wcov_src))  deallocate(wm%wcov_src)
        if (allocated(wm%wcov_dst))  deallocate(wm%wcov_dst)
        if (allocated(wm%q_src))     deallocate(wm%q_src)
        if (allocated(wm%r_dst))     deallocate(wm%r_dst)
        wm%n_src = 0; wm%n_dst = 0; wm%n_links = 0
    end subroutine weight_map_free

    logical function weight_map_has_coverage(wm)
        ! Does the map carry the coverage payload (wfull/wcov, q/r)?
        type(weight_map_t), intent(in) :: wm
        weight_map_has_coverage = allocated(wm%wfull_src) .and. allocated(wm%wfull_dst) .and. &
                                  allocated(wm%q_src) .and. allocated(wm%r_dst)
    end function weight_map_has_coverage

    subroutine weight_map_coverage(wm, set_cov)
        ! Derive q_src / r_dst from wcov_* and the current link sums. With
        ! set_cov=.true. (the generator, before any pruning), wcov_* are first
        ! set to the current link sums. A no-op for a map without wfull_*.
        type(weight_map_t), intent(inout) :: wm
        logical, optional,  intent(in)    :: set_cov
        real(dp), allocatable :: s_src(:), s_dst(:)
        integer :: k, j

        if (wm%kind /= MAP_WEIGHT) return
        if (.not. (allocated(wm%wfull_src) .and. allocated(wm%wfull_dst))) return

        allocate(s_src(wm%n_src), s_dst(wm%n_dst))
        s_src = 0.0_dp; s_dst = 0.0_dp
        do k = 1, wm%n_dst
            do j = wm%dst_off(k), wm%dst_off(k+1) - 1
                s_dst(k)          = s_dst(k) + wm%w(j)
                s_src(wm%src(j))  = s_src(wm%src(j)) + wm%w(j)
            end do
        end do
        if (present(set_cov)) then
            if (set_cov) then
                wm%wcov_src = s_src
                wm%wcov_dst = s_dst
            end if
        end if
        if (.not. (allocated(wm%wcov_src) .and. allocated(wm%wcov_dst))) return

        if (allocated(wm%q_src)) deallocate(wm%q_src)
        if (allocated(wm%r_dst)) deallocate(wm%r_dst)
        allocate(wm%q_src(wm%n_src), wm%r_dst(wm%n_dst))
        wm%q_src = 1.0_dp
        wm%r_dst = 1.0_dp
        where (s_src > 0.0_dp) wm%q_src = wm%wcov_src/s_src
        where (s_dst > 0.0_dp) wm%r_dst = wm%wcov_dst/s_dst
    end subroutine weight_map_coverage

    subroutine weight_map_index(wm)
        ! Build the CSR offset array dst_off from dst(:). Requires links to be
        ! grouped by destination in ascending order (as map_init produces them).
        type(weight_map_t), intent(inout) :: wm
        integer :: j, k
        if (.not. allocated(wm%dst_off)) allocate(wm%dst_off(wm%n_dst+1))
        wm%dst_off = wm%n_links + 1
        wm%dst_off(1) = 1
        ! For each target k, find first link with dst >= k
        k = 1
        do j = 1, wm%n_links
            do while (k <= wm%dst(j))
                wm%dst_off(k) = j
                k = k + 1
            end do
        end do
        do while (k <= wm%n_dst+1)
            wm%dst_off(k) = wm%n_links + 1
            k = k + 1
        end do
    end subroutine weight_map_index

    subroutine weight_map_prune(wm, cutoff)
        ! Drop the boundary-sliver links a separable/tensor generator produces:
        ! for each target, remove links whose weight is below `cutoff` times
        ! that target's largest link weight, then compact and reindex. apply
        ! already renormalizes per target over the survivors, so dropping
        ! sub-cutoff links changes results negligibly while shrinking the map
        ! (memory, the on-disk cache, and apply work). The largest link of a
        ! target always survives (a threshold relative to the link sum would
        ! drop every link of a target with more than 1/cutoff near-equal links,
        ! e.g. a 1 km source on a 2 deg cell). The coverage removed is kept in
        ! q_src / r_dst (weight_map_coverage). No-op unless it is a MAP_WEIGHT
        ! map with cutoff > 0 (MAP_DISTANCE weights are built at apply time).
        type(weight_map_t), intent(inout) :: wm
        real(dp),           intent(in)    :: cutoff
        integer  :: k, j, j1, j2, nkeep
        real(dp) :: wmax, thresh
        integer,  allocatable :: src2(:), dst2(:), off2(:)
        real(dp), allocatable :: w2(:)

        if (cutoff <= 0.0_dp)         return
        if (wm%kind /= MAP_WEIGHT)    return
        if (wm%n_links == 0)          return

        allocate(src2(wm%n_links), dst2(wm%n_links), w2(wm%n_links))
        allocate(off2(wm%n_dst+1))
        off2(1) = 1
        nkeep = 0
        do k = 1, wm%n_dst
            j1 = wm%dst_off(k); j2 = wm%dst_off(k+1) - 1
            wmax = 0.0_dp
            do j = j1, j2
                wmax = max(wmax, wm%w(j))
            end do
            thresh = cutoff * wmax
            do j = j1, j2
                if (wm%w(j) >= thresh) then
                    nkeep = nkeep + 1
                    src2(nkeep) = wm%src(j)
                    dst2(nkeep) = wm%dst(j)
                    w2(nkeep)   = wm%w(j)
                end if
            end do
            off2(k+1) = nkeep + 1
        end do

        deallocate(wm%src, wm%dst, wm%w)
        allocate(wm%src(nkeep), wm%dst(nkeep), wm%w(nkeep))
        wm%src = src2(1:nkeep)
        wm%dst = dst2(1:nkeep)
        wm%w   = w2(1:nkeep)
        wm%dst_off = off2
        wm%n_links = nkeep
        call weight_map_coverage(wm)
    end subroutine weight_map_prune

    subroutine weight_map_apply(wm, var1, var2, kernel, stat, missing_value, radius, mask2, &
                                norm, frac2, area1, area2)
        ! Apply the map: var2(n_dst) = aggregation of var1(n_src) over each
        ! target's links. Two roles, deliberately separate:
        !   kernel - the interpolation kernel ("nn"/"shepard"/"quadrant"), used
        !            only for MAP_DISTANCE to build per-link weights from the
        !            stored neighbor distances/quadrants. Ignored for MAP_WEIGHT
        !            (those weights are already baked, e.g. conservative/cdo).
        !   stat   - the aggregation/reduction over the weighted links
        !            ("mean" [default]/"count"/"stdev"). Same meaning for both
        !            map kinds, so a given `stat` is always valid.
        ! Currently-missing sources are skipped (and, for MAP_DISTANCE, weights
        ! are recomputed over the survivors).
        !
        ! Normalisation and coverage (MAP_WEIGHT maps with the coverage payload;
        ! docs/dev/coupler-global-design.md, 1.3). V_k: links of k with w > 0
        ! and a valid source value; phi_jk = q_j*w_jk/wfull_src_j.
        !   norm  - "frac" (default; the weighted mean over V_k, today's path)
        !           "dst"  r_k*sum_V w v / wfull_dst_k, or with area1+area2
        !                  sum_V phi_jk*area1_j*v_j / area2_k (exact in that metric)
        !           "sum"  sum_V phi_jk*alpha_j*v_j, alpha = area1 (else 1)
        !           norm /= "frac" needs stat="mean".
        !   frac2 - covered fraction min(1, r_k*sum_V w / wfull_dst_k); 0 where
        !           nothing valid is linked.
        !   area1, area2 - model areas of the source / target cells.
        ! mask2 = (sum_V w > 0) for every norm.
        type(weight_map_t), intent(in)  :: wm
        real(dp),           intent(in)  :: var1(:)
        real(dp),           intent(out) :: var2(:)
        character(len=*), optional, intent(in)  :: kernel
        character(len=*), optional, intent(in)  :: stat
        real(dp), optional, intent(in)  :: missing_value
        real(dp), optional, intent(in)  :: radius
        logical,  optional, intent(out) :: mask2(:)
        character(len=*), optional, intent(in)  :: norm
        real(dp), optional, intent(out) :: frac2(:)
        real(dp), optional, intent(in)  :: area1(:)
        real(dp), optional, intent(in)  :: area2(:)

        character(len=32) :: knl, sta, nrm
        real(dp) :: miss, rad
        integer  :: k, j, j1, j2, nlink, nmax
        real(dp), allocatable :: vloc(:), wloc(:), dloc(:)
        integer,  allocatable :: qloc(:)
        real(dp) :: val, wsum, vsum, w
        logical  :: ok

        knl = "nn"
        if (present(kernel)) knl = trim(kernel)
        sta = "mean"
        if (present(stat)) sta = trim(stat)
        miss = mv_dp
        if (present(missing_value)) miss = missing_value
        rad = huge(1.0_dp)
        if (present(radius)) rad = radius
        nrm = "frac"
        if (present(norm)) nrm = trim(norm)

        if (trim(nrm) /= "frac" .or. present(frac2)) then
            if (wm%kind /= MAP_WEIGHT .or. .not. weight_map_has_coverage(wm)) &
                call fesm_error("weight_map_apply", "norm /= 'frac' and frac2 need a conservative map "// &
                                "with coverage data (coords generator).")
        end if

        if (present(frac2)) call coverage_frac(wm, var1, miss, frac2)

        if (trim(nrm) /= "frac") then
            if (trim(sta) /= "mean") &
                call fesm_error("weight_map_apply", "norm='"//trim(nrm)//"' needs stat='mean' (stat='"// &
                                trim(sta)//"').")
            call apply_norm(wm, var1, var2, miss, nrm, mask2, area1, area2)
            return
        end if

        var2 = miss
        if (present(mask2)) mask2 = .false.

        ! Fast path: baked weights + mean (the per-timestep coupling case, e.g.
        ! conservative/cdo maps applied every model step). Reduce the stored
        ! weights in place with zero per-target heap traffic -- gathering the
        ! source values into a scratch array first buys nothing for a plain
        ! weighted mean and dominates the cost at high target-point counts.
        if (wm%kind == MAP_WEIGHT .and. trim(sta) == "mean") then
            ! Disjoint var2(k)/mask2(k) writes and thread-local accumulators, so
            ! the target loop parallelizes directly with no shared scratch.
            !$omp parallel do default(shared) schedule(guided) &
            !$omp   private(k, j, j1, j2, wsum, vsum, w, val)
            do k = 1, wm%n_dst
                j1 = wm%dst_off(k); j2 = wm%dst_off(k+1) - 1
                if (j2 < j1) cycle
                wsum = 0.0_dp; vsum = 0.0_dp
                do j = j1, j2
                    w = wm%w(j)
                    if (w > 0.0_dp) then
                        val = var1(wm%src(j))
                        if (val /= miss) then
                            wsum = wsum + w
                            vsum = vsum + w*val
                        end if
                    end if
                end do
                if (wsum > 0.0_dp) then
                    var2(k) = vsum/wsum
                    if (present(mask2)) mask2(k) = .true.
                end if
            end do
            !$omp end parallel do
            return
        end if

        ! General path (count/stdev, or distance maps whose weights are built
        ! per apply from the chosen kernel). The per-target scratch is sized to
        ! the largest link block and allocated ONCE, then reused across targets
        ! via (1:nlink) sections -- never allocated inside the target loop.
        nmax = 0
        do k = 1, wm%n_dst
            nmax = max(nmax, wm%dst_off(k+1) - wm%dst_off(k))
        end do
        if (nmax < 1) return

        ! Each thread allocates its own scratch once (private allocatables,
        ! sized to nmax) and reuses it across its share of the target loop;
        ! reduce/distance_weights are pure and var2(k)/mask2(k) writes disjoint.
        !$omp parallel default(shared) &
        !$omp   private(k, j, j1, j2, nlink, vloc, wloc, dloc, qloc, val, ok)
        allocate(vloc(nmax), wloc(nmax))
        if (wm%kind /= MAP_WEIGHT) allocate(dloc(nmax), qloc(nmax))
        !$omp do schedule(guided)
        do k = 1, wm%n_dst
            j1 = wm%dst_off(k); j2 = wm%dst_off(k+1) - 1
            nlink = j2 - j1 + 1
            if (nlink < 1) cycle

            do j = 1, nlink
                vloc(j) = var1(wm%src(j1+j-1))
            end do

            if (wm%kind == MAP_WEIGHT) then
                wloc(1:nlink) = wm%w(j1:j2)
                call reduce(sta, vloc(1:nlink), wloc(1:nlink), miss, val, ok)
            else
                dloc(1:nlink) = wm%dist(j1:j2)
                qloc(1:nlink) = wm%quadrant(j1:j2)
                call distance_weights(knl, vloc(1:nlink), dloc(1:nlink), qloc(1:nlink), &
                                      miss, rad, wloc(1:nlink))
                call reduce(sta, vloc(1:nlink), wloc(1:nlink), miss, val, ok)
            end if

            if (ok) then
                var2(k) = val
                if (present(mask2)) mask2(k) = .true.
            end if
        end do
        !$omp end do
        if (allocated(dloc)) deallocate(dloc, qloc)
        deallocate(vloc, wloc)
        !$omp end parallel
    end subroutine weight_map_apply

    subroutine coverage_frac(wm, var1, miss, frac2)
        ! frac2_k = min(1, r_k*sum_{V_k} w / wfull_dst_k), 0 where nothing valid.
        type(weight_map_t), intent(in)  :: wm
        real(dp),           intent(in)  :: var1(:), miss
        real(dp),           intent(out) :: frac2(:)
        integer  :: k, j
        real(dp) :: wsum

        !$omp parallel do default(shared) schedule(guided) private(k, j, wsum)
        do k = 1, wm%n_dst
            wsum = 0.0_dp
            do j = wm%dst_off(k), wm%dst_off(k+1) - 1
                if (wm%w(j) > 0.0_dp) then
                    if (var1(wm%src(j)) /= miss) wsum = wsum + wm%w(j)
                end if
            end do
            frac2(k) = 0.0_dp
            if (wsum > 0.0_dp .and. wm%wfull_dst(k) > 0.0_dp) &
                frac2(k) = min(1.0_dp, wm%r_dst(k)*wsum/wm%wfull_dst(k))
        end do
        !$omp end parallel do
    end subroutine coverage_frac

    subroutine apply_norm(wm, var1, var2, miss, nrm, mask2, area1, area2)
        ! norm="dst" / "sum" apply (see weight_map_apply).
        type(weight_map_t), intent(in)  :: wm
        real(dp),           intent(in)  :: var1(:), miss
        real(dp),           intent(out) :: var2(:)
        character(len=*),   intent(in)  :: nrm
        logical,  optional, intent(out) :: mask2(:)
        real(dp), optional, intent(in)  :: area1(:), area2(:)
        integer  :: k, j, js
        real(dp) :: wsum, vsum, phi, alpha
        logical  :: with_areas

        select case (trim(nrm))
            case ("dst")
                if (present(area1) .neqv. present(area2)) &
                    call fesm_error("weight_map_apply", "norm='dst' needs both area1 and area2, or neither.")
                with_areas = present(area1)
            case ("sum")
                with_areas = present(area1)
            case default
                call fesm_error("weight_map_apply", "unknown norm '"//trim(nrm)//"'.", &
                                "expected 'frac', 'dst' or 'sum'.")
        end select

        var2 = miss
        if (present(mask2)) mask2 = .false.

        !$omp parallel do default(shared) schedule(guided) private(k, j, js, wsum, vsum, phi, alpha)
        do k = 1, wm%n_dst
            wsum = 0.0_dp; vsum = 0.0_dp
            do j = wm%dst_off(k), wm%dst_off(k+1) - 1
                if (wm%w(j) <= 0.0_dp) cycle
                js = wm%src(j)
                if (var1(js) == miss) cycle
                wsum = wsum + wm%w(j)
                if (trim(nrm) == "dst" .and. .not. with_areas) then
                    vsum = vsum + wm%w(j)*var1(js)
                else
                    phi = 0.0_dp
                    if (wm%wfull_src(js) > 0.0_dp) phi = wm%q_src(js)*wm%w(j)/wm%wfull_src(js)
                    alpha = 1.0_dp
                    if (with_areas) alpha = area1(js)
                    vsum = vsum + phi*alpha*var1(js)
                end if
            end do
            if (wsum > 0.0_dp) then
                if (trim(nrm) == "sum") then
                    var2(k) = vsum
                else if (with_areas) then
                    var2(k) = vsum/area2(k)
                else
                    var2(k) = 0.0_dp
                    if (wm%wfull_dst(k) > 0.0_dp) var2(k) = wm%r_dst(k)*vsum/wm%wfull_dst(k)
                end if
                if (present(mask2)) mask2(k) = .true.
            end if
        end do
        !$omp end parallel do
    end subroutine apply_norm

    subroutine distance_weights(method, v, dist, quad, miss, rad, w)
        ! Compute per-link weights for a MAP_DISTANCE target, over the valid
        ! (non-missing, within-radius) neighbors, for the given method.
        character(len=*), intent(in)  :: method
        real(dp),         intent(in)  :: v(:), dist(:)
        integer,          intent(in)  :: quad(:)
        real(dp),         intent(in)  :: miss, rad
        real(dp),         intent(out) :: w(:)
        integer  :: n, j, jq, jbest
        real(dp) :: dbest
        logical  :: valid

        n = size(v)
        w = 0.0_dp

        select case (trim(method))

            case ("nn", "nearest")
                ! single nearest valid neighbor
                jbest = 0; dbest = huge(1.0_dp)
                do j = 1, n
                    if (v(j) /= miss .and. dist(j) <= rad .and. dist(j) < dbest) then
                        dbest = dist(j); jbest = j
                    end if
                end do
                if (jbest > 0) w(jbest) = 1.0_dp

            case ("radius", "shepard")
                ! inverse-distance over all valid neighbors within radius
                do j = 1, n
                    if (v(j) /= miss .and. dist(j) <= rad) then
                        w(j) = 1.0_dp / max(dist(j), 1.0e-12_dp)**2
                    end if
                end do

            case ("quadrant")
                ! nearest valid neighbor in each of the 4 quadrants, IDW-weighted
                do jq = 1, 4
                    jbest = 0; dbest = huge(1.0_dp)
                    do j = 1, n
                        if (quad(j) == jq .and. v(j) /= miss .and. dist(j) <= rad &
                            .and. dist(j) < dbest) then
                            dbest = dist(j); jbest = j
                        end if
                    end do
                    if (jbest > 0) w(jbest) = 1.0_dp / max(dbest, 1.0e-12_dp)**2
                end do

            case default
                write(*,*) "weight_map:: distance_weights: unknown method '"//trim(method)//"'"
                stop

        end select
    end subroutine distance_weights

    subroutine reduce(method, v, w, miss, val, ok)
        ! Reduce values v with weights w (skipping missing) into a single value.
        character(len=*), intent(in)  :: method
        real(dp),         intent(in)  :: v(:), w(:)
        real(dp),         intent(in)  :: miss
        real(dp),         intent(out) :: val
        logical,          intent(out) :: ok
        integer  :: n, j
        real(dp) :: wsum, vsum, mean, varsum

        n = size(v)
        val = miss; ok = .false.

        select case (trim(method))

            case ("mean")
                wsum = 0.0_dp; vsum = 0.0_dp
                do j = 1, n
                    if (v(j) /= miss .and. w(j) > 0.0_dp) then
                        wsum = wsum + w(j)
                        vsum = vsum + w(j)*v(j)
                    end if
                end do
                if (wsum > 0.0_dp) then
                    val = vsum / wsum; ok = .true.
                end if

            case ("count")
                ! weighted mode (most total weight per distinct value)
                call weighted_mode(v, w, miss, val, ok)

            case ("stdev")
                wsum = 0.0_dp; vsum = 0.0_dp
                do j = 1, n
                    if (v(j) /= miss .and. w(j) > 0.0_dp) then
                        wsum = wsum + w(j); vsum = vsum + w(j)*v(j)
                    end if
                end do
                if (wsum > 0.0_dp) then
                    mean = vsum / wsum
                    varsum = 0.0_dp
                    do j = 1, n
                        if (v(j) /= miss .and. w(j) > 0.0_dp) &
                            varsum = varsum + w(j)*(v(j)-mean)**2
                    end do
                    val = sqrt(varsum / wsum); ok = .true.
                end if

            case default
                write(*,*) "weight_map:: reduce: unknown method '"//trim(method)//"'"
                stop

        end select
    end subroutine reduce

    subroutine weighted_mode(v, w, miss, val, ok)
        real(dp), intent(in)  :: v(:), w(:), miss
        real(dp), intent(out) :: val
        logical,  intent(out) :: ok
        integer  :: n, i, j
        real(dp) :: tot, best
        n = size(v)
        val = miss; ok = .false.; best = -1.0_dp
        do i = 1, n
            if (v(i) == miss .or. w(i) <= 0.0_dp) cycle
            tot = 0.0_dp
            do j = 1, n
                if (v(j) == v(i) .and. w(j) > 0.0_dp) tot = tot + w(j)
            end do
            if (tot > best) then
                best = tot; val = v(i); ok = .true.
            end if
        end do
    end subroutine weighted_mode

end module weight_map
