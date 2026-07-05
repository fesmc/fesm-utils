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

    use precision,       only: dp
    use constants, only: mv_dp

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
    end type

    public :: weight_map_t
    public :: weight_map_alloc, weight_map_free, weight_map_index, weight_map_apply

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
        wm%n_src = 0; wm%n_dst = 0; wm%n_links = 0
    end subroutine weight_map_free

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

    subroutine weight_map_apply(wm, var1, var2, kernel, stat, missing_value, radius, mask2)
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
        type(weight_map_t), intent(in)  :: wm
        real(dp),           intent(in)  :: var1(:)
        real(dp),           intent(out) :: var2(:)
        character(len=*), optional, intent(in)  :: kernel
        character(len=*), optional, intent(in)  :: stat
        real(dp), optional, intent(in)  :: missing_value
        real(dp), optional, intent(in)  :: radius
        logical,  optional, intent(out) :: mask2(:)

        character(len=32) :: knl, sta
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
