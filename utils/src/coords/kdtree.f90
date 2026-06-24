module kdtree
    ! Array-backed balanced k-d tree for fast nearest-neighbor and radius
    ! queries. Self-contained, no external dependencies.
    !
    ! Points are passed as pts(ndim, n). For geographic neighbor search the
    ! caller embeds lon/lat on the unit sphere in 3D (ndim=3); for a shared
    ! projected/Cartesian system the projected x,y are used directly (ndim=2).
    ! Chord distance in the embedding is monotonic in the true metric, so the
    ! k nearest in the tree are the k nearest in reality.
    !
    ! Build:  median split on the widest dimension (quickselect), O(n log n).
    ! Query:  recursive descent with branch pruning.

    use precision, only: dp

    implicit none
    private

    type kdtree_t
        integer :: n    = 0
        integer :: ndim = 0
        real(dp), allocatable :: pts(:,:)      ! (ndim, n) copy of the points
        integer,  allocatable :: idx(:)        ! permutation of 1..n (tree layout)
        integer,  allocatable :: splitdim(:)   ! split dimension at each tree node
    end type

    public :: kdtree_t
    public :: kdtree_build, kdtree_knn, kdtree_radius, kdtree_free

contains

    subroutine kdtree_build(tree, pts)
        type(kdtree_t), intent(out) :: tree
        real(dp),       intent(in)  :: pts(:,:)   ! (ndim, n)
        integer :: i

        tree%ndim = size(pts,1)
        tree%n    = size(pts,2)
        allocate(tree%pts(tree%ndim, tree%n))
        tree%pts = pts
        allocate(tree%idx(tree%n))
        allocate(tree%splitdim(tree%n))
        do i = 1, tree%n
            tree%idx(i) = i
        end do
        if (tree%n > 0) call build_range(tree, 1, tree%n)
    end subroutine kdtree_build

    subroutine kdtree_free(tree)
        type(kdtree_t), intent(inout) :: tree
        if (allocated(tree%pts))      deallocate(tree%pts)
        if (allocated(tree%idx))      deallocate(tree%idx)
        if (allocated(tree%splitdim)) deallocate(tree%splitdim)
        tree%n = 0; tree%ndim = 0
    end subroutine kdtree_free

    recursive subroutine build_range(tree, lo, hi)
        type(kdtree_t), intent(inout) :: tree
        integer,        intent(in)    :: lo, hi
        integer :: mid, d
        if (lo > hi) return
        d   = widest_dim(tree, lo, hi)
        mid = (lo + hi) / 2
        call select_kth(tree, lo, hi, mid, d)   ! partition so idx(mid) is the median along d
        tree%splitdim(mid) = d
        call build_range(tree, lo,    mid-1)
        call build_range(tree, mid+1, hi)
    end subroutine build_range

    integer function widest_dim(tree, lo, hi) result(d)
        type(kdtree_t), intent(in) :: tree
        integer,        intent(in) :: lo, hi
        real(dp) :: lo_v(tree%ndim), hi_v(tree%ndim), v
        integer  :: i, j
        lo_v =  huge(1.0_dp)
        hi_v = -huge(1.0_dp)
        do i = lo, hi
            do j = 1, tree%ndim
                v = tree%pts(j, tree%idx(i))
                if (v < lo_v(j)) lo_v(j) = v
                if (v > hi_v(j)) hi_v(j) = v
            end do
        end do
        d = 1
        do j = 2, tree%ndim
            if (hi_v(j) - lo_v(j) > hi_v(d) - lo_v(d)) d = j
        end do
    end function widest_dim

    subroutine select_kth(tree, lo, hi, k, d)
        ! Hoare quickselect: rearrange tree%idx(lo:hi) so that position k holds
        ! the element that would be at k if sorted by coordinate d, with all
        ! smaller coordinates to its left and larger to its right.
        type(kdtree_t), intent(inout) :: tree
        integer,        intent(in)    :: lo, hi, k, d
        integer  :: l, r, i, j, tmp
        real(dp) :: pivot
        l = lo; r = hi
        do while (l < r)
            pivot = tree%pts(d, tree%idx((l+r)/2))
            i = l; j = r
            do
                do while (tree%pts(d, tree%idx(i)) < pivot)
                    i = i + 1
                end do
                do while (tree%pts(d, tree%idx(j)) > pivot)
                    j = j - 1
                end do
                if (i <= j) then
                    tmp = tree%idx(i); tree%idx(i) = tree%idx(j); tree%idx(j) = tmp
                    i = i + 1; j = j - 1
                end if
                if (i > j) exit
            end do
            if (k <= j) then
                r = j
            else if (k >= i) then
                l = i
            else
                exit
            end if
        end do
    end subroutine select_kth

    subroutine kdtree_knn(tree, q, k, idx_out, d2_out, nfound)
        ! k nearest neighbors to query point q. idx_out/d2_out (size >= k) are
        ! returned in ascending squared-distance order; nfound <= k.
        type(kdtree_t), intent(in)  :: tree
        real(dp),       intent(in)  :: q(:)
        integer,        intent(in)  :: k
        integer,        intent(out) :: idx_out(:)
        real(dp),       intent(out) :: d2_out(:)
        integer,        intent(out) :: nfound
        nfound = 0
        idx_out = 0
        d2_out  = huge(1.0_dp)
        if (tree%n > 0 .and. k > 0) call knn_range(tree, 1, tree%n, q, k, idx_out, d2_out, nfound)
    end subroutine kdtree_knn

    recursive subroutine knn_range(tree, lo, hi, q, k, idx_out, d2_out, nfound)
        type(kdtree_t), intent(in)    :: tree
        integer,        intent(in)    :: lo, hi, k
        real(dp),       intent(in)    :: q(:)
        integer,        intent(inout) :: idx_out(:)
        real(dp),       intent(inout) :: d2_out(:)
        integer,        intent(inout) :: nfound
        integer  :: mid, dsp, p
        real(dp) :: dd, diff, worst
        if (lo > hi) return
        mid = (lo + hi) / 2
        dsp = tree%splitdim(mid)
        p   = tree%idx(mid)
        dd  = sum((q - tree%pts(:,p))**2)
        call knn_insert(idx_out, d2_out, nfound, k, p, dd)
        diff  = q(dsp) - tree%pts(dsp, p)
        if (diff < 0.0_dp) then
            call knn_range(tree, lo, mid-1, q, k, idx_out, d2_out, nfound)
            worst = merge(d2_out(k), huge(1.0_dp), nfound >= k)
            if (diff*diff < worst) call knn_range(tree, mid+1, hi, q, k, idx_out, d2_out, nfound)
        else
            call knn_range(tree, mid+1, hi, q, k, idx_out, d2_out, nfound)
            worst = merge(d2_out(k), huge(1.0_dp), nfound >= k)
            if (diff*diff < worst) call knn_range(tree, lo, mid-1, q, k, idx_out, d2_out, nfound)
        end if
    end subroutine knn_range

    subroutine knn_insert(idx_out, d2_out, nfound, k, p, dd)
        ! Insert candidate (p, dd) into the ascending-d2 bounded list of size k.
        integer,  intent(inout) :: idx_out(:)
        real(dp), intent(inout) :: d2_out(:)
        integer,  intent(inout) :: nfound
        integer,  intent(in)    :: k, p
        real(dp), intent(in)    :: dd
        integer :: i, pos
        if (nfound >= k) then
            if (dd >= d2_out(k)) return
        end if
        ! find insertion position
        pos = min(nfound, k-1) + 1
        do i = 1, min(nfound, k)
            if (dd < d2_out(i)) then
                pos = i
                exit
            end if
        end do
        ! shift down (drop the last if list already full)
        do i = min(nfound+1, k), pos+1, -1
            d2_out(i)  = d2_out(i-1)
            idx_out(i) = idx_out(i-1)
        end do
        d2_out(pos)  = dd
        idx_out(pos) = p
        if (nfound < k) nfound = nfound + 1
    end subroutine knn_insert

    subroutine kdtree_radius(tree, q, r, idx_out, d2_out, nfound)
        ! All neighbors within radius r of q (squared distance <= r^2).
        ! Output capped at size(idx_out); nfound is the number returned.
        type(kdtree_t), intent(in)  :: tree
        real(dp),       intent(in)  :: q(:)
        real(dp),       intent(in)  :: r
        integer,        intent(out) :: idx_out(:)
        real(dp),       intent(out) :: d2_out(:)
        integer,        intent(out) :: nfound
        nfound = 0
        if (tree%n > 0) call radius_range(tree, 1, tree%n, q, r*r, idx_out, d2_out, nfound)
    end subroutine kdtree_radius

    recursive subroutine radius_range(tree, lo, hi, q, r2, idx_out, d2_out, nfound)
        type(kdtree_t), intent(in)    :: tree
        integer,        intent(in)    :: lo, hi
        real(dp),       intent(in)    :: q(:), r2
        integer,        intent(inout) :: idx_out(:)
        real(dp),       intent(inout) :: d2_out(:)
        integer,        intent(inout) :: nfound
        integer  :: mid, dsp, p, cap
        real(dp) :: dd, diff
        if (lo > hi) return
        cap = size(idx_out)
        mid = (lo + hi) / 2
        dsp = tree%splitdim(mid)
        p   = tree%idx(mid)
        dd  = sum((q - tree%pts(:,p))**2)
        if (dd <= r2 .and. nfound < cap) then
            nfound = nfound + 1
            idx_out(nfound) = p
            d2_out(nfound)  = dd
        end if
        diff = q(dsp) - tree%pts(dsp, p)
        if (diff < 0.0_dp) then
            call radius_range(tree, lo, mid-1, q, r2, idx_out, d2_out, nfound)
            if (diff*diff <= r2) call radius_range(tree, mid+1, hi, q, r2, idx_out, d2_out, nfound)
        else
            call radius_range(tree, mid+1, hi, q, r2, idx_out, d2_out, nfound)
            if (diff*diff <= r2) call radius_range(tree, lo, mid-1, q, r2, idx_out, d2_out, nfound)
        end if
    end subroutine radius_range

end module kdtree
