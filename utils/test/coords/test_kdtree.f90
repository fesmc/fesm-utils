program test_kdtree
    ! Validate kdtree knn + radius queries against brute force on deterministic
    ! pseudo-random 3D points.

    use kdtree
    use precision, only: dp

    implicit none

    integer, parameter :: n = 3000, ndim = 3, k = 12
    real(dp) :: pts(ndim, n)
    type(kdtree_t) :: tree
    integer(8) :: s
    integer :: i, j, iq, nf, nfb
    real(dp) :: q(ndim), r
    integer  :: kidx(k), bidx(k)
    real(dp) :: kd2(k), bd2(k)
    integer  :: ridx(n)
    real(dp) :: rd2(n)
    integer  :: fails

    ! Deterministic pseudo-random points in [0,1)^3 via a simple LCG
    s = 123456789_8
    do i = 1, n
        do j = 1, ndim
            s = mod(1103515245_8*s + 12345_8, 2147483648_8)
            pts(j,i) = real(s, dp) / 2147483648.0_dp
        end do
    end do

    call kdtree_build(tree, pts)

    fails = 0

    ! --- knn vs brute force ---
    do iq = 1, 40
        q = pts(:, mod(iq*101, n) + 1)
        ! perturb so the query is not exactly a data point
        q = q + 1.0e-3_dp
        call kdtree_knn(tree, q, k, kidx, kd2, nf)
        call brute_knn(q)
        if (nf /= k) then
            write(*,*) "FAIL: knn nfound /= k", nf; fails = fails + 1
        end if
        do i = 1, k
            if (abs(kd2(i) - bd2(i)) > 1.0e-12_dp) then
                write(*,*) "FAIL: knn d2 mismatch at q", iq, " pos", i, kd2(i), bd2(i)
                fails = fails + 1
            end if
            if (kidx(i) /= bidx(i)) then
                write(*,*) "FAIL: knn idx mismatch at q", iq, " pos", i, kidx(i), bidx(i)
                fails = fails + 1
            end if
        end do
    end do

    ! --- radius vs brute force ---
    r = 0.08_dp
    do iq = 1, 40
        q = pts(:, mod(iq*131, n) + 1) + 1.0e-3_dp
        call kdtree_radius(tree, q, r, ridx, rd2, nf)
        nfb = 0
        do i = 1, n
            if (sum((q - pts(:,i))**2) <= r*r) nfb = nfb + 1
        end do
        if (nf /= nfb) then
            write(*,*) "FAIL: radius count mismatch at q", iq, nf, nfb; fails = fails + 1
        end if
        do i = 1, nf
            if (rd2(i) > r*r + 1.0e-12_dp) then
                write(*,*) "FAIL: radius returned point outside r"; fails = fails + 1
            end if
        end do
    end do

    call kdtree_free(tree)

    if (fails > 0) then
        write(*,*) "test_kdtree FAILED with", fails, "errors"
        stop 1
    end if
    write(*,*) "PASS: test_kdtree (knn + radius match brute force over 80 queries)"

contains

    subroutine brute_knn(qq)
        real(dp), intent(in) :: qq(ndim)
        real(dp) :: d2all(n)
        logical  :: used(n)
        integer  :: ii, jj, best
        real(dp) :: bestv
        do ii = 1, n
            d2all(ii) = sum((qq - pts(:,ii))**2)
        end do
        used = .false.
        do ii = 1, k
            best = 0; bestv = huge(1.0_dp)
            do jj = 1, n
                if (.not. used(jj) .and. d2all(jj) < bestv) then
                    bestv = d2all(jj); best = jj
                end if
            end do
            used(best) = .true.
            bidx(ii) = best
            bd2(ii)  = bestv
        end do
    end subroutine brute_knn

end program test_kdtree
