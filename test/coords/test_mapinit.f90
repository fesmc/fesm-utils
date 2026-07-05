program test_mapinit
    ! Validate that the k-d-tree neighbor search in map_init returns the same
    ! neighbor sets and distances as a brute-force geodesic search.

    use coords

    implicit none

    type(points_class) :: src, tgt
    type(map_class)    :: map
    integer, parameter :: nx = 12, ny = 12, ns = nx*ny, nt = 6, k = 8
    real(dp) :: slon(ns), slat(ns), tlon(nt), tlat(nt)
    real(dp) :: a, f
    integer  :: i, j, c, ix, iy, iq, j1, j2, nl
    integer  :: bidx(k)
    real(dp) :: bd2(k)
    integer  :: fails

    ! Source: 12x12 lon/lat lattice over Greenland-ish region
    c = 0
    do iy = 1, ny
        do ix = 1, nx
            c = c + 1
            slon(c) = -60.0_dp + real(ix-1,dp)*3.0_dp
            slat(c) =  60.0_dp + real(iy-1,dp)*1.5_dp
        end do
    end do
    ! Targets: a handful of points inside the region
    tlon = [-50.0_dp, -40.0_dp, -30.0_dp, -55.0_dp, -35.0_dp, -45.0_dp]
    tlat = [ 65.0_dp,  70.0_dp,  72.0_dp,  68.0_dp,  75.0_dp,  63.0_dp]

    call points_init(src, name="src", mtype="latlon", units="degrees", x=slon, y=slat)
    call points_init(tgt, name="tgt", mtype="latlon", units="degrees", x=tlon, y=tlat)

    call map_init(map, src, tgt, max_neighbors=k)

    a = src%cs%planet%a; f = src%cs%planet%f
    fails = 0

    do i = 1, nt
        ! neighbor links for target i
        j1 = map%wm%dst_off(i); j2 = map%wm%dst_off(i+1) - 1
        nl = j2 - j1 + 1
        if (nl /= k) then
            write(*,*) "FAIL: target", i, "got", nl, "neighbors, expected", k
            fails = fails + 1
        end if

        ! brute-force k nearest by geodesic distance
        call brute_knn(tlon(i), tlat(i))

        ! The meaningful invariant is the multiset of neighbor distances: the
        ! map's sorted distances must equal the brute-force k-nearest distances.
        ! Specific indices among exactly-equidistant sources are arbitrary, so
        ! we compare distances (ascending), not index order.
        do c = 1, min(nl, k)
            if (abs(map%wm%dist(j1+c-1) - sqrt(bd2(c))) > 1.0e-6_dp) then
                write(*,*) "FAIL: target", i, "neighbor", c, "dist", &
                           map%wm%dist(j1+c-1), "vs brute", sqrt(bd2(c))
                fails = fails + 1
            end if
        end do
    end do

    if (fails > 0) then
        write(*,*) "test_mapinit FAILED with", fails, "errors"; stop 1
    end if
    write(*,*) "PASS: test_mapinit (kdtree neighbors == brute-force geodesic over", nt, "targets)"

contains

    subroutine brute_knn(qlon, qlat)
        real(dp), intent(in) :: qlon, qlat
        real(dp) :: d(ns)
        logical  :: used(ns)
        integer  :: ii, jj, best
        real(dp) :: bv, dd
        do ii = 1, ns
            dd = planet_distance(a, f, qlon, qlat, slon(ii), slat(ii))
            d(ii) = dd*dd          ! store squared for comparison convenience
        end do
        used = .false.
        do ii = 1, k
            best = 0; bv = huge(1.0_dp)
            do jj = 1, ns
                if (.not. used(jj) .and. d(jj) < bv) then
                    bv = d(jj); best = jj
                end if
            end do
            used(best) = .true.
            bidx(ii) = best
            bd2(ii)  = bv
        end do
    end subroutine brute_knn

end program test_mapinit
