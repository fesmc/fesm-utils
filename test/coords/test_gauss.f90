program test_gauss
    ! Gaussian-latitude generator — check against known Legendre-root latitudes.

    use gaussian_latitudes
    use precision, only: dp

    implicit none

    real(dp) :: lat2(2), w2(2), lat4(4)

    ! N=2: roots of P2 at x = +-1/sqrt(3) -> lat = +-35.2643897 deg; weights = 1 each
    call gaussian_latitudes_calc(2, lat2, w2)
    write(*,*) "N=2 lats:", lat2
    write(*,*) "N=2 wts :", w2, " sum=", sum(w2)
    if (abs(lat2(1) + 35.2643897_dp) > 1.d-4 .or. abs(lat2(2) - 35.2643897_dp) > 1.d-4) then
        write(*,*) "FAIL: N=2 Gaussian latitudes wrong"; stop 1
    end if
    if (abs(sum(w2) - 2.0_dp) > 1.d-12) then
        write(*,*) "FAIL: weights do not sum to 2"; stop 1
    end if

    ! N=4: latitudes approx +-19.8758, +-59.4400 deg (ascending S->N)
    call gaussian_latitudes_calc(4, lat4)
    write(*,*) "N=4 lats:", lat4
    if (abs(lat4(1) + 59.44_dp) > 1.d-2 .or. abs(lat4(2) + 19.876_dp) > 1.d-2 .or. &
        abs(lat4(3) - 19.876_dp) > 1.d-2 .or. abs(lat4(4) - 59.44_dp) > 1.d-2) then
        write(*,*) "FAIL: N=4 Gaussian latitudes wrong"; stop 1
    end if

    write(*,*) "PASS: test_gauss"

end program test_gauss
