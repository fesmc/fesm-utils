program test_geodinverse
    ! Geodesic inverse distance (Karney/GeographicLib) — checks planet_distance
    ! against an analytically exact case: along the equator the geodesic is the
    ! equator itself, so the distance from (0E,0N) to (90E,0N) is a*pi/2.

    use planet

    implicit none

    double precision, parameter :: pi = 4.d0*atan(1.d0)
    double precision :: a, f, dist, expected, err

    ! WGS84 values
    a = 6378137d0
    f = 1d0/298.257223563d0

    ! planet_distance takes (a,f, lon1,lat1, lon2,lat2) — longitude first.
    ! Equator arc 0E->90E at lat 0: exact length = a * (pi/2).
    dist     = planet_distance(a, f, 0.d0, 0.d0, 90.d0, 0.d0)
    expected = a * pi / 2.d0
    err      = abs(dist - expected)

    write(*,*) "equator arc 0E->90E at lat 0"
    write(*,*) "dist     =", dist,     "m"
    write(*,*) "expected =", expected, "m  (a*pi/2)"
    write(*,*) "abs err  =", err,      "m"

    if (err > 1.d-2) then
        write(*,*) "FAIL: geodesic distance off by more than 1 cm"
        stop 1
    end if
    write(*,*) "PASS: test_geodinverse"

end program test_geodinverse
