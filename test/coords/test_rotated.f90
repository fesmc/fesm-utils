program test_rotated
    ! Rotated-pole projection: round-trip, pole anchor, and points_init wiring.
    ! Uses the CORDEX EUR-11 rotated north pole (lon=-162.0, lat=39.25).

    use coords

    implicit none

    type(projection_class) :: proj
    type(planet_class)     :: pl
    type(points_class)     :: pts
    real(dp) :: glon, glat, rlon, rlat, glon2, glat2
    real(dp), parameter :: lon_p = -162.0_dp, lat_p = 39.25_dp

    call planet_init(pl, "Spherical Earth")
    call projection_init(proj, "rotated_latitude_longitude", pl, lambda=lon_p, phi=lat_p)

    ! Anchor: the rotated north pole (rlon=0,rlat=90) is the geographic pole position
    call oblimap_projection_inverse(0.d0, 90.d0, glon, glat, proj)
    write(*,*) "rotated N-pole -> geographic:", glon, glat, " (expect", lon_p, lat_p, ")"
    if (abs(glat - lat_p) > 1.d-6 .or. abs(glon - lon_p) > 1.d-6) then
        write(*,*) "FAIL: rotated-pole anchor wrong"; stop 1
    end if

    ! Round-trip: geographic -> rotated -> geographic
    call oblimap_projection(10.d0, 50.d0, rlon, rlat, proj)
    call oblimap_projection_inverse(rlon, rlat, glon2, glat2, proj)
    write(*,*) "geo(10,50) -> rot(", rlon, rlat, ") -> geo(", glon2, glat2, ")"
    if (abs(glon2 - 10.d0) > 1.d-9 .or. abs(glat2 - 50.d0) > 1.d-9) then
        write(*,*) "FAIL: rotated-pole round-trip not identity"; stop 1
    end if

    ! points_init integration: x/y are rotated lon/lat; lon/lat come out geographic
    call points_init(pts, name="rot-pole", mtype="rotated_latitude_longitude", units="degrees", &
                     lambda=lon_p, phi=lat_p, lon180=.TRUE., x=[0.d0], y=[90.d0])
    write(*,*) "points_init rotated pole -> lon,lat:", pts%lon(1), pts%lat(1)
    if (abs(pts%lat(1) - lat_p) > 1.d-4 .or. abs(pts%lon(1) - lon_p) > 1.d-4) then
        write(*,*) "FAIL: points_init rotated-pole wiring wrong"; stop 1
    end if

    write(*,*) "PASS: test_rotated"

end program test_rotated
