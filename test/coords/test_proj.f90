program test_proj
    ! Projection check — define a polar-stereographic point set (Bamber corners)
    ! and inverse-project x/y to lon/lat via the coords/oblimap machinery.

    use coords

    implicit none

    type(points_class) :: pts
    integer :: i

    call points_init(pts, name="Bamber-corners", mtype="polar_stereographic", units="kilometers", &
                     lambda=-39.d0, phi=90.d0, alpha=19.d0, lon180=.TRUE., &
                     x=[-800.d0,-800.d0,700.d0,700.d0], y=[-3400.d0,-600.d0,-600.d0,-3400.d0])

    write(*,"(a)") "        x         y    =>        lat       lon"
    do i = 1, size(pts%x,1)
        write(*,"(2f10.3,a,2f10.3)") pts%x(i), pts%y(i), " => ", pts%lat(i), pts%lon(i)
    end do

    ! Sanity: all corners should land in the northern hemisphere over Greenland.
    if (any(pts%lat < 55.d0) .or. any(pts%lat > 85.d0)) then
        write(*,*) "FAIL: projected latitudes outside expected Greenland range"
        stop 1
    end if
    ! Transverse Mercator (UTM zone 18S, WGS84): forward and inverse against
    ! pyproj (EPSG:4326 -> EPSG:32718) reference values [m].
    call check_utm18s()

    write(*,*) "PASS: test_proj"

contains

    subroutine check_utm18s()
        integer, parameter :: n = 7
        real(dp), parameter :: lon_ref(n) = [-75.0d0, -73.5d0, -73.9318d0, -73.2159d0, -72.0d0, -78.0d0, -75.0d0]
        real(dp), parameter :: lat_ref(n) = [  0.0d0, -46.7d0, -46.5568d0, -46.8647d0, -50.0d0, -10.0d0, -60.0d0]
        real(dp), parameter :: x_ref(n) = [500000.0000d0, 614674.3413d0, 581879.3585d0, 635977.8524d0, &
                                           714984.2367d0, 171071.2639d0, 500000.0000d0]
        real(dp), parameter :: y_ref(n) = [10000000.0000d0, 4827080.3547d0, 4843530.9793d0, 4808325.9553d0, &
                                           4457055.9814d0, 8893091.1458d0, 3348588.8096d0]
        type(points_class) :: p_ll, p_xy
        real(dp) :: err_x, err_y, err_lon, err_lat

        ! Forward: lon/lat -> x/y [m]
        call points_init(p_ll, name="utm18s-ll", mtype="transverse_mercator", units="meters", planet="WGS84", &
                         lambda=-75.d0, phi=0.d0, k0=0.9996d0, x_e=500000.d0, y_n=10000000.d0, &
                         x=lon_ref, y=lat_ref, latlon=.TRUE., lon180=.TRUE.)
        err_x = maxval(abs(p_ll%x - x_ref))
        err_y = maxval(abs(p_ll%y - y_ref))

        ! Inverse: x/y [km] -> lon/lat
        call points_init(p_xy, name="utm18s-xy", mtype="transverse_mercator", units="kilometers", planet="WGS84", &
                         lambda=-75.d0, phi=0.d0, k0=0.9996d0, x_e=500000.d0, y_n=10000000.d0, &
                         x=x_ref*1.d-3, y=y_ref*1.d-3, lon180=.TRUE.)
        err_lon = maxval(abs(p_xy%lon - lon_ref))
        err_lat = maxval(abs(p_xy%lat - lat_ref))

        write(*,"(a,2es10.2,a,2es10.2,a)") " transverse_mercator vs pyproj: |dx|,|dy| = ", err_x, err_y, &
                                            " m, |dlon|,|dlat| = ", err_lon, err_lat, " deg"
        if (err_x > 1.d-3 .or. err_y > 1.d-3 .or. err_lon > 1.d-8 .or. err_lat > 1.d-8) then
            write(*,*) "FAIL: transverse_mercator does not match pyproj (UTM 18S)"
            stop 1
        end if
    end subroutine check_utm18s

end program test_proj
