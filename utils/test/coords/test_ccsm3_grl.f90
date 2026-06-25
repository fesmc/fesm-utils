program test_ccsm3_grl
    ! End-to-end regression: global CCSM3 lat-lon -> regional Greenland oblique-
    ! stereographic, the canonical use case. Map a latitude field and check it
    ! reproduces the target-grid latitudes (a smooth field a good interpolant
    ! must recover closely). Exercises grid_init (both types), cross-projection
    ! map_init (geodesic path), and map_field.

    use coords

    implicit none

    type(grid_class) :: gCCSM3, gGRL
    type(map_class)  :: map
    real(dp), allocatable :: lat_src(:,:), lat_out(:,:)
    logical,  allocatable :: m2(:,:)
    real(dp) :: emax_bilin, emax_shep
    integer  :: fails

    ! Global 2-degree lat-lon source
    call grid_init(gCCSM3, name="CCSM3-T42", mtype="latlon", units="degrees", &
                   x0=0.0_dp, dx=2.0_dp, nx=180, y0=-90.0_dp, dy=2.0_dp, ny=91)

    ! Regional oblique-stereographic target over Greenland
    call grid_init(gGRL, name="GRL-20KM", mtype="stereographic", units="kilometers", &
                   dx=20.0_dp, nx=76, dy=20.0_dp, ny=151, &
                   lambda=-40.0_dp, phi=72.0_dp, alpha=7.5_dp)

    allocate(lat_src(180,91), lat_out(76,151), m2(76,151))
    lat_src = gCCSM3%lat

    call map_init(map, gCCSM3, gGRL, max_neighbors=10)

    fails = 0

    ! Bilinear
    call map_field(map, "lat", lat_src, lat_out, method="bilinear", mask2=m2)
    emax_bilin = maxval(abs(lat_out - gGRL%lat), mask=m2)
    write(*,*) "filled (bilinear): ", count(m2), " of", size(m2)
    write(*,*) "bilinear max |lat_out - lat_grl| =", emax_bilin, "deg"

    ! Shepard (IDW)
    call map_field(map, "lat", lat_src, lat_out, method="shepard", mask2=m2)
    emax_shep = maxval(abs(lat_out - gGRL%lat), mask=m2)
    write(*,*) "shepard  max |lat_out - lat_grl| =", emax_shep, "deg"

    if (count(m2) < size(m2)/2) then
        write(*,*) "FAIL: fewer than half the target points were filled"; fails = fails + 1
    end if
    ! A 2-degree source interpolated to a smooth latitude field should be well
    ! under 1 degree everywhere.
    if (emax_bilin > 1.0_dp) then
        write(*,*) "FAIL: bilinear latitude error too large"; fails = fails + 1
    end if
    if (emax_shep > 1.0_dp) then
        write(*,*) "FAIL: shepard latitude error too large"; fails = fails + 1
    end if

    if (fails > 0) stop 1
    write(*,*) "PASS: test_ccsm3_grl"

end program test_ccsm3_grl
