program test_map_init_names
    ! Check the map_init grid-names interface: building a map from grid *names*
    ! (reading grid_<name>.txt from a folder and reconstructing the grids
    ! internally) must produce exactly the same map as building it from the
    ! grid objects directly.

    use coords

    implicit none

    type(grid_class) :: gll, gant
    type(map_class)  :: m_obj, m_nam
    real(dp), allocatable :: vs(:,:), vt_obj(:,:), vt_nam(:,:)
    logical,  allocatable :: m2_obj(:,:), m2_nam(:,:)
    character(len=*), parameter :: fldr = "maps_names_test"
    integer :: i, j, fails, ndiff
    real(dp) :: dmax

    call execute_command_line("rm -rf "//fldr)
    call execute_command_line("mkdir -p "//fldr)

    fails = 0

    ! Source: 5deg global lon/lat; target: Antarctic polar-stereographic
    call grid_init(gll, name="ll5", mtype="latlon", units="degrees", planet="WGS84", &
                   x0=0.0_dp, dx=5.0_dp, nx=72, y0=-90.0_dp, dy=5.0_dp, ny=37)
    call grid_init(gant, name="ANT-128KM", mtype="polar_stereographic", units="kilometers", planet="WGS84", &
                   x0=-3040.0_dp, dx=128.0_dp, nx=48, y0=-3040.0_dp, dy=128.0_dp, ny=48, &
                   lambda=0.0_dp, phi=-71.0_dp, alpha=19.0_dp)

    ! Publish the grid descriptions so the names interface can find them
    call grid_cdo_write_desc_short(gll,  fldr)
    call grid_cdo_write_desc_short(gant, fldr)

    ! Reference map: from the grid objects directly
    call map_init(m_obj, gll, gant, max_neighbors=10, method="nn", fldr=fldr, load=.false.)

    ! Map under test: from the grid names alone
    call map_init(m_nam, "ll5", "ANT-128KM", max_neighbors=10, method="nn", fldr=fldr, load=.false.)

    ! Apply both maps to the same source field and compare the targets
    allocate(vs(72,37), vt_obj(48,48), vt_nam(48,48), m2_obj(48,48), m2_nam(48,48))
    do j = 1, 37
        do i = 1, 72
            vs(i,j) = gll%lat(i,j)
        end do
    end do

    call map_field(m_obj, "lat", vs, vt_obj, mask2=m2_obj)
    call map_field(m_nam, "lat", vs, vt_nam, mask2=m2_nam)

    dmax  = 0.0_dp
    ndiff = 0
    do j = 1, 48
        do i = 1, 48
            if (m2_obj(i,j) .neqv. m2_nam(i,j)) ndiff = ndiff + 1
            if (m2_obj(i,j) .and. m2_nam(i,j)) dmax = max(dmax, abs(vt_obj(i,j) - vt_nam(i,j)))
        end do
    end do

    write(*,"(a,es10.2)") " names-vs-objects max field diff =", dmax
    if (ndiff /= 0) then
        write(*,"(a,i0,a)") " FAIL: mask differs at ", ndiff, " target cells"
        fails = fails + 1
    end if
    if (dmax > 1.0e-12_dp) then
        write(*,*) "FAIL: mapped field differs between names and objects interfaces"
        fails = fails + 1
    end if

    call execute_command_line("rm -rf "//fldr)

    write(*,*)
    if (fails == 0) then
        write(*,*) "PASS: test_map_init_names (names interface matches objects interface)"
    else
        write(*,"(a,i0,a)") " FAIL: test_map_init_names (", fails, " check(s) failed)"
        stop 1
    end if

end program test_map_init_names
