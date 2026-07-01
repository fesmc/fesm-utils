program test_coupler
    ! Exercise the coupler layer over the coords mapping core:
    !   * bilin downscale (coarse -> fine) is exact for a linear field,
    !   * the map cache builds each (src,dst,method) edge once and reuses it,
    !   * method defaults to "con", and is part of the cache key,
    !   * remap sizes the destination array on demand (2D and 3D),
    !   * coupler_prime pre-builds without adding duplicate maps.

    use coords
    use coupler

    implicit none

    type(coupler_class)   :: cpl
    type(grid_class)      :: gs, gt
    real(dp), allocatable :: fs(:,:), ft(:,:), fback(:,:)
    real(dp), allocatable :: fs3(:,:,:), ft3(:,:,:)
    integer  :: i, j, k, fails
    real(dp) :: emax
    character(len=*), parameter :: fldr = "maps_coupler_test"

    fails = 0
    call execute_command_line("rm -rf "//fldr)

    ! Coarse (100 km) and fine (50 km) grids over the same [0,1000] km extent.
    call grid_init(gs, name="coupler-coarse", mtype="cartesian", units="kilometers", &
                   x0=0.0_dp, dx=100.0_dp, nx=11, y0=0.0_dp, dy=100.0_dp, ny=11)
    call grid_init(gt, name="coupler-fine",   mtype="cartesian", units="kilometers", &
                   x0=0.0_dp, dx=50.0_dp,  nx=21, y0=0.0_dp, dy=50.0_dp,  ny=21)

    call coupler_init(cpl, map_fldr=fldr)
    call coupler_add_grid(cpl, "coupler-coarse", gs)
    call coupler_add_grid(cpl, "coupler-fine",   gt)

    ! Linear source field on the coarse grid.
    allocate(fs(11,11))
    do j = 1, 11
        do i = 1, 11
            fs(i,j) = flin(gs%x(i,j), gs%y(i,j))
        end do
    end do

    ! --- Test 1: bilin downscale is exact for a linear field (interior) ------
    call remap(cpl, fs, "coupler-coarse", ft, "coupler-fine", method="bilin")
    if (.not. allocated(ft)) then
        write(*,*) "FAIL: remap did not allocate the destination"; fails = fails + 1
    else if (size(ft,1) /= 21 .or. size(ft,2) /= 21) then
        write(*,*) "FAIL: destination sized ", shape(ft), " expected 21 x 21"; fails = fails + 1
    else
        emax = 0.0_dp
        do j = 2, 20
            do i = 2, 20
                emax = max(emax, abs(ft(i,j) - flin(gt%x(i,j), gt%y(i,j))))
            end do
        end do
        write(*,*) "bilin downscale interior max err =", emax
        if (emax > 1.0e-6_dp) then
            write(*,*) "FAIL: bilin downscale not exact for linear field"; fails = fails + 1
        end if
    end if

    ! --- Test 2: the map cache builds one edge and reuses it -----------------
    if (cpl%nmaps /= 1) then
        write(*,*) "FAIL: expected 1 cached map, got", cpl%nmaps; fails = fails + 1
    end if
    call remap(cpl, fs, "coupler-coarse", ft, "coupler-fine", method="bilin")
    if (cpl%nmaps /= 1) then
        write(*,*) "FAIL: repeated remap rebuilt the map (nmaps=", cpl%nmaps, ")"; fails = fails + 1
    end if

    ! --- Test 3: method defaults to "con", and is a distinct cache edge ------
    call remap(cpl, ft, "coupler-fine", fback, "coupler-coarse")   ! no method -> con
    if (cpl%nmaps /= 2) then
        write(*,*) "FAIL: default-method remap did not add a con edge (nmaps=", cpl%nmaps, ")"; fails = fails + 1
    end if
    if (trim(cpl%maps(2)%method) /= "con") then
        write(*,*) "FAIL: default method was '", trim(cpl%maps(2)%method), "', expected 'con'"; fails = fails + 1
    end if
    if (size(fback,1) /= 11 .or. size(fback,2) /= 11) then
        write(*,*) "FAIL: aggregated field sized ", shape(fback), " expected 11 x 11"; fails = fails + 1
    else
        ! Conservative aggregation of a linear field ~ value at coarse centres (interior).
        emax = 0.0_dp
        do j = 2, 10
            do i = 2, 10
                emax = max(emax, abs(fback(i,j) - flin(gs%x(i,j), gs%y(i,j))))
            end do
        end do
        write(*,*) "con aggregate interior max err  =", emax
        if (emax > 1.0e-3_dp) then
            write(*,*) "FAIL: con aggregation of linear field off by", emax; fails = fails + 1
        end if
    end if

    ! --- Test 4: coupler_prime does not duplicate an existing edge -----------
    call coupler_prime(cpl, "coupler-coarse", "coupler-fine", "bilin")
    if (cpl%nmaps /= 2) then
        write(*,*) "FAIL: prime duplicated an existing edge (nmaps=", cpl%nmaps, ")"; fails = fails + 1
    end if

    ! --- Test 5: 3D remap matches per-level 2D remap ------------------------
    allocate(fs3(11,11,3))
    do k = 1, 3
        fs3(:,:,k) = fs + real(k, dp)          ! three offset linear levels
    end do
    call remap(cpl, fs3, "coupler-coarse", ft3, "coupler-fine", method="bilin")
    if (size(ft3,1) /= 21 .or. size(ft3,2) /= 21 .or. size(ft3,3) /= 3) then
        write(*,*) "FAIL: 3D destination sized ", shape(ft3), " expected 21 x 21 x 3"; fails = fails + 1
    else
        emax = 0.0_dp
        do k = 1, 3
            call remap(cpl, fs3(:,:,k), "coupler-coarse", ft, "coupler-fine", method="bilin")
            emax = max(emax, maxval(abs(ft3(:,:,k) - ft)))
        end do
        write(*,*) "3D vs per-level 2D max diff      =", emax
        if (emax > 1.0e-12_dp) then
            write(*,*) "FAIL: 3D remap disagrees with per-level 2D remap"; fails = fails + 1
        end if
    end if

    call execute_command_line("rm -rf "//fldr)

    if (fails > 0) stop 1
    write(*,*) "PASS: test_coupler"

contains

    pure function flin(x, y) result(f)
        real(dp), intent(in) :: x, y
        real(dp) :: f
        f = 2.0_dp + 0.003_dp*x + 0.005_dp*y
    end function flin

end program test_coupler
