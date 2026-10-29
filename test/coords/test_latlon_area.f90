program test_latlon_area
    ! Lat-lon cell areas and map metric (coords), sphere vs the default WGS84:
    !  (1) planet="Spherical Earth": area = R^2 dlon (sin(lat2)-sin(lat1)) to 0.25*dphi^2
    !      (the great-circle-edge corner polygon), sum(area) = 4 pi R^2;
    !  (2) the conservative lat-lon map conserves sum(f*area) in grid%area to the
    !      corner-polygon area error of the target (sphere);
    !  (3) grid_init_gaussian with no planet given: area = R^2 dlon w_j (its stated
    !      contract, the SHTns quadrature weight) -- fails while the default is WGS84,
    !      whose authalic zone areas differ from the sphere by up to 0.9% (+0.9% poles,
    !      -0.45% equator).
    use coords
    implicit none
    integer :: fails
    fails = 0
    call check_area("Spherical Earth", 0.5_dp)
    call check_area("Spherical Earth", 2.0_dp)
    call check_map("Spherical Earth")
    call check_gauss_default()
    if (fails > 0) then
        write(*,*) "FAILED:", fails
        error stop 1
    end if
    write(*,*) "OK"
contains
    subroutine mkgrid(g, nm, planet, dl)
        type(grid_class), intent(inout) :: g
        character(len=*), intent(in) :: nm, planet
        real(dp), intent(in) :: dl
        real(dp), allocatable :: x(:), y(:)
        integer :: i, nx, ny
        nx = nint(360.0_dp/dl); ny = nint(180.0_dp/dl)
        allocate(x(nx), y(ny))
        do i = 1, nx; x(i) = (real(i,dp)-0.5_dp)*dl; end do
        do i = 1, ny; y(i) = -90.0_dp + (real(i,dp)-0.5_dp)*dl; end do
        call grid_init(g, name=nm, mtype="latlon", units="degrees", planet=planet, x=x, y=y)
    end subroutine mkgrid

    subroutine check_area(planet, dl)
        character(len=*), intent(in) :: planet
        real(dp), intent(in) :: dl
        type(grid_class) :: g
        real(dp) :: R, ex, emax, tol
        integer :: j, ny
        call mkgrid(g, "area-"//planet, planet, dl)
        R = g%cs%planet%a; ny = g%G%ny; emax = 0.0_dp
        do j = 1, ny
            ex = R*R*dl*degrees_to_radians*(sin(g%G%y_bnds(j)*degrees_to_radians) - sin(g%G%y_bnds(j-1)*degrees_to_radians))
            emax = max(emax, abs(g%area(1,j)/ex - 1.0_dp))
        end do
        tol = 0.25_dp*(dl*degrees_to_radians)**2
        write(*,'(a,a,f5.2,a,es10.2,a,es10.2)') trim(planet), " dl=", dl, " max rel area err =", emax, "  tol", tol
        if (emax > tol) then; write(*,*) "FAIL: area"; fails = fails + 1; end if
        ex = abs(sum(g%area)/(4.0_dp*pi*R*R) - 1.0_dp)
        if (ex > 1.0e-12_dp) then; write(*,*) "FAIL: total area", ex; fails = fails + 1; end if
    end subroutine check_area

    subroutine check_gauss_default()
        type(grid_class) :: g
        real(dp) :: lat(48), w(48), emax
        integer :: j
        call grid_init_gaussian(g, "gauss-default", nlon=96, nlat=48)
        call gaussian_latitudes_calc(48, lat, w)
        emax = 0.0_dp
        do j = 1, 48
            emax = max(emax, abs(g%area(1,j)/(g%cs%planet%a**2*(2.0_dp*pi/96.0_dp)*w(j)) - 1.0_dp))
        end do
        write(*,'(a,es10.2)') "grid_init_gaussian default planet: max |area/(R2 dlon w) - 1| =", emax
        if (emax > 1.0e-12_dp) then; write(*,*) "FAIL: Gauss default area"; fails = fails + 1; end if
    end subroutine check_gauss_default

    subroutine check_map(planet)
        character(len=*), intent(in) :: planet
        type(grid_class) :: gs, gt
        type(map_class)  :: map
        real(dp), allocatable :: f(:,:), ft(:,:)
        logical, allocatable :: m2(:,:)
        real(dp) :: r
        integer :: i, j
        call mkgrid(gs, "src-"//planet, planet, 0.5_dp)
        call mkgrid(gt, "tgt-"//planet, planet, 3.0_dp)
        allocate(f(gs%G%nx,gs%G%ny), ft(gt%G%nx,gt%G%ny), m2(gt%G%nx,gt%G%ny))
        do j = 1, gs%G%ny; do i = 1, gs%G%nx
            f(i,j) = 1.0_dp + sin(gs%G%y(j)*degrees_to_radians)**2 + 0.3_dp*cos(gs%G%x(i)*degrees_to_radians)
        end do; end do
        call map_init_conservative(map, gs, gt)
        call map_field(map, "f", f, ft, stat="mean", mask2=m2)
        r = sum(ft*gt%area)/sum(f*gs%area) - 1.0_dp
        write(*,'(a,a,es10.2)') trim(planet), " map: sum(f*A) rel. change =", r
        if (abs(r) > 0.25_dp*(3.0_dp*degrees_to_radians)**2) then; write(*,*) "FAIL: map not conservative in grid%area"; fails = fails + 1; end if
    end subroutine check_map
end program
