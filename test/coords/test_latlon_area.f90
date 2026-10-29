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
    call check_gauss_radius()
    call check_radius_carried()
    call check_radius_signature()
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

    subroutine check_gauss_radius()
        ! radius override: Gaussian grid with R = 6371000 (VILMA2, CLIMBER-X)
        type(grid_class) :: g, g0, gd
        real(dp), parameter :: R = 6371000.0_dp
        real(dp) :: lat(48), w(48), emax, dlon
        integer :: j
        logical :: same
        character(len=256) :: why
        call grid_init_gaussian(g, "gauss-R", nlon=96, nlat=48, radius=R)
        call grid_init_gaussian(g0, "gauss-R0", nlon=96, nlat=48)
        call gaussian_latitudes_calc(48, lat, w)
        dlon = 2.0_dp*pi/96.0_dp
        emax = 0.0_dp
        do j = 1, 48
            emax = max(emax, abs(g%area(1,j)/(R*R*dlon*w(j)) - 1.0_dp))
        end do
        write(*,'(a,es10.2)') "grid_init_gaussian radius=6371000: max |area/(R2 dlon w) - 1| =", emax
        if (emax > 1.0e-12_dp) then; write(*,*) "FAIL: radius area"; fails = fails + 1; end if
        if (abs(sum(g%area)/(4.0_dp*pi*R*R) - 1.0_dp) > 1.0e-12_dp) then
            write(*,*) "FAIL: radius total area"; fails = fails + 1
        end if
        if (g%cs%planet%a /= R .or. g%cs%planet%R /= R .or. g%cs%planet%f /= 0.0_dp .or. .not. g%cs%planet%is_sphere) then
            write(*,*) "FAIL: planet a/R/f"; fails = fails + 1
        end if
        ! default unchanged
        if (g0%cs%planet%a /= 6.371221E6_dp) then; write(*,*) "FAIL: default radius changed"; fails = fails + 1; end if
        ! the cache signature distinguishes the radius
        same = grid_signature_match(grid_signature_str(g), grid_signature_str(g0), why)
        if (same) then; write(*,*) "FAIL: signature ignores radius"; fails = fails + 1; end if
        ! a grid derived from it keeps the radius; the planet-name form works for grid_init
        call grid_init(gd, g, "derived", nx=10, ny=5, x0=0.0_dp, dx=36.0_dp, y0=-80.0_dp, dy=32.0_dp)
        if (gd%cs%planet%a /= R) then; write(*,*) "FAIL: derived grid lost the radius"; fails = fails + 1; end if
        call grid_init(gd, "named", mtype="latlon", units="degrees", planet="Spherical Earth", radius=R, &
                       x0=0.0_dp, dx=10.0_dp, nx=36, y0=-85.0_dp, dy=10.0_dp, ny=18)
        if (gd%cs%planet%a /= R) then; write(*,*) "FAIL: grid_init radius"; fails = fails + 1; end if
        if (abs(sum(gd%area)/(4.0_dp*pi*R*R) - 1.0_dp) > 1.0e-12_dp) then
            write(*,*) "FAIL: grid_init radius total area"; fails = fails + 1
        end if
    end subroutine check_gauss_radius

    subroutine check_radius_carried()
        ! The radius is a field of the planet: the name stays "Spherical Earth", and every
        ! path that copies a planet keeps it (grid from grid / points, points from grid / points,
        ! the cdo grid description, the map namelist); WGS84 copies stay WGS84.
        type(grid_class)   :: g, g2, g3
        type(points_class) :: p, p2
        real(dp), parameter :: R = 6371000.0_dp
        character(len=*), parameter :: fldr = "logs/check"
        character(len=*), parameter :: nml = "logs/check/test_latlon_area_radius.nml"
        integer :: u

        call grid_init(g, "ps-R", mtype="polar_stereographic", units="km", planet="Spherical Earth", radius=R, &
                       lambda=-40.0_dp, phi=71.0_dp, x0=-300.0_dp, dx=20.0_dp, nx=31, y0=-300.0_dp, dy=20.0_dp, ny=31)
        if (trim(g%cs%planet%name) /= "Spherical Earth") then
            write(*,*) "FAIL: planet name carries the radius: ", trim(g%cs%planet%name); fails = fails + 1
        end if
        if (g%cs%planet%a /= R .or. g%cs%planet%R /= R .or. g%cs%proj%R /= R) then
            write(*,*) "FAIL: radius not in planet / projection"; fails = fails + 1
        end if

        call grid_init(g2, g, "from-grid", nx=11, ny=11, x0=-100.0_dp, dx=20.0_dp, y0=-100.0_dp, dy=20.0_dp)
        call expect_radius("grid_init from grid", g2%cs%planet, R)

        call points_init(p, g, "pts-from-grid", x=[0.0_dp, 50.0_dp], y=[0.0_dp, 20.0_dp])
        call expect_radius("points_init from grid", p%cs%planet, R)
        call points_init(p2, p, "pts-from-pts", x=[0.0_dp, 50.0_dp], y=[0.0_dp, 20.0_dp])
        call expect_radius("points_init from points", p2%cs%planet, R)
        call grid_init(g3, p, "grid-from-pts", nx=11, ny=11, x0=-100.0_dp, dx=20.0_dp, y0=-100.0_dp, dy=20.0_dp)
        call expect_radius("grid_init from points", g3%cs%planet, R)

        ! cdo grid description round trip
        call execute_command_line("mkdir -p "//fldr)
        call grid_cdo_write_desc_short(g2, fldr)
        call grid_cdo_read_desc(g3, "from-grid", fldr)
        call expect_radius("cdo grid description round trip", g3%cs%planet, R)

        ! map namelist: planet_radius [m]
        open(newunit=u, file=nml, status="replace")
        write(u,"(a)") "&map"
        write(u,"(a)") "  name = 'nml-R', mtype = 'latlon', units = 'degrees', planet = 'Spherical Earth',"
        write(u,"(a)") "  planet_radius = 6371000.0, lon180 = .false."
        write(u,"(a)") "/"
        close(u)
        call grid_init(g3, nml, x0=0.0_dp, dx=10.0_dp, nx=36, y0=-85.0_dp, dy=10.0_dp, ny=18)
        call expect_radius("grid_init from namelist (planet_radius)", g3%cs%planet, R)
        if (abs(sum(g3%area)/(4.0_dp*pi*R*R) - 1.0_dp) > 1.0e-12_dp) then
            write(*,*) "FAIL: namelist radius total area"; fails = fails + 1
        end if

        ! an ellipsoid (WGS84) stays an ellipsoid through the copy
        call grid_init(g, "wgs", mtype="latlon", units="degrees", x0=0.0_dp, dx=10.0_dp, nx=36, y0=-85.0_dp, &
                       dy=10.0_dp, ny=18)
        call grid_init(g2, g, "wgs-derived", nx=10, ny=5, x0=0.0_dp, dx=36.0_dp, y0=-80.0_dp, dy=32.0_dp)
        if (g2%cs%planet%is_sphere .or. g2%cs%planet%a /= g%cs%planet%a .or. g2%cs%planet%f /= g%cs%planet%f) then
            write(*,*) "FAIL: WGS84 copy changed"; fails = fails + 1
        end if
    end subroutine check_radius_carried

    subroutine expect_radius(what, pl, R)
        character(len=*),   intent(in) :: what
        type(planet_class), intent(in) :: pl
        real(dp),           intent(in) :: R
        if (pl%a /= R .or. pl%R /= R .or. pl%f /= 0.0_dp .or. .not. pl%is_sphere &
            .or. trim(pl%name) /= "Spherical Earth") then
            write(*,*) "FAIL: radius lost: "//what; fails = fails + 1
        end if
    end subroutine expect_radius

    subroutine check_radius_signature()
        ! The cache signature of a grid with a radius is unchanged against the records
        ! computed at 3ba56a9 (radius carried in the planet name), so existing caches stay valid.
        type(grid_class) :: g, g2
        real(dp), parameter :: R = 6371000.0_dp
        character(len=*), parameter :: sig_gauss = "coords-sig-1 mtype=gaussian units=degrees nx=96 ny=48 lon180=F "// &
            "r=6.3710000000000000E+006,0.0000000000000000E+000,1.0000000000000000E+000,0.0000000000000000E+000,"// &
            "0.0000000000000000E+000,0.0000000000000000E+000,0.0000000000000000E+000,0.0000000000000000E+000,"// &
            "1.0000000000000000E+000,0.0000000000000000E+000,3.5625000000000000E+002,1.7100000000000000E+004,"// &
            "-8.7159094555862850E+001,8.7159094555862850E+001,1.7053025658242404E-013,-1.8750000000000000E+000,"// &
            "3.5812500000000000E+002,1.7278125000000000E+004,-9.0000000000000000E+001,9.0000000000000000E+001,"// &
            "0.0000000000000000E+000"
        character(len=*), parameter :: sig_latlon = "coords-sig-1 mtype=latlon units=degrees nx=36 ny=18 lon180=F "// &
            "r=6.3710000000000000E+006,0.0000000000000000E+000,1.0000000000000000E+000,0.0000000000000000E+000,"// &
            "0.0000000000000000E+000,0.0000000000000000E+000,0.0000000000000000E+000,0.0000000000000000E+000,"// &
            "1.0000000000000000E+000,0.0000000000000000E+000,3.5000000000000000E+002,6.3000000000000000E+003,"// &
            "-8.5000000000000000E+001,8.5000000000000000E+001,0.0000000000000000E+000,-5.0000000000000000E+000,"// &
            "3.5500000000000000E+002,6.4750000000000000E+003,-9.0000000000000000E+001,9.0000000000000000E+001,"// &
            "0.0000000000000000E+000"
        character(len=*), parameter :: sig_ps = "coords-sig-1 mtype=polar_stereographic units=km nx=11 ny=11 lon180=F "// &
            "r=6.3710000000000000E+006,0.0000000000000000E+000,1.0000000000000000E+003,-4.0000000000000000E+001,"// &
            "7.1000000000000000E+001,1.9000000000000000E+001,0.0000000000000000E+000,0.0000000000000000E+000,"// &
            "1.0000000000000000E+000,-1.0000000000000000E+002,1.0000000000000000E+002,0.0000000000000000E+000,"// &
            "-1.0000000000000000E+002,1.0000000000000000E+002,0.0000000000000000E+000,-1.1000000000000000E+002,"// &
            "1.1000000000000000E+002,0.0000000000000000E+000,-1.1000000000000000E+002,1.1000000000000000E+002,"// &
            "0.0000000000000000E+000"

        call grid_init_gaussian(g, "sig-gauss", nlon=96, nlat=48, radius=R)
        call expect_sig("gaussian radius", g, sig_gauss)
        call grid_init(g, "sig-latlon", mtype="latlon", units="degrees", planet="Spherical Earth", radius=R, &
                       x0=0.0_dp, dx=10.0_dp, nx=36, y0=-85.0_dp, dy=10.0_dp, ny=18)
        call expect_sig("latlon radius", g, sig_latlon)
        call grid_init(g, "sig-ps", mtype="polar_stereographic", units="km", planet="Spherical Earth", radius=R, &
                       lambda=-40.0_dp, phi=71.0_dp, x0=-300.0_dp, dx=20.0_dp, nx=31, y0=-300.0_dp, dy=20.0_dp, ny=31)
        call grid_init(g2, g, "sig-ps-derived", nx=11, ny=11, x0=-100.0_dp, dx=20.0_dp, y0=-100.0_dp, dy=20.0_dp)
        call expect_sig("derived polar_stereographic radius", g2, sig_ps)
    end subroutine check_radius_signature

    subroutine expect_sig(what, g, sig0)
        character(len=*), intent(in) :: what, sig0
        type(grid_class), intent(in) :: g
        character(len=256) :: why
        if (.not. grid_signature_match(sig0, grid_signature_str(g), why)) then
            write(*,*) "FAIL: signature differs from 3ba56a9: "//what//": "//trim(why); fails = fails + 1
        end if
    end subroutine expect_sig

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
