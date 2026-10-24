program test_gauss_bounds
    ! Gaussian grid bounds from the Gauss weights (docs/dev/coupler-global-design.md, 3):
    ! grid_init_gaussian(nlat=96, nlon=192, "Spherical Earth"):
    !   * bounds monotonic, +-90 at the ends, each node inside its band,
    !   * area = R^2*dlon*w_j to 1e-13, sum(area) = 4*pi*R^2 to 1e-13,
    !   * lat-lon -> Gauss "con" of 1 is 1,
    !   * sum_k dlon*w_j*v_k equals the band-metric source integral to 1e-13,
    !   * grid_cdo write/read keeps y_bnds (and the explicit-bounds flag),
    !   * a polar-cell spherical clip (projected source over the pole onto the
    !     Gauss grid, whose polar cells have two corners at the pole) returns
    !     finite weights, every source cell is covered once (sum_k phi_jk = 1)
    !     and the polar-row cells inside the source are covered once.
    ! Also: explicit bounds are validated (a centre outside its cell stops;
    ! child process, exit code checked).

    use coords

    implicit none

    character(len=*), parameter :: fldr = "maps_gauss_bounds_test"
    type(grid_class) :: gg, gr, gll, gps
    type(map_class)  :: map
    real(dp) :: lat(96), wts(96), R, dlon, emax, s_src, s_tgt, a_exp
    real(dp), allocatable :: v1(:,:), v2(:,:)
    character(len=512) :: self, arg
    integer  :: i, j, fails, rc

    call get_command_argument(0, self)
    arg = ""
    if (command_argument_count() > 0) call get_command_argument(1, arg)
    if (trim(arg) == "bad-bounds") then
        call grid_init(gr, name="bad", mtype="latlon", units="degrees", planet="Spherical Earth", &
                       x0=0.5_dp, dx=1.0_dp, nx=4, y=[0.5_dp, 1.5_dp, 2.5_dp], &
                       y_bnds=[0.0_dp, 1.0_dp, 1.2_dp, 3.0_dp])        ! centre 1.5 outside [1.0,1.2]
        stop 0
    end if

    fails = 0
    call execute_command_line("rm -rf "//fldr//" ; mkdir -p "//fldr)

    call grid_init_gaussian(gg, "gauss-t64", nlon=192, nlat=96, planet="Spherical Earth")
    call gaussian_latitudes_calc(96, lat, wts)
    R    = gg%cs%planet%a
    dlon = 360.0_dp/192.0_dp

    ! --- bounds --------------------------------------------------------------------
    if (.not. gg%G%explicit_bnds) then
        write(*,*) "FAIL: explicit_bnds not set"; fails = fails + 1
    end if
    if (gg%G%y_bnds(0) /= -90.0_dp .or. gg%G%y_bnds(96) /= 90.0_dp) then
        write(*,*) "FAIL: end bounds", gg%G%y_bnds(0), gg%G%y_bnds(96); fails = fails + 1
    end if
    do j = 1, 96
        if (gg%G%y_bnds(j) <= gg%G%y_bnds(j-1) .or. gg%G%y(j) <= gg%G%y_bnds(j-1) .or. &
            gg%G%y(j) >= gg%G%y_bnds(j)) then
            write(*,*) "FAIL: band", j, gg%G%y_bnds(j-1), gg%G%y(j), gg%G%y_bnds(j); fails = fails + 1
            exit
        end if
    end do

    ! --- areas: R^2 dlon w_j ----------------------------------------------------------
    emax = 0.0_dp
    do j = 1, 96
        a_exp = R*R*dlon*degrees_to_radians*wts(j)
        emax = max(emax, maxval(abs(gg%area(:,j) - a_exp))/a_exp)
    end do
    write(*,*) "area vs R^2*dlon*w_j: max rel err =", emax, &
               "  sum(area)/(4 pi R^2) - 1 =", sum(gg%area)/(4.0_dp*pi*R*R) - 1.0_dp
    if (emax > 1.0e-13_dp) then
        write(*,*) "FAIL: Gauss cell area /= quadrature weight"; fails = fails + 1
    end if
    if (abs(sum(gg%area)/(4.0_dp*pi*R*R) - 1.0_dp) > 1.0e-13_dp) then
        write(*,*) "FAIL: sum of Gauss cell areas /= 4 pi R^2"; fails = fails + 1
    end if

    ! --- lat-lon -> Gauss con: constants, band-metric integral ---------------------
    call grid_init(gll, name="gauss-src-ll1", mtype="latlon", units="degrees", planet="Spherical Earth", &
                   x0=0.5_dp, dx=1.0_dp, nx=360, y0=-89.5_dp, dy=1.0_dp, ny=180)
    call map_init(map, gll, gg, method="con", fldr=fldr)
    allocate(v1(360,180), v2(192,96))
    v1 = 1.0_dp
    call map_field(map, "one", v1, v2)
    write(*,*) "con of 1: max |v2 - 1| =", maxval(abs(v2 - 1.0_dp))
    if (maxval(abs(v2 - 1.0_dp)) > 1.0e-13_dp) then
        write(*,*) "FAIL: con of 1 /= 1"; fails = fails + 1
    end if
    do j = 1, 180
        do i = 1, 360
            v1(i,j) = 2.0_dp + sin(gll%lat(i,j)*degrees_to_radians*3.0_dp) + cos(gll%lon(i,j)*degrees_to_radians)
        end do
    end do
    call map_field(map, "v", v1, v2)
    s_src = 0.0_dp
    do j = 1, 180
        do i = 1, 360
            s_src = s_src + v1(i,j)*abs(gll%G%x_bnds(i) - gll%G%x_bnds(i-1))*degrees_to_radians* &
                    abs(sin(gll%G%y_bnds(j)*degrees_to_radians) - sin(gll%G%y_bnds(j-1)*degrees_to_radians))
        end do
    end do
    s_tgt = 0.0_dp
    do j = 1, 96
        s_tgt = s_tgt + sum(v2(:,j))*dlon*degrees_to_radians*wts(j)
    end do
    write(*,*) "band-metric integral source / Gauss quadrature, rel diff =", s_src, s_tgt, &
               abs(s_tgt - s_src)/abs(s_src)
    if (abs(s_tgt - s_src) > 1.0e-13_dp*abs(s_src)) then
        write(*,*) "FAIL: Gauss quadrature integral /= band-metric source integral"; fails = fails + 1
    end if

    ! --- grid_cdo round trip keeps the bounds -----------------------------------------
    call grid_cdo_write_desc_short(gg, fldr)
    call grid_cdo_read_desc(gr, "gauss-t64", fldr)
    emax = maxval(abs(gr%G%y_bnds - gg%G%y_bnds))
    write(*,*) "grid_cdo round trip: max |y_bnds diff| =", emax, " explicit:", gr%G%explicit_bnds
    if (.not. gr%G%explicit_bnds .or. emax > 1.0e-10_dp .or. &
        maxval(abs(gr%area - gg%area))/maxval(gg%area) > 1.0e-10_dp) then
        write(*,*) "FAIL: grid_cdo round trip lost the bounds"; fails = fails + 1
    end if
    if (.not. grid_signature_match(grid_signature_str(gg), grid_signature_str(gr))) then
        write(*,*) "FAIL: round-trip signature differs"; fails = fails + 1
    end if

    ! --- polar cells: projected source over the pole onto the Gauss grid ------------
    call grid_init(gps, name="gauss-src-ps", mtype="polar_stereographic", units="kilometers", &
                   planet="Spherical Earth", lon180=.false., x0=-600.0_dp, dx=40.0_dp, nx=31, &
                   y0=-600.0_dp, dy=40.0_dp, ny=31, lambda=0.0_dp, phi=-90.0_dp)
    call map_init(map, gps, gg, method="con", fldr=fldr)
    write(*,*) "polar clip: links =", map%wm%n_links, " polar-row links =", &
               count(map%wm%dst <= 192), " finite:", all(abs(map%wm%w) <= huge(1.0_dp))
    if (map%wm%n_links == 0 .or. count(map%wm%dst <= 192) == 0 .or. &
        .not. all(abs(map%wm%w) <= huge(1.0_dp)) .or. .not. all(map%wm%w == map%wm%w) .or. &
        .not. all(map%wm%wfull_dst == map%wm%wfull_dst)) then
        write(*,*) "FAIL: polar-cell clip gives no / non-finite weights"; fails = fails + 1
    end if
    ! The source (+-600 km around the pole) covers the two southern polar rows
    ! (to ~86.7 S) and lies inside the global target.
    emax = maxval(abs(map%wm%wcov_src/map%wm%wfull_src - 1.0_dp))
    a_exp = maxval(abs(map%wm%wcov_dst(1:2*192)/map%wm%wfull_dst(1:2*192) - 1.0_dp))
    write(*,*) "polar clip: max |sum_k phi_jk - 1| =", emax, " polar rows max |coverage - 1| =", a_exp
    if (emax > 1.0e-9_dp .or. a_exp > 1.0e-9_dp) then
        write(*,*) "FAIL: polar-cell clip does not cover the source / the polar rows once"; fails = fails + 1
    end if

    ! --- explicit bounds are validated (child) ---------------------------------------
    call execute_command_line(trim(self)//" bad-bounds > "//fldr//"/child.log 2>&1", exitstat=rc)
    if (rc == 0) then
        write(*,*) "FAIL: a centre outside its explicit cell was accepted"; fails = fails + 1
    end if

    call execute_command_line("rm -rf "//fldr)
    if (fails > 0) then
        write(*,*) "FAIL: test_gauss_bounds (", fails, " failure(s))"
        stop 1
    end if
    write(*,*) "PASS: test_gauss_bounds"

end program test_gauss_bounds
