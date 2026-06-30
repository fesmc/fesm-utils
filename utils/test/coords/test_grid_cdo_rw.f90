program test_grid_cdo_rw
    ! Round-trip check for the cdo grid description writer/reader:
    !   grid_cdo_write_desc_short  ->  grid_<name>.txt  ->  grid_cdo_read_desc
    ! A grid reconstructed from its description file must match the original
    ! (axes + coordinate system) for every grid family the writer emits:
    ! latlon, polar_stereographic, oblique stereographic, cartesian, gaussian.

    use coords

    implicit none

    type(grid_class) :: g, g2
    real(dp), allocatable :: glat(:), gwts(:)
    character(len=*), parameter :: fldr = "maps_cdo_rw_test"
    integer :: fails

    call execute_command_line("rm -rf "//fldr)
    call execute_command_line("mkdir -p "//fldr)

    fails = 0

    ! --- latlon (global 5deg) ---
    call grid_init(g, name="ll5", mtype="latlon", units="degrees", planet="Spherical Earth", &
                   x0=0.0_dp, dx=5.0_dp, nx=72, y0=-90.0_dp, dy=5.0_dp, ny=37)
    call roundtrip(g, g2, fldr, "latlon", fails)

    ! --- polar_stereographic (Antarctic 64km) ---
    call grid_init(g, name="ANT-64KM", mtype="polar_stereographic", units="kilometers", planet="WGS84", &
                   x0=-3040.0_dp, dx=64.0_dp, nx=96, y0=-3040.0_dp, dy=64.0_dp, ny=96, &
                   lambda=0.0_dp, phi=-71.0_dp, alpha=19.0_dp)
    call roundtrip(g, g2, fldr, "polar_stereographic", fails)

    ! --- oblique stereographic (Greenland-like, non-trivial alpha) ---
    call grid_init(g, name="GRL-32KM", mtype="stereographic", units="kilometers", planet="WGS84", &
                   x0=-800.0_dp, dx=32.0_dp, nx=50, y0=-3400.0_dp, dy=32.0_dp, ny=90, &
                   lambda=-39.0_dp, phi=72.0_dp, alpha=8.4_dp)
    call roundtrip(g, g2, fldr, "stereographic", fails)

    ! --- cartesian ---
    call grid_init(g, name="cart", mtype="cartesian", units="kilometers", planet="WGS84", &
                   x0=0.0_dp, dx=10.0_dp, nx=20, y0=0.0_dp, dy=10.0_dp, ny=15)
    call roundtrip(g, g2, fldr, "cartesian", fails)

    ! --- gaussian (explicit latitudes via yvals) ---
    allocate(glat(48), gwts(48))
    call gaussian_latitudes_calc(48, glat, gwts)
    call grid_init(g, name="gauss48", mtype="gaussian", units="degrees", planet="Spherical Earth", &
                   x0=0.0_dp, dx=3.75_dp, nx=96, y=glat)
    call roundtrip(g, g2, fldr, "gaussian", fails)

    call execute_command_line("rm -rf "//fldr)

    write(*,*)
    if (fails == 0) then
        write(*,*) "PASS: test_grid_cdo_rw (all grid families round-trip)"
    else
        write(*,"(a,i0,a)") " FAIL: test_grid_cdo_rw (", fails, " mismatch(es))"
        stop 1
    end if

contains

    subroutine roundtrip(g0, gr, fldr, label, fails)
        type(grid_class), intent(in)    :: g0
        type(grid_class), intent(inout) :: gr
        character(len=*), intent(in)    :: fldr, label
        integer,          intent(inout) :: fails

        real(dp), parameter :: tol = 1.0e-8_dp
        integer :: nbad

        call grid_cdo_write_desc_short(g0, fldr)
        call grid_cdo_read_desc(gr, trim(g0%name), fldr)

        nbad = 0

        ! Axes
        call eq_i(g0%G%nx, gr%G%nx, "nx", nbad)
        call eq_i(g0%G%ny, gr%G%ny, "ny", nbad)
        call eq_r(g0%G%dx, gr%G%dx, tol, "dx", nbad)
        call eq_r(g0%G%dy, gr%G%dy, tol, "dy", nbad)
        call eq_vec(g0%G%x, gr%G%x, tol, "x-axis", nbad)
        call eq_vec(g0%G%y, gr%G%y, tol, "y-axis", nbad)

        ! Coordinate system
        if (trim(g0%cs%mtype) /= trim(gr%cs%mtype)) then
            write(*,*) "  mismatch mtype: '"//trim(g0%cs%mtype)//"' vs '"//trim(gr%cs%mtype)//"'"; nbad = nbad + 1
        end if
        if (trim(g0%cs%units) /= trim(gr%cs%units)) then
            write(*,*) "  mismatch units: '"//trim(g0%cs%units)//"' vs '"//trim(gr%cs%units)//"'"; nbad = nbad + 1
        end if
        if (g0%cs%is_cartesian  .neqv. gr%cs%is_cartesian)  then; write(*,*) "  mismatch is_cartesian";  nbad = nbad + 1; end if
        if (g0%cs%is_projection .neqv. gr%cs%is_projection) then; write(*,*) "  mismatch is_projection"; nbad = nbad + 1; end if
        if (g0%cs%is_lon180     .neqv. gr%cs%is_lon180)     then; write(*,*) "  mismatch is_lon180";     nbad = nbad + 1; end if
        call eq_r(g0%cs%planet%a, gr%cs%planet%a, tol, "planet%a", nbad)
        call eq_r(g0%cs%planet%f, gr%cs%planet%f, tol, "planet%f", nbad)
        call eq_r(g0%cs%xy_conv,  gr%cs%xy_conv,  tol, "xy_conv",  nbad)

        ! Projection parameters (only meaningful for projected grids)
        if (g0%cs%is_projection) then
            call eq_r(g0%cs%proj%lambda, gr%cs%proj%lambda, tol, "proj%lambda", nbad)
            call eq_r(g0%cs%proj%phi,    gr%cs%proj%phi,    tol, "proj%phi",    nbad)
            call eq_r(g0%cs%proj%alpha,  gr%cs%proj%alpha,  tol, "proj%alpha",  nbad)
        end if

        ! Derived lon/lat fields (the strongest end-to-end check). Skipped for a
        ! pure cartesian grid, which has no geographic lon/lat (left undefined).
        if (g0%cs%is_projection .or. .not. g0%cs%is_cartesian) then
            call eq_mat(g0%lon, gr%lon, 1.0e-6_dp, "lon", nbad)
            call eq_mat(g0%lat, gr%lat, 1.0e-6_dp, "lat", nbad)
        end if

        if (nbad == 0) then
            write(*,"(a)") " ok: "//label
        else
            write(*,"(a,i0,a)") " BAD: "//label//" (", nbad, " field mismatch(es))"
            fails = fails + 1
        end if

    end subroutine roundtrip

    subroutine eq_i(a, b, name, nbad)
        integer, intent(in) :: a, b
        character(len=*), intent(in) :: name
        integer, intent(inout) :: nbad
        if (a /= b) then
            write(*,"(a,i0,a,i0)") "  mismatch "//name//": ", a, " vs ", b
            nbad = nbad + 1
        end if
    end subroutine eq_i

    subroutine eq_r(a, b, tol, name, nbad)
        real(dp), intent(in) :: a, b, tol
        character(len=*), intent(in) :: name
        integer, intent(inout) :: nbad
        if (abs(a - b) > tol) then
            write(*,"(a,es14.6,a,es14.6)") "  mismatch "//name//": ", a, " vs ", b
            nbad = nbad + 1
        end if
    end subroutine eq_r

    subroutine eq_vec(a, b, tol, name, nbad)
        real(dp), intent(in) :: a(:), b(:)
        real(dp), intent(in) :: tol
        character(len=*), intent(in) :: name
        integer, intent(inout) :: nbad
        if (size(a) /= size(b)) then
            write(*,*) "  mismatch "//name//" size"; nbad = nbad + 1; return
        end if
        if (maxval(abs(a - b)) > tol) then
            write(*,"(a,es14.6)") "  mismatch "//name//" max diff = ", maxval(abs(a - b))
            nbad = nbad + 1
        end if
    end subroutine eq_vec

    subroutine eq_mat(a, b, tol, name, nbad)
        real(dp), intent(in) :: a(:,:), b(:,:)
        real(dp), intent(in) :: tol
        character(len=*), intent(in) :: name
        integer, intent(inout) :: nbad
        if (size(a,1) /= size(b,1) .or. size(a,2) /= size(b,2)) then
            write(*,*) "  mismatch "//name//" shape"; nbad = nbad + 1; return
        end if
        if (maxval(abs(a - b)) > tol) then
            write(*,"(a,es14.6)") "  mismatch "//name//" max diff = ", maxval(abs(a - b))
            nbad = nbad + 1
        end if
    end subroutine eq_mat

end program test_grid_cdo_rw
