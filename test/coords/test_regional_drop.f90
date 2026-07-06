program test_regional_drop
    ! Regression test for the coords bilinear/nn out-of-span clamp bound.
    !
    ! locate_cell clamps a target just beyond the source span (within one
    ! boundary cell) onto the nearest source edge -- this fills a small polar cap
    ! (see test_pole_fill). But a target FARTHER than one cell outside the span is
    ! off the source domain and must be left UNMAPPED, matching cdo genbil.
    !
    ! Without that bound, a regional source (e.g. an NH-only ice grid, ~40-90N)
    ! mapped onto a global target clamped every southern/Antarctic target onto the
    ! NH source edge. In climber-x's geo module that fed a spurious "no ice" mask
    ! onto Antarctica (mask_ice_geo=0 -> geo%hires%h_ice zeroed), deleting the
    ! Antarctic ice sheet and blowing up the ocean with CFL instabilities.
    use coords
    implicit none
    type(grid_class) :: gsrc, gtgt
    type(map_class)  :: m
    integer :: k, i, j, n_ant, n_ant_cov, n_edge, n_edge_cov, fails
    logical, allocatable :: covered(:)
    real(dp) :: la
    character(len=*), parameter :: fldr = "maps_regional_test"

    call execute_command_line("rm -rf "//fldr//"; mkdir -p "//fldr)
    fails = 0

    ! Regional NH source: lat centers 42..88 (dy=4, ny=12), global lon.
    call grid_init(gsrc, name="nhsrc", mtype="latlon", units="degrees", planet="WGS84", &
                   x0=-177.5_dp, dx=5.0_dp, nx=72, y0=42.0_dp, dy=4.0_dp, ny=12)
    ! Global 1-deg target.
    call grid_init(gtgt, name="glob1", mtype="latlon", units="degrees", planet="WGS84", &
                   x0=-179.5_dp, dx=1.0_dp, nx=360, y0=-89.5_dp, dy=1.0_dp, ny=180)

    call map_init(m, gsrc, gtgt, method="bil", gen="coords", &
                  fldr=fldr, load=.false., clean=.false.)

    allocate(covered(m%wm%n_dst)); covered = .false.
    do k = 1, m%wm%n_links
        covered(m%wm%dst(k)) = .true.
    end do

    n_ant = 0; n_ant_cov = 0; n_edge = 0; n_edge_cov = 0
    do j = 1, 180
        do i = 1, 360
            la = gtgt%lat(i,j); k = (j-1)*360 + i
            if (la < -60.0_dp) then
                n_ant = n_ant + 1
                if (covered(k)) n_ant_cov = n_ant_cov + 1
            end if
            ! Just beyond the last source row (88N): within one 4-deg cell.
            if (la > 88.0_dp .and. la < 91.0_dp) then
                n_edge = n_edge + 1
                if (covered(k)) n_edge_cov = n_edge_cov + 1
            end if
        end do
    end do

    write(*,"(a,i0,a,i0)") " Antarctic (lat<-60) targets covered  = ", n_ant_cov, " / ", n_ant
    write(*,"(a,i0,a,i0)") " near-edge (88<lat<91) targets covered = ", n_edge_cov, " / ", n_edge

    if (n_ant_cov /= 0) then
        write(*,"(a,i0,a)") " FAIL: ", n_ant_cov, " Antarctic target(s) clamped onto the NH source edge"
        fails = fails + 1
    end if
    if (n_edge_cov /= n_edge) then
        write(*,"(a,i0,a,i0,a)") " FAIL: only ", n_edge_cov, " / ", n_edge, " near-edge overshoot targets filled"
        fails = fails + 1
    end if

    call execute_command_line("rm -rf "//fldr)

    write(*,*)
    if (fails == 0) then
        write(*,*) "PASS: test_regional_drop (far targets dropped, near-edge overshoot filled)"
    else
        write(*,"(a,i0,a)") " FAIL: test_regional_drop (", fails, " check(s) failed)"
        stop 1
    end if
end program test_regional_drop
