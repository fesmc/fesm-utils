program test_pole_fill
    ! Regression test for the coords bilinear generator at the pole.
    !
    ! A bilinear map from a lat-lon source spanning only +/-87.5 (like CMN-5x5)
    ! to a target whose rows reach beyond that span (as a polar-stereographic
    ! ice grid does at the pole) must cover EVERY target cell. Previously,
    ! map_init_structured dropped targets poleward of the last source row,
    ! leaving a hole that fed uninitialised garbage into the model and crashed
    ! it with atmosphere instabilities. locate_cell now clamps out-of-span
    ! targets onto the nearest source boundary cell (nearest-neighbour edge
    ! fill, matching cdo genbil), so the hole is gone.
    use coords
    implicit none
    type(grid_class) :: gsrc, gtgt
    type(map_class)  :: m
    character(len=*), parameter :: fldr = "maps_pole_test"
    integer :: ndst, k, uncovered, i, j, nbeyond, fails
    logical, allocatable :: covered(:)

    call execute_command_line("rm -rf "//fldr)
    call execute_command_line("mkdir -p "//fldr)
    fails = 0

    ! Source: lat centers -87.5..87.5 (ny=36), lon -177.5..177.5 (nx=72)
    call grid_init(gsrc, name="src5", mtype="latlon", units="degrees", planet="WGS84", &
                   x0=-177.5_dp, dx=5.0_dp, nx=72, y0=-87.5_dp, dy=5.0_dp, ny=36)

    ! Target: lat centers -89..89 (ny=90) -> first and last rows lie OUTSIDE
    ! the source span, the same situation as an ice grid's polar cap.
    call grid_init(gtgt, name="tgt2", mtype="latlon", units="degrees", planet="WGS84", &
                   x0=-179.0_dp, dx=2.0_dp, nx=180, y0=-89.0_dp, dy=2.0_dp, ny=90)

    call map_init(m, gsrc, gtgt, method="bil", gen="coords", &
                  fldr=fldr, load=.false., clean=.false.)

    ndst = m%wm%n_dst
    allocate(covered(ndst)); covered = .false.
    do k = 1, m%wm%n_links
        covered(m%wm%dst(k)) = .true.
    end do
    uncovered = count(.not. covered)

    ! Number of target cells beyond the source span (the two polar rows).
    nbeyond = 0
    do j = 1, 90
        do i = 1, 180
            if (abs(gtgt%lat(i,j)) > 87.5_dp) nbeyond = nbeyond + 1
        end do
    end do

    write(*,"(a,i0)") " target cells (n_dst)        = ", ndst
    write(*,"(a,i0)") " target cells beyond +/-87.5 = ", nbeyond
    write(*,"(a,i0)") " UNCOVERED target cells      = ", uncovered

    if (nbeyond == 0) then
        write(*,*) "FAIL: test grids do not exercise the beyond-span (pole) path"
        fails = fails + 1
    end if
    if (uncovered /= 0) then
        write(*,"(a,i0,a)") " FAIL: ", uncovered, " target cells left unmapped (pole hole)"
        fails = fails + 1
    end if

    call execute_command_line("rm -rf "//fldr)

    write(*,*)
    if (fails == 0) then
        write(*,*) "PASS: test_pole_fill (every target cell mapped, incl. polar cap)"
    else
        write(*,"(a,i0,a)") " FAIL: test_pole_fill (", fails, " check(s) failed)"
        stop 1
    end if
end program test_pole_fill
