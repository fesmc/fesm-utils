program test_grid_cf_read
    ! grid_cdo_read_desc must reconstruct a grid from a cdo-native description
    ! file (CF projection keys, no fesm-utils '#' header) just as well as from a
    ! full fesm-utils file. For each reference grid we:
    !   1. write it with grid_cdo_write_desc_short (# header + cdo/CF body),
    !   2. read it back (uses the # header) and check it round-trips,
    !   3. strip the # header (-> pure CF file) and read again (CF fallback path),
    !      checking it still reproduces the reference grid.
    !
    ! Reference grids use lon180=.false. so they match the CF-read default (the
    ! CF description carries no lon180 flag).

    use coords

    implicit none

    type(grid_class) :: g_ps, g_st
    integer  :: fails
    character(len=*), parameter :: fldr = "maps_cf_test"

    fails = 0
    call execute_command_line("rm -rf "//fldr//" ; mkdir -p "//fldr)

    ! --- polar_stereographic (ANT-like) ---
    call grid_init(g_ps, name="cf-ps", mtype="polar_stereographic", units="kilometers", &
                   lon180=.false., x0=-3040.0_dp, dx=16.0_dp, nx=381, &
                   y0=-3040.0_dp, dy=16.0_dp, ny=381, lambda=0.0_dp, phi=-71.0_dp)
    call roundtrip(g_ps, "cf-ps", fldr, fails)

    ! --- (oblique) stereographic (GRL-like), exercises the stereographic keys ---
    call grid_init(g_st, name="cf-st", mtype="stereographic", units="kilometers", &
                   lon180=.false., x0=-720.0_dp, dx=40.0_dp, nx=37, &
                   y0=-3450.0_dp, dy=40.0_dp, ny=61, lambda=-39.0_dp, phi=90.0_dp, alpha=71.0_dp)
    call roundtrip(g_st, "cf-st", fldr, fails)

    call execute_command_line("rm -rf "//fldr)

    if (fails > 0) stop 1
    write(*,*) "PASS: test_grid_cf_read"

contains

    subroutine roundtrip(gref, name, fldr, fails)
        type(grid_class), intent(in)    :: gref
        character(len=*), intent(in)    :: name, fldr
        integer,          intent(inout) :: fails
        type(grid_class) :: g_full, g_cf
        character(len=512) :: f_full, f_cf

        f_full = trim(fldr)//"/grid_"//trim(name)//".txt"
        f_cf   = trim(fldr)//"/grid_"//trim(name)//"-cf.txt"

        ! Write full fesm-utils file, then read it back (# header path).
        call grid_cdo_write_desc_short(gref, fldr)
        call grid_cdo_read_desc(g_full, trim(name), fldr)
        call check(gref, g_full, trim(name)//":header", fails)

        ! Strip the # header to a cdo-native file, then read it (CF fallback path).
        call execute_command_line("grep -v '^#' "//trim(f_full)//" > "//trim(f_cf))
        call grid_cdo_read_desc(g_cf, trim(name)//"-cf", fldr)
        call check(gref, g_cf, trim(name)//":cf", fails)
    end subroutine roundtrip

    subroutine check(ga, gb, tag, fails)
        type(grid_class), intent(in)    :: ga, gb
        character(len=*), intent(in)    :: tag
        integer,          intent(inout) :: fails
        real(dp) :: ex, ey, elon, elat

        if (gb%G%nx /= ga%G%nx .or. gb%G%ny /= ga%G%ny) then
            write(*,*) "FAIL ["//tag//"]: dims ", gb%G%nx, gb%G%ny, " expected ", ga%G%nx, ga%G%ny
            fails = fails + 1
            return
        end if
        if (trim(gb%cs%mtype) /= trim(ga%cs%mtype)) then
            write(*,*) "FAIL ["//tag//"]: mtype '"//trim(gb%cs%mtype)//"' expected '"//trim(ga%cs%mtype)//"'"
            fails = fails + 1
        end if

        ex   = maxval(abs(gb%x   - ga%x))
        ey   = maxval(abs(gb%y   - ga%y))
        elon = maxval(abs(gb%lon - ga%lon))
        elat = maxval(abs(gb%lat - ga%lat))
        write(*,"(a,4es12.3)") " "//tag//": max|dx,dy,dlon,dlat| =", ex, ey, elon, elat
        if (ex > 1.0e-6_dp .or. ey > 1.0e-6_dp .or. elon > 1.0e-4_dp .or. elat > 1.0e-4_dp) then
            write(*,*) "FAIL ["//tag//"]: grid coordinates disagree"
            fails = fails + 1
        end if
    end subroutine check

end program test_grid_cf_read
