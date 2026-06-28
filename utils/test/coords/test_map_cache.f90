program test_map_cache
    ! map_init caches each generated map to a SCRIP(-superset) file under `fldr`
    ! and, with load=.true., reloads it instead of regenerating. A reloaded map
    ! must reproduce the freshly generated one exactly, for both a distance map
    ! (MAP_DISTANCE store) and a conservative map (MAP_WEIGHT store).

    use coords

    implicit none

    type(grid_class) :: gs, gt
    type(map_class)  :: m_gen, m_load
    real(dp), allocatable :: vs(:,:), vt_gen(:,:), vt_load(:,:)
    character(len=*), parameter :: fldr = "maps_cache_test"
    integer  :: i, j, fails
    real(dp) :: emax

    call grid_init(gs, name="cache-src", mtype="cartesian", units="kilometers", &
                   x0=0.0_dp,   dx=100.0_dp, nx=10, y0=0.0_dp,   dy=100.0_dp, ny=10)
    call grid_init(gt, name="cache-tgt", mtype="cartesian", units="kilometers", &
                   x0=150.0_dp, dx=100.0_dp, nx=7,  y0=150.0_dp, dy=100.0_dp, ny=7)

    allocate(vs(10,10), vt_gen(7,7), vt_load(7,7))
    do j = 1, 10
        do i = 1, 10
            vs(i,j) = 2.0_dp + 0.003_dp*gs%x(i,j) + 0.005_dp*gs%y(i,j)
        end do
    end do

    fails = 0
    call execute_command_line("rm -rf "//fldr)

    ! --- distance map (shepard): generate+save, then load, compare ---
    call map_init(m_gen,  gs, gt, method="shepard", max_neighbors=8, fldr=fldr, load=.false.)
    call map_init(m_load, gs, gt, method="shepard", max_neighbors=8, fldr=fldr, load=.true.)
    call map_field(m_gen,  "v", vs, vt_gen,  stat="mean")
    call map_field(m_load, "v", vs, vt_load, stat="mean")
    emax = maxval(abs(vt_gen - vt_load))
    write(*,*) "distance cache round-trip max diff     =", emax
    if (emax > 1.0e-12_dp) then
        write(*,*) "FAIL: distance map cache did not round-trip"; fails = fails + 1
    end if

    ! --- conservative map: generate+save, then load, compare ---
    call map_init(m_gen,  gs, gt, method="con", fldr=fldr, load=.false.)
    call map_init(m_load, gs, gt, method="con", fldr=fldr, load=.true.)
    call map_field(m_gen,  "v", vs, vt_gen,  stat="mean")
    call map_field(m_load, "v", vs, vt_load, stat="mean")
    emax = maxval(abs(vt_gen - vt_load))
    write(*,*) "conservative cache round-trip max diff =", emax
    if (emax > 1.0e-12_dp) then
        write(*,*) "FAIL: conservative map cache did not round-trip"; fails = fails + 1
    end if

    call execute_command_line("rm -rf "//fldr)

    if (fails > 0) stop 1
    write(*,*) "PASS: test_map_cache"

end program test_map_cache
