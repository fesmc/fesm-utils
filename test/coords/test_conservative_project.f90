program test_conservative_project
    ! Conservative map from a global lon/lat grid onto a PROJECTED (polar-stereographic) target: the path of
    ! conservative_planar that first marks the source cells that are candidates of some target cell (arrays cidx0,
    ! cd20) and then runs the OpenMP region with the private arrays cidx, cd2. With ifx an allocatable that was
    ! allocated and freed before the region and is `private` in it was "already allocated" on entry (forrtl 151);
    ! test_map_con_proj covers the geometry but needs cdo and skips without it. Here: no cdo, two maps in a row
    ! (the second enters the path again), a constant field preserved by both, every covered target cell.
    use coords

    implicit none

    type(grid_class) :: gs1, gs2, gt
    type(map_class)  :: m1, m2
    integer, parameter :: nxt = 45, nyt = 45
    character(len=*), parameter :: fldr = "maps_con_project_test"
    real(dp), allocatable :: vs(:,:), vt(:,:)
    logical,  allocatable :: mk(:,:)
    real(dp) :: e1, e2
    integer :: fails

    call execute_command_line("rm -rf "//fldr)

    call grid_init(gs1, name="llp4", mtype="latlon", units="degrees", &
                   x0=-178.0_dp, dx=4.0_dp, nx=90, y0=-88.0_dp, dy=4.0_dp, ny=45)
    call grid_init(gs2, name="llp3", mtype="latlon", units="degrees", &
                   x0=-178.5_dp, dx=3.0_dp, nx=120, y0=-88.5_dp, dy=3.0_dp, ny=60)
    call grid_init(gt, name="NHp-128KM", mtype="polar_stereographic", units="kilometers", &
                   x0=-2816.0_dp, dx=128.0_dp, nx=nxt, y0=-2816.0_dp, dy=128.0_dp, ny=nyt, &
                   lambda=-45.0_dp, phi=70.0_dp, alpha=1.0e-3_dp, lon180=.true.)

    fails = 0
    allocate(vt(nxt,nyt), mk(nxt,nyt))

    call map_init(m1, gs1, gt, method="con", gen="coords", fldr=fldr, load=.false.)
    allocate(vs(gs1%G%nx, gs1%G%ny)); vs = 5.0_dp
    call map_field(m1, "c1", vs, vt, stat="mean", mask2=mk)
    e1 = maxval(abs(vt - 5.0_dp), mask=mk)
    deallocate(vs)

    call map_init(m2, gs2, gt, method="con", gen="coords", fldr=fldr, load=.false.)
    allocate(vs(gs2%G%nx, gs2%G%ny)); vs = 5.0_dp
    call map_field(m2, "c2", vs, vt, stat="mean", mask2=mk)
    e2 = maxval(abs(vt - 5.0_dp), mask=mk)

    if (count(mk) < nxt*nyt/2) then
        write(*,*) "FAIL: too few covered target cells:", count(mk); fails = fails + 1
    end if
    if (e1 > 1.0e-9_dp .or. e2 > 1.0e-9_dp) then
        write(*,*) "FAIL: constant field not preserved:", e1, e2; fails = fails + 1
    end if

    call execute_command_line("rm -rf "//fldr)
    if (fails > 0) error stop 1
    write(*,*) "test_conservative_project: OK (constant error", e1, e2, ", covered cells", count(mk), ")"
end program test_conservative_project
