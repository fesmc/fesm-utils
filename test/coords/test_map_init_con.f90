program test_map_init_con
    ! Exercise the unified conservative seam: map_init(..., method="con") must
    ! select the in-package analytic generator by default (gen="coords") and
    ! produce exactly the same map as the direct map_init_conservative call.
    ! Uses the same exact 2x refinement as test_conservative (each target cell
    ! lies in one source cell, so the area-weighted mean transfers the parent
    ! value exactly).

    use coords

    implicit none

    type(grid_class) :: gs, gt
    type(map_class)  :: map_seam, map_direct
    real(dp) :: vs(4,4), vt(8,8), vt2(8,8)
    integer  :: i, j, fails
    real(dp) :: emax

    call grid_init(gs, name="src", mtype="cartesian", units="kilometers", &
                   x0=50.0_dp, dx=100.0_dp, nx=4, y0=50.0_dp, dy=100.0_dp, ny=4)
    call grid_init(gt, name="tgt", mtype="cartesian", units="kilometers", &
                   x0=25.0_dp, dx=50.0_dp,  nx=8, y0=25.0_dp, dy=50.0_dp,  ny=8)

    do j = 1, 4
        do i = 1, 4
            vs(i,j) = real(i, dp) + 10.0_dp*real(j, dp)
        end do
    end do

    ! New generic seam (gen defaults to "coords")
    call map_init(map_seam, gs, gt, method="con", load=.false.)
    call map_field(map_seam, "v", vs, vt, stat="mean")

    ! Direct analytic call for comparison
    call map_init_conservative(map_direct, gs, gt)
    call map_field(map_direct, "v", vs, vt2, stat="mean")

    fails = 0

    ! Seam must reproduce the direct analytic result bit-for-bit
    if (maxval(abs(vt - vt2)) > 0.0_dp) then
        write(*,*) "FAIL: map_init(method='con') differs from map_init_conservative"
        fails = fails + 1
    end if

    ! And it must transfer parent values exactly
    emax = 0.0_dp
    do j = 1, 8
        do i = 1, 8
            emax = max(emax, abs(vt(i,j) - vs((i+1)/2, (j+1)/2)))
        end do
    end do
    write(*,*) "max |target - parent source value| =", emax
    if (emax > 1.0e-9_dp) then
        write(*,*) "FAIL: conservative seam did not transfer parent values exactly"
        fails = fails + 1
    end if

    ! Explicit gen="coords" must match the default
    call map_init(map_seam, gs, gt, method="con", gen="coords", load=.false.)
    call map_field(map_seam, "v", vs, vt2, stat="mean")
    if (maxval(abs(vt - vt2)) > 0.0_dp) then
        write(*,*) "FAIL: gen='coords' differs from the default"
        fails = fails + 1
    end if

    if (fails > 0) stop 1
    write(*,*) "PASS: test_map_init_con"

end program test_map_init_con
