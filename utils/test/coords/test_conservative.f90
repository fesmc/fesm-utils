program test_conservative
    ! Planar conservative remapping on an exact 2x grid refinement (coincident
    ! domains). Each target cell lies in exactly one source cell, so the area-
    ! weighted mean must transfer the parent value exactly and conserve the
    ! area integral.

    use coords

    implicit none

    type(grid_class) :: gs, gt
    type(map_class)  :: map
    real(dp) :: vs(4,4), vt(8,8)
    logical  :: m2(8,8)
    integer  :: i, j, fails
    real(dp) :: int_s, int_t, emax
    real(dp), parameter :: As = 100.0_dp*100.0_dp, At = 50.0_dp*50.0_dp

    ! Source 4x4 cells partition [0,400]^2 ; target 8x8 cells partition the same.
    call grid_init(gs, name="src", mtype="cartesian", units="kilometers", &
                   x0=50.0_dp, dx=100.0_dp, nx=4, y0=50.0_dp, dy=100.0_dp, ny=4)
    call grid_init(gt, name="tgt", mtype="cartesian", units="kilometers", &
                   x0=25.0_dp, dx=50.0_dp,  nx=8, y0=25.0_dp, dy=50.0_dp,  ny=8)

    do j = 1, 4
        do i = 1, 4
            vs(i,j) = real(i, dp) + 10.0_dp*real(j, dp)
        end do
    end do

    call map_init_conservative(map, gs, gt)
    call map_field(map, "v", vs, vt, method="mean", mask2=m2)

    fails = 0
    if (.not. all(m2)) then
        write(*,*) "FAIL: not all target cells covered"; fails = fails + 1
    end if

    ! exact parent-value transfer (target (it,jt) lies in source ((it+1)/2,(jt+1)/2))
    emax = 0.0_dp
    do j = 1, 8
        do i = 1, 8
            emax = max(emax, abs(vt(i,j) - vs((i+1)/2, (j+1)/2)))
        end do
    end do
    write(*,*) "max |target - parent source value| =", emax
    if (emax > 1.0e-9_dp) then
        write(*,*) "FAIL: conservative did not transfer parent values exactly"; fails = fails + 1
    end if

    ! area-integral conservation
    int_s = sum(vs) * As
    int_t = sum(vt) * At
    write(*,*) "source integral =", int_s, " target integral =", int_t
    if (abs(int_s - int_t) > 1.0e-6_dp * abs(int_s)) then
        write(*,*) "FAIL: area integral not conserved"; fails = fails + 1
    end if

    if (fails > 0) stop 1
    write(*,*) "PASS: test_conservative"

end program test_conservative
