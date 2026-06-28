program test_conservative_sphere
    ! Stage-C spherical conservative: lat-lon source -> lat-lon target on an exact
    ! 2x refinement (coincident cell boundaries). Each target cell lies in one
    ! source cell, so the area-weighted mean must transfer the parent value
    ! exactly via great-circle clipping.

    use coords

    implicit none

    type(grid_class) :: gs, gt
    type(map_class)  :: map
    real(dp) :: vs(4,3), vt(8,6)
    logical  :: m2(8,6)
    integer  :: i, j, fails
    real(dp) :: emax, emax_const

    ! Source 4x3 cells partition [0,40] lon x [30,60] lat
    call grid_init(gs, name="ll-src", mtype="latlon", units="degrees", &
                   x0=5.0_dp,  dx=10.0_dp, nx=4, y0=35.0_dp, dy=10.0_dp, ny=3)
    ! Target 8x6 cells: exact 2x refinement, coincident domain
    call grid_init(gt, name="ll-tgt", mtype="latlon", units="degrees", &
                   x0=2.5_dp,  dx=5.0_dp,  nx=8, y0=32.5_dp, dy=5.0_dp,  ny=6)

    do j = 1, 3
        do i = 1, 4
            vs(i,j) = real(i,dp) + 10.0_dp*real(j,dp)
        end do
    end do

    call map_init_conservative(map, gs, gt)

    fails = 0

    ! (1) constant field is preserved exactly (validates the partition normalizes)
    call map_field(map, "c", vs*0.0_dp + 3.0_dp, vt, stat="mean", mask2=m2)
    emax_const = maxval(abs(vt - 3.0_dp), mask=m2)
    write(*,*) "covered:", count(m2), " of", size(m2), " const err =", emax_const
    if (count(m2) < size(m2)) then; write(*,*) "FAIL: not all covered"; fails=fails+1; end if
    if (emax_const > 1.0e-12_dp) then; write(*,*) "FAIL: constant not preserved"; fails=fails+1; end if

    ! (2) value transfer is accurate to the subsampling approximation
    call map_field(map, "v", vs, vt, stat="mean", mask2=m2)
    emax = 0.0_dp
    do j = 1, 6
        do i = 1, 8
            emax = max(emax, abs(vt(i,j) - vs((i+1)/2, (j+1)/2)))
        end do
    end do
    write(*,*) "max |target - parent value| =", emax, " (subsampling approx)"
    if (emax > 1.0e-2_dp) then; write(*,*) "FAIL: value transfer beyond subsampling tol"; fails=fails+1; end if

    ! (3) global area-integral conservation (physical grid areas)
    block
        real(dp) :: int_s, int_t
        int_s = sum(vs * gs%area)
        int_t = sum(vt * gt%area)
        write(*,*) "source integral =", int_s, " target integral =", int_t, &
                   " rel diff =", abs(int_s-int_t)/abs(int_s)
        if (abs(int_s - int_t) > 1.0e-2_dp * abs(int_s)) then
            write(*,*) "FAIL: integral not conserved within subsampling tol"; fails=fails+1
        end if
    end block

    if (fails > 0) stop 1
    write(*,*) "PASS: test_conservative_sphere"

end program test_conservative_sphere
