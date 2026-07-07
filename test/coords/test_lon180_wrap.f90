program test_lon180_wrap
    ! Regression: a regular lon/lat source built from 0..360 centers with
    ! lon180=.true. stores its longitude centers remapped in place (0..175,
    ! -180..-5) -- non-monotonic across the seam. The structured nn/bilinear
    ! locate and the separable conservative edges both require a monotonic
    ! longitude axis; without the unwrap they mislocate every (bil/nn) or corrupt
    ! the seam column (con), which silently wrecked e.g. the VILMA Gauss-grid
    ! coupling maps in climber-x.
    !
    ! Invariant (no cdo needed): the SAME physical grid described in 0..360
    ! (lon180=.false.) and in -180..180 (lon180=.true.) must produce an identical
    ! map. Checks nn, bilinear and conservative.

    use coords

    implicit none

    type(grid_class) :: g360, g180, gt
    integer, parameter :: nxs=72, nys=36, nxt=72, nyt=36
    real(dp), parameter :: d2r = 3.14159265358979_dp/180.0_dp
    real(dp), parameter :: rtol = 1.0e-9_dp
    real(dp) :: lons(nxs), lats(nys)
    real(dp) :: vs(nxs,nys)
    integer  :: i, j, fails

    ! source longitude centers in 0..360 (as VILMA/Gaussian grids are stored)
    do i = 1, nxs
        lons(i) = 2.5_dp + real(i-1,dp)*5.0_dp     ! 2.5 .. 357.5
    end do
    do j = 1, nys
        lats(j) = -87.5_dp + real(j-1,dp)*5.0_dp
    end do

    ! Same centers, two conventions. g180 remaps lon>180 by -360 in place -> its
    ! stored G%x is non-monotonic across the seam.
    call grid_init(g360, name="src360", mtype="latlon", units="degrees", &
                   x=lons, y=lats, lon180=.false.)
    call grid_init(g180, name="src180", mtype="latlon", units="degrees", &
                   x=lons, y=lats, lon180=.true.)

    ! target: a -180..180 regular grid overlapping the source
    call grid_init(gt, name="tgt180", mtype="latlon", units="degrees", &
                   x0=-177.5_dp, dx=5.0_dp, nx=nxt, y0=-87.5_dp, dy=5.0_dp, ny=nyt)

    ! 360-periodic-in-lon source field (identical array for both conventions)
    do j = 1, nys
        do i = 1, nxs
            vs(i,j) = 1000.0_dp*cos(3.0_dp*lons(i)*d2r)*cos(5.0_dp*lats(j)*d2r)
        end do
    end do

    call execute_command_line("rm -rf maps_lon180")

    fails = 0
    call check_same("nn",  g360, g180, gt, vs, fails)
    call check_same("bil", g360, g180, gt, vs, fails)
    call check_same("con", g360, g180, gt, vs, fails)

    call execute_command_line("rm -rf maps_lon180")

    if (fails > 0) then
        write(*,*) "test_lon180_wrap FAILED with", fails, "errors"; stop 1
    end if
    write(*,*) "PASS: test_lon180_wrap"

contains

    subroutine check_same(method, ga, gb, gt, vs, fails)
        character(len=*), intent(in)    :: method
        type(grid_class), intent(in)    :: ga, gb, gt
        real(dp),         intent(in)    :: vs(:,:)
        integer,          intent(inout) :: fails
        type(map_class) :: ma, mb
        real(dp) :: va(gt%G%nx,gt%G%ny), vb(gt%G%nx,gt%G%ny)
        logical  :: mka(gt%G%nx,gt%G%ny), mkb(gt%G%nx,gt%G%ny)
        real(dp) :: dmax
        integer  :: nmiss

        call map_init(ma, ga, gt, method=method, gen="coords", fldr="maps_lon180", load=.false.)
        call map_init(mb, gb, gt, method=method, gen="coords", fldr="maps_lon180", load=.false.)
        call map_field(ma, "v", vs, va, stat="mean", mask2=mka)
        call map_field(mb, "v", vs, vb, stat="mean", mask2=mkb)

        if (count(mka .neqv. mkb) /= 0) then
            write(*,*) "  FAIL "//method//": coverage differs between conventions"; fails = fails + 1
        end if
        dmax  = maxval(abs(va - vb), mask=(mka .and. mkb))
        nmiss = count(mka)
        write(*,"(a,a,a,es10.2,a,i5)") "  ", method, ": max|v360 - v180| =", dmax, "  covered=", nmiss
        if (dmax > rtol) then
            write(*,*) "  FAIL "//method//": lon180 convention changes the map"; fails = fails + 1
        end if
        if (nmiss == 0) then
            write(*,*) "  FAIL "//method//": no target cells covered"; fails = fails + 1
        end if
    end subroutine check_same

end program test_lon180_wrap
