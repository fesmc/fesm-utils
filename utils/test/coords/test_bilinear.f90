program test_bilinear
    ! Bilinear interpolation must reproduce a linear field exactly on a regular
    ! Cartesian grid. Source 10x10 grid, target interior grid, f = 2 + bx + cy.

    use coords

    implicit none

    type(grid_class) :: gsrc, gtgt
    type(map_class)  :: map
    real(dp), allocatable :: vsrc(:,:), vtgt(:,:), vexp(:,:)
    logical,  allocatable :: m2(:,:)
    real(dp), parameter :: b = 0.003_dp, c = 0.005_dp, a0 = 2.0_dp
    real(dp) :: emax
    integer  :: i, j, fails

    call grid_init(gsrc, name="src-cart", mtype="cartesian", units="kilometers", &
                   x0=0.0_dp,   dx=100.0_dp, nx=10, y0=0.0_dp,   dy=100.0_dp, ny=10)
    call grid_init(gtgt, name="tgt-cart", mtype="cartesian", units="kilometers", &
                   x0=150.0_dp, dx=100.0_dp, nx=7,  y0=150.0_dp, dy=100.0_dp, ny=7)

    allocate(vsrc(10,10), vtgt(7,7), vexp(7,7), m2(7,7))

    do j = 1, 10
        do i = 1, 10
            vsrc(i,j) = a0 + b*gsrc%x(i,j) + c*gsrc%y(i,j)
        end do
    end do
    do j = 1, 7
        do i = 1, 7
            vexp(i,j) = a0 + b*gtgt%x(i,j) + c*gtgt%y(i,j)
        end do
    end do

    call map_init(map, gsrc, gtgt, max_neighbors=8, load=.false.)
    call map_field(map, "lin", vsrc, vtgt, method="bilinear", mask2=m2)

    emax = maxval(abs(vtgt - vexp))
    write(*,*) "bilinear max abs error =", emax
    write(*,*) "all targets filled? ", all(m2)

    fails = 0
    if (.not. all(m2)) then
        write(*,*) "FAIL: not all interior targets filled"; fails = fails + 1
    end if
    if (emax > 1.0e-9_dp) then
        write(*,*) "FAIL: bilinear did not reproduce the linear field exactly"; fails = fails + 1
    end if

    if (fails > 0) stop 1
    write(*,*) "PASS: test_bilinear"

end program test_bilinear
