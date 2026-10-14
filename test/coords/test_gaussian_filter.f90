program test_gaussian_filter
    ! gaussian_filter: filter_gaussian_fast (no mask, mask), filter_gaussian
    ! without dx, hgrad boundaries on a non-square grid.

    use gaussian_filter, only: filter_gaussian, filter_gaussian_fast, hgrad
    use precision, only: dp

    implicit none

    integer, parameter :: nx = 40, ny = 30
    real(dp), parameter :: dx = 1.0_dp, sigma = 3.0_dp
    real(dp) :: var(nx,ny), var0(nx,ny), var_dx(nx,ny)
    real(4)  :: u(nx,ny), du(nx,ny)
    logical  :: mask(nx,ny)
    integer  :: i

    ! 1. filter_gaussian_fast: constant field stays constant (incl. edges and corners)
    var = 5.0_dp
    call filter_gaussian_fast(var, sigma, dx)
    write(*,*) "fast constant: min/max = ", minval(var), maxval(var)
    if (maxval(abs(var - 5.0_dp)) > 1.d-5) then
        write(*,*) "FAIL: constant field changed"; stop 1
    end if

    ! 2. filter_gaussian_fast: spike in the interior is smoothed, sum conserved
    var0 = 0.0_dp
    var0(nx/2,ny/2) = 1.0_dp
    var = var0
    call filter_gaussian_fast(var, sigma, dx)
    write(*,*) "fast spike: max, sum = ", maxval(var), sum(var)
    if (maxval(var) > 0.5_dp .or. var(nx/2+1,ny/2) <= 0.0_dp) then
        write(*,*) "FAIL: spike not smoothed"; stop 1
    end if
    if (abs(sum(var) - 1.0_dp) > 1.d-5) then
        write(*,*) "FAIL: spike sum not conserved"; stop 1
    end if

    ! 3. filter_gaussian_fast with mask: masked points are not used and stay unchanged
    mask = .TRUE.
    mask(1:nx/2,:) = .FALSE.
    var = 1.0_dp
    where (.not. mask) var = 100.0_dp
    call filter_gaussian_fast(var, sigma, dx, mask=mask)
    write(*,*) "fast mask: valid min/max = ", minval(var, mask=mask), maxval(var, mask=mask)
    write(*,*) "fast mask: masked min/max = ", minval(var, mask=.not. mask), maxval(var, mask=.not. mask)
    if (maxval(abs(var - 1.0_dp), mask=mask) > 1.d-5) then
        write(*,*) "FAIL: masked points used"; stop 1
    end if
    if (any(var /= 100.0_dp .and. .not. mask)) then
        write(*,*) "FAIL: masked points changed"; stop 1
    end if

    ! 4. filter_gaussian without dx: sigma in grid points, same as dx=1
    var = var0
    call filter_gaussian(var, sigma)
    var_dx = var0
    call filter_gaussian(var_dx, sigma, dx=dx)
    write(*,*) "filter_gaussian spike: max, max diff (dx absent/present) = ", &
               maxval(var), maxval(abs(var - var_dx))
    if (maxval(var) > 0.5_dp .or. maxval(abs(var - var_dx)) > 1.d-6) then
        write(*,*) "FAIL: filter_gaussian without dx"; stop 1
    end if

    ! 5. hgrad of u = x on a non-square grid is 1 everywhere, incl. boundaries
    do i = 1, nx
        u(i,:) = real(i)
    end do
    du = hgrad(u, 1.0)
    write(*,*) "hgrad(x): min/max = ", minval(du), maxval(du)
    if (maxval(abs(du - 1.0)) > 1.e-5) then
        write(*,*) "FAIL: hgrad boundary"; stop 1
    end if

    write(*,*) "PASS: test_gaussian_filter"

end program test_gaussian_filter
