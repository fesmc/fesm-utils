program test_gaussian_fast
    ! filter_gaussian_fast: no mask (constant field, spike) and with mask.

    use gaussian_filter, only: filter_gaussian_fast
    use precision, only: dp

    implicit none

    integer, parameter :: nx = 40, ny = 30
    real(dp), parameter :: dx = 1.0_dp, sigma = 3.0_dp
    real(dp) :: var(nx,ny), var0(nx,ny)
    logical  :: mask(nx,ny)

    ! 1. Constant field stays constant (incl. edges and corners)
    var = 5.0_dp
    call filter_gaussian_fast(var, sigma, dx)
    write(*,*) "constant: min/max = ", minval(var), maxval(var)
    if (maxval(abs(var - 5.0_dp)) > 1.d-5) then
        write(*,*) "FAIL: constant field changed"; stop 1
    end if

    ! 2. Spike in the interior is smoothed, sum conserved
    var0 = 0.0_dp
    var0(nx/2,ny/2) = 1.0_dp
    var = var0
    call filter_gaussian_fast(var, sigma, dx)
    write(*,*) "spike: max, sum = ", maxval(var), sum(var)
    if (maxval(var) > 0.5_dp .or. var(nx/2+1,ny/2) <= 0.0_dp) then
        write(*,*) "FAIL: spike not smoothed"; stop 1
    end if
    if (abs(sum(var) - 1.0_dp) > 1.d-5) then
        write(*,*) "FAIL: spike sum not conserved"; stop 1
    end if

    ! 3. Masked points are not used: valid points (=1) stay 1
    mask = .TRUE.
    mask(1:nx/2,:) = .FALSE.
    var = 1.0_dp
    where (.not. mask) var = 100.0_dp
    call filter_gaussian_fast(var, sigma, dx, mask=mask)
    write(*,*) "mask: valid min/max = ", minval(var, mask=mask), maxval(var, mask=mask)
    if (maxval(abs(var - 1.0_dp), mask=mask) > 1.d-5) then
        write(*,*) "FAIL: masked points used"; stop 1
    end if

    write(*,*) "PASS: test_gaussian_fast"

end program test_gaussian_fast
