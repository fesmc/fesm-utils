module gaussian_latitudes
    ! Generate Gaussian grid latitudes and quadrature weights.
    !
    ! A global Gaussian grid places its latitudes at the roots of the Legendre
    ! polynomial of degree nlat (Gauss-Legendre nodes), with longitudes evenly
    ! spaced. This module computes those latitudes (in degrees) and the
    ! associated quadrature weights, so a Gaussian grid can be defined from a
    ! resolution alone rather than read from a file.

    use precision, only: dp
    use constants, only: pi, radians_to_degrees

    implicit none

    private
    public :: gaussian_latitudes_calc

contains

    subroutine gaussian_latitudes_calc(nlat, lat, wts)
        ! Return nlat Gaussian latitudes (degrees, ascending south -> north)
        ! and, optionally, the quadrature weights (which sum to 2).

        implicit none

        integer,  intent(in)  :: nlat
        real(dp), intent(out) :: lat(nlat)
        real(dp), intent(out), optional :: wts(nlat)

        real(dp) :: x(nlat), w(nlat)
        real(dp) :: z, z1, p1, p2, p3, pp
        real(dp), parameter :: eps = 1.0e-15_dp
        integer  :: i, j, m, iter
        integer, parameter :: max_iter = 100

        m = (nlat + 1) / 2   ! roots are symmetric about the equator

        do i = 1, m
            ! Initial guess for the i-th root (Newton-Raphson on P_nlat)
            z = cos(pi * (real(i,dp) - 0.25_dp) / (real(nlat,dp) + 0.5_dp))

            do iter = 1, max_iter
                p1 = 1.0_dp
                p2 = 0.0_dp
                ! Recurrence for the Legendre polynomial at z
                do j = 1, nlat
                    p3 = p2
                    p2 = p1
                    p1 = ((2.0_dp*real(j,dp) - 1.0_dp)*z*p2 - (real(j,dp) - 1.0_dp)*p3) / real(j,dp)
                end do
                ! Derivative pp via the standard relation
                pp = real(nlat,dp) * (z*p1 - p2) / (z*z - 1.0_dp)
                z1 = z
                z  = z1 - p1/pp
                if (abs(z - z1) <= eps) exit
            end do

            ! Symmetric nodes (x ascending: -1 .. +1) and weights
            x(i)            = -z
            x(nlat + 1 - i) =  z
            w(i)            = 2.0_dp / ((1.0_dp - z*z) * pp*pp)
            w(nlat + 1 - i) = w(i)
        end do

        ! Convert nodes to latitude in degrees (south -> north)
        do i = 1, nlat
            lat(i) = asin(x(i)) * radians_to_degrees
        end do

        if (present(wts)) wts = w

        return

    end subroutine gaussian_latitudes_calc

end module gaussian_latitudes
