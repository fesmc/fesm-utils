module gaussian_quadrature

    use, intrinsic :: iso_fortran_env, only : input_unit, output_unit, error_unit    
    use precision

    implicit none
    

contains

    subroutine get_gauss_interp_points(f, u, i, j, grid)
        
        implicit none
        
        real(wp), intent(OUT) :: f(4)       ! Values at 4 Gaussian quadrature points
        real(wp), intent(in)  :: u(:,:)     ! Variable to be interpolated
        integer,  intent(in)  :: i          ! x-index of current cell
        integer,  intent(in)  :: j          ! y-index of current cell
        character(len=*), intent(IN) :: grid    ! "aa", "ab", "acx", "acy"

        ! Local variables
        integer :: nx, ny, q

        ! Quadrature point locations in reference element (-1,1)
        real(kind=8), parameter :: sqrt3_inv = 1.0d0 / sqrt(3.0d0)
        real(kind=8) :: xi(4), eta(4)

        real(kind=8) :: N(4,4)              ! Shape functions N at Gaussian quadrature points
        real(kind=8) :: u1, u2, u3, u4      ! Values of u at the four corners of the cell

        nx = size(u,1)
        ny = size(u,2)

        ! Define quadrature points
        xi  = [ -sqrt3_inv, sqrt3_inv, sqrt3_inv, -sqrt3_inv ]
        eta = [ -sqrt3_inv, -sqrt3_inv, sqrt3_inv, sqrt3_inv ]

        ! Define shape functions N(q,corner)
        do q = 1, 4
            N(q,1) = 0.25d0 * (1 - xi(q)) * (1 - eta(q))  ! N1
            N(q,2) = 0.25d0 * (1 + xi(q)) * (1 - eta(q))  ! N2
            N(q,3) = 0.25d0 * (1 + xi(q)) * (1 + eta(q))  ! N3
            N(q,4) = 0.25d0 * (1 - xi(q)) * (1 + eta(q))  ! N4
        end do

        ! Compute values of u at the four cell corners

        select case(trim(grid))

            case("aa")
                u1 = 0.25d0 * (u(i-1, j-1) + u(i, j-1) + u(i-1, j) + u(i, j))  ! Bottom-left
                u2 = 0.25d0 * (u(i, j-1) + u(i+1, j-1) + u(i, j) + u(i+1, j))  ! Bottom-right
                u3 = 0.25d0 * (u(i, j) + u(i+1, j) + u(i, j+1) + u(i+1, j+1))  ! Top-right
                u4 = 0.25d0 * (u(i-1, j) + u(i, j) + u(i-1, j+1) + u(i, j+1))  ! Top-left
            case("acx")
                u1 = 0.5d0 * (u(i-1, j-1) + u(i-1, j))      ! Bottom-left
                u2 = 0.5d0 * (u(i, j-1) + u(i, j))          ! Bottom-right
                u3 = 0.5d0 * (u(i, j) + u(i, j+1))          ! Top-right
                u4 = 0.5d0 * (u(i-1, j) + u(i-1, j+1))      ! Top-left
            case("acy")
                u1 = 0.5d0 * (u(i-1, j-1) + u(i, j-1))      ! Bottom-left
                u2 = 0.5d0 * (u(i, j-1) + u(i+1, j-1))      ! Bottom-right
                u3 = 0.5d0 * (u(i, j) + u(i+1, j))          ! Top-right
                u4 = 0.5d0 * (u(i-1, j) + u(i, j))          ! Top-left
            case DEFAULT
                write(error_unit,*) "get_gauss_interp_points:: Error: grid not recognized."
                write(error_unit,*) "grid = ", trim(grid)
                stop
        end select

        ! Compute function values at quadrature points
        do q = 1, 4
            f(q) = N(q,1) * u1 + N(q,2) * u2 + N(q,3) * u3 + N(q,4) * u4
        end do

        return
        
    end subroutine get_gauss_interp_points


end module gaussian_quadrature