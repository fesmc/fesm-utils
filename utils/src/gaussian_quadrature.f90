module gaussian_quadrature

    use, intrinsic :: iso_fortran_env, only : input_unit, output_unit, error_unit    
    use precision

    implicit none
    
    type gq2D_class
        real(8) :: xi(4)
        real(8) :: eta(4)
        real(8) :: N(4,4)
    end type

contains

    subroutine gq2D_init(gq)

        implicit none

        type(gq2D_class), intent(OUT) :: gq

        ! Local variables
        integer :: q

        ! Quadrature point locations in reference element (-1,1)
        real(8), parameter :: sqrt3_inv = 1.0d0 / sqrt(3.0d0)
        
        ! Define quadrature points
        gq%xi  = [ -sqrt3_inv, sqrt3_inv, sqrt3_inv, -sqrt3_inv ]
        gq%eta = [ -sqrt3_inv, -sqrt3_inv, sqrt3_inv, sqrt3_inv ]

        ! Define shape functions phi(q,corner)
        do q = 1, 4
            gq%N(q,1) = 0.25d0 * (1 - gq%xi(q)) * (1 - gq%eta(q))  ! N1
            gq%N(q,2) = 0.25d0 * (1 + gq%xi(q)) * (1 - gq%eta(q))  ! N2
            gq%N(q,3) = 0.25d0 * (1 + gq%xi(q)) * (1 + gq%eta(q))  ! N3
            gq%N(q,4) = 0.25d0 * (1 - gq%xi(q)) * (1 + gq%eta(q))  ! N4
        end do


        return

    end subroutine gq2D_init

    subroutine gq2D_get_points(f, gq, u, i, j, grid_type)
        
        implicit none
        
        real(wp), intent(OUT) :: f(4)           ! Values at 4 Gaussian quadrature points
        type(gq2D_class), intent(IN) :: gq      ! Gaussian Quadrature 2D object
        real(wp), intent(in)  :: u(:,:)         ! Variable to be interpolated
        integer,  intent(in)  :: i              ! x-index of current cell
        integer,  intent(in)  :: j              ! y-index of current cell
        character(len=*), intent(IN) :: grid_type   ! "aa", "ab", "acx", "acy"

        ! Local variables
        integer :: nx, ny, q
        real(8) :: u1, u2, u3, u4               ! Values of u at the four corners of the cell

        nx = size(u,1)
        ny = size(u,2)

        ! Compute values of u at the four cell corners

        select case(trim(grid_type))

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
                write(error_unit,*) "gq2D_get_points:: Error: grid_type not recognized."
                write(error_unit,*) "grid_type = ", trim(grid_type)
                stop
        end select

        ! Compute function values at quadrature points
        do q = 1, 4
            f(q) = gq%N(q,1) * u1 + gq%N(q,2) * u2 + gq%N(q,3) * u3 + gq%N(q,4) * u4
        end do

        return
        
    end subroutine gq2D_get_points


end module gaussian_quadrature