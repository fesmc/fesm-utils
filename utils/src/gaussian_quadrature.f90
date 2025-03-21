module gaussian_quadrature

    use, intrinsic :: iso_fortran_env, only : input_unit, output_unit, error_unit    
    use precision

    implicit none
    
    type gq2D_class
        real(8) :: xi(4)
        real(8) :: eta(4)
        real(8) :: N(4,4)
        real(8) :: dNdx(4,4)
        real(8) :: dNdy(4,4)
    end type

    private
    public :: gq2D_class
    public :: gq2D_init
    public :: gq2D_get_points
    
contains

    subroutine gq2D_init(gq)

        implicit none

        type(gq2D_class), intent(OUT) :: gq

        ! Local variables
        integer :: q, n

        ! Quadrature point locations in reference element (-1,1)
        real(8), parameter :: sqrt3_inv = 1.0d0 / sqrt(3.0d0)

        ! Define quadrature points
        gq%xi  = [ -sqrt3_inv, sqrt3_inv, sqrt3_inv, -sqrt3_inv ]
        gq%eta = [ -sqrt3_inv, -sqrt3_inv, sqrt3_inv, sqrt3_inv ]

        ! Define shape functions N(corner, q) and derivatives
        ! where q represents each quadrature point
        do q = 1, 4

            gq%N(1,q) = 0.25d0 * (1 - gq%xi(q)) * (1 - gq%eta(q))  ! N1
            gq%N(2,q) = 0.25d0 * (1 + gq%xi(q)) * (1 - gq%eta(q))  ! N2
            gq%N(3,q) = 0.25d0 * (1 + gq%xi(q)) * (1 + gq%eta(q))  ! N3
            gq%N(4,q) = 0.25d0 * (1 - gq%xi(q)) * (1 + gq%eta(q))  ! N4

            gq%dNdx(1,q) = -(1.d0 - gq%eta(q)) / 4.d0 
            gq%dNdx(2,q) =  (1.d0 - gq%eta(q)) / 4.d0 
            gq%dNdx(3,q) =  (1.d0 + gq%eta(q)) / 4.d0 
            gq%dNdx(4,q) = -(1.d0 + gq%eta(q)) / 4.d0

            gq%dNdy(1,q) = -(1.d0 - gq%xi(q)) / 4.d0 
            gq%dNdy(2,q) = -(1.d0 + gq%xi(q)) / 4.d0 
            gq%dNdy(3,q) =  (1.d0 + gq%xi(q)) / 4.d0 
            gq%dNdy(4,q) =  (1.d0 - gq%xi(q)) / 4.d0 

            if (.TRUE.) then
                write(*,*) " "
                write(*,*) "Quad point, q =", q
                write(*,*) "n, N, dNdx, dNdy:"
                do n = 1, 4
                    write(*,*) n, gq%N(n, q), gq%dNdx(n, q), gq%dNdy(n, q)
                end do
                write(*,*) "sum(N)", sum(gq%N(:, q))           ! Verified that sum = 1
                write(*,*) "sum(dN/dx)", sum(gq%dNdx(:, q))    ! Verified that sum = 0 (within roundoff)
                write(*,*) "sum(dN/dy)", sum(gq%dNdy(:, q))    ! Verified that sum = 0 (within roundoff)
            endif

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
        real(8) :: ux1, ux2, ux3, ux4           ! Derivatives du/dx at the four corners of the cell
        real(8) :: uy1, uy2, uy3, uy4           ! Derivatives du/dy at the four corners of the cell
        
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

if (.TRUE.) then
        ! Compute function values at quadrature points
        do q = 1, 4
            f(q) = gq%N(1,q) * u1 + gq%N(2,q) * u2 + gq%N(3,q) * u3 + gq%N(4,q) * u4
        end do

else
        ! Compute function values at quadrature points with Jacobian transformation
        do q = 1, 4
            f(q) = gq%N(1,q) * u1 + gq%N(2,q) * u2 + gq%N(3,q) * u3 + gq%N(4,q) * u4 + &
                gq%dNdx(1,q) * ux1 + gq%dNdx(2,q) * ux2 + gq%dNdx(3,q) * ux3 + gq%dNdx(4,q) * ux4 + &
                gq%dNdy(1,q) * uy1 + gq%dNdy(2,q) * uy2 + gq%dNdy(3,q) * uy3 + gq%dNdy(4,q) * uy4
        end do

end if

        return
        
    end subroutine gq2D_get_points


end module gaussian_quadrature