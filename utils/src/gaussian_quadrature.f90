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

    type gq3D_class
        real(8) :: xi(8)
        real(8) :: eta(8)
        real(8) :: zeta(8)
        real(8) :: N(8,8)
        real(8) :: dNdx(8,8)
        real(8) :: dNdy(8,8)
        real(8) :: dNdz(8,8)
    end type


    private
    public :: gq2D_class
    public :: gq2D_init
    public :: gq2D_to_nodes

    public :: gq3D_class
    public :: gq3D_init
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

            gq%N(1,q) = (1.d0 - gq%xi(q)) * (1.d0 - gq%eta(q)) / 4.d0  ! N1
            gq%N(2,q) = (1.d0 + gq%xi(q)) * (1.d0 - gq%eta(q)) / 4.d0  ! N2
            gq%N(3,q) = (1.d0 + gq%xi(q)) * (1.d0 + gq%eta(q)) / 4.d0  ! N3
            gq%N(4,q) = (1.d0 - gq%xi(q)) * (1.d0 + gq%eta(q)) / 4.d0  ! N4

            gq%dNdx(1,q) = -(1.d0 - gq%eta(q)) / 4.d0 
            gq%dNdx(2,q) =  (1.d0 - gq%eta(q)) / 4.d0 
            gq%dNdx(3,q) =  (1.d0 + gq%eta(q)) / 4.d0 
            gq%dNdx(4,q) = -(1.d0 + gq%eta(q)) / 4.d0

            gq%dNdy(1,q) = -(1.d0 - gq%xi(q)) / 4.d0 
            gq%dNdy(2,q) = -(1.d0 + gq%xi(q)) / 4.d0 
            gq%dNdy(3,q) =  (1.d0 + gq%xi(q)) / 4.d0 
            gq%dNdy(4,q) =  (1.d0 - gq%xi(q)) / 4.d0 

            if (.FALSE.) then
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

        subroutine gq3D_init(gq)

        implicit none

        type(gq3D_class), intent(OUT) :: gq

        ! Local variables
        integer :: q, n

        ! Quadrature point locations in reference element (-1,1)
        real(8), parameter :: sqrt3_inv = 1.0d0 / sqrt(3.0d0)

        ! Define quadrature points
        gq%xi   = [ -sqrt3_inv,  sqrt3_inv,  sqrt3_inv, -sqrt3_inv, -sqrt3_inv,  sqrt3_inv,  sqrt3_inv, -sqrt3_inv ]
        gq%eta  = [ -sqrt3_inv, -sqrt3_inv,  sqrt3_inv,  sqrt3_inv, -sqrt3_inv, -sqrt3_inv,  sqrt3_inv,  sqrt3_inv ]
        gq%zeta = [ -sqrt3_inv, -sqrt3_inv, -sqrt3_inv, -sqrt3_inv,  sqrt3_inv,  sqrt3_inv,  sqrt3_inv,  sqrt3_inv ]

        ! Define shape functions N(corner, q) and derivatives
        ! where q represents each quadrature point
        do q = 1, 8

            gq%N(1,q) = (1.d0 - gq%xi(q)) * (1.d0 - gq%eta(q)) * (1.d0 - gq%zeta(q)) / 8.d0  ! N1
            gq%N(2,q) = (1.d0 + gq%xi(q)) * (1.d0 - gq%eta(q)) * (1.d0 - gq%zeta(q)) / 8.d0  ! N2
            gq%N(3,q) = (1.d0 + gq%xi(q)) * (1.d0 + gq%eta(q)) * (1.d0 - gq%zeta(q)) / 8.d0  ! N3
            gq%N(4,q) = (1.d0 - gq%xi(q)) * (1.d0 + gq%eta(q)) * (1.d0 - gq%zeta(q)) / 8.d0  ! N4
            gq%N(5,q) = (1.d0 - gq%xi(q)) * (1.d0 - gq%eta(q)) * (1.d0 + gq%zeta(q)) / 8.d0  ! N5
            gq%N(6,q) = (1.d0 + gq%xi(q)) * (1.d0 - gq%eta(q)) * (1.d0 + gq%zeta(q)) / 8.d0  ! N6
            gq%N(7,q) = (1.d0 + gq%xi(q)) * (1.d0 + gq%eta(q)) * (1.d0 + gq%zeta(q)) / 8.d0  ! N7
            gq%N(8,q) = (1.d0 - gq%xi(q)) * (1.d0 + gq%eta(q)) * (1.d0 + gq%zeta(q)) / 8.d0  ! N8

            gq%dNdx(1,q) = -(1.d0 - gq%eta(q)) * (1.d0 - gq%zeta(q)) / 8.d0 
            gq%dNdx(2,q) =  (1.d0 - gq%eta(q)) * (1.d0 - gq%zeta(q)) / 8.d0 
            gq%dNdx(3,q) =  (1.d0 + gq%eta(q)) * (1.d0 - gq%zeta(q)) / 8.d0 
            gq%dNdx(4,q) = -(1.d0 + gq%eta(q)) * (1.d0 - gq%zeta(q)) / 8.d0 
            gq%dNdx(5,q) = -(1.d0 - gq%eta(q)) * (1.d0 + gq%zeta(q)) / 8.d0 
            gq%dNdx(6,q) =  (1.d0 - gq%eta(q)) * (1.d0 + gq%zeta(q)) / 8.d0 
            gq%dNdx(7,q) =  (1.d0 + gq%eta(q)) * (1.d0 + gq%zeta(q)) / 8.d0 
            gq%dNdx(8,q) = -(1.d0 + gq%eta(q)) * (1.d0 + gq%zeta(q)) / 8.d0 

            gq%dNdy(1,q) = -(1.d0 - gq%xi(q)) * (1.d0 - gq%zeta(q)) / 8.d0 
            gq%dNdy(2,q) = -(1.d0 + gq%xi(q)) * (1.d0 - gq%zeta(q)) / 8.d0 
            gq%dNdy(3,q) =  (1.d0 + gq%xi(q)) * (1.d0 - gq%zeta(q)) / 8.d0 
            gq%dNdy(4,q) =  (1.d0 - gq%xi(q)) * (1.d0 - gq%zeta(q)) / 8.d0 
            gq%dNdy(5,q) = -(1.d0 - gq%xi(q)) * (1.d0 + gq%zeta(q)) / 8.d0 
            gq%dNdy(6,q) = -(1.d0 + gq%xi(q)) * (1.d0 + gq%zeta(q)) / 8.d0 
            gq%dNdy(7,q) =  (1.d0 + gq%xi(q)) * (1.d0 + gq%zeta(q)) / 8.d0 
            gq%dNdy(8,q) =  (1.d0 - gq%xi(q)) * (1.d0 + gq%zeta(q)) / 8.d0 

            gq%dNdz(1,q) = -(1.d0 - gq%xi(q)) * (1.d0 - gq%eta(q)) / 8.d0 
            gq%dNdz(2,q) = -(1.d0 + gq%xi(q)) * (1.d0 - gq%eta(q)) / 8.d0 
            gq%dNdz(3,q) = -(1.d0 + gq%xi(q)) * (1.d0 + gq%eta(q)) / 8.d0 
            gq%dNdz(4,q) = -(1.d0 - gq%xi(q)) * (1.d0 + gq%eta(q)) / 8.d0 
            gq%dNdz(5,q) =  (1.d0 - gq%xi(q)) * (1.d0 - gq%eta(q)) / 8.d0 
            gq%dNdz(6,q) =  (1.d0 + gq%xi(q)) * (1.d0 - gq%eta(q)) / 8.d0 
            gq%dNdz(7,q) =  (1.d0 + gq%xi(q)) * (1.d0 + gq%eta(q)) / 8.d0 
            gq%dNdz(8,q) =  (1.d0 - gq%xi(q)) * (1.d0 + gq%eta(q)) / 8.d0 

            if (.TRUE.) then
                write(*,*) " "
                write(*,*) "Quad point, q =", q
                write(*,*) "n, phi_3d, dphi_dxr_3d, dphi_dyr_3d, dphi_dzr_3d:"
                do n = 1, 8
                    write(*,*) n, gq%N(n,q), gq%dNdx(n,q), gq%dNdy(n,q), gq%dNdz(n,q)
                enddo
                write(*,*) " "
                write(*,*) "sum(N)", sum(gq%N(:,q))  ! verified that sum = 1
                write(*,*) "sum(dN/dx)", sum(gq%dNdx(:,q))  ! verified that sum = 0 (within roundoff)
                write(*,*) "sum(dN/dy)", sum(gq%dNdy(:,q))  ! verified that sum = 0 (within roundoff)
                write(*,*) "sum(dN/dz)", sum(gq%dNdz(:,q))  ! verified that sum = 0 (within roundoff)
            end if

        end do

        return

    end subroutine gq3D_init

    subroutine gq2D_to_nodes(f, gq, u, i, j, grid_type)
        
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
                write(error_unit,*) "gq2D_to_nodes:: Error: grid_type not recognized."
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
        
    end subroutine gq2D_to_nodes

    subroutine gq3D_to_nodes(f, gq, u, i, j, k, grid_type)
        
        implicit none
        
        real(wp), intent(OUT) :: f(8)           ! Values at 8 Gaussian quadrature points
        type(gq3D_class), intent(IN) :: gq      ! Gaussian Quadrature 3D object
        real(wp), intent(in)  :: u(:,:,:)       ! Variable to be interpolated
        integer,  intent(in)  :: i              ! x-index of current cell
        integer,  intent(in)  :: j              ! y-index of current cell
        integer,  intent(in)  :: k              ! z-index of current cell
        character(len=*), intent(IN) :: grid_type   ! "aa", "ab", "acx", "acy", "acz"

        ! Local variables
        integer :: nx, ny, nz, q
        real(8) :: u1, u2, u3, u4, u5, u6, u7, u8    ! Values of u at the eight corners of the cell
        real(8) :: ux(8), uy(8), uz(8)               ! Derivatives at the eight corners
        
        nx = size(u,1)
        ny = size(u,2)
        nz = size(u,3)

        ! Compute values of u at the eight cell corners
        select case(trim(grid_type))
            case("aa")
                u1 = u(i-1, j-1, k-1)
                u2 = u(i, j-1, k-1)
                u3 = u(i, j, k-1)
                u4 = u(i-1, j, k-1)
                u5 = u(i-1, j-1, k)
                u6 = u(i, j-1, k)
                u7 = u(i, j, k)
                u8 = u(i-1, j, k)
            case DEFAULT
                write(error_unit,*) "gq3D_to_nodes:: Error: grid_type not recognized."
                write(error_unit,*) "grid_type = ", trim(grid_type)
                stop
        end select

if (.TRUE.) then
        ! Compute function values at quadrature points
        do q = 1, 8
            f(q) = gq%N(1,q) * u1 + gq%N(2,q) * u2 + gq%N(3,q) * u3 + gq%N(4,q) * u4 + &
                   gq%N(5,q) * u5 + gq%N(6,q) * u6 + gq%N(7,q) * u7 + gq%N(8,q) * u8
        end do
else
        ! Compute function values at quadrature points with Jacobian transformation
        ! do q = 1, 8
        !     f(q) = gq%N(1,q) * u1 + gq%N(2,q) * u2 + gq%N(3,q) * u3 + gq%N(4,q) * u4 + &
        !            gq%N(5,q) * u5 + gq%N(6,q) * u6 + gq%N(7,q) * u7 + gq%N(8,q) * u8 + &
        !            gq%dNdx(1,q) * ux1 + gq%dNdx(2,q) * ux2 + gq%dNdx(3,q) * ux3 + gq%dNdx(4,q) * ux4 + &
        !            gq%dNdx(5,q) * ux5 + gq%dNdx(6,q) * ux6 + gq%dNdx(7,q) * ux7 + gq%dNdx(8,q) * ux8 + &
        !            gq%dNdy(1,q) * uy1 + gq%dNdy(2,q) * uy2 + gq%dNdy(3,q) * uy3 + gq%dNdy(4,q) * uy4 + &
        !            gq%dNdy(5,q) * uy5 + gq%dNdy(6,q) * uy6 + gq%dNdy(7,q) * uy7 + gq%dNdy(8,q) * uy8 + &
        !            gq%dNdz(1,q) * uz1 + gq%dNdz(2,q) * uz2 + gq%dNdz(3,q) * uz3 + gq%dNdz(4,q) * uz4 + &
        !            gq%dNdz(5,q) * uz5 + gq%dNdz(6,q) * uz6 + gq%dNdz(7,q) * uz7 + gq%dNdz(8,q) * uz8
        ! end do
end if

        return
        
    end subroutine gq3D_to_nodes

end module gaussian_quadrature