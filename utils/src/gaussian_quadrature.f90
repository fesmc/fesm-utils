module gaussian_quadrature

    use, intrinsic :: iso_fortran_env, only : input_unit, output_unit, error_unit    
    use precision

    implicit none
    
    logical :: verbose_init = .true. 
    logical :: verbose_jac  = .false.
    logical, parameter :: check_symmetry = .true.       ! if true, then check symmetry of assembled matrix

    integer :: itest, jtest, ktest                      ! coordinates of diagnostic point
    integer, parameter :: rtest = -9999                 ! task number for processor containing diagnostic point
    integer, parameter :: this_rank = 0

    !----------------------------------------------------------------
    ! Finite element properties
    ! Assume 3D hexahedral elements.
    !----------------------------------------------------------------

    real(dp), parameter :: rsqrt3 = 1.d0/sqrt(3.d0)     ! for quadrature points
    
    ! 3D
    integer, parameter :: nNodesPerElement_3d = 8       ! 8 nodes for hexahedral elements
    integer, parameter :: nQuadPoints_3d = 8            ! number of quadrature points per hexahedral element
                                                        ! These live at +- 1/sqrt(3) for reference hexahedron
    integer, parameter :: nNodeNeighbors_3d = 27        ! number of nearest node neighbors in 3D (including the node itself)

    ! 2D
    integer, parameter :: nNodesPerElement_2d = 4       ! 4 nodes for faces of hexahedral elements
    integer, parameter :: nQuadPoints_2d = 4            ! number of quadrature points per element face
                                                        ! These live at +- 1/sqrt(3) for reference square
    integer, parameter :: nNodeNeighbors_2d = 9         ! number of nearest node neighbors in 2D (including the node itself)

    !----------------------------------------------------------------
    ! Arrays used for finite-element calculations
    !
    ! Most integals are done over 3D hexahedral elements.
    ! Surface integrals are done over 2D faces of these elements. 
    !----------------------------------------------------------------

    real(dp), dimension(nNodesPerElement_3d, nQuadPoints_3d) ::   & 
       phi_3d,         &    ! trilinear basis function, evaluated at quad pts
       dphi_dxr_3d,    &    ! dphi/dx for reference hexehedral element, evaluated at quad pts
       dphi_dyr_3d,    &    ! dphi/dy for reference hexahedral element, evaluated at quad pts
       dphi_dzr_3d          ! dphi/dy for reference hexahedral element, evaluated at quad pts

    real(dp), dimension(nNodesPerElement_3d) ::   & 
       phi_3d_ctr,         &! trilinear basis function, evaluated at cell ctr
       dphi_dxr_3d_ctr,    &! dphi/dx for reference hexahedral element, evaluated at cell ctr
       dphi_dyr_3d_ctr,    &! dphi/dy for reference hexahedral element, evaluated at cell ctr
       dphi_dzr_3d_ctr      ! dphi/dz for reference hexahedral element, evaluated at cell ctr

    real(dp), dimension(nQuadPoints_3d) ::  &
       xqp_3d, yqp_3d, zqp_3d,  &! quad pt coordinates in reference element
       wqp_3d                    ! quad pt weights

    real(dp), dimension(nNodesPerElement_2d, nQuadPoints_2d) ::   & 
       phi_2d,         &    ! bilinear basis function, evaluated at quad pts
       dphi_dxr_2d,    &    ! dphi/dx for reference rectangular element, evaluated at quad pts
       dphi_dyr_2d          ! dphi/dy for reference rectangular element, evaluated at quad pts

    real(dp), dimension(nNodesPerElement_2d) ::   & 
       phi_2d_ctr,         &! bilinear basis function, evaluated at cell ctr
       dphi_dxr_2d_ctr,    &! dphi/dx for reference rectangular element, evaluated at cell ctr
       dphi_dyr_2d_ctr      ! dphi/dy for reference rectangular element, evaluated at cell ctr

    real(dp), dimension(nQuadPoints_2d) ::  &
       xqp_2d, yqp_2d, &    ! quad pt coordinates in reference square
       wqp_2d               ! quad pt weights

    integer, dimension(nNodesPerElement_3d, nNodesPerElement_3d) ::  &
       ishift, jshift, kshift   ! matrices describing relative indices of nodes in an element

    integer, dimension(-1:1,-1:1,-1:1) :: &
       indxA_3d              ! maps relative (x,y,z) coordinates to an index between 1 and 27
                             ! index order is (i,j,k)

    integer, dimension(-1:1,-1:1) :: &
       indxA_2d              ! maps relative (x,y) coordinates to an index between 1 and 9
                             ! index order is (i,j)

    real(dp), dimension(3,3) ::  &
       identity3             ! 3 x 3 identity matrix

    real(dp) :: vol0    ! volume scale (m^3), used to scale 3D matrix values

    !WHL - debug for efvs
    real(dp), dimension(nNodesPerElement_3d, nQuadPoints_2d) ::   & 
       phi_3d_vav,         &! vertical avg of phi_3d
       dphi_dxr_3d_vav,    &! vertical avg of dphi_dxr_3d
       dphi_dyr_3d_vav,    &! vertical avg of dphi_dyr_3d
       dphi_dzr_3d_vav      ! vertical avg of dphi_dzr_3d

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


    ! private
    ! public :: gq2D_class
    ! public :: gq2D_init
    ! public :: gq2D_to_nodes

    ! public :: gq3D_class
    ! public :: gq3D_init

    ! public :: gaussian_quadrature_init
    ! public :: gaussian_quadrature_2D_to_nodes
    public

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

! == Using CISM2.1 methods and variables ==

    subroutine gaussian_quadrature_2D_to_nodes(var_qp,var,dx,dy,i,j,grid_type)

        implicit none

        real(wp), intent(OUT) :: var_qp(:)
        real(wp), intent(IN)  :: var(:,:)
        ! real(wp), intent(IN)  :: xx(:,:)
        ! real(wp), intent(IN)  :: yy(:,:)
        real(wp), intent(IN)  :: dx
        real(wp), intent(IN)  :: dy
        integer,  intent(IN)  :: i
        integer,  intent(IN)  :: j
        character(len=*), intent(IN) :: grid_type

        ! Local variables

        integer :: q, n, p 
        real(dp) :: x(nNodesPerElement_2d)
        real(dp) :: y(nNodesPerElement_2d)
        real(dp) :: v(nNodesPerElement_2d)

        real(dp), dimension(nQuadPoints_2d) :: detJ                 ! determinant of J

        ! derivatives of basis function, evaluated at quad pts
        ! set dphi_dz = 0 for 2D problem
        real(dp) :: dphi_dx_2d(nNodesPerElement_2d)
        real(dp) :: dphi_dy_2d(nNodesPerElement_2d)
        real(dp) :: dphi_dz_2d(nNodesPerElement_2d)
                                            
        ! Step 1: determine x and y values of input array values

        ! Map xc,yc onto x,y vectors. Account for whether input
        ! variable is on aa, acx or acy grid. Account for boundary
        ! conditions (periodic, etc)

        ! Set x and y for each node

        !     4-----3       y
        !     |     |       ^
        !     |     |       |
        !     1-----2       ---> x

        ! Compute values of u at the four cell corners
        
        ! First get x and y values (xx and yy defined on aa-nodes)

        x(1) = -dx/2.d0
        x(2) =  dx/2.d0
        x(3) =  dx/2.d0
        x(4) = -dx/2.d0

        y(1) = -dy/2.d0
        y(2) = -dy/2.d0
        y(3) =  dy/2.d0
        y(4) =  dy/2.d0
        
        select case(trim(grid_type))

            case("aa")

                v(1) = 0.25d0 * (var(i-1, j-1) + var(i, j-1) + var(i-1, j) + var(i, j))  ! Bottom-left
                v(2) = 0.25d0 * (var(i, j-1) + var(i+1, j-1) + var(i, j) + var(i+1, j))  ! Bottom-right
                v(3) = 0.25d0 * (var(i, j) + var(i+1, j) + var(i, j+1) + var(i+1, j+1))  ! Top-right
                v(4) = 0.25d0 * (var(i-1, j) + var(i, j) + var(i-1, j+1) + var(i, j+1))  ! Top-left
            
            case("ab")

                v(1) = var(i-1,j-1)         ! Bottom-left
                v(2) = var(i,j-1)           ! Bottom-right
                v(3) = var(i,j)             ! Top-right
                v(4) = var(i-1,j)           ! Top-left

            case("acx")
                
                v(1) = 0.5d0 * (var(i-1, j-1) + var(i-1, j))      ! Bottom-left
                v(2) = 0.5d0 * (var(i, j-1) + var(i, j))          ! Bottom-right
                v(3) = 0.5d0 * (var(i, j) + var(i, j+1))          ! Top-right
                v(4) = 0.5d0 * (var(i-1, j) + var(i-1, j+1))      ! Top-left
            
            case("acy")

                v(1) = 0.5d0 * (var(i-1, j-1) + var(i, j-1))      ! Bottom-left
                v(2) = 0.5d0 * (var(i, j-1) + var(i+1, j-1))      ! Bottom-right
                v(3) = 0.5d0 * (var(i, j) + var(i+1, j))          ! Top-right
                v(4) = 0.5d0 * (var(i-1, j) + var(i, j))          ! Top-left

            case DEFAULT

                write(error_unit,*) "gaussian_quadrature_2D_to_nodes:: Error: grid_type not recognized."
                write(error_unit,*) "grid_type = ", trim(grid_type)
                stop
        
        end select

        ! Loop over quadrature points for this element
    
        do p = 1, nQuadPoints_2d

            ! Compute basis function derivatives and det(J) for this quadrature point
            ! For now, pass in i, j, k, p for debugging
            !TODO - Modify this subroutine so that the output derivatives are optional?

            call get_basis_function_derivatives_2d(x(:),             y(:),               &
                                                    dphi_dxr_2d(:,p), dphi_dyr_2d(:,p),   &
                                                    dphi_dx_2d(:),    dphi_dy_2d(:),      &
                                                    detJ(p),                                 &
                                                    itest, jtest, rtest,                  &
                                                    i, j, p)

            dphi_dz_2d(:) = 0.d0

            ! Evaluate var at this quadrature point, taking a phi-weighted sum over neighboring vertices.
            var_qp(p) = 0.d0
            do n = 1, nNodesPerElement_2d
                var_qp(p) = var_qp(p) + phi_2d(n,p) * v(n)
            end do
        
        end do

        return

    end subroutine gaussian_quadrature_2D_to_nodes

    subroutine gaussian_quadrature_3D_to_nodes(var_qp,var,dx,dy,dz0,dz1,i,j,k,grid_type)

        implicit none

        real(wp), intent(OUT) :: var_qp(:)
        real(wp), intent(IN)  :: var(:,:,:)
        ! real(wp), intent(IN)  :: xx(:,:)
        ! real(wp), intent(IN)  :: yy(:,:)
        real(wp), intent(IN)  :: dx
        real(wp), intent(IN)  :: dy
        real(wp), intent(IN)  :: dz0
        real(wp), intent(IN)  :: dz1
        integer,  intent(IN)  :: i
        integer,  intent(IN)  :: j
        integer,  intent(IN)  :: k
        character(len=*), intent(IN) :: grid_type

        ! Local variables

        integer :: q, n, p 
        real(dp) :: x(nNodesPerElement_3d)
        real(dp) :: y(nNodesPerElement_3d)
        real(dp) :: z(nNodesPerElement_3d)
        real(dp) :: v(nNodesPerElement_3d)

        real(dp), dimension(nQuadPoints_3d) :: detJ                 ! determinant of J

        ! derivatives of basis function, evaluated at quad pts
        real(dp) :: dphi_dx_3d(nNodesPerElement_3d)
        real(dp) :: dphi_dy_3d(nNodesPerElement_3d)
        real(dp) :: dphi_dz_3d(nNodesPerElement_3d)
                                            
        ! Step 1: determine x and y values of input array values

        ! Map xc,yc onto x,y vectors. Account for whether input
        ! variable is on aa, acx or acy grid. Account for boundary
        ! conditions (periodic, etc)

        ! Set x and y for each node

        !     4-----3       y
        !     |     |       ^
        !     |     |       |
        !     1-----2       ---> x

        ! Compute values of u at the four cell corners
        
        ! First get x and y values (xx and yy defined on aa-nodes)

        x(1) = -dx/2.d0
        x(2) =  dx/2.d0
        x(3) =  dx/2.d0
        x(4) = -dx/2.d0
        x(5) = -dx/2.d0
        x(6) =  dx/2.d0
        x(7) =  dx/2.d0
        x(8) = -dx/2.d0

        y(1) = -dy/2.d0
        y(2) = -dy/2.d0
        y(3) =  dy/2.d0
        y(4) =  dy/2.d0
        y(5) = -dy/2.d0
        y(6) = -dy/2.d0
        y(7) =  dy/2.d0
        y(8) =  dy/2.d0
        
        z(1) = -dz0/2.d0
        z(2) = -dz0/2.d0
        z(3) = -dz0/2.d0
        z(4) = -dz0/2.d0
        z(5) =  dz1/2.d0
        z(6) =  dz1/2.d0
        z(7) =  dz1/2.d0
        z(8) =  dz1/2.d0
        
        ! select case(trim(grid_type))

        !     case("aa-aa")

        !         v(1) = 0.25d0 * (var(i-1, j-1) + var(i, j-1) + var(i-1, j) + var(i, j))  ! Bottom-left
        !         v(2) = 0.25d0 * (var(i, j-1) + var(i+1, j-1) + var(i, j) + var(i+1, j))  ! Bottom-right
        !         v(3) = 0.25d0 * (var(i, j) + var(i+1, j) + var(i, j+1) + var(i+1, j+1))  ! Top-right
        !         v(4) = 0.25d0 * (var(i-1, j) + var(i, j) + var(i-1, j+1) + var(i, j+1))  ! Top-left
            
        !     case("ab")

        !         v(1) = var(i-1,j-1)         ! Bottom-left
        !         v(2) = var(i,j-1)           ! Bottom-right
        !         v(3) = var(i,j)             ! Top-right
        !         v(4) = var(i-1,j)           ! Top-left

        !     case("acx")
                
        !         v(1) = 0.5d0 * (var(i-1, j-1) + var(i-1, j))      ! Bottom-left
        !         v(2) = 0.5d0 * (var(i, j-1) + var(i, j))          ! Bottom-right
        !         v(3) = 0.5d0 * (var(i, j) + var(i, j+1))          ! Top-right
        !         v(4) = 0.5d0 * (var(i-1, j) + var(i-1, j+1))      ! Top-left
            
        !     case("acy")

        !         v(1) = 0.5d0 * (var(i-1, j-1) + var(i, j-1))      ! Bottom-left
        !         v(2) = 0.5d0 * (var(i, j-1) + var(i+1, j-1))      ! Bottom-right
        !         v(3) = 0.5d0 * (var(i, j) + var(i+1, j))          ! Top-right
        !         v(4) = 0.5d0 * (var(i-1, j) + var(i, j))          ! Top-left

        !     case DEFAULT

        !         write(error_unit,*) "gaussian_quadrature_3D_to_nodes:: Error: grid_type not recognized."
        !         write(error_unit,*) "grid_type = ", trim(grid_type)
        !         stop
        
        ! end select

        ! Loop over quadrature points for this element
    
        do p = 1, nQuadPoints_3d

            ! Compute basis function derivatives and det(J) for this quadrature point
            ! For now, pass in i, j, k, p for debugging
            !TODO - Modify this subroutine so that the output derivatives are optional?

            call get_basis_function_derivatives_3d(x(:),             y(:),          z(:),               &
                                                    dphi_dxr_3d(:,p), dphi_dyr_3d(:,p), dphi_dzr_3d(:,p),   &
                                                    dphi_dx_3d(:),    dphi_dy_3d(:),    dphi_dz_3d(:),  &
                                                    detJ(p),                                 &
                                                    itest, jtest, rtest,                  &
                                                    i, j, k, p)

            ! Evaluate var at this quadrature point, taking a phi-weighted sum over neighboring vertices.
            var_qp(p) = 0.d0
            do n = 1, nNodesPerElement_3d
                var_qp(p) = var_qp(p) + phi_3d(n,p) * v(n)
            end do
        
        end do

        return

    end subroutine gaussian_quadrature_3D_to_nodes

! == Directly from CISM2.1 ==

  subroutine gaussian_quadrature_init()

    !----------------------------------------------------------------
    ! Initial calculations for glissade higher-order solver.
    !----------------------------------------------------------------

    integer :: i, j, k, m, n, p
    integer :: pplus
    real(dp) :: xctr, yctr, zctr
    real(dp) :: sumx, sumy, sumz

    !----------------------------------------------------------------
    ! Initialize some time-independent finite element arrays
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Trilinear basis set for reference hexahedron, x=(-1,1), y=(-1,1), z=(-1,1)             
    ! Indexing is counter-clockwise from SW corner, with 1-4 on lower surface
    !  and 5-8 on upper surface
    ! The code uses "phi_3d" to denote these basis functions. 
    !
    ! N1 = (1-x)*(1-y)*(1-z)/8             N4----N3
    ! N2 = (1+x)*(1-y)*(1-z)/8             |     |    Lower layer        
    ! N3 = (1+x)*(1+y)*(1-z)/8             |     |
    ! N4 = (1-x)*(1+y)*(1-z)/8             N1----N2

    ! N5 = (1-x)*(1-y)*(1+z)/8             N8----N7
    ! N6 = (1+x)*(1-y)*(1+z)/8             |     |    Upper layer
    ! N7 = (1+x)*(1+y)*(1+z)/8             |     |
    ! N8 = (1-x)*(1+y)*(1+z)/8             N5----N6
    !----------------------------------------------------------------
   
    ! Set coordinates and weights of quadrature points for reference hexahedral element.
    ! Numbering is counter-clockwise from southwest, lower face (1-4) followed by
    !  upper face (5-8).

    xqp_3d(1) = -rsqrt3; yqp_3d(1) = -rsqrt3; zqp_3d(1) = -rsqrt3
    wqp_3d(1) =  1.d0

    xqp_3d(2) =  rsqrt3; yqp_3d(2) = -rsqrt3; zqp_3d(2) = -rsqrt3
    wqp_3d(2) =  1.d0

    xqp_3d(3) =  rsqrt3; yqp_3d(3) =  rsqrt3; zqp_3d(3) = -rsqrt3
    wqp_3d(3) =  1.d0

    xqp_3d(4) = -rsqrt3; yqp_3d(4) =  rsqrt3; zqp_3d(4) = -rsqrt3
    wqp_3d(4) =  1.d0

    xqp_3d(5) = -rsqrt3; yqp_3d(5) = -rsqrt3; zqp_3d(5) =  rsqrt3
    wqp_3d(5) =  1.d0

    xqp_3d(6) =  rsqrt3; yqp_3d(6) = -rsqrt3; zqp_3d(6) =  rsqrt3
    wqp_3d(6) =  1.d0

    xqp_3d(7) =  rsqrt3; yqp_3d(7) =  rsqrt3; zqp_3d(7) =  rsqrt3
    wqp_3d(7) =  1.d0

    xqp_3d(8) = -rsqrt3; yqp_3d(8) =  rsqrt3; zqp_3d(8) =  rsqrt3
    wqp_3d(8) =  1.d0

    if (verbose_init) then
       print*, ' '
       print*, 'Hexahedral elements, quad points, x, y, z:'
       sumx = 0.d0; sumy = 0.d0; sumz = 0.d0
       do p = 1, nQuadPoints_3d
          print*, p, xqp_3d(p), yqp_3d(p), zqp_3d(p)
          sumx = sumx + xqp_3d(p); sumy = sumy + yqp_3d(p); sumz = sumz + zqp_3d(p)
       enddo
       print*, ' '
       print*, 'sums:', sumx, sumy, sumz
    endif

    ! Evaluate trilinear basis functions and their derivatives at each quad pt

    do p = 1, nQuadPoints_3d

       phi_3d(1,p) = (1.d0 - xqp_3d(p)) * (1.d0 - yqp_3d(p)) * (1.d0 - zqp_3d(p)) / 8.d0
       phi_3d(2,p) = (1.d0 + xqp_3d(p)) * (1.d0 - yqp_3d(p)) * (1.d0 - zqp_3d(p)) / 8.d0
       phi_3d(3,p) = (1.d0 + xqp_3d(p)) * (1.d0 + yqp_3d(p)) * (1.d0 - zqp_3d(p)) / 8.d0
       phi_3d(4,p) = (1.d0 - xqp_3d(p)) * (1.d0 + yqp_3d(p)) * (1.d0 - zqp_3d(p)) / 8.d0
       phi_3d(5,p) = (1.d0 - xqp_3d(p)) * (1.d0 - yqp_3d(p)) * (1.d0 + zqp_3d(p)) / 8.d0
       phi_3d(6,p) = (1.d0 + xqp_3d(p)) * (1.d0 - yqp_3d(p)) * (1.d0 + zqp_3d(p)) / 8.d0
       phi_3d(7,p) = (1.d0 + xqp_3d(p)) * (1.d0 + yqp_3d(p)) * (1.d0 + zqp_3d(p)) / 8.d0
       phi_3d(8,p) = (1.d0 - xqp_3d(p)) * (1.d0 + yqp_3d(p)) * (1.d0 + zqp_3d(p)) / 8.d0

       dphi_dxr_3d(1,p) = -(1.d0 - yqp_3d(p)) * (1.d0 - zqp_3d(p)) / 8.d0 
       dphi_dxr_3d(2,p) =  (1.d0 - yqp_3d(p)) * (1.d0 - zqp_3d(p)) / 8.d0 
       dphi_dxr_3d(3,p) =  (1.d0 + yqp_3d(p)) * (1.d0 - zqp_3d(p)) / 8.d0 
       dphi_dxr_3d(4,p) = -(1.d0 + yqp_3d(p)) * (1.d0 - zqp_3d(p)) / 8.d0
       dphi_dxr_3d(5,p) = -(1.d0 - yqp_3d(p)) * (1.d0 + zqp_3d(p)) / 8.d0 
       dphi_dxr_3d(6,p) =  (1.d0 - yqp_3d(p)) * (1.d0 + zqp_3d(p)) / 8.d0 
       dphi_dxr_3d(7,p) =  (1.d0 + yqp_3d(p)) * (1.d0 + zqp_3d(p)) / 8.d0 
       dphi_dxr_3d(8,p) = -(1.d0 + yqp_3d(p)) * (1.d0 + zqp_3d(p)) / 8.d0

       dphi_dyr_3d(1,p) = -(1.d0 - xqp_3d(p)) * (1.d0 - zqp_3d(p)) / 8.d0 
       dphi_dyr_3d(2,p) = -(1.d0 + xqp_3d(p)) * (1.d0 - zqp_3d(p)) / 8.d0 
       dphi_dyr_3d(3,p) =  (1.d0 + xqp_3d(p)) * (1.d0 - zqp_3d(p)) / 8.d0 
       dphi_dyr_3d(4,p) =  (1.d0 - xqp_3d(p)) * (1.d0 - zqp_3d(p)) / 8.d0 
       dphi_dyr_3d(5,p) = -(1.d0 - xqp_3d(p)) * (1.d0 + zqp_3d(p)) / 8.d0 
       dphi_dyr_3d(6,p) = -(1.d0 + xqp_3d(p)) * (1.d0 + zqp_3d(p)) / 8.d0 
       dphi_dyr_3d(7,p) =  (1.d0 + xqp_3d(p)) * (1.d0 + zqp_3d(p)) / 8.d0 
       dphi_dyr_3d(8,p) =  (1.d0 - xqp_3d(p)) * (1.d0 + zqp_3d(p)) / 8.d0 

       dphi_dzr_3d(1,p) = -(1.d0 - xqp_3d(p)) * (1.d0 - yqp_3d(p)) / 8.d0 
       dphi_dzr_3d(2,p) = -(1.d0 + xqp_3d(p)) * (1.d0 - yqp_3d(p)) / 8.d0 
       dphi_dzr_3d(3,p) = -(1.d0 + xqp_3d(p)) * (1.d0 + yqp_3d(p)) / 8.d0 
       dphi_dzr_3d(4,p) = -(1.d0 - xqp_3d(p)) * (1.d0 + yqp_3d(p)) / 8.d0 
       dphi_dzr_3d(5,p) =  (1.d0 - xqp_3d(p)) * (1.d0 - yqp_3d(p)) / 8.d0 
       dphi_dzr_3d(6,p) =  (1.d0 + xqp_3d(p)) * (1.d0 - yqp_3d(p)) / 8.d0 
       dphi_dzr_3d(7,p) =  (1.d0 + xqp_3d(p)) * (1.d0 + yqp_3d(p)) / 8.d0 
       dphi_dzr_3d(8,p) =  (1.d0 - xqp_3d(p)) * (1.d0 + yqp_3d(p)) / 8.d0 

       if (verbose_init) then
          print*, ' '
          print*, 'Quad point, p =', p
          print*, 'n, phi_3d, dphi_dxr_3d, dphi_dyr_3d, dphi_dzr_3d:'
          do n = 1, 8
             print*, n, phi_3d(n,p), dphi_dxr_3d(n,p), dphi_dyr_3d(n,p), dphi_dzr_3d(n,p)
          enddo
          print*, ' '
          print*, 'sum(phi_3d)', sum(phi_3d(:,p))  ! verified that sum = 1
          print*, 'sum(dphi/dx)', sum(dphi_dxr_3d(:,p))  ! verified that sum = 0 (within roundoff)
          print*, 'sum(dphi/dy)', sum(dphi_dyr_3d(:,p))  ! verified that sum = 0 (within roundoff)
          print*, 'sum(dphi/dz)', sum(dphi_dzr_3d(:,p))  ! verified that sum = 0 (within roundoff)
       endif

    enddo   ! nQuadPoints_3d

    ! Evaluate trilinear basis functions and their derivatives at cell center
    ! Full formulas are not really needed at (x,y,z) = (0,0,0), but are included for completeness

    xctr = 0.d0
    yctr = 0.d0
    zctr = 0.d0

    phi_3d_ctr(1) = (1.d0 - xctr) * (1.d0 - yctr) * (1.d0 - zctr) / 8.d0
    phi_3d_ctr(2) = (1.d0 + xctr) * (1.d0 - yctr) * (1.d0 - zctr) / 8.d0
    phi_3d_ctr(3) = (1.d0 + xctr) * (1.d0 + yctr) * (1.d0 - zctr) / 8.d0
    phi_3d_ctr(4) = (1.d0 - xctr) * (1.d0 + yctr) * (1.d0 - zctr) / 8.d0
    phi_3d_ctr(5) = (1.d0 - xctr) * (1.d0 - yctr) * (1.d0 + zctr) / 8.d0
    phi_3d_ctr(6) = (1.d0 + xctr) * (1.d0 - yctr) * (1.d0 + zctr) / 8.d0
    phi_3d_ctr(7) = (1.d0 + xctr) * (1.d0 + yctr) * (1.d0 + zctr) / 8.d0
    phi_3d_ctr(8) = (1.d0 - xctr) * (1.d0 + yctr) * (1.d0 + zctr) / 8.d0
    
    dphi_dxr_3d_ctr(1) = -(1.d0 - yctr) * (1.d0 - zctr) / 8.d0 
    dphi_dxr_3d_ctr(2) =  (1.d0 - yctr) * (1.d0 - zctr) / 8.d0 
    dphi_dxr_3d_ctr(3) =  (1.d0 + yctr) * (1.d0 - zctr) / 8.d0 
    dphi_dxr_3d_ctr(4) = -(1.d0 + yctr) * (1.d0 - zctr) / 8.d0
    dphi_dxr_3d_ctr(5) = -(1.d0 - yctr) * (1.d0 + zctr) / 8.d0 
    dphi_dxr_3d_ctr(6) =  (1.d0 - yctr) * (1.d0 + zctr) / 8.d0 
    dphi_dxr_3d_ctr(7) =  (1.d0 + yctr) * (1.d0 + zctr) / 8.d0 
    dphi_dxr_3d_ctr(8) = -(1.d0 + yctr) * (1.d0 + zctr) / 8.d0
    
    dphi_dyr_3d_ctr(1) = -(1.d0 - xctr) * (1.d0 - zctr) / 8.d0 
    dphi_dyr_3d_ctr(2) = -(1.d0 + xctr) * (1.d0 - zctr) / 8.d0 
    dphi_dyr_3d_ctr(3) =  (1.d0 + xctr) * (1.d0 - zctr) / 8.d0 
    dphi_dyr_3d_ctr(4) =  (1.d0 - xctr) * (1.d0 - zctr) / 8.d0 
    dphi_dyr_3d_ctr(5) = -(1.d0 - xctr) * (1.d0 + zctr) / 8.d0 
    dphi_dyr_3d_ctr(6) = -(1.d0 + xctr) * (1.d0 + zctr) / 8.d0 
    dphi_dyr_3d_ctr(7) =  (1.d0 + xctr) * (1.d0 + zctr) / 8.d0 
    dphi_dyr_3d_ctr(8) =  (1.d0 - xctr) * (1.d0 + zctr) / 8.d0 
    
    dphi_dzr_3d_ctr(1) = -(1.d0 - xctr) * (1.d0 - yctr) / 8.d0 
    dphi_dzr_3d_ctr(2) = -(1.d0 + xctr) * (1.d0 - yctr) / 8.d0 
    dphi_dzr_3d_ctr(3) = -(1.d0 + xctr) * (1.d0 + yctr) / 8.d0 
    dphi_dzr_3d_ctr(4) = -(1.d0 - xctr) * (1.d0 + yctr) / 8.d0 
    dphi_dzr_3d_ctr(5) =  (1.d0 - xctr) * (1.d0 - yctr) / 8.d0 
    dphi_dzr_3d_ctr(6) =  (1.d0 + xctr) * (1.d0 - yctr) / 8.d0 
    dphi_dzr_3d_ctr(7) =  (1.d0 + xctr) * (1.d0 + yctr) / 8.d0 
    dphi_dzr_3d_ctr(8) =  (1.d0 - xctr) * (1.d0 + yctr) / 8.d0 

    ! Identity matrix
    identity3(1,:) = (/ 1.d0, 0.d0, 0.d0 /)
    identity3(2,:) = (/ 0.d0, 1.d0, 0.d0 /)
    identity3(3,:) = (/ 0.d0, 0.d0, 1.d0 /)

    ! Initialize some matrices that describe how the i, j and k indices of each node
    ! in each element are related to one another.

    ! The ishift matrix describes how the i indices of the 8 nodes are related to one another.
    ! E.g, if ishift (1,2) = 1, this means that node 2 has an i index
    ! one greater than the i index of node 1.

    ishift(1,:) = (/ 0,  1,  1,  0,  0,  1,  1,  0/)   
    ishift(2,:) = (/-1,  0,  0, -1, -1,  0,  0, -1/)   
    ishift(3,:) = ishift(2,:)
    ishift(4,:) = ishift(1,:)
    ishift(5,:) = ishift(1,:)
    ishift(6,:) = ishift(2,:)
    ishift(7,:) = ishift(2,:)
    ishift(8,:) = ishift(1,:)

    ! The jshift matrix describes how the j indices of the 8 nodes are related to one another.
    ! E.g, if jshift (1,4) = 1, this means that node 4 has a j index
    ! one greater than the j index of node 1.

    jshift(1,:) = (/ 0,  0,  1,  1,  0,  0,  1,  1/)   
    jshift(2,:) = jshift(1,:)
    jshift(3,:) = (/-1, -1,  0,  0, -1, -1,  0,  0/)   
    jshift(4,:) = jshift(3,:)
    jshift(5,:) = jshift(1,:)
    jshift(6,:) = jshift(1,:)
    jshift(7,:) = jshift(3,:)
    jshift(8,:) = jshift(3,:)

    ! The kshift matrix describes how the k indices of the 8 nodes are related to one another.
    ! E.g, if kshift (1,5) = -1, this means that node 5 has a k index
    ! one less than the k index of node 1.  (Assume that k increases downward.)

    kshift(1,:) = (/ 0,  0,  0,  0, -1, -1, -1, -1/)   
    kshift(2,:) = kshift(1,:)
    kshift(3,:) = kshift(1,:)
    kshift(4,:) = kshift(1,:)
    kshift(5,:) = (/ 1,  1,  1,  1,  0,  0,  0,  0/)
    kshift(6,:) = kshift(5,:)
    kshift(7,:) = kshift(5,:)
    kshift(8,:) = kshift(5,:)

    if (verbose_init) then
       print*, ' '
       print*, 'ishift:'
       do n = 1, 8
          write (6,'(8i4)') ishift(n,:)
       enddo
       print*, ' '
       print*, 'jshift:'
       do n = 1, 8
          write (6,'(8i4)') jshift(n,:)
       enddo
       print*, ' '
       print*, 'kshift:'
       do n = 1, 8
          write (6,'(8i4)') kshift(n,:)
       enddo
    endif

    !----------------------------------------------------------------
    ! Bilinear basis set for reference square, x=(-1,1), y=(-1,1)             
    ! Indexing is counter-clockwise from SW corner
    ! The code uses "phi_2d" to denote these basis functions. 
    !
    ! N1 = (1-x)*(1-y)/4             N4----N3
    ! N2 = (1+x)*(1-y)/4             |     |
    ! N3 = (1+x)*(1+y)/4             |     |
    ! N4 = (1-x)*(1+y)/4             N1----N2
    !----------------------------------------------------------------

    ! Set coordinates and weights of quadrature points for reference square.
    ! Numbering is counter-clockwise from southwest

    xqp_2d(1) = -rsqrt3; yqp_2d(1) = -rsqrt3
    wqp_2d(1) =  1.d0

    xqp_2d(2) =  rsqrt3; yqp_2d(2) = -rsqrt3
    wqp_2d(2) =  1.d0

    xqp_2d(3) =  rsqrt3; yqp_2d(3) =  rsqrt3
    wqp_2d(3) =  1.d0

    xqp_2d(4) = -rsqrt3; yqp_2d(4) =  rsqrt3
    wqp_2d(4) =  1.d0

    if (verbose_init) then
       print*, ' '
       print*, ' '
       print*, 'Quadrilateral elements, quad points, x, y:'
       sumx = 0.d0; sumy = 0.d0; sumz = 0.d0
       do p = 1, nQuadPoints_2d
          print*, p, xqp_2d(p), yqp_2d(p)
          sumx = sumx + xqp_2d(p); sumy = sumy + yqp_2d(p)
       enddo
       print*, ' '
       print*, 'sumx, sumy:', sumx, sumy
    endif

    ! Evaluate bilinear basis functions and their derivatives at each quad pt

    do p = 1, nQuadPoints_2d

       phi_2d(1,p) = (1.d0 - xqp_2d(p)) * (1.d0 - yqp_2d(p)) / 4.d0 
       phi_2d(2,p) = (1.d0 + xqp_2d(p)) * (1.d0 - yqp_2d(p)) / 4.d0
       phi_2d(3,p) = (1.d0 + xqp_2d(p)) * (1.d0 + yqp_2d(p)) / 4.d0 
       phi_2d(4,p) = (1.d0 - xqp_2d(p)) * (1.d0 + yqp_2d(p)) / 4.d0

       dphi_dxr_2d(1,p) = -(1.d0 - yqp_2d(p)) / 4.d0 
       dphi_dxr_2d(2,p) =  (1.d0 - yqp_2d(p)) / 4.d0 
       dphi_dxr_2d(3,p) =  (1.d0 + yqp_2d(p)) / 4.d0 
       dphi_dxr_2d(4,p) = -(1.d0 + yqp_2d(p)) / 4.d0

       dphi_dyr_2d(1,p) = -(1.d0 - xqp_2d(p)) / 4.d0 
       dphi_dyr_2d(2,p) = -(1.d0 + xqp_2d(p)) / 4.d0 
       dphi_dyr_2d(3,p) =  (1.d0 + xqp_2d(p)) / 4.d0 
       dphi_dyr_2d(4,p) =  (1.d0 - xqp_2d(p)) / 4.d0 

       if (verbose_init) then
          print*, ' '
          print*, 'Quad point, p =', p
          print*, 'n, phi_2d, dphi_dxr_2d, dphi_dyr_2d:'
          do n = 1, 4
             print*, n, phi_2d(n,p), dphi_dxr_2d(n,p), dphi_dyr_2d(n,p)
          enddo
          print*, 'sum(phi_2d)', sum(phi_2d(:,p))        ! verified that sum = 1
          print*, 'sum(dphi/dx_2d)', sum(dphi_dxr_2d(:,p))  ! verified that sum = 0 (within roundoff)
          print*, 'sum(dphi/dy_2d)', sum(dphi_dyr_2d(:,p))  ! verified that sum = 0 (within roundoff)
       endif

    enddo   ! nQuadPoints_2d

    ! Evaluate bilinear basis functions and their derivatives at cell center
    ! Full formulas are not really needed at (x,y) = (0,0), but are included for completeness

    xctr = 0.d0
    yctr = 0.d0

    phi_2d_ctr(1) = (1.d0 - xctr) * (1.d0 - yctr) / 4.d0 
    phi_2d_ctr(2) = (1.d0 + xctr) * (1.d0 - yctr) / 4.d0
    phi_2d_ctr(3) = (1.d0 + xctr) * (1.d0 + yctr) / 4.d0 
    phi_2d_ctr(4) = (1.d0 - xctr) * (1.d0 + yctr) / 4.d0
    
    dphi_dxr_2d_ctr(1) = -(1.d0 - yctr) / 4.d0 
    dphi_dxr_2d_ctr(2) =  (1.d0 - yctr) / 4.d0 
    dphi_dxr_2d_ctr(3) =  (1.d0 + yctr) / 4.d0 
    dphi_dxr_2d_ctr(4) = -(1.d0 + yctr) / 4.d0

    dphi_dyr_2d_ctr(1) = -(1.d0 - xctr) / 4.d0 
    dphi_dyr_2d_ctr(2) = -(1.d0 + xctr) / 4.d0 
    dphi_dyr_2d_ctr(3) =  (1.d0 + xctr) / 4.d0 
    dphi_dyr_2d_ctr(4) =  (1.d0 - xctr) / 4.d0 

    !----------------------------------------------------------------
    ! Compute indxA_3d; maps displacements i,j,k = (-1,0,1) onto an index from 1 to 27
    ! Numbering starts in SW corner of layers k-1, finishes in NE corner of layer k+1
    ! Diagonal term has index 14
    !----------------------------------------------------------------

    ! Layer k-1:           Layer k:            Layer k+1:
    !
    !   7    8    9          16   17   18        25   26   27 
    !   4    5    6          13   14   15        22   23   24
    !   1    2    3          10   11   12        19   20   21                                                                                               

    m = 0
    do k = -1,1
       do j = -1,1
          do i = -1,1
             m = m + 1
             indxA_3d(i,j,k) = m
          enddo
       enddo
    enddo

    !----------------------------------------------------------------
    ! Compute indxA_2d; maps displacements i,j = (-1,0,1) onto an index from 1 to 9
    ! Same as indxA_3d, but for a single layer
    !----------------------------------------------------------------

    m = 0
    do j = -1,1
       do i = -1,1
          m = m + 1
          indxA_2d(i,j) = m
       enddo
    enddo

    !WHL - debug for efvs

    ! Evaluate vertical averages of dphi_dxr_3d, dphi_dyr_3d and dphi_dzr_3d at each 2d quad pts.
    ! Using these instead of the full 3d basis functions can result in similar accuracy with
    !  only half as many QP computations.

    do p = 1, nQuadPoints_2d
       pplus = p + nQuadPoints_3d/2  ! p + 4 for hexahedra
       do n = 1, nNodesPerElement_3d
          phi_3d_vav(n,p) = 0.5d0 * (phi_3d(n,p) + phi_3d(n,pplus))
          dphi_dxr_3d_vav(n,p) = 0.5d0 * (dphi_dxr_3d(n,p) + dphi_dxr_3d(n,pplus))
          dphi_dyr_3d_vav(n,p) = 0.5d0 * (dphi_dyr_3d(n,p) + dphi_dyr_3d(n,pplus))
          dphi_dzr_3d_vav(n,p) = 0.5d0 * (dphi_dzr_3d(n,p) + dphi_dzr_3d(n,pplus))
       enddo
    enddo

  end subroutine gaussian_quadrature_init

!****************************************************************************

  subroutine get_basis_function_derivatives_3d(xNode,       yNode,       zNode,       &
                                               dphi_dxr_3d, dphi_dyr_3d, dphi_dzr_3d, &
                                               dphi_dx_3d,  dphi_dy_3d,  dphi_dz_3d,  &
                                               detJ,                                  &
                                               itest, jtest, rtest,                   &
                                               i, j, k, p)

    !------------------------------------------------------------------
    ! Evaluate the x, y and z derivatives of the element basis functions
    ! at a particular quadrature point.
    !
    ! Also determine the Jacobian of the transformation between the
    ! reference element and the true element.
    ! 
    ! This subroutine should work for any 3D element with any number of nodes.
    !------------------------------------------------------------------
 
    real(dp), dimension(nNodesPerElement_3d), intent(in) :: &
       xNode, yNode, zNode,          &! nodal coordinates
       dphi_dxr_3d, dphi_dyr_3d, dphi_dzr_3d   ! derivatives of basis functions at quad pt
                                               !  wrt x, y and z in reference element

    real(dp), dimension(nNodesPerElement_3d), intent(out) :: &
       dphi_dx_3d, dphi_dy_3d, dphi_dz_3d      ! derivatives of basis functions at quad pt
                                               !  wrt x, y and z in true Cartesian coordinates  

    real(dp), intent(out) :: &
         detJ      ! determinant of Jacobian matrix

    real(dp), dimension(3,3) ::  &
         Jac,      &! Jacobian matrix
         Jinv,     &! inverse Jacobian matrix
         cofactor   ! matrix of cofactors

    integer, intent(in) :: &
       itest, jtest, rtest              ! coordinates of diagnostic point

    integer, intent(in) :: i, j, k, p   ! indices passed in for debugging

    integer :: n, row, col

    logical, parameter :: Jac_bug_check = .false.   ! set to true for debugging
    real(dp), dimension(3,3) :: prod     ! Jac * Jinv (should be identity matrix)

    !------------------------------------------------------------------
    ! Compute the Jacobian for the transformation from the reference
    ! coordinates to the true coordinates:
    !
    !                 |                                                                          |
    !                 | sum_n{dphi_n/dxr * xn}   sum_n{dphi_n/dxr * yn}   sum_n{dphi_n/dxr * zn} |
    !   J(xr,yr,zr) = |                                                                          |
    !                 | sum_n{dphi_n/dyr * xn}   sum_n{dphi_n/dyr * yn}   sum_n{dphi_n/dyr * zn} |
    !                 |                                                                          |
    !                 | sum_n{dphi_n/dzr * xn}   sum_n{dphi_n/dzr * yn}   sum_n{dphi_n/dzr * zn} |
    !                 !                                                                          |
    !
    ! where (xn,yn,zn) are the true Cartesian nodal coordinates,
    !       (xr,yr,zr) are the coordinates of the quad point in the reference element,
    !       and sum_n denotes a sum over nodes.
    !------------------------------------------------------------------

    if (verbose_Jac .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
       print*, ' '
       print*, 'In get_basis_function_derivatives_3d: i, j, k, p =', i, j, k, p
    endif

    Jac(:,:) = 0.d0

    do n = 1, nNodesPerElement_3d
       Jac(1,1) = Jac(1,1) + dphi_dxr_3d(n) * xNode(n)
       Jac(1,2) = Jac(1,2) + dphi_dxr_3d(n) * yNode(n)
       Jac(1,3) = Jac(1,3) + dphi_dxr_3d(n) * zNode(n)
       Jac(2,1) = Jac(2,1) + dphi_dyr_3d(n) * xNode(n)
       Jac(2,2) = Jac(2,2) + dphi_dyr_3d(n) * yNode(n)
       Jac(2,3) = Jac(2,3) + dphi_dyr_3d(n) * zNode(n)
       Jac(3,1) = Jac(3,1) + dphi_dzr_3d(n) * xNode(n)
       Jac(3,2) = Jac(3,2) + dphi_dzr_3d(n) * yNode(n)
       Jac(3,3) = Jac(3,3) + dphi_dzr_3d(n) * zNode(n)
    enddo

    !------------------------------------------------------------------
    ! Compute the determinant and inverse of J
    !------------------------------------------------------------------

    cofactor(1,1) =   Jac(2,2)*Jac(3,3) - Jac(2,3)*Jac(3,2)
    cofactor(1,2) = -(Jac(2,1)*Jac(3,3) - Jac(2,3)*Jac(3,1))
    cofactor(1,3) =   Jac(2,1)*Jac(3,2) - Jac(2,2)*Jac(3,1)
    cofactor(2,1) = -(Jac(1,2)*Jac(3,3) - Jac(1,3)*Jac(3,2))
    cofactor(2,2) =   Jac(1,1)*Jac(3,3) - Jac(1,3)*Jac(3,1)
    cofactor(2,3) = -(Jac(1,1)*Jac(3,2) - Jac(1,2)*Jac(3,1))
    cofactor(3,1) =   Jac(1,2)*Jac(2,3) - Jac(1,3)*Jac(2,2)
    cofactor(3,2) = -(Jac(1,1)*Jac(2,3) - Jac(1,3)*Jac(2,1))
    cofactor(3,3) =   Jac(1,1)*Jac(2,2) - Jac(1,2)*Jac(2,1)

    detJ = Jac(1,1)*cofactor(1,1) + Jac(1,2)*cofactor(1,2) + Jac(1,3)*cofactor(1,3)

    if (verbose_Jac .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
       print*, ' '
       print*, 'detJ1:', Jac(1,1)*cofactor(1,1) + Jac(1,2)*cofactor(1,2) + Jac(1,3)*cofactor(1,3)
       print*, 'detJ2:', Jac(2,1)*cofactor(2,1) + Jac(2,2)*cofactor(2,2) + Jac(2,3)*cofactor(2,3)
       print*, 'detJ3:', Jac(3,1)*cofactor(3,1) + Jac(3,2)*cofactor(3,2) + Jac(3,3)*cofactor(3,3)
    endif

    if (abs(detJ) > 0.d0) then
       do col = 1, 3
          do row = 1, 3
             Jinv(row,col) = cofactor(col,row)
          enddo
       enddo
       Jinv(:,:) = Jinv(:,:) / detJ
    else
       print*, 'stopping, det J = 0'
       print*, 'i, j, k, p:', i, j, k, p
       print*, 'Jacobian matrix:'
       print*, Jac(1,:)
       print*, Jac(2,:)
       print*, Jac(3,:) 
       !call write_log('Jacobian matrix is singular', GM_FATAL)
       write(error_unit,*) 'Jacobian matrix is singular'
       stop
    endif

    if (verbose_Jac .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
       print*, ' '
       print*, 'Jacobian calc, p =', p
       print*, 'det J =', detJ
       print*, ' '
       print*, 'Jacobian matrix:'
       print*, Jac(1,:)
       print*, Jac(2,:)
       print*, Jac(3,:)
       print*, ' '
       print*, 'cofactor matrix:'
       print*, cofactor(1,:)
       print*, cofactor(2,:)
       print*, cofactor(3,:)
       print*, ' '
       print*, 'Inverse matrix:'
       print*, Jinv(1,:)
       print*, Jinv(2,:)
       print*, Jinv(3,:)
       print*, ' '
       prod = matmul(Jac, Jinv)
       print*, 'Jac*Jinv:'
       print*, prod(1,:)
       print*, prod(2,:)
       print*, prod(3,:)
    endif

    ! Optional bug check: Verify that J * Jinv = I

    if (Jac_bug_check) then
       prod = matmul(Jac,Jinv)
       do col = 1, 3
          do row = 1, 3
             if (abs(prod(row,col) - identity3(row,col)) > 1.d-11) then
                print*, 'stopping, Jac * Jinv /= identity'
                print*, 'i, j, k, p:', i, j, k, p
                print*, 'Jac*Jinv:'
                print*, prod(1,:)
                print*, prod(2,:)
                print*, prod(3,:)
                !call write_log('Jacobian matrix was not correctly inverted', GM_FATAL)
                write(error_unit,*) 'Jacobian matrix was not correctly inverted'
                stop
             endif
          enddo
       enddo
    endif  ! Jac_bug_check

    !------------------------------------------------------------------
    ! Compute the contribution of this quadrature point to dphi/dx and dphi/dy
    ! for each basis function.
    !
    !   | dphi_n/dx |          | dphi_n/dxr |
    !   |           |          |            | 
    !   | dphi_n/dy | = Jinv * | dphi_n/dyr |
    !   |           |          |            |
    !   | dphi_n/dz |          | dphi_n/dzr |
    !
    !------------------------------------------------------------------

    dphi_dx_3d(:) = 0.d0
    dphi_dy_3d(:) = 0.d0
    dphi_dz_3d(:) = 0.d0

    do n = 1, nNodesPerElement_3d
       dphi_dx_3d(n) = Jinv(1,1)*dphi_dxr_3d(n)  &
                     + Jinv(1,2)*dphi_dyr_3d(n)  &
                     + Jinv(1,3)*dphi_dzr_3d(n)
       dphi_dy_3d(n) = Jinv(2,1)*dphi_dxr_3d(n)  &
                     + Jinv(2,2)*dphi_dyr_3d(n)  &
                     + Jinv(2,3)*dphi_dzr_3d(n)
       dphi_dz_3d(n) = Jinv(3,1)*dphi_dxr_3d(n)  &
                     + Jinv(3,2)*dphi_dyr_3d(n)  &
                     + Jinv(3,3)*dphi_dzr_3d(n)
    enddo

    if (Jac_bug_check) then

       ! Check that the sum of dphi_dx, etc. is close to zero  

       if (abs( sum(dphi_dx_3d)/maxval(dphi_dx_3d) ) > 1.d-11) then
          print*, 'stopping, sum over basis functions of dphi_dx > 0'
          print*, 'dphi_dx_3d =', dphi_dx_3d(:)
          print*, 'sum =', sum(dphi_dx_3d)
          print*, 'i, j, k, p =', i, j, k, p
          !call write_log('Sum over basis functions of dphi_dx /= 0', GM_FATAL)
          write(error_unit,*) 'Sum over basis functions of dphi_dx /= 0'
          stop
       endif

       if (abs( sum(dphi_dy_3d)/maxval(dphi_dy_3d) ) > 1.d-11) then
          print*, 'stopping, sum over basis functions of dphi_dy > 0'
          print*, 'dphi_dy_3d =', dphi_dy_3d(:)
          print*, 'sum =', sum(dphi_dy_3d)
          print*, 'i, j, k, p =', i, j, k, p
          !call write_log('Sum over basis functions of dphi_dy /= 0', GM_FATAL)
          write(error_unit,*) 'Sum over basis functions of dphi_dy /= 0'
          stop
       endif

       if (abs( sum(dphi_dz_3d)/maxval(dphi_dz_3d) ) > 1.d-11) then
          print*, 'stopping, sum over basis functions of dphi_dz > 0'
          print*, 'dphi_dz_3d =', dphi_dz_3d(:)
          print*, 'sum =', sum(dphi_dz_3d)
          print*, 'i, j, k, p =', i, j, k, p
          !call write_log('Sum over basis functions of dphi_dz /= 0', GM_FATAL)
          write(error_unit,*) 'Sum over basis functions of dphi_dz /= 0'
          stop
       endif

    endif  ! Jac_bug_check

  end subroutine get_basis_function_derivatives_3d

!****************************************************************************

  subroutine get_basis_function_derivatives_2d(xNode,       yNode,         &
                                               dphi_dxr_2d, dphi_dyr_2d,   &
                                               dphi_dx_2d,  dphi_dy_2d,    &
                                               detJ,                       &
                                               itest, jtest, rtest,        &
                                               i, j, p)

    !------------------------------------------------------------------
    ! Evaluate the x and y derivatives of 2D element basis functions
    ! at a particular quadrature point.
    !
    ! Also determine the Jacobian of the transformation between the
    ! reference element and the true element.
    ! 
    ! This subroutine should work for any 2D element with any number of nodes.
    !------------------------------------------------------------------

    real(dp), dimension(nNodesPerElement_2d), intent(in) :: &
       xNode, yNode,                   &! nodal coordinates
       dphi_dxr_2d, dphi_dyr_2d         ! derivatives of basis functions at quad pt
                                        !  wrt x and y in reference element

    real(dp), dimension(nNodesPerElement_2d), intent(out) :: &
       dphi_dx_2d, dphi_dy_2d           ! derivatives of basis functions at quad pt
                                        !  wrt x and y in true Cartesian coordinates  

    real(dp), intent(out) :: &
                detJ      ! determinant of Jacobian matrix

    real(dp), dimension(2,2) ::  &
                Jac,      &! Jacobian matrix
                Jinv       ! inverse Jacobian matrix

    integer, intent(in) :: &
       itest, jtest, rtest              ! coordinates of diagnostic point

    integer, intent(in) :: i, j, p

    integer :: n, row, col

    logical, parameter :: Jac_bug_check = .false.   ! set to true for debugging
    real(dp), dimension(2,2) :: prod     ! Jac * Jinv (should be identity matrix)

    !------------------------------------------------------------------
    ! Compute the Jacobian for the transformation from the reference
    ! coordinates to the true coordinates:
    !
    !              |                                                  |
    !              | sum_n{dphi_n/dxr * xn}   sum_n{dphi_n/dxr * yn}  |
    !   J(xr,yr) = |                                                  |
    !              | sum_n{dphi_n/dyr * xn}   sum_n{dphi_n/dyr * yn}  |
    !              |                                                  |
    !
    ! where (xn,yn) are the true Cartesian nodal coordinates,
    !       (xr,yr) are the coordinates of the quad point in the reference element,
    !       and sum_n denotes a sum over nodes.
    !------------------------------------------------------------------

    Jac(:,:) = 0.d0

    if (verbose_Jac .and. this_rank==rtest .and. i==itest .and. j==jtest) then
       print*, ' '
       print*, 'In get_basis_function_derivatives_2d: i, j, p =', i, j, p
    endif

    do n = 1, nNodesPerElement_2d
       if (verbose_Jac .and. this_rank==rtest .and. i==itest .and. j==jtest) then
          print*, ' '
          print*, 'n, x, y:', n, xNode(n), yNode(n)
          print*, 'dphi_dxr_2d, dphi_dyr_2d:', dphi_dxr_2d(n), dphi_dyr_2d(n)
       endif
       Jac(1,1) = Jac(1,1) + dphi_dxr_2d(n) * xNode(n)
       Jac(1,2) = Jac(1,2) + dphi_dxr_2d(n) * yNode(n)
       Jac(2,1) = Jac(2,1) + dphi_dyr_2d(n) * xNode(n)
       Jac(2,2) = Jac(2,2) + dphi_dyr_2d(n) * yNode(n)
    enddo

    !------------------------------------------------------------------
    ! Compute the determinant and inverse of J
    !------------------------------------------------------------------

    detJ = Jac(1,1)*Jac(2,2) - Jac(1,2)*Jac(2,1)

    if (abs(detJ) > 0.d0) then
       Jinv(1,1) =  Jac(2,2)/detJ
       Jinv(1,2) = -Jac(1,2)/detJ
       Jinv(2,1) = -Jac(2,1)/detJ
       Jinv(2,2) =  Jac(1,1)/detJ
    else
       print*, 'stopping, det J = 0'
       print*, 'i, j, p:', i, j, p
       print*, 'Jacobian matrix:'
       print*, Jac(1,:)
       print*, Jac(2,:)
       !call write_log('Jacobian matrix is singular', GM_FATAL)
       write(error_unit,*) 'Jacobian matrix is singular'
       stop
    endif

    if (verbose_Jac .and. this_rank==rtest .and. i==itest .and. j==jtest) then
       print*, ' '
       print*, 'Jacobian calc, p =', p
       print*, 'det J =', detJ
       print*, ' '
       print*, 'Jacobian matrix:'
       print*, Jac(1,:)
       print*, Jac(2,:)
       print*, ' '
       print*, 'Inverse matrix:'
       print*, Jinv(1,:)
       print*, Jinv(2,:)
       print*, ' '
       prod = matmul(Jac, Jinv)
       print*, 'Jac*Jinv:'
       print*, prod(1,:)
       print*, prod(2,:)
    endif

    ! Optional bug check - Verify that J * Jinv = I

    if (Jac_bug_check) then
       prod = matmul(Jac,Jinv)
       do col = 1, 2
          do row = 1, 2
             if (abs(prod(row,col) - identity3(row,col)) > 1.d-12) then
                print*, 'stopping, Jac * Jinv /= identity'
                print*, 'i, j, p:', i, j, p
                print*, 'Jac*Jinv:'
                print*, prod(1,:)
                print*, prod(2,:)
                !call write_log('Jacobian matrix was not correctly inverted', GM_FATAL)
                write(error_unit,*) 'Jacobian matrix was not correctly inverted'
                stop
             endif
          enddo
       enddo
    endif

    !------------------------------------------------------------------
    ! Compute the contribution of this quadrature point to dphi/dx and dphi/dy
    ! for each basis function.
    !
    !   | dphi_n/dx |          | dphi_n/dxr |
    !   |           | = Jinv * |            |
    !   | dphi_n/dy |          | dphi_n/dyr |
    !
    !------------------------------------------------------------------

    dphi_dx_2d(:) = 0.d0
    dphi_dy_2d(:) = 0.d0

    do n = 1, nNodesPerElement_2d
       dphi_dx_2d(n) = dphi_dx_2d(n) + Jinv(1,1)*dphi_dxr_2d(n)  &
                                     + Jinv(1,2)*dphi_dyr_2d(n)
       dphi_dy_2d(n) = dphi_dy_2d(n) + Jinv(2,1)*dphi_dxr_2d(n)  &
                                     + Jinv(2,2)*dphi_dyr_2d(n)
    enddo

    if (Jac_bug_check) then

       ! Check that the sum of dphi_dx, etc. is close to zero  
       if (abs( sum(dphi_dx_2d)/maxval(dphi_dx_2d) ) > 1.d-11) then
          print*, 'stopping, sum over basis functions of dphi_dx > 0'
          print*, 'dphi_dx_2d =', dphi_dx_2d(:)
          print*, 'i, j, p =', i, j, p
          !call write_log('Sum over basis functions of dphi_dx /= 0', GM_FATAL)
          write(error_unit,*) 'Sum over basis functions of dphi_dx /= 0'
          stop
       endif

       if (abs( sum(dphi_dy_2d)/maxval(dphi_dy_2d) ) > 1.d-11) then
          print*, 'stopping, sum over basis functions of dphi_dy > 0'
          print*, 'dphi_dy =', dphi_dy_2d(:)
          print*, 'i, j, p =', i, j, p
          !call write_log('Sum over basis functions of dphi_dy /= 0', GM_FATAL)
          write(error_unit,*) 'Sum over basis functions of dphi_dy /= 0'
          stop
       endif

    endif

  end subroutine get_basis_function_derivatives_2d

!****************************************************************************

end module gaussian_quadrature