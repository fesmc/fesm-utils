module gaussian_quadrature

    use, intrinsic :: iso_fortran_env, only : input_unit, output_unit, error_unit    
    use precision

    implicit none
    
    logical :: verbose_init = .true. 
    logical, parameter :: check_symmetry = .true.   ! if true, then check symmetry of assembled matrix

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


    private
    public :: gq2D_class
    public :: gq2D_init
    public :: gq2D_to_nodes

    public :: gq3D_class
    public :: gq3D_init

    public :: gaussian_quadrature_init

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


! == Directly from CISM2.1 ==

  subroutine gaussian_quadrature_init

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

end module gaussian_quadrature