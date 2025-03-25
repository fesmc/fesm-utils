module gaussian_quadrature

    use, intrinsic :: iso_fortran_env, only : input_unit, output_unit, error_unit    
    use precision

    implicit none
    
    logical :: verbose_init_global = .false. 
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
        integer :: n_qp
        integer :: n_nodes
        real(8) :: xr(4)
        real(8) :: yr(4)
        real(8) :: wt(4)
        real(8) :: wt_tot
        real(8) :: N(4,4)
        real(8) :: dNdxr(4,4)
        real(8) :: dNdyr(4,4)
        
        real(8) :: dNdx(4)  ! At nodes
        real(8) :: dNdy(4)  ! At nodes
        real(8) :: detJ(4)  ! At nodes
        real(8) :: v(4)     ! variable value at quadrature points
    end type

    type gq3D_class
        integer :: n_qp
        integer :: n_nodes
        real(8) :: xr(8)
        real(8) :: yr(8)
        real(8) :: zr(8)
        real(8) :: wt(8)
        real(8) :: wt_tot
        real(8) :: N(8,8)
        real(8) :: dNdxr(8,8)
        real(8) :: dNdyr(8,8)
        real(8) :: dNdzr(8,8)
        real(8) :: vol0
        
        real(8) :: dNdx(8)  ! At nodes
        real(8) :: dNdy(8)  ! At nodes
        real(8) :: dNdz(8)  ! At nodes
        real(8) :: detJ(8)  ! At nodes
        real(8) :: v(8)     ! variable value at quadrature points
    end type


    private
    public :: gq2D_class
    public :: gq2D_init
    public :: gq2D_to_nodes

    public :: gq3D_class
    public :: gq3D_init
    public :: gq3D_to_nodes
    
contains

    subroutine gq2D_init(gq,verbose)

        implicit none

        type(gq2D_class), intent(OUT) :: gq
        logical, intent(IN), optional :: verbose 

        ! Local variables
        integer :: p, n
        logical :: verbose_init

        ! Quadrature point locations in reference element (-1,1)
        real(8), parameter :: sqrt3_inv = 1.0d0 / sqrt(3.0d0)

        if (present(verbose)) then
            verbose_init = verbose
        else
            verbose_init = verbose_init_global
        end if

        ! Define how many quadrature points and nodes (corners) of cell
        gq%n_qp    = 4
        gq%n_nodes = 4

        ! Define quadrature points
        gq%xr = [ -sqrt3_inv, sqrt3_inv, sqrt3_inv, -sqrt3_inv ]
        gq%yr = [ -sqrt3_inv, -sqrt3_inv, sqrt3_inv, sqrt3_inv ]

        ! Define weights
        gq%wt = [1.0d0, 1.0d0, 1.0d0, 1.0d0]
        gq%wt_tot = 4.0d0           ! Surface area of square [-1:1,-1:1]=> 2x2 => 4 

        ! Define shape functions N(n, p) and derivatives
        ! where p represents each quadrature point
        do p = 1, gq%n_qp

            gq%N(1,p) = (1.d0 - gq%xr(p)) * (1.d0 - gq%yr(p)) / 4.d0  ! N1
            gq%N(2,p) = (1.d0 + gq%xr(p)) * (1.d0 - gq%yr(p)) / 4.d0  ! N2
            gq%N(3,p) = (1.d0 + gq%xr(p)) * (1.d0 + gq%yr(p)) / 4.d0  ! N3
            gq%N(4,p) = (1.d0 - gq%xr(p)) * (1.d0 + gq%yr(p)) / 4.d0  ! N4

            gq%dNdxr(1,p) = -(1.d0 - gq%yr(p)) / 4.d0 
            gq%dNdxr(2,p) =  (1.d0 - gq%yr(p)) / 4.d0 
            gq%dNdxr(3,p) =  (1.d0 + gq%yr(p)) / 4.d0 
            gq%dNdxr(4,p) = -(1.d0 + gq%yr(p)) / 4.d0

            gq%dNdyr(1,p) = -(1.d0 - gq%xr(p)) / 4.d0 
            gq%dNdyr(2,p) = -(1.d0 + gq%xr(p)) / 4.d0 
            gq%dNdyr(3,p) =  (1.d0 + gq%xr(p)) / 4.d0 
            gq%dNdyr(4,p) =  (1.d0 - gq%xr(p)) / 4.d0 

            if (verbose_init) then
                write(*,*) " "
                write(*,*) "Quad point, p =", p
                write(*,*) "n, N, dNdxr, dNdyr:"
                do n = 1, gq%n_nodes
                    write(*,*) n, gq%N(n, p), gq%dNdxr(n, p), gq%dNdyr(n, p)
                end do
                write(*,*) "sum(N)", sum(gq%N(:, p))            ! Verified that sum = 1
                write(*,*) "sum(dN/dxr)", sum(gq%dNdxr(:, p))   ! Verified that sum = 0 (within roundoff)
                write(*,*) "sum(dN/dyr)", sum(gq%dNdyr(:, p))   ! Verified that sum = 0 (within roundoff)
            endif

            ! Set node calculation values to zero to start
            gq%dNdx = 0.0
            gq%dNdy = 0.0
            gq%detJ = 0.0
            gq%v = 0.0

        end do

        return

    end subroutine gq2D_init

    subroutine gq3D_init(gq,verbose)

        implicit none

        type(gq3D_class), intent(OUT) :: gq
        logical, intent(IN), optional :: verbose

        ! Local variables
        integer :: p, n
        logical :: verbose_init

        ! Quadrature point locations in reference element (-1,1)
        real(8), parameter :: sqrt3_inv = 1.0d0 / sqrt(3.0d0)

        if (present(verbose)) then
            verbose_init = verbose
        else
            verbose_init = verbose_init_global
        end if

        ! Define how many quadrature points and nodes (corners) of cell
        gq%n_qp    = 8
        gq%n_nodes = 8

        ! Define quadrature points
        gq%xr = [ -sqrt3_inv,  sqrt3_inv,  sqrt3_inv, -sqrt3_inv, -sqrt3_inv,  sqrt3_inv,  sqrt3_inv, -sqrt3_inv ]
        gq%yr = [ -sqrt3_inv, -sqrt3_inv,  sqrt3_inv,  sqrt3_inv, -sqrt3_inv, -sqrt3_inv,  sqrt3_inv,  sqrt3_inv ]
        gq%zr = [ -sqrt3_inv, -sqrt3_inv, -sqrt3_inv, -sqrt3_inv,  sqrt3_inv,  sqrt3_inv,  sqrt3_inv,  sqrt3_inv ]

        ! Define weights
        gq%wt = [1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0]
        gq%wt_tot = 8.0d0       ! Volume of square [-1:1,-1:1, -1:1]=> 2x2x2 => 8
        
        ! Set volume scale
        ! This is not strictly necessary, but dividing by this scale gives matrix coefficients 
        !  that are ~1.
        gq%vol0  = 1.0d9    ! volume scale (m^3)

        ! Define shape functions N(n, p) and derivatives
        ! where p represents each quadrature point
        do p = 1, gq%n_qp

            gq%N(1,p) = (1.d0 - gq%xr(p)) * (1.d0 - gq%yr(p)) * (1.d0 - gq%zr(p)) / 8.d0  ! N1
            gq%N(2,p) = (1.d0 + gq%xr(p)) * (1.d0 - gq%yr(p)) * (1.d0 - gq%zr(p)) / 8.d0  ! N2
            gq%N(3,p) = (1.d0 + gq%xr(p)) * (1.d0 + gq%yr(p)) * (1.d0 - gq%zr(p)) / 8.d0  ! N3
            gq%N(4,p) = (1.d0 - gq%xr(p)) * (1.d0 + gq%yr(p)) * (1.d0 - gq%zr(p)) / 8.d0  ! N4
            gq%N(5,p) = (1.d0 - gq%xr(p)) * (1.d0 - gq%yr(p)) * (1.d0 + gq%zr(p)) / 8.d0  ! N5
            gq%N(6,p) = (1.d0 + gq%xr(p)) * (1.d0 - gq%yr(p)) * (1.d0 + gq%zr(p)) / 8.d0  ! N6
            gq%N(7,p) = (1.d0 + gq%xr(p)) * (1.d0 + gq%yr(p)) * (1.d0 + gq%zr(p)) / 8.d0  ! N7
            gq%N(8,p) = (1.d0 - gq%xr(p)) * (1.d0 + gq%yr(p)) * (1.d0 + gq%zr(p)) / 8.d0  ! N8

            gq%dNdxr(1,p) = -(1.d0 - gq%yr(p)) * (1.d0 - gq%zr(p)) / 8.d0 
            gq%dNdxr(2,p) =  (1.d0 - gq%yr(p)) * (1.d0 - gq%zr(p)) / 8.d0 
            gq%dNdxr(3,p) =  (1.d0 + gq%yr(p)) * (1.d0 - gq%zr(p)) / 8.d0 
            gq%dNdxr(4,p) = -(1.d0 + gq%yr(p)) * (1.d0 - gq%zr(p)) / 8.d0 
            gq%dNdxr(5,p) = -(1.d0 - gq%yr(p)) * (1.d0 + gq%zr(p)) / 8.d0 
            gq%dNdxr(6,p) =  (1.d0 - gq%yr(p)) * (1.d0 + gq%zr(p)) / 8.d0 
            gq%dNdxr(7,p) =  (1.d0 + gq%yr(p)) * (1.d0 + gq%zr(p)) / 8.d0 
            gq%dNdxr(8,p) = -(1.d0 + gq%yr(p)) * (1.d0 + gq%zr(p)) / 8.d0 

            gq%dNdyr(1,p) = -(1.d0 - gq%xr(p)) * (1.d0 - gq%zr(p)) / 8.d0 
            gq%dNdyr(2,p) = -(1.d0 + gq%xr(p)) * (1.d0 - gq%zr(p)) / 8.d0 
            gq%dNdyr(3,p) =  (1.d0 + gq%xr(p)) * (1.d0 - gq%zr(p)) / 8.d0 
            gq%dNdyr(4,p) =  (1.d0 - gq%xr(p)) * (1.d0 - gq%zr(p)) / 8.d0 
            gq%dNdyr(5,p) = -(1.d0 - gq%xr(p)) * (1.d0 + gq%zr(p)) / 8.d0 
            gq%dNdyr(6,p) = -(1.d0 + gq%xr(p)) * (1.d0 + gq%zr(p)) / 8.d0 
            gq%dNdyr(7,p) =  (1.d0 + gq%xr(p)) * (1.d0 + gq%zr(p)) / 8.d0 
            gq%dNdyr(8,p) =  (1.d0 - gq%xr(p)) * (1.d0 + gq%zr(p)) / 8.d0 

            gq%dNdzr(1,p) = -(1.d0 - gq%xr(p)) * (1.d0 - gq%yr(p)) / 8.d0 
            gq%dNdzr(2,p) = -(1.d0 + gq%xr(p)) * (1.d0 - gq%yr(p)) / 8.d0 
            gq%dNdzr(3,p) = -(1.d0 + gq%xr(p)) * (1.d0 + gq%yr(p)) / 8.d0 
            gq%dNdzr(4,p) = -(1.d0 - gq%xr(p)) * (1.d0 + gq%yr(p)) / 8.d0 
            gq%dNdzr(5,p) =  (1.d0 - gq%xr(p)) * (1.d0 - gq%yr(p)) / 8.d0 
            gq%dNdzr(6,p) =  (1.d0 + gq%xr(p)) * (1.d0 - gq%yr(p)) / 8.d0 
            gq%dNdzr(7,p) =  (1.d0 + gq%xr(p)) * (1.d0 + gq%yr(p)) / 8.d0 
            gq%dNdzr(8,p) =  (1.d0 - gq%xr(p)) * (1.d0 + gq%yr(p)) / 8.d0 

            if (verbose_init) then
                write(*,*) " "
                write(*,*) "Quad point, p =", p
                write(*,*) "n, N_3d, dNdxr, dNdyr, dNdzr:"
                do n = 1, gq%n_nodes
                    write(*,*) n, gq%N(n,p), gq%dNdxr(n,p), gq%dNdyr(n,p), gq%dNdzr(n,p)
                enddo
                write(*,*) " "
                write(*,*) "sum(N)", sum(gq%N(:,p))             ! verified that sum = 1
                write(*,*) "sum(dN/dxr)", sum(gq%dNdxr(:,p))    ! verified that sum = 0 (within roundoff)
                write(*,*) "sum(dN/dyr)", sum(gq%dNdyr(:,p))    ! verified that sum = 0 (within roundoff)
                write(*,*) "sum(dN/dzr)", sum(gq%dNdzr(:,p))    ! verified that sum = 0 (within roundoff)
            end if

            ! Set node calculation values to zero to start
            gq%dNdx = 0.0
            gq%dNdy = 0.0
            gq%dNdz = 0.0
            gq%detJ = 0.0
            gq%v = 0.0
            
        end do

        return

    end subroutine gq3D_init

    subroutine gq2D_to_nodes(gq, v_qp, var, dx, dy, grid_type, i, j, im1, ip1, jm1, jp1)
        
        implicit none
        
        type(gq2D_class), intent(INOUT) :: gq       ! Gaussian Quadrature 2D object
        real(wp),         intent(INOUT) :: v_qp(4)  ! Variable values at quadrature points
        real(wp), intent(in)  :: var(:,:)           ! Variable to be interpolated
        real(wp), intent(IN)  :: dx 
        real(wp), intent(IN)  :: dy
        character(len=*), intent(IN) :: grid_type   ! "aa", "ab", "acx", "acy"
        integer,  intent(IN)  :: i, j               ! [x,y] indices of current cell
        integer,  intent(IN)  :: im1, ip1, jm1, jp1 ! Neighbor indices of current cell
        
        ! Local variables
        integer :: nx, ny, p, n
        real(8) :: x(4)                 ! Real x-coordinates at the four corners of the cell
        real(8) :: y(4)                 ! Real y-coordinates at the four corners of the cell
        real(8) :: v(4)                 ! Values of u at the four corners of the cell
        real(8) :: vx(4)                ! Derivatives du/dx at the four corners of the cell
        real(8) :: vy(4)                ! Derivatives du/dy at the four corners of the cell
        
        nx = size(var,1)
        ny = size(var,2)

        ! Step 1: determine x and y values of input array values

        ! Map xc,yc onto x,y vectors. Account for whether input
        ! variable is on aa, acx or acy grid. Account for boundary
        ! conditions (periodic, etc).

        !----------------------------------------------------------------
        ! Bilinear basis set for reference square, x=(-1,1), y=(-1,1)             
        ! Indexrng is counter-clockwise from SW corner
        ! The code uses "gq%N" to denote these basis functions. 
        !
        ! N1 = (1-x)*(1-y)/4             N4----N3
        ! N2 = (1+x)*(1-y)/4             |     |
        ! N3 = (1+x)*(1+y)/4             |     |
        ! N4 = (1-x)*(1+y)/4             N1----N2
        !----------------------------------------------------------------

        ! First get x and y values (xx and yy defined on aa-nodes)
        ! Note: these could be actual coordinates passed in as arguments,
        ! as done in e.g. CISM. However, we can assume the horizontal grid
        ! is regular with constant spacing and define relative x/y coordinates
        ! here.

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
                v(1) = 0.25d0 * (var(im1, jm1) + var(i, jm1) + var(im1, j) + var(i, j))  ! Bottom-left
                v(2) = 0.25d0 * (var(i, jm1) + var(ip1, jm1) + var(i, j) + var(ip1, j))  ! Bottom-right
                v(3) = 0.25d0 * (var(i, j) + var(ip1, j) + var(i, jp1) + var(ip1, jp1))  ! Top-right
                v(4) = 0.25d0 * (var(im1, j) + var(i, j) + var(im1, jp1) + var(i, jp1))  ! Top-left
            case("ab")
                v(1) = var(im1,jm1)         ! Bottom-left
                v(2) = var(i,jm1)           ! Bottom-right
                v(3) = var(i,j)             ! Top-right
                v(4) = var(im1,j)           ! Top-left
            case("acx")
                v(1) = 0.5d0 * (var(im1, jm1) + var(im1, j))      ! Bottom-left
                v(2) = 0.5d0 * (var(i, jm1) + var(i, j))          ! Bottom-right
                v(3) = 0.5d0 * (var(i, j) + var(i, jp1))          ! Top-right
                v(4) = 0.5d0 * (var(im1, j) + var(im1, jp1))      ! Top-left
            case("acy")
                v(1) = 0.5d0 * (var(im1, jm1) + var(i, jm1))      ! Bottom-left
                v(2) = 0.5d0 * (var(i, jm1) + var(ip1, jm1))      ! Bottom-right
                v(3) = 0.5d0 * (var(i, j) + var(ip1, j))          ! Top-right
                v(4) = 0.5d0 * (var(im1, j) + var(i, j))          ! Top-left
            case DEFAULT
                write(error_unit,*) "gq2D_to_nodes:: Error: grid_type not recognized."
                write(error_unit,*) "grid_type = ", trim(grid_type)
                stop
        end select

        ! Step 2: Calculate values at each quadrature point

        ! Loop over quadrature points for this element
        do p = 1, gq%n_qp

            ! Compute basis function derivatives and det(J) for this quadrature point
            ! For now, pass in i, j, k, p for debugging
            !TODO - Modify this subroutine so that the output derivatives are optional?

            ! call get_basis_function_derivatives_2d(x(:),             y(:),          &
            !                                         gq%dNdxr(:,p), gq%dNdyr(:,p),   &
            !                                         gq%dNdx(:),    gq%dNdy(:),      &
            !                                         gq%detJ(p),                     &
            !                                         itest, jtest, rtest,            &
            !                                         i, j, p)

            ! Evaluate var at this quadrature point, taking a N-weighted sum over neighboring vertices.
            gq%v(p) = 0.d0
            do n = 1, gq%n_nodes
                gq%v(p) = gq%v(p) + gq%N(n,p) * v(n)
            end do
        
        end do

        ! Store in output variable too
        v_qp = gq%v

        return
        
    end subroutine gq2D_to_nodes

    subroutine gq3D_to_nodes(gq, v_qp, var, dx, dy, dz0, dz1, grid_type, i, j, k, im1, ip1, jm1, jp1, km1, kp1)
        
        implicit none
        
        type(gq3D_class), intent(INOUT) :: gq       ! Gaussian Quadrature 3D object
        real(wp),         intent(INOUT) :: v_qp(8)  ! Variable values at quadrature points
        real(wp), intent(in)  :: var(:,:,:)         ! Variable to be interpolated
        real(wp), intent(IN)  :: dx                 ! Horizontal grid spacing (const)
        real(wp), intent(IN)  :: dy                 ! Horizontal grid spacing (const)
        real(wp), intent(IN)  :: dz0                ! dz to cell below
        real(wp), intent(IN)  :: dz1                ! dz to cell above
        character(len=*), intent(IN) :: grid_type   ! "aa", "ab", "acx", "acy", "acz"
        integer,  intent(IN)  :: i, j, k            ! [x,y] indices of current cell
        integer,  intent(IN)  :: im1, ip1, jm1, jp1 ! Neighbor indices of current cell
        integer,  intent(IN)  :: km1, kp1

        ! Local variables
        integer :: nx, ny, nz, p, n
        real(8) :: x(8)                         ! Real x-coordinates at the four corners of the cell
        real(8) :: y(8)                         ! Real y-coordinates at the four corners of the cell
        real(8) :: z(8)                         ! Real z-coordinates at the four corners of the cell
        real(8) :: v(8)                         ! Values of var at the eight corners of the cell
        real(8) :: vx(8), vy(8), vz(8)          ! Derivatives at the eight corners
        
        nx = size(var,1)
        ny = size(var,2)
        nz = size(var,3)

        ! Step 1: determine x, y and z values of input array values

        ! Map xc,yc,zc onto x,y,z vectors. Account for whether input
        ! variable is on aa, acx or acy grid. Account for boundary
        ! conditions (periodic, etc).

        !----------------------------------------------------------------
        ! Trilinear basis set for reference hexahedron, x=(-1,1), y=(-1,1), z=(-1,1)             
        ! Indexrng is counter-clockwise from SW corner, with 1-4 on lower surface
        !  and 5-8 on upper surface
        ! The code uses "gq%N" to denote these basis functions. 
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
        
        ! First get x and y values (xx and yy defined on aa-nodes)
        ! Note: these could be actual coordinates passed in as arguments,
        ! as done in e.g. CISM. However, we can assume the horizontal grid
        ! is regular with constant spacing and define relative x/y coordinates
        ! here. 
        ! The vertical axis z may contain different spacing above and below this cell,
        ! which is why a dz0 (distance to cell below) and dz1 (distance to cell above)
        ! is used.

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
        
        ! Compute values of u at the eight cell corners
        select case(trim(grid_type))
            case("aa")
                v(1) = 0.25d0 * ( 0.5d0 * (var(im1, jm1, km1) + var(im1, jm1, k))       &
                                + 0.5d0 * (var(im1, j, km1)   + var(im1, j, k))         &
                                + 0.5d0 * (var(i, j, km1)     + var(i, j, k))           &
                                + 0.5d0 * (var(i, jm1, km1)   + var(i, jm1, k)) )

                v(2) = 0.25d0 * ( 0.5d0 * (var(i, jm1, km1)   + var(i, jm1, k))         &
                                + 0.5d0 * (var(i, j, km1)     + var(i, j, k))           &
                                + 0.5d0 * (var(ip1, j, km1)   + var(ip1, j, k))         &
                                + 0.5d0 * (var(ip1, jm1, km1) + var(ip1, jm1, k)) )

                v(3) = 0.25d0 * ( 0.5d0 * (var(i, j, km1)     + var(i, j, k))           &
                                + 0.5d0 * (var(i, jp1, km1)   + var(i, jp1, k))         &
                                + 0.5d0 * (var(ip1, jp1, km1) + var(ip1, jp1, k))       &
                                + 0.5d0 * (var(ip1, j, km1)   + var(ip1, j, k)) )

                v(4) = 0.25d0 * ( 0.5d0 * (var(im1, j, km1)   + var(im1, j, k))         &
                                + 0.5d0 * (var(im1, jp1, km1) + var(im1, jp1, k))       &
                                + 0.5d0 * (var(i, jp1, km1)   + var(i, jp1, k))         &
                                + 0.5d0 * (var(i, j, km1)     + var(i, j, k)) )

                v(5) = 0.25d0 * ( 0.5d0 * (var(im1, jm1, k)   + var(im1, jm1, kp1))     &
                                + 0.5d0 * (var(im1, j, k)     + var(im1, j, kp1))       &
                                + 0.5d0 * (var(i, j, k)       + var(i, j, kp1))         &
                                + 0.5d0 * (var(i, jm1, k)     + var(i, jm1, kp1)) )

                v(6) = 0.25d0 * ( 0.5d0 * (var(i, jm1, k)     + var(i, jm1, kp1))       &
                                + 0.5d0 * (var(i, j, k)       + var(i, j, kp1))         &
                                + 0.5d0 * (var(ip1, j, k)     + var(ip1, j, kp1))       &
                                + 0.5d0 * (var(ip1, jm1, k)   + var(ip1, jm1, kp1)) )

                v(7) = 0.25d0 * ( 0.5d0 * (var(i, j, k)       + var(i, j, kp1))         &
                                + 0.5d0 * (var(i, jp1, k)     + var(i, jp1, kp1))       &
                                + 0.5d0 * (var(ip1, jp1, k)   + var(ip1, jp1, kp1))     &
                                + 0.5d0 * (var(ip1, j, k)     + var(ip1, j, kp1)) )

                v(8) = 0.25d0 * ( 0.5d0 * (var(im1, j, k)     + var(im1, j, kp1))       &
                                + 0.5d0 * (var(im1, jp1, k)   + var(im1, jp1, kp1))     &
                                + 0.5d0 * (var(i, jp1, k)     + var(i, jp1, kp1))       &
                                + 0.5d0 * (var(i, j, k)       + var(i, j, kp1)) )
            case("ab")
                v(1) = var(im1, jm1, km1)
                v(2) = var(i, jm1, km1)
                v(3) = var(i, j, km1)
                v(4) = var(im1, j, km1)
                v(5) = var(im1, jm1, k)
                v(6) = var(i, jm1, k)
                v(7) = var(i, j, k)
                v(8) = var(im1, j, k)
            case("acx")
                v(1) = 0.25d0 * ( var(im1, jm1, km1)+ var(im1, jm1, k)     &
                                + var(im1, j, km1)  + var(im1, j, k) )
                v(2) = 0.25d0 * ( var(i, jm1, km1)  + var(i, jm1, k)         &
                                + var(i, j, km1)    + var(i, j, k) )
                v(3) = 0.25d0 * ( var(i, j, km1)    + var(i, j, k)           &
                                + var(i, jp1, km1)  + var(i, jp1, k) )
                v(4) = 0.25d0 * ( var(im1, j, km1)  + var(im1, j, k)       &
                                + var(im1, jp1, km1)+ var(im1, jp1, k) )

                v(5) = 0.25d0 * ( var(im1, jm1, k)  + var(im1, jm1, kp1)     &
                                + var(im1, j, k)    + var(im1, j, kp1) )
                v(6) = 0.25d0 * ( var(i, jm1, k)    + var(i, jm1, kp1)         &
                                + var(i, j, k)      + var(i, j, kp1) )
                v(7) = 0.25d0 * ( var(i, j, k)      + var(i, j, kp1)           &
                                + var(i, jp1, k)    + var(i, jp1, kp1) )
                v(8) = 0.25d0 * ( var(im1, j, k)    + var(im1, j, kp1)       &
                                + var(im1, jp1, k)  + var(im1, jp1, kp1) )
            case("acy")
                v(1) = 0.25d0 * ( var(im1, jm1, km1)+ var(i, jm1, km1)     &
                                + var(im1, jm1, k)  + var(i, jm1, k) )
                v(2) = 0.25d0 * ( var(i, jm1, km1)  + var(ip1, jm1, km1)     &
                                + var(i, jm1, k)    + var(ip1, jm1, k) )
                v(3) = 0.25d0 * ( var(i, j, km1)    + var(ip1, j, km1)         &
                                + var(i, j, k)      + var(ip1, j, k) )
                v(4) = 0.25d0 * ( var(im1, j, km1)  + var(i, j, km1)         &
                                + var(im1, j, k)    + var(i, j, k) )

                v(5) = 0.25d0 * ( var(im1, jm1, k)  + var(i, jm1, k)         &
                                + var(im1, jm1, kp1)+ var(i, jm1, kp1) )
                v(6) = 0.25d0 * ( var(i, jm1, k)    + var(ip1, jm1, k)         &
                                + var(i, jm1, kp1)  + var(ip1, jm1, kp1) )
                v(7) = 0.25d0 * ( var(i, j, k)      + var(ip1, j, k)             &
                                + var(i, j, kp1)    + var(ip1, j, kp1) )
                v(8) = 0.25d0 * ( var(im1, j, k)    + var(i, j, k)             &
                                + var(im1, j, kp1)  + var(i, j, kp1) )
            case("acz")
                v(1) = 0.25d0 * ( var(im1, jm1, km1)+ var(i, jm1, km1)   &
                                + var(i, j, km1)    + var(im1, j, km1) )
                v(2) = 0.25d0 * ( var(i, jm1, km1)  + var(ip1, jm1, km1)     &
                                + var(ip1, j, km1)  + var(i, j, km1) )
                v(3) = 0.25d0 * ( var(i, j, km1)    + var(ip1, j, km1)       &
                                + var(ip1, jp1, km1)+ var(i, jp1, km1) )
                v(4) = 0.25d0 * ( var(im1, j, km1)  + var(i, j, km1)     &
                                + var(i, jp1, km1)  + var(im1, jp1, km1) )
                                            
                v(5) = 0.25d0 * ( var(im1, jm1, k)  + var(i, jm1, k)   &
                                + var(i, j, k)      + var(im1, j, k) )
                v(6) = 0.25d0 * ( var(i, jm1, k)    + var(ip1, jm1, k)     &
                                + var(ip1, j, k)    + var(i, j, k) )
                v(7) = 0.25d0 * ( var(i, j, k)      + var(ip1, j, k)       &
                                + var(ip1, jp1, k)  + var(i, jp1, k) )
                v(8) = 0.25d0 * ( var(im1, j, k)    + var(i, j, k)     &
                                + var(i, jp1, k)    + var(im1, jp1, k) )

            case DEFAULT
                write(error_unit,*) "gq3D_to_nodes:: Error: grid_type not recognized."
                write(error_unit,*) "grid_type = ", trim(grid_type)
                stop
        end select

        ! Loop over quadrature points for this element
        do p = 1, gq%n_qp

            ! Compute basis function derivatives and det(J) for this quadrature point
            ! For now, pass in i, j, k, p for debugging
            !TODO - Modify this subroutine so that the output derivatives are optional?

            ! call get_basis_function_derivatives_3d(x(:),             y(:),          z(:),           &
            !                                         gq%dNdxr(:,p), gq%dNdyr(:,p), gq%dNdzr(:,p),    &
            !                                         gq%dNdx(:),    gq%dNdy(:),    gq%dNdz(:),       &
            !                                         gq%detJ(p),                                     &
            !                                         itest, jtest, rtest,                            &
            !                                         i, j, k, p)

            ! Evaluate var at this quadrature point, taking a N-weighted sum over neighboring vertices.
            gq%v(p) = 0.d0
            do n = 1, gq%n_nodes
                gq%v(p) = gq%v(p) + gq%N(n,p) * v(n)
            end do
        
        end do
        
        ! Store in output variable too
        v_qp = gq%v
        
        return
        
    end subroutine gq3D_to_nodes

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

end module gaussian_quadrature