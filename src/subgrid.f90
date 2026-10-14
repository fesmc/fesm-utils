module subgrid

    use, intrinsic :: iso_fortran_env, only : input_unit, output_unit, error_unit
    use precision

    implicit none

    private
    public :: calc_subgrid_array            ! generic (sp/dp)
    public :: calc_subgrid_array_mask       ! generic (sp/dp), with mask
    public :: calc_subgrid_array_cell       ! generic (sp/dp)
    public :: calc_subgrid_array_quad       ! generic (sp/dp)
    public :: calc_subgrid_array_dp         ! double-precision worker
    public :: calc_subgrid_array_mask_dp    ! double-precision worker (mask)
    public :: calc_subgrid_array_cell_dp    ! double-precision worker
    public :: calc_subgrid_array_quad_dp    ! double-precision worker

    ! The subgrid interpolation is always performed in double precision (the _dp
    ! procedures). The _sp procedures are thin wrappers that promote the 3x3
    ! neighbourhood of cell (i,j) (all the workers read), call the _dp worker on
    ! it and demote the result. Promoting the whole field instead would cost
    ! O(nx*ny) per cell.
    ! Exposing each operation as a generic interface lets callers built in
    ! either precision (e.g. Yelmo with wp=sp or wp=dp) use the generic name and
    ! have it resolve to the matching specific by argument kind.

    interface calc_subgrid_array
        module procedure calc_subgrid_array_sp
        module procedure calc_subgrid_array_dp
    end interface calc_subgrid_array

    interface calc_subgrid_array_mask
        module procedure calc_subgrid_array_mask_sp
        module procedure calc_subgrid_array_mask_dp
    end interface calc_subgrid_array_mask

    interface calc_subgrid_array_cell
        module procedure calc_subgrid_array_cell_sp
        module procedure calc_subgrid_array_cell_dp
    end interface calc_subgrid_array_cell

    interface calc_subgrid_array_quad
        module procedure calc_subgrid_array_quad_sp
        module procedure calc_subgrid_array_quad_dp
    end interface calc_subgrid_array_quad

contains

    ! ===================================================================
    ! Single-precision wrappers: promote -> _dp worker -> demote
    ! ===================================================================

    subroutine calc_subgrid_array_sp(vint,v,nxi,i,j,im1,ip1,jm1,jp1)
        ! Single-precision wrapper for calc_subgrid_array_dp.

        implicit none

        real(sp), intent(INOUT) :: vint(:,:)
        real(sp), intent(IN)  :: v(:,:)
        integer,  intent(IN)  :: nxi                    ! Number of interpolation points
        integer,  intent(IN)  :: i, j                   ! Indices of current cell
        integer,  intent(IN)  :: im1, ip1, jm1, jp1     ! Indices of neighbors

        ! Local variables
        real(dp) :: v_nb(3,3)
        real(dp) :: vint_dble(nxi,nxi)

        ! Neighbourhood of (i,j), at local indices 1:3 with the cell at (2,2)
        v_nb = neighborhood_dp(v,i,j,im1,ip1,jm1,jp1)

        call calc_subgrid_array_dp(vint_dble,v_nb,nxi,2,2,1,3,1,3)

        vint = real(vint_dble,sp)

        return

    end subroutine calc_subgrid_array_sp

    subroutine calc_subgrid_array_mask_sp(vint,v,mask,nxi,i,j,im1,ip1,jm1,jp1)
        ! Single-precision wrapper for calc_subgrid_array_mask_dp.

        implicit none

        real(sp), intent(INOUT) :: vint(:,:)
        real(sp), intent(IN)  :: v(:,:)
        logical,  intent(IN)  :: mask(:,:)
        integer,  intent(IN)  :: nxi                    ! Number of interpolation points
        integer,  intent(IN)  :: i, j                   ! Indices of current cell
        integer,  intent(IN)  :: im1, ip1, jm1, jp1     ! Indices of neighbors

        ! Local variables
        real(dp) :: v_nb(3,3)
        logical  :: mask_nb(3,3)
        real(dp) :: vint_dble(nxi,nxi)

        ! Neighbourhood of (i,j), at local indices 1:3 with the cell at (2,2)
        v_nb = neighborhood_dp(v,i,j,im1,ip1,jm1,jp1)
        mask_nb = neighborhood_mask(mask,i,j,im1,ip1,jm1,jp1)

        call calc_subgrid_array_mask_dp(vint_dble,v_nb,mask_nb,nxi,2,2,1,3,1,3)

        vint = real(vint_dble,sp)

        return

    end subroutine calc_subgrid_array_mask_sp

    subroutine calc_subgrid_array_cell_sp(vint,v1,v2,v3,v4,nxi)
        ! Single-precision wrapper for calc_subgrid_array_cell_dp.

        implicit none

        real(sp), intent(INOUT) :: vint(:,:)
        real(sp), intent(IN)  :: v1, v2, v3, v4
        integer,  intent(IN)  :: nxi                    ! Number of interpolation points

        ! Local variables
        real(dp), allocatable :: vint_dble(:,:)

        allocate(vint_dble(nxi,nxi))

        call calc_subgrid_array_cell_dp(vint_dble,real(v1,dp),real(v2,dp), &
                                                real(v3,dp),real(v4,dp),nxi)

        vint = real(vint_dble,sp)

        return

    end subroutine calc_subgrid_array_cell_sp

    subroutine calc_subgrid_array_quad_sp(vint,v,nxi,i,j,im1,ip1,jm1,jp1)
        ! Single-precision wrapper for calc_subgrid_array_quad_dp.

        implicit none

        real(sp), intent(INOUT) :: vint(:,:)
        real(sp), intent(IN)  :: v(:,:)
        integer,  intent(IN)  :: nxi                    ! Number of interpolation points
        integer,  intent(IN)  :: i, j                   ! Indices of current cell
        integer,  intent(IN)  :: im1, ip1, jm1, jp1     ! Indices of neighbors

        ! Local variables
        real(dp) :: v_nb(3,3)
        real(dp) :: vint_dble(nxi,nxi)

        ! Neighbourhood of (i,j), at local indices 1:3 with the cell at (2,2)
        v_nb = neighborhood_dp(v,i,j,im1,ip1,jm1,jp1)

        call calc_subgrid_array_quad_dp(vint_dble,v_nb,nxi,2,2,1,3,1,3)

        vint = real(vint_dble,sp)

        return

    end subroutine calc_subgrid_array_quad_sp

    ! ===================================================================
    ! Double-precision workers: the actual computation
    ! ===================================================================

    subroutine calc_subgrid_array_dp(vint,v,nxi,i,j,im1,ip1,jm1,jp1)

        implicit none

        real(dp), intent(INOUT) :: vint(:,:)
        real(dp), intent(IN)  :: v(:,:)
        integer,  intent(IN)  :: nxi                    ! Number of interpolation points
        integer,  intent(IN)  :: i, j                   ! Indices of current cell
        integer,  intent(IN)  :: im1, ip1, jm1, jp1     ! Indices of neighbors

        ! Local variables
        real(dp) :: v1, v2, v3, v4

        if (nxi .eq. 1) then
            ! Case of no interpolation, just set subgrid array equal to current value

            vint = v(i,j)

        else
            ! Subgrid interpolation necessary

            ! First calculate corner values of current cell (ab-nodes)
            v1 = 0.25_dp*(v(i,j) + v(ip1,j) + v(ip1,jp1) + v(i,jp1))
            v2 = 0.25_dp*(v(i,j) + v(im1,j) + v(im1,jp1) + v(i,jp1))
            v3 = 0.25_dp*(v(i,j) + v(im1,j) + v(im1,jm1) + v(i,jm1))
            v4 = 0.25_dp*(v(i,j) + v(ip1,j) + v(ip1,jm1) + v(i,jm1))

            ! Next calculate the subgrid array of values for this cell
            call calc_subgrid_array_cell_dp(vint,v1,v2,v3,v4,nxi)

        end if

        return

    end subroutine calc_subgrid_array_dp

    subroutine calc_subgrid_array_mask_dp(vint,v,mask,nxi,i,j,im1,ip1,jm1,jp1)

        implicit none

        real(dp), intent(INOUT) :: vint(:,:)
        real(dp), intent(IN)    :: v(:,:)
        logical,  intent(IN)    :: mask(:,:)              ! mask array
        integer,  intent(IN)    :: nxi                    ! Number of interpolation points
        integer,  intent(IN)    :: i, j                   ! Indices of current cell
        integer,  intent(IN)    :: im1, ip1, jm1, jp1     ! Indices of neighbors

        ! Local variables
        real(dp) :: v1, v2, v3, v4
        real(dp) :: sumval
        integer  :: count
        logical  :: use_mask

        use_mask = .TRUE.       !present(mask)

        if (nxi .eq. 1) then
            ! Case of no interpolation, just set subgrid array equal to current value
            vint = v(i,j)

        else
            ! Subgrid interpolation necessary
            ! Each subgrid corner is based on 4 surrounding v-values,
            ! but we only include those allowed by mask.

            ! v1: (i,j), (ip1,j), (ip1,jp1), (i,jp1)
            sumval = 0.0_dp
            count  = 0
            if (mask(i,j)) then
                sumval = sumval + v(i,j); count = count + 1
            end if
            if (mask(ip1,j)) then
                sumval = sumval + v(ip1,j); count = count + 1
            end if
            if (mask(ip1,jp1)) then
                sumval = sumval + v(ip1,jp1); count = count + 1
            end if
            if (mask(i,jp1)) then
                sumval = sumval + v(i,jp1); count = count + 1
            end if
            if (count > 0) then
                v1 = sumval / real(count,dp)
            else
                v1 = v(i,j)  ! fallback
            end if

            ! v2: (i,j), (im1,j), (im1,jp1), (i,jp1)
            sumval = 0.0_dp; count = 0
            if (mask(i,j)) then
                sumval = sumval + v(i,j); count = count + 1
            end if
            if (mask(im1,j)) then
                sumval = sumval + v(im1,j); count = count + 1
            end if
            if (mask(im1,jp1)) then
                sumval = sumval + v(im1,jp1); count = count + 1
            end if
            if (mask(i,jp1)) then
                sumval = sumval + v(i,jp1); count = count + 1
            end if
            if (count > 0) then
                v2 = sumval / real(count,dp)
            else
                v2 = v(i,j)
            end if

            ! v3: (i,j), (im1,j), (im1,jm1), (i,jm1)
            sumval = 0.0_dp; count = 0
            if (mask(i,j)) then
                sumval = sumval + v(i,j); count = count + 1
            end if
            if (mask(im1,j)) then
                sumval = sumval + v(im1,j); count = count + 1
            end if
            if (mask(im1,jm1)) then
                sumval = sumval + v(im1,jm1); count = count + 1
            end if
            if (mask(i,jm1)) then
                sumval = sumval + v(i,jm1); count = count + 1
            end if
            if (count > 0) then
                v3 = sumval / real(count,dp)
            else
                v3 = v(i,j)
            end if

            ! v4: (i,j), (ip1,j), (ip1,jm1), (i,jm1)
            sumval = 0.0_dp; count = 0
            if (mask(i,j)) then
                sumval = sumval + v(i,j); count = count + 1
            end if
            if (mask(ip1,j)) then
                sumval = sumval + v(ip1,j); count = count + 1
            end if
            if (mask(ip1,jm1)) then
                sumval = sumval + v(ip1,jm1); count = count + 1
            end if
            if (mask(i,jm1)) then
                sumval = sumval + v(i,jm1); count = count + 1
            end if
            if (count > 0) then
                v4 = sumval / real(count,dp)
            else
                v4 = v(i,j)
            end if

            ! Next calculate the subgrid array of values for this cell
            call calc_subgrid_array_cell_dp(vint,v1,v2,v3,v4,nxi)

        end if

        return

    end subroutine calc_subgrid_array_mask_dp

    subroutine calc_subgrid_array_cell_dp(vint,v1,v2,v3,v4,nxi)
        ! Given the four corners of a cell in quadrants 1,2,3,4,
        ! calculate the subgrid values via linear interpolation
        ! Assumes vint is a square array of dimensions nxi: vint[nxi,nxi]
        ! Convention:
        !
        !    v2---v1
        !    |     |
        !    |     |
        !    v3---v4

        implicit none

        real(dp), intent(INOUT) :: vint(:,:)
        real(dp), intent(IN)  :: v1, v2, v3, v4
        integer,  intent(IN)  :: nxi                    ! Number of interpolation points

        ! Local variables
        integer :: i, j
        real(dp) :: x(nxi), y(nxi)

        if (nxi .eq. 1) then
            ! Make sure interpolation point represents the center of the subgrid array
            x(1) = 0.5_dp
            y(1) = 0.5_dp
        else
            ! Populate x,y axes for interpolation points (between 0 and 1)
            do i = 1, nxi
                x(i) = 0.0_dp + real(i-1,dp)/real(nxi-1,dp)
            end do
            y = x
        end if

        ! Calculate interpolated value
        vint = 0.0_dp
        do i = 1, nxi
        do j = 1, nxi

            vint(i,j) = interp_bilin_pt(v1,v2,v3,v4,x(i),y(j))

        end do
        end do

        return

    end subroutine calc_subgrid_array_cell_dp

    subroutine calc_subgrid_array_quad_dp(vint,v,nxi,i,j,im1,ip1,jm1,jp1)
        ! Subgrid values of cell (i,j) from the bilinear interpolation of v
        ! between cell centres. Each quadrant of the cell is bilinear between
        ! the cell centre v(i,j), the two face midpoints (mean of v(i,j) and
        ! the neighbour across the face) and the cell corner (mean of the four
        ! cells around it), as the grounded-fraction quadrants of Leguy et al.
        ! (2021). Unlike calc_subgrid_array, which interpolates between the
        ! corner means only, the field passes through v(i,j) at the centre.
        ! The nxi*nxi points are the centres of an nxi x nxi partition of the
        ! cell, so each point stands for the same area (nxi = 1: v(i,j)).
        ! vint(i1,j1): i1 along x, j1 along y.

        implicit none

        real(dp), intent(INOUT) :: vint(:,:)
        real(dp), intent(IN)  :: v(:,:)
        integer,  intent(IN)  :: nxi                    ! Number of interpolation points per side
        integer,  intent(IN)  :: i, j                   ! Indices of current cell
        integer,  intent(IN)  :: im1, ip1, jm1, jp1     ! Indices of neighbors

        ! Local variables
        integer  :: i1, j1, ii, jj
        real(dp) :: x, y, s, t
        real(dp) :: v_x, v_y, v_xy

        do j1 = 1, nxi

            ! Quadrant in y and distance from the centre (0) to the face (1)
            y  = (real(j1,dp)-0.5_dp)/real(nxi,dp)
            jj = merge(jp1,jm1,y .ge. 0.5_dp)
            t  = abs(2.0_dp*y-1.0_dp)

            do i1 = 1, nxi

                x  = (real(i1,dp)-0.5_dp)/real(nxi,dp)
                ii = merge(ip1,im1,x .ge. 0.5_dp)
                s  = abs(2.0_dp*x-1.0_dp)

                ! Face midpoints and corner of the quadrant
                v_x  = 0.5_dp *(v(i,j)+v(ii,j))
                v_y  = 0.5_dp *(v(i,j)+v(i,jj))
                v_xy = 0.25_dp*(v(i,j)+v(ii,j)+v(i,jj)+v(ii,jj))

                vint(i1,j1) = (1.0_dp-s)*(1.0_dp-t)*v(i,j) + s*(1.0_dp-t)*v_x &
                            + (1.0_dp-s)*t*v_y + s*t*v_xy

            end do
        end do

        return

    end subroutine calc_subgrid_array_quad_dp

    pure function neighborhood_dp(v,i,j,im1,ip1,jm1,jp1) result(v_nb)
        ! 3x3 neighbourhood of cell (i,j) of a single-precision field, in
        ! double precision, with the cell at (2,2). Gathered by index, so
        ! wrapped (periodic) neighbour indices are handled.

        implicit none

        real(sp), intent(IN) :: v(:,:)
        integer,  intent(IN) :: i, j, im1, ip1, jm1, jp1
        real(dp) :: v_nb(3,3)

        v_nb(1,1) = real(v(im1,jm1),dp)
        v_nb(2,1) = real(v(i,  jm1),dp)
        v_nb(3,1) = real(v(ip1,jm1),dp)
        v_nb(1,2) = real(v(im1,j),  dp)
        v_nb(2,2) = real(v(i,  j),  dp)
        v_nb(3,2) = real(v(ip1,j),  dp)
        v_nb(1,3) = real(v(im1,jp1),dp)
        v_nb(2,3) = real(v(i,  jp1),dp)
        v_nb(3,3) = real(v(ip1,jp1),dp)

    end function neighborhood_dp

    pure function neighborhood_mask(mask,i,j,im1,ip1,jm1,jp1) result(mask_nb)
        ! 3x3 neighbourhood of cell (i,j) of a mask, with the cell at (2,2)

        implicit none

        logical, intent(IN) :: mask(:,:)
        integer, intent(IN) :: i, j, im1, ip1, jm1, jp1
        logical :: mask_nb(3,3)

        mask_nb(1,1) = mask(im1,jm1)
        mask_nb(2,1) = mask(i,  jm1)
        mask_nb(3,1) = mask(ip1,jm1)
        mask_nb(1,2) = mask(im1,j)
        mask_nb(2,2) = mask(i,  j)
        mask_nb(3,2) = mask(ip1,j)
        mask_nb(1,3) = mask(im1,jp1)
        mask_nb(2,3) = mask(i,  jp1)
        mask_nb(3,3) = mask(ip1,jp1)

    end function neighborhood_mask

    function interp_bilin_pt(z1,z2,z3,z4,xout,yout) result(zout)
        ! Interpolate a point given four neighbors at corners of square (0:1,0:1)
        ! z2    z1
        !    x,y
        ! z3    z4
        !

        implicit none

        real(dp), intent(IN) :: z1, z2, z3, z4
        real(dp), intent(IN) :: xout, yout
        real(dp) :: zout

        ! Local variables
        real(dp) :: x0, x1, y0, y1
        real(dp) :: alpha1, alpha2, p0, p1

        x0 = 0.0_dp
        x1 = 1.0_dp
        y0 = 0.0_dp
        y1 = 1.0_dp

        alpha1  = (xout - x0) / (x1-x0)
        p0      = z3 + alpha1*(z4-z3)
        p1      = z2 + alpha1*(z1-z2)

        alpha2  = (yout - y0) / (y1-y0)
        zout    = p0 + alpha2*(p1-p0)

        return

    end function interp_bilin_pt

end module subgrid
