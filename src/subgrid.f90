module subgrid

    use, intrinsic :: iso_fortran_env, only : input_unit, output_unit, error_unit
    use precision

    implicit none

    private
    public :: calc_subgrid_array            ! generic (sp/dp)
    public :: calc_subgrid_array_mask       ! generic (sp/dp), with mask
    public :: calc_subgrid_array_cell       ! generic (sp/dp)
    public :: calc_subgrid_array_dp         ! double-precision worker
    public :: calc_subgrid_array_mask_dp    ! double-precision worker (mask)
    public :: calc_subgrid_array_cell_dp    ! double-precision worker

    ! The subgrid interpolation is always performed in double precision (the _dp
    ! procedures). The _sp procedures are thin wrappers that promote their
    ! single-precision arguments, call the _dp worker and demote the result.
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

contains

    ! ===================================================================
    ! Single-precision wrappers: promote -> _dp worker -> demote
    ! ===================================================================

    subroutine calc_subgrid_array_sp(vint,v,nxi,i,j,im1,ip1,jm1,jp1)

        implicit none

        real(sp), intent(INOUT) :: vint(:,:)
        real(sp), intent(IN)  :: v(:,:)
        integer,  intent(IN)  :: nxi                    ! Number of interpolation points
        integer,  intent(IN)  :: i, j                   ! Indices of current cell
        integer,  intent(IN)  :: im1, ip1, jm1, jp1     ! Indices of neighbors

        ! Local variables
        real(dp), allocatable :: vint_dble(:,:)

        allocate(vint_dble(nxi,nxi))

        call calc_subgrid_array_dp(vint_dble,real(v,dp),nxi,i,j,im1,ip1,jm1,jp1)

        vint = real(vint_dble,sp)

        return

    end subroutine calc_subgrid_array_sp

    subroutine calc_subgrid_array_mask_sp(vint,v,mask,nxi,i,j,im1,ip1,jm1,jp1)

        implicit none

        real(sp), intent(INOUT) :: vint(:,:)
        real(sp), intent(IN)  :: v(:,:)
        logical,  intent(IN)  :: mask(:,:)
        integer,  intent(IN)  :: nxi                    ! Number of interpolation points
        integer,  intent(IN)  :: i, j                   ! Indices of current cell
        integer,  intent(IN)  :: im1, ip1, jm1, jp1     ! Indices of neighbors

        ! Local variables
        real(dp), allocatable :: vint_dble(:,:)

        allocate(vint_dble(nxi,nxi))

        call calc_subgrid_array_mask_dp(vint_dble,real(v,dp),mask,nxi,i,j,im1,ip1,jm1,jp1)

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
