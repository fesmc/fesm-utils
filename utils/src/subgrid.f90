module subgrid

    use, intrinsic :: iso_fortran_env, only : input_unit, output_unit, error_unit    
    use precision

    implicit none
    
    private
    public :: calc_subgrid_array
    public :: calc_subgrid_array_cell

contains

    subroutine calc_subgrid_array(vint,v,nxi,i,j,im1,ip1,jm1,jp1)

        implicit none

        real(wp), intent(INOUT) :: vint(:,:)  
        real(wp), intent(IN)  :: v(:,:)
        integer,  intent(IN)  :: nxi                    ! Number of interpolation points 
        integer,  intent(IN)  :: i, j                   ! Indices of current cell
        integer,  intent(IN)  :: im1, ip1, jm1, jp1     ! Indices of neighbors

        ! Local variables
        real(wp) :: v1, v2, v3, v4

        if (nxi .eq. 1) then
            ! Case of no interpolation, just set subgrid array equal to current value

            vint = v(i,j)
        
        else
            ! Subgrid interpolation necessary

            ! First calculate corner values of current cell (ab-nodes)
            v1 = 0.25_wp*(v(i,j) + v(ip1,j) + v(ip1,jp1) + v(i,jp1))
            v2 = 0.25_wp*(v(i,j) + v(im1,j) + v(im1,jp1) + v(i,jp1))
            v3 = 0.25_wp*(v(i,j) + v(im1,j) + v(im1,jm1) + v(i,jm1))
            v4 = 0.25_wp*(v(i,j) + v(ip1,j) + v(ip1,jm1) + v(i,jm1))
                
            ! Next calculate the subgrid array of values for this cell
            call calc_subgrid_array_cell(vint,v1,v2,v3,v4,nxi)

        end if

        return
        
    end subroutine calc_subgrid_array

    subroutine calc_subgrid_array_cell(vint,v1,v2,v3,v4,nxi)
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

        real(wp), intent(INOUT) :: vint(:,:)  
        real(wp), intent(IN)  :: v1, v2, v3, v4
        integer,  intent(IN)  :: nxi                    ! Number of interpolation points 

        ! Local variables 
        integer :: i, j 
        real(wp) :: x(nxi), y(nxi) 

        if (nxi .eq. 1) then
            ! Make sure interpolation point represents the center of the subgrid array
            x(1) = 0.5
            y(1) = 0.5
        else
            ! Populate x,y axes for interpolation points (between 0 and 1)
            do i = 1, nxi 
                x(i) = 0.0 + real(i-1)/real(nxi-1)
            end do 
            y = x 
        end if

        ! Calculate interpolated value      
        vint = 0.0 
        do i = 1, nxi 
        do j = 1, nxi 

            vint(i,j) = interp_bilin_pt(v1,v2,v3,v4,x(i),y(j))

        end do 
        end do 

        return 

    end subroutine calc_subgrid_array_cell

    function interp_bilin_pt(z1,z2,z3,z4,xout,yout) result(zout)
        ! Interpolate a point given four neighbors at corners of square (0:1,0:1)
        ! z2    z1
        !    x,y
        ! z3    z4 
        ! 

        implicit none 

        real(wp), intent(IN) :: z1, z2, z3, z4 
        real(wp), intent(IN) :: xout, yout 
        real(wp) :: zout 

        ! Local variables 
        real(wp) :: x0, x1, y0, y1 
        real(wp) :: alpha1, alpha2, p0, p1 

        x0 = 0.0 
        x1 = 1.0 
        y0 = 0.0 
        y1 = 1.0 

        alpha1  = (xout - x0) / (x1-x0)
        p0      = z3 + alpha1*(z4-z3)
        p1      = z2 + alpha1*(z1-z2)
            
        alpha2  = (yout - y0) / (y1-y0)
        zout    = p0 + alpha2*(p1-p0)

        return 

    end function interp_bilin_pt

end module subgrid