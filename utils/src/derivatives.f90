module derivatives

    use, intrinsic :: iso_fortran_env, only : input_unit, output_unit, error_unit    
    use precision

    implicit none
    
    private
    public :: calc_dvdx_2D
    public :: calc_dvdy_2D
    public :: calc_dvdx_1D
    public :: calc_dvdy_1D

contains

    subroutine calc_dvdx_2D(dvdx,v,dx,mask,bc)
        
        implicit none

        real(wp), intent(OUT) :: dvdx(:,:)
        real(wp), intent(IN)  :: v(:,:)
        real(wp), intent(IN)  :: dx
        logical,  intent(IN)  :: mask(:,:)
        character(len=*), intent(IN) :: bc

        ! Local variables
        integer :: nx, ny, i, j, im1, ip1, jm1, jp1

        real(dp) :: wu(3), wc(3), wd(3)
        real(dp) :: f(3)

        nx = size(dvdx,1)
        ny = size(dvdx,2)

        ! Define stencil weights
        wu = [-3.0d0,  4.0d0, -1.0d0] / (2.0d0*dx)      ! upstream
        wc = [-1.0d0,  0.0d0,  1.0d0] / (2.0d0*dx)      ! centered
        wd = [ 3.0d0, -4.0d0,  1.0d0] / (2.0d0*dx)      ! downstream

        do j = 1, ny
        do i = 1, nx
            if (mask(i,j)) then
                ! Calculate derivative at this point

                if (i .eq. 1) then
                    ! Left-boundary
                    if (trim(bc) .eq. "periodic") then
                        f = [v(nx,j),v(i,j),v(i+1,j)]
                        dvdx(i,j) = sum(wc*f)
                    else
                        f = [v(i,j),v(i+1,j),v(i+2,j)]
                        dvdx(i,j) = sum(wd*f)
                    end if
                else if (i .eq. nx) then
                    ! Right-boundary
                    if (trim(bc) .eq. "periodic") then
                        f = [v(i-1,j),v(i,j),v(1,j)]
                        dvdx(i,j) = sum(wc*f)
                    else
                        f = [v(i-2,j),v(i-1,j),v(i,j)]
                        dvdx(i,j) = sum(wu*f)
                    end if
                else
                    ! Inner point
                    f = [v(i-1,j),v(i,j),v(i+1,j)]
                    dvdx(i,j) = sum(wc*f)
                end if

            end if
        end do
        end do

        return
        
    end subroutine calc_dvdx_2D

    subroutine calc_dvdy_2D(dvdy,v,dy,mask,bc)
        
        implicit none

        real(wp), intent(OUT) :: dvdy(:,:)
        real(wp), intent(IN)  :: v(:,:)
        real(wp), intent(IN)  :: dy
        logical,  intent(IN)  :: mask(:,:)
        character(len=*), intent(IN) :: bc

        ! Local variables
        integer :: nx, ny, i, j, im1, ip1, jm1, jp1

        real(dp) :: wu(3), wc(3), wd(3)
        real(dp) :: f(3)

        nx = size(dvdy,1)
        ny = size(dvdy,2)

        ! Define stencil weights
        wu = [-3.0d0,  4.0d0, -1.0d0] / (2.0d0*dy)      ! upstream
        wc = [-1.0d0,  0.0d0,  1.0d0] / (2.0d0*dy)      ! centered
        wd = [ 3.0d0, -4.0d0,  1.0d0] / (2.0d0*dy)      ! downstream

        do j = 1, ny
        do i = 1, nx
            if (mask(i,j)) then
                ! Calculate derivative at this point

                if (j .eq. 1) then
                    ! Lower-boundary
                    if (trim(bc) .eq. "periodic") then
                        f = [v(i,ny),v(i,j),v(i,j+1)]
                        dvdy(i,j) = sum(wc*f)
                    else
                        f = [v(i,j),v(i,j+1),v(i,j+2)]
                        dvdy(i,j) = sum(wd*f)
                    end if
                else if (i .eq. ny) then
                    ! Upper-boundary
                    if (trim(bc) .eq. "periodic") then
                        f = [v(i,j-1),v(i,j),v(i,1)]
                        dvdy(i,j) = sum(wc*f)
                    else
                        f = [v(i,j-2),v(i,j-1),v(i,j)]
                        dvdy(i,j) = sum(wu*f)
                    end if
                else
                    ! Inner point
                    f = [v(i,j-1),v(i,j),v(i,j+1)]
                    dvdy(i,j) = sum(wc*f)
                end if

            end if
        end do
        end do

        return
        
    end subroutine calc_dvdy_2D

    subroutine calc_dvdx_1D(dvdx,v,dx,mask,bc)
        
        implicit none

        real(wp), intent(OUT) :: dvdx(:)
        real(wp), intent(IN)  :: v(:)
        real(wp), intent(IN)  :: dx
        logical,  intent(IN)  :: mask(:)
        character(len=*), intent(IN) :: bc

        ! Local variables
        integer :: nx, i, im1, ip1

        real(dp) :: wu(3), wc(3), wd(3)
        real(dp) :: f(3)

        nx = size(dvdx,1)

        ! Define stencil weights
        wu = [-3.0d0,  4.0d0, -1.0d0] / (2.0d0*dx)      ! upstream
        wc = [-1.0d0,  0.0d0,  1.0d0] / (2.0d0*dx)      ! centered
        wd = [ 3.0d0, -4.0d0,  1.0d0] / (2.0d0*dx)      ! downstream

        do i = 1, nx
            if (mask(i)) then
                ! Calculate derivative at this point

                if (i .eq. 1) then
                    ! Left-boundary
                    if (trim(bc) .eq. "periodic") then
                        f = [v(nx),v(i),v(i+1)]
                        dvdx(i) = sum(wc*f)
                    else
                        f = [v(i),v(i+1),v(i+2)]
                        dvdx(i) = sum(wd*f)
                    end if
                else if (i .eq. nx) then
                    ! Right-boundary
                    if (trim(bc) .eq. "periodic") then
                        f = [v(i-1),v(i),v(1)]
                        dvdx(i) = sum(wc*f)
                    else
                        f = [v(i-2),v(i-1),v(i)]
                        dvdx(i) = sum(wu*f)
                    end if
                else
                    ! Inner point
                    f = [v(i-1),v(i),v(i+1)]
                    dvdx(i) = sum(wc*f)
                end if

            end if
        end do

        return
        
    end subroutine calc_dvdx_1D

    subroutine calc_dvdy_1D(dvdy,v,dy,mask,bc)
        
        implicit none

        real(wp), intent(OUT) :: dvdy(:)
        real(wp), intent(IN)  :: v(:)
        real(wp), intent(IN)  :: dy
        logical,  intent(IN)  :: mask(:)
        character(len=*), intent(IN) :: bc

        ! Just call the dvdx routine - in 1D dvdy=dvdx
        call calc_dvdx_1D(dvdy,v,dy,mask,bc)
        
        return
        
    end subroutine calc_dvdy_1D

end module derivatives