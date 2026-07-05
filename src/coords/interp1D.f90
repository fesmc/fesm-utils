
module interp1D

    use precision, only: dp, sp, wp
    use constants, only: mv_dp, pi
    implicit none

    ! interp_linear and interp_spline are generic over sp and dp. The dp routines
    ! are authoritative (the numerics live in one place); the sp variants promote
    ! their arguments to dp, call the dp routine, and demote the result back to sp.
    interface interp_linear
        module procedure interp_linear_pt_sp,  interp_linear_pt_dp
        module procedure interp_linear_vec_sp, interp_linear_vec_dp
    end interface

    interface interp_spline
        module procedure interp_spline_sp, interp_spline_dp
    end interface

    private
    public :: interp_align
    public :: interp_linear
    public :: interp1D_bins
    public :: interp_spline

contains

    function interp_align(x,y,xout,missing_value,tol) result(yout)
        ! Fill in a vector of times with the nearest input
        ! times to within a certain tolerance
        implicit none

        real(wp), dimension(:), intent(IN) :: x, y
        real(wp), dimension(:), intent(IN) :: xout
        real(wp), intent(IN), optional :: missing_value, tol
        real(wp), dimension(size(xout)) :: yout
        real(wp) :: tolerance
        integer :: i, n, nout

        ! Length of output vector
        nout = size(xout)

        ! Define tolerance for filling points
!         tolerance = 0.1d0
        tolerance = minval(xout(2:nout)-xout(1:nout-1),dim=1)*0.49d0
        if (present(tol)) tolerance = tol

        ! Define missing values and fill vector intially
        yout = mv_dp
        if (present(missing_value)) yout = missing_value

        do i = 1, nout
            n = minloc(abs(x-xout(i)),dim=1)
            if (abs(x(n)-xout(i)) <= tolerance) yout(i) = y(n)
        end do

        return

    end function interp_align

    function interp_linear_pt_sp(x,y,xout) result(yout)
        ! sp wrapper around interp_linear_pt_dp
        implicit none
        real(sp), dimension(:), intent(IN) :: x, y
        real(sp),               intent(IN) :: xout
        real(sp) :: yout
        yout = real(interp_linear_pt_dp(real(x,dp),real(y,dp),real(xout,dp)), sp)
    end function interp_linear_pt_sp

    function interp_linear_pt_dp(x,y,xout) result(yout)
        ! Interpolate y from ordered x to ordered xout positions

        implicit none

        real(dp), dimension(:), intent(IN) :: x, y
        real(dp), intent(IN) :: xout
        real(dp) :: yout
        integer :: i, j, n, nout

        n    = size(x)

        if (xout .lt. x(1)) then
            yout = y(1)
        else if (xout .gt. x(n)) then
            yout = y(n)
        else
            do j = 1, n
                if (x(j) .ge. xout) exit
            end do

            if (j .eq. 1) then
                yout = y(1)
            else if (j .eq. n+1) then
                yout = y(n)
            else
                yout = interp_linear_internal(x(j-1:j),y(j-1:j),xout)
            end if
        end if

        return

      end function interp_linear_pt_dp

    function interp_linear_vec_sp(x,y,xout) result(yout)
        ! sp wrapper around interp_linear_vec_dp
        implicit none
        real(sp), dimension(:), intent(IN) :: x, y
        real(sp), dimension(:), intent(IN) :: xout
        real(sp), dimension(size(xout)) :: yout
        yout = real(interp_linear_vec_dp(real(x,dp),real(y,dp),real(xout,dp)), sp)
    end function interp_linear_vec_sp

    function interp_linear_vec_dp(x,y,xout) result(yout)
        ! Interpolate y from ordered x to ordered xout positions

        implicit none

        real(dp), dimension(:), intent(IN) :: x, y
        real(dp), dimension(:), intent(IN) :: xout
        real(dp), dimension(size(xout)) :: yout
        integer :: i, j, n, nout

        n    = size(x)
        nout = size(xout)

        do i = 1, nout
            if (xout(i) .lt. x(1)) then
                yout(i) = y(1)
            else if (xout(i) .gt. x(n)) then
                yout(i) = y(n)
            else
                do j = 1, n
                    if (x(j) .ge. xout(i)) exit
                end do

                if (j .eq. 1) then
                    yout(i) = y(1)
                else if (j .eq. n+1) then
                    yout(i) = y(n)
                else
                    yout(i) = interp_linear_internal(x(j-1:j),y(j-1:j),xout(i))
                end if
            end if
        end do

        return

      end function interp_linear_vec_dp

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !   Subroutine :  interp_linear_internal
    !   Author     :  Alex Robinson
    !   Purpose    :  Interpolates for the y value at the desired x value,
    !                 given x and y values around the desired point.
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    function interp_linear_internal(x,y,xout) result(yout)

        implicit none

        real(dp), intent(IN)  :: x(2), y(2), xout
        real(dp) :: yout
        real(dp) :: alph

        if ( xout .lt. x(1) .or. xout .gt. x(2) ) then
            write(*,*) "interp1: xout < x0 or xout > x1 !"
            write(*,*) "xout = ",xout
            write(*,*) "x0   = ",x(1)
            write(*,*) "x1   = ",x(2)
            stop
        end if

        alph = (xout - x(1)) / (x(2) - x(1))
        yout = y(1) + alph*(y(2) - y(1))

        return

    end function interp_linear_internal


    subroutine interp1D_bins(y,x,y0,x0)

        implicit none

        real(wp), intent(OUT) :: y(:)
        real(wp), intent(IN)  :: x(:)
        real(wp), intent(IN)  :: y0(:)
        real(wp), intent(IN)  :: x0(:)

        ! Local variables
        integer :: k, k0, k1
        integer :: nx, nx0
        real(wp) :: dx, dx0, xl0, xr0, xl, xr
        real(wp) :: dx_now, wt_now

        real(wp), allocatable :: wts(:)

        nx  = size(x,1)
        nx0 = size(x0,1)

        ! Consistency check
        if (x(1) .ne. x0(1) .or. x(nx) .ne. x0(nx0)) then
            write(*,*) "interp1D_bins:: Error: for this routine, boundary values of axes must be the same."
            stop
        end if

        allocate(wts(nx))

        ! Initially set target vector and weights to zero
        y   = 0.0_wp
        wts = 0.0_wp

        ! Store boundary values, which
        ! do not require interpolation
        y(1)  = y0(1)
        y(nx) = y0(nx0)

        wts(1)  = 1.0_wp
        wts(nx) = 1.0_wp

        ! Loop over internal values of
        ! original vector and
        ! distribute contributions to each
        ! bin in target vector
        do k0 = 2, nx0-1

            ! Define left and right boundaries
            ! of current original cells
            if (k0 .eq. 2) then
                xl0 = x0(1)
            else
                xl0 = 0.5_wp*(x0(k0-1)+x0(k0))
            end if

            if (k0 .eq. nx0-1) then
                xr0 = x0(nx0)
            else
                xr0 = 0.5_wp*(x0(k0)+x0(k0+1))
            end if

            ! Add contribution if any from
            ! current original cell to any
            ! target cells.

            do k = 2, nx-1

                ! Define left and right boundaries
                ! of target cell
                if (k .eq. 2) then
                    xl = x(1)
                else
                    xl = 0.5_wp*(x(k-1)+x(k))
                end if

                if (k .eq. nx-1) then
                    xr = x(nx)
                else
                    xr = 0.5_wp*(x(k)+x(k+1))
                end if

                ! Total target cell width
                dx = xr-xl

                if (xl0 .ge. xr .or. xr0 .le. xl) then
                    ! No part of original cell appears in target cell, set
                    ! weight to zero
                    dx_now = 0.0_wp
                else
                    ! Determine width of original cell that appears in target cell
                    dx_now = min(xr,xr0) - max(xl,xl0)
                end if

                wt_now = dx_now / dx

                y(k)   = y(k) + wt_now*y0(k0)

                wts(k) = wts(k) + wt_now

            end do

        end do

        return

    end subroutine interp1D_bins

    function interp_spline_sp(x,y,xout) result(yout)
        ! sp wrapper around interp_spline_dp
        implicit none
        real(sp), dimension(:), intent(IN) :: x, y
        real(sp), dimension(:), intent(IN) :: xout
        real(sp), dimension(size(xout)) :: yout
        yout = real(interp_spline_dp(real(x,dp),real(y,dp),real(xout,dp)), sp)
    end function interp_spline_sp

    function interp_spline_dp(x,y,xout) result(yout)

        implicit none

        real(dp), dimension(:), intent(IN) :: x, y
        real(dp), dimension(:), intent(IN) :: xout
        real(dp), dimension(size(xout)) :: yout
        real(dp), dimension(:), allocatable :: b, c, d
        real(dp) :: uh, dx, yh
        integer :: i, n, nout

        n    = size(x)
        nout = size(xout)

        ! Get spline coefficients b, c, d
        allocate(b(n),c(n),d(n))
        call spline (x, y, b, c, d, n)

        do i = 1, nout
            if (xout(i) .lt. x(1)) then
                dx = x(1)-xout(i)
                uh = x(1)+dx
                yh = ispline(uh,x,y,b,c,d,n)
                yout(i) = y(1) + (y(1)-yh)
            else if (xout(i) .gt. x(n)) then
                dx = xout(i)-x(n)
                uh = x(n)-dx
                yh = ispline(uh,x,y,b,c,d,n)
                yout(i) = y(n) + (y(n)-yh)
            else
                yout(i) = ispline(xout(i), x, y, b, c, d, n)
            end if
        end do

        return

    end function interp_spline_dp

subroutine spline (x, y, b, c, d, n)
!======================================================================
!  SOURCE: http://ww2.odu.edu/~agodunov/computing/programs/book2/Ch01/spline.f90
!  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n
!  for cubic spline interpolation
!  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
!  for  x(i) <= x <= x(i+1)
!  Alex G: January 2010
!----------------------------------------------------------------------
!  input..
!  x = the arrays of data abscissas (in strictly increasing order)
!  y = the arrays of data ordinates
!  n = size of the arrays xi() and yi() (n>=2)
!  output..
!  b, c, d  = arrays of spline coefficients
!  comments ...
!  spline.f90 program is based on fortran version of program spline.f
!  the accompanying function fspline can be used for interpolation
!======================================================================
implicit none
integer n
real(dp), dimension(:) :: x, y, b, c, d
integer i, j, gap
real(dp) h

gap = n-1
! check input
if ( n < 2 ) return
if ( n < 3 ) then
  b(1) = (y(2)-y(1))/(x(2)-x(1))   ! linear interpolation
  c(1) = 0.
  d(1) = 0.
  b(2) = b(1)
  c(2) = 0.
  d(2) = 0.
  return
end if
!
! step 1: preparation
!
d(1) = x(2) - x(1)
c(2) = (y(2) - y(1))/d(1)
do i = 2, gap
  d(i) = x(i+1) - x(i)
  b(i) = 2.0*(d(i-1) + d(i))
  c(i+1) = (y(i+1) - y(i))/d(i)
  c(i) = c(i+1) - c(i)
end do
!
! step 2: end conditions
!
b(1) = -d(1)
b(n) = -d(n-1)
c(1) = 0.0
c(n) = 0.0
if(n /= 3) then
  c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
  c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
  c(1) = c(1)*d(1)**2/(x(4)-x(1))
  c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
end if
!
! step 3: forward elimination
!
do i = 2, n
  h = d(i-1)/b(i-1)
  b(i) = b(i) - h*d(i-1)
  c(i) = c(i) - h*c(i-1)
end do
!
! step 4: back substitution
!
c(n) = c(n)/b(n)
do j = 1, gap
  i = n-j
  c(i) = (c(i) - d(i)*c(i+1))/b(i)
end do
!
! step 5: compute spline coefficients
!
b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2.0*c(n))
do i = 1, gap
  b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0*c(i))
  d(i) = (c(i+1) - c(i))/d(i)
  c(i) = 3.*c(i)
end do
c(n) = 3.0*c(n)
d(n) = d(n-1)
end subroutine spline

function ispline(u, x, y, b, c, d, n)
!======================================================================
! SOURCE: http://ww2.odu.edu/~agodunov/computing/programs/book2/Ch01/spline.f90
! function ispline evaluates the cubic spline interpolation at point z
! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
! where  x(i) <= u <= x(i+1)
!----------------------------------------------------------------------
! input..
! u       = the abscissa at which the spline is to be evaluated
! x, y    = the arrays of given data points
! b, c, d = arrays of spline coefficients computed by spline
! n       = the number of data points
! output:
! ispline = interpolated value at point u
!=======================================================================
implicit none
real(dp) ispline
integer n
real(dp) :: u
real(dp), dimension(:) :: x, y, b, c, d
integer i, j, k
real(dp) dx

! if u is ouside the x() interval take a boundary value (left or right)
if(u <= x(1)) then
  ispline = y(1)
  return
end if
if(u >= x(n)) then
  ispline = y(n)
  return
end if

!*
!  binary search for for i, such that x(i) <= u <= x(i+1)
!*
i = 1
j = n+1
do while (j > i+1)
  k = (i+j)/2
  if(u < x(k)) then
    j=k
    else
    i=k
   end if
end do
!*
!  evaluate spline interpolation
!*
dx = u - x(i)
ispline = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
end function ispline


end module interp1D
