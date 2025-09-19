module distances

    use, intrinsic :: iso_fortran_env, only : input_unit, output_unit, error_unit    
    use precision

    implicit none

    private
    public :: compute_distance_to_mask

contains

subroutine compute_distance_to_mask(dist, mask)
    
    implicit none
    
    real(wp), intent(OUT) :: dist(:,:)
    integer,  intent(IN)  :: mask(:,:)  ! -1: inside, 0: at location, 1: outside
    
    ! Loca variables
    integer  :: i, j, nx, ny
    
    ! Chamfer weights
    real(dp), parameter :: w1 = 1.0d0           ! horizontal/vertical
    real(dp), parameter :: w2 = sqrt(2.0d0)     ! diagonal
    real(dp), parameter :: w3 = sqrt(5.0d0)     ! knight move (optional, improves accuracy)

    real(dp), parameter :: dist_max = 1.d9      ! Maximum allowed distance

    nx = size(dist,1)
    ny = size(dist,2)

    ! Initialize: distance is 0 at mask values, big number elsewhere
    do j = 1, ny
    do i = 1, nx
        if (mask(i,j) == 0) then
            dist(i,j) = 0.0d0
        else
            dist(i,j) = dist_max
        end if
    end do
    end do

    ! Forward pass
    do j = 2, ny
    do i = 2, nx
        dist(i,j) = min(dist(i,j), dist(i-1,j) + w1)
        dist(i,j) = min(dist(i,j), dist(i,j-1) + w1)
        dist(i,j) = min(dist(i,j), dist(i-1,j-1) + w2)
        if (i < nx .and. j > 1) dist(i,j) = min(dist(i,j), dist(i+1,j-1) + w2)
    end do
    end do

    ! Backward pass
    do j = ny-1, 1, -1
    do i = nx-1, 1, -1
        dist(i,j) = min(dist(i,j), dist(i+1,j) + w1)
        dist(i,j) = min(dist(i,j), dist(i,j+1) + w1)
        dist(i,j) = min(dist(i,j), dist(i+1,j+1) + w2)
        if (i > 1 .and. j < ny) dist(i,j) = min(dist(i,j), dist(i-1,j+1) + w2)
    end do
    end do

    ! Optional: apply sign based on original mask
    do j = 1, ny
    do i = 1, nx
        if (mask(i,j) < 0) dist(i,j) = -dist(i,j)
    end do
    end do

end subroutine compute_distance_to_mask

end module distances