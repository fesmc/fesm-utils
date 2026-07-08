module distances

    use, intrinsic :: iso_fortran_env, only : input_unit, output_unit, error_unit    
    use precision

    implicit none

    private
    public :: compute_distance_to_mask

contains

subroutine compute_distance_to_mask(dist, mask, periodic_x, periodic_y)

    implicit none

    real(wp), intent(OUT) :: dist(:,:)
    integer,  intent(IN)  :: mask(:,:)  ! -1: inside, 0: at location, 1: outside
    logical,  intent(IN), optional :: periodic_x   ! Wrap x-axis (i) neighbors
    logical,  intent(IN), optional :: periodic_y   ! Wrap y-axis (j) neighbors

    ! Loca variables
    integer  :: i, j, nx, ny, iter
    integer  :: im1, ip1, jm1, jp1
    logical  :: perx, pery
    real(wp), allocatable :: dist_ref(:,:)

    ! Chamfer weights
    real(dp), parameter :: w1 = 1.0d0           ! horizontal/vertical
    real(dp), parameter :: w2 = sqrt(2.0d0)     ! diagonal
    real(dp), parameter :: w3 = sqrt(5.0d0)     ! knight move (optional, improves accuracy)

    real(dp), parameter :: dist_max = 1.d9      ! Maximum allowed distance
    integer,  parameter :: iter_max = 100       ! Safety cap on periodic iterations

    ! Resolve optional periodic flags (default: non-periodic)
    perx = .FALSE.
    pery = .FALSE.
    if (present(periodic_x)) perx = periodic_x
    if (present(periodic_y)) pery = periodic_y

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

    if (.not. perx .and. .not. pery) then
        ! Non-periodic path: single forward+backward chamfer sweep.
        ! This block is kept bit-identical to the original implementation.

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

    else
        ! Periodic path (x and/or y): a single two-pass sweep cannot propagate
        ! distance across a periodic seam, so iterate {forward; backward} until
        ! the dist array is unchanged, with a safety cap of iter_max.
        ! Neighbor indices are wrapped on periodic axes; on a non-periodic axis
        ! an out-of-range neighbor (marked 0) is skipped, matching the guarded
        ! boundary behavior of the non-periodic sweep.

        allocate(dist_ref(nx,ny))

        do iter = 1, iter_max

            dist_ref = dist

            ! Forward pass over the full grid, with wrapping
            do j = 1, ny
            do i = 1, nx

                call get_neighbors(im1,ip1,jm1,jp1,i,j,nx,ny,perx,pery)

                if (im1 /= 0)                dist(i,j) = min(dist(i,j), dist(im1,j)   + w1)
                if (jm1 /= 0)                dist(i,j) = min(dist(i,j), dist(i,jm1)   + w1)
                if (im1 /= 0 .and. jm1 /= 0) dist(i,j) = min(dist(i,j), dist(im1,jm1) + w2)
                if (ip1 /= 0 .and. jm1 /= 0) dist(i,j) = min(dist(i,j), dist(ip1,jm1) + w2)

            end do
            end do

            ! Backward pass over the full grid, with wrapping
            do j = ny, 1, -1
            do i = nx, 1, -1

                call get_neighbors(im1,ip1,jm1,jp1,i,j,nx,ny,perx,pery)

                if (ip1 /= 0)                dist(i,j) = min(dist(i,j), dist(ip1,j)   + w1)
                if (jp1 /= 0)                dist(i,j) = min(dist(i,j), dist(i,jp1)   + w1)
                if (ip1 /= 0 .and. jp1 /= 0) dist(i,j) = min(dist(i,j), dist(ip1,jp1) + w2)
                if (im1 /= 0 .and. jp1 /= 0) dist(i,j) = min(dist(i,j), dist(im1,jp1) + w2)

            end do
            end do

            ! Converged when nothing changed this iteration
            if (all(dist == dist_ref)) exit

        end do

        deallocate(dist_ref)

    end if

    ! Apply sign based on original mask (once, after convergence)
    do j = 1, ny
    do i = 1, nx
        if (mask(i,j) < 0) dist(i,j) = -dist(i,j)
    end do
    end do

end subroutine compute_distance_to_mask

subroutine get_neighbors(im1,ip1,jm1,jp1,i,j,nx,ny,perx,pery)
    ! Return neighbor indices with wrapping on periodic axes.
    ! A neighbor that falls outside a non-periodic axis is returned as 0
    ! (an invalid marker), so callers can skip it.

    implicit none

    integer, intent(OUT) :: im1, ip1, jm1, jp1
    integer, intent(IN)  :: i, j, nx, ny
    logical, intent(IN)  :: perx, pery

    if (i == 1) then
        if (perx) then ; im1 = nx ; else ; im1 = 0 ; end if
    else
        im1 = i-1
    end if

    if (i == nx) then
        if (perx) then ; ip1 = 1  ; else ; ip1 = 0 ; end if
    else
        ip1 = i+1
    end if

    if (j == 1) then
        if (pery) then ; jm1 = ny ; else ; jm1 = 0 ; end if
    else
        jm1 = j-1
    end if

    if (j == ny) then
        if (pery) then ; jp1 = 1  ; else ; jp1 = 0 ; end if
    else
        jp1 = j+1
    end if

    return

end subroutine get_neighbors

end module distances