module topoprep
    ! Cleaning of a reference topography on its grid, for the hub (htopo) to
    ! hand to ice models as given: removal of englacial lakes and a limit on
    ! the bed gradient. Ported from Yelmo (remove_englacial_lakes of
    ! physics/topography.f90, adjust_topography_gradients with smooth_gauss_2D
    ! and gauss_values of yelmo_tools.f90) with the same operations in the same
    ! order, so that both give the same numbers in single precision (wp = sp).

    use precision
    use constants,  only : pi
    use fesm_log,   only : fesm_error, to_str

    implicit none
    private

    public :: topoprep_remove_lakes
    public :: topoprep_limit_gradients

contains

    subroutine topoprep_remove_lakes(H_ice, z_bed, z_srf, z_sl, rho_ice, rho_sw)
        ! Diagnose where ice should be grounded, but a gap exists. In these
        ! locations, increase the ice thickness to fill in the englacial lake.
        ! Also where the ice base is below the bedrock, reduce the thickness.
        ! The surface z_srf is kept, as the reference.
        real(wp), intent(inout) :: H_ice(:,:)   ! [m] ice thickness
        real(wp), intent(in)    :: z_bed(:,:)   ! [m] bedrock elevation
        real(wp), intent(in)    :: z_srf(:,:)   ! [m] surface elevation
        real(wp), intent(in)    :: z_sl         ! [m] sea level
        real(wp), intent(in)    :: rho_ice      ! [kg m-3]
        real(wp), intent(in)    :: rho_sw       ! [kg m-3]

        integer  :: i, j, nx, ny
        real(wp) :: H_grnd_now, rho_sw_ice

        rho_sw_ice = rho_sw/rho_ice

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        do j = 1, ny
        do i = 1, nx

            ! Ice thickness overburden (grounded if > 0)
            H_grnd_now = H_ice(i,j) - rho_sw_ice*max(z_sl-z_bed(i,j),0.0_wp)

            if (H_grnd_now .gt. 0.0_wp .and. z_srf(i,j) - H_ice(i,j) .gt. z_bed(i,j)) then
                ! Grounded but a gap exists: a lake, increase the thickness
                H_ice(i,j) = z_srf(i,j) - z_bed(i,j)

            else if (z_srf(i,j) - H_ice(i,j) .lt. z_bed(i,j)) then
                ! Erroneously thick, reduce it
                H_ice(i,j) = z_srf(i,j) - z_bed(i,j)

            end if

        end do
        end do

    end subroutine topoprep_remove_lakes

    subroutine topoprep_limit_gradients(z_bed, H_ice, grad_lim, dx, n_iter)
        ! Smooth the bedrock and the ice thickness so that the limit on the bed
        ! gradient is not exceeded, only where the gradient is too large and its
        ! neighbour on the staggered node (up to 50 passes of a Gaussian filter
        ! with sigma = 2 dx).
        !
        ! Yelmo passes its boundary treatment, which never matters here: the
        ! gradients are only taken at i = 3..nx-3, j = 3..ny-3, where all
        ! treatments use the neighbour i+1, j+1 (and the clamping at the
        ! borders it also sets lies outside that range).
        real(wp), intent(inout) :: z_bed(:,:)   ! [m] bedrock elevation
        real(wp), intent(inout) :: H_ice(:,:)   ! [m] ice thickness
        real(wp), intent(in)    :: grad_lim     ! [m m-1] limit on |dz_bed/dx|
        real(wp), intent(in)    :: dx           ! [m] grid spacing (dy = dx)
        integer,  intent(out), optional :: n_iter   ! passes of smoothing

        integer  :: i, j, q, nx, ny
        real(wp) :: dzbdx, dzbdy
        logical, allocatable :: mask_apply(:,:)

        integer,  parameter :: iter_max = 50
        real(wp), parameter :: grad_max = 100.0_wp   ! gradients are clipped as in Yelmo's calc_gradient_ac*

        nx = size(z_bed,1)
        ny = size(z_bed,2)

        allocate(mask_apply(nx,ny))

        if (present(n_iter)) n_iter = 0

        do q = 1, iter_max

            ! Where the gradients (on the staggered nodes) are too large
            mask_apply = .false.
            do j = 3, ny-3
            do i = 3, nx-3

                dzbdx = clip((z_bed(i+1,j)-z_bed(i,j))/dx, grad_max)
                dzbdy = clip((z_bed(i,j+1)-z_bed(i,j))/dx, grad_max)

                if (abs(dzbdx) .ge. grad_lim) then
                    mask_apply(i,j)   = .true.
                    mask_apply(i+1,j) = .true.
                end if

                if (abs(dzbdy) .ge. grad_lim) then
                    mask_apply(i,j)   = .true.
                    mask_apply(i,j+1) = .true.
                end if

            end do
            end do

            if (count(mask_apply) .eq. 0) exit

            ! Smooth z_bed at these locations, and H_ice so that it avoids spurious patterns
            call smooth_gauss_2D(z_bed, dx, 2.0_wp, mask_apply)
            call smooth_gauss_2D(H_ice, dx, 2.0_wp, mask_apply)

            if (present(n_iter)) n_iter = q

        end do

    end subroutine topoprep_limit_gradients

    elemental function clip(var, var_lim) result(var_c)
        ! var limited to [-var_lim, var_lim]
        real(wp), intent(in) :: var, var_lim
        real(wp) :: var_c

        var_c = var
        if (var .lt. -var_lim) then
            var_c = -var_lim
        else if (var .gt. var_lim) then
            var_c = var_lim
        end if

    end function clip

    subroutine smooth_gauss_2D(var, dx, f_sigma, mask_apply)
        ! Gaussian smoothing with sigma = f_sigma*dx (f_sigma >= 1) at the
        ! points of mask_apply, from all points of the window (mirrored at the
        ! borders). As Yelmo's smooth_gauss_2D with every point in mask_use,
        ! where the normalised filter is the same at all points.
        real(wp), intent(inout) :: var(:,:)
        real(wp), intent(in)    :: dx, f_sigma
        logical,  intent(in)    :: mask_apply(:,:)

        integer  :: i, j, nx, ny, n, n2
        real(wp) :: sigma
        real(wp), allocatable :: filter(:,:), var_old(:,:)

        nx = size(var,1)
        ny = size(var,2)

        if (f_sigma .lt. 1.0_wp) then
            call fesm_error("smooth_gauss_2D", "f_sigma must be >= 1.", "f_sigma = "//to_str(f_sigma))
        end if

        sigma = dx*f_sigma

        ! Half-width of the filter: 3 sigma
        n2 = 3*ceiling(f_sigma)
        n  = 2*n2+1

        if (nx .lt. n2 .or. ny .lt. n2) then
            call fesm_error("smooth_gauss_2D", "the grid is smaller than the filter half-width.", &
                            "nx, ny, half-width = "//to_str(nx)//" "//to_str(ny)//" "//to_str(n2))
        end if

        allocate(var_old(nx+2*n2,ny+2*n2))
        allocate(filter(n,n))

        filter = gauss_values(dx, dx, sigma, n)
        filter = filter/sum(filter)

        var_old = 0.0_wp
        var_old(n2+1:n2+nx,n2+1:n2+ny) = var

        ! Edge halos by mirror reflection about each border
        var_old(1:n2,n2+1:n2+ny)            = var(n2:1:-1,:)
        var_old(n2+nx+1:nx+2*n2,n2+1:n2+ny) = var(nx:nx-n2+1:-1,:)
        var_old(n2+1:n2+nx,1:n2)            = var(:,n2:1:-1)
        var_old(n2+1:n2+nx,n2+ny+1:ny+2*n2) = var(:,ny:ny-n2+1:-1)

        ! Corner halos by double mirror reflection
        var_old(1:n2,1:n2)                       = var(n2:1:-1,n2:1:-1)
        var_old(n2+nx+1:nx+2*n2,1:n2)            = var(nx:nx-n2+1:-1,n2:1:-1)
        var_old(1:n2,n2+ny+1:ny+2*n2)            = var(n2:1:-1,ny:ny-n2+1:-1)
        var_old(n2+nx+1:nx+2*n2,n2+ny+1:ny+2*n2) = var(nx:nx-n2+1:-1,ny:ny-n2+1:-1)

        do j = 1, ny
        do i = 1, nx
            if (mask_apply(i,j)) var(i,j) = sum(var_old(i:i+2*n2,j:j+2*n2)*filter)
        end do
        end do

    end subroutine smooth_gauss_2D

    function gauss_values(dx, dy, sigma, n) result(filt)
        ! 2D Gaussian smoothing kernel (n odd)
        ! https://en.wikipedia.org/wiki/Gaussian_blur
        real(wp), intent(in) :: dx, dy, sigma
        integer,  intent(in) :: n
        real(wp) :: filt(n,n)

        integer  :: n2, i, j, i1, j1
        real(wp) :: x, y
        real(wp), parameter :: pi_wp = real(pi,wp)

        n2 = (n-1)/2

        do j = -n2, n2
        do i = -n2, n2
            x = i*dx
            y = j*dy

            i1 = i+1+n2
            j1 = j+1+n2
            filt(i1,j1) = 1.0/(2.0*pi_wp*sigma**2)*exp(-(x**2+y**2)/(2*sigma**2))
        end do
        end do

        ! Normalize to ensure sum to 1
        filt = filt / sum(filt)

    end function gauss_values

end module topoprep
