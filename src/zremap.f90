module zremap
    ! Elevation remapping of a 2D field f (e.g. SMB, surface temperature).
    !
    ! From a reference state (f_ref, z_ref, mask_ref) a local profile f(z) is
    ! built around each table point, and the field is then remapped to the
    ! current elevation: f_now(i,j) = f(z_now(i,j)). This follows the lookup
    ! tables of Goelzer et al. (2020, TC), but the neighbourhood of a point is a
    ! Gaussian kernel of width sigma instead of a drainage basin (cf. the
    ! search radius of Helsen et al., 2012, TC), so no basins are needed.
    !
    ! Profile at a table point (kernel regression in x,y and z):
    !   * each valid reference point contributes, with the horizontal kernel
    !     weight, to the two elevation bins around z_ref (hat weights);
    !   * the bin sums are smoothed along z with a Gaussian of width sigma_z, so
    !     that sparse bins borrow from their neighbours;
    !   * a bin with data gives a node (z_bar, f_bar): the weighted mean
    !     elevation and value of the (smoothed) sums (exact for f linear in z);
    !   * the profile is the piecewise-linear interpolation of the nodes on the
    !     bins zb, extrapolated with gradient dfdz_lo below the lowest node and
    !     dfdz_hi above the highest one.
    ! The support W of a table point is the kernel weight of all its data (0-1:
    ! ~1 inside the data, ~0.5 at their edge, decaying over ~sigma beyond).
    ! A grid point interpolates bilinearly between the profiles of the table
    ! points around it, each weighted with its W, and blends the result with the
    ! fallback f_ref + dfdz*(z_now-z_ref), with dfdz_lo below z_ref and dfdz_hi
    ! above it, by alpha = min(1, W/w_full), where W
    ! is the bilinear support at the grid point. The remapped field is thus
    ! continuous and fades to the fallback away from the data.
    !
    ! Tables are stored every `stride` grid points (and at the last row/column).
    ! The grid is regular with spacing dx, and fields are on that same grid.
    !
    ! Namelist group (e.g. "zremap"):
    !   sigma    [m]      width of the horizontal Gaussian kernel (cut at 3 sigma)
    !   z_min    [m]      lowest elevation bin
    !   z_max    [m]      highest elevation bin
    !   dz       [m]      bin spacing
    !   sigma_z  [m]      width of the Gaussian smoothing along z (0 = none)
    !   w_full   [1]      support from which the profiles are fully used (e.g. 0.25)
    !   stride   [-]      table spacing in grid points
    !   dfdz_lo  [f m-1]  gradient below the data (and below z_ref in the fallback)
    !   dfdz_hi  [f m-1]  gradient above the data (and above z_ref in the fallback)

    use precision
    use constants, only : mv
    use nml
    use ncio

    implicit none
    private

    type zremap_par_class
        real(wp) :: sigma
        real(wp) :: z_min, z_max, dz
        real(wp) :: sigma_z
        real(wp) :: w_full
        integer  :: stride
        real(wp) :: dfdz_lo, dfdz_hi
    end type

    type zremap_class
        type(zremap_par_class) :: par
        integer  :: nx, ny              ! grid size
        real(wp) :: dx                  ! [m] grid spacing
        integer  :: nxo, nyo, nb        ! table points and bins
        integer  :: nk                  ! kernel half-width [points]
        real(wp), allocatable :: kern(:)        ! (-nk:nk) 1D kernel, sums to 1
        integer  :: nkz                 ! smoothing half-width along z [bins]
        real(wp), allocatable :: kernz(:)       ! (-nkz:nkz) smoothing kernel along z
        integer,  allocatable :: io(:), jo(:)   ! grid indices of the table points
        integer,  allocatable :: ki(:), kj(:)   ! table point at or left/below of each grid point
        real(wp), allocatable :: wi(:), wj(:)   ! bilinear weight of the next table point
        real(wp), allocatable :: zb(:)          ! [m] bin elevations
        real(wp), allocatable :: tab(:,:,:)     ! (nb,nxo,nyo) profile f(zb)
        real(wp), allocatable :: wt(:,:,:)      ! (nb,nxo,nyo) kernel weight of the data in each bin
        real(wp), allocatable :: wsup(:,:)      ! (nxo,nyo) support W (kernel weight of all data)
        real(wp), allocatable :: z_lo(:,:)      ! (nxo,nyo) [m] elevation range of the data (mean
        real(wp), allocatable :: z_hi(:,:)      ! (nxo,nyo) [m] elevation of the lowest/highest bin)
        real(wp), allocatable :: f_ref(:,:)     ! reference field
        real(wp), allocatable :: z_ref(:,:)     ! [m] reference elevation
        real(wp), allocatable :: w(:,:)         ! [1] support W at the grid points
        real(wp), allocatable :: dfdz(:,:)      ! [f m-1] gradient of the remapped field at z_ref
    end type

    public :: zremap_par_class, zremap_class
    public :: zremap_init, zremap_update, zremap_map, zremap_write, zremap_end

contains

    subroutine zremap_init(zr, filename, group, nx, ny, dx)
        ! Load the parameters and set up the kernel, the table points and the bins.
        type(zremap_class), intent(out) :: zr
        character(len=*),   intent(in)  :: filename     ! parameter file
        character(len=*),   intent(in)  :: group        ! namelist group, e.g. "zremap"
        integer,            intent(in)  :: nx, ny
        real(wp),           intent(in)  :: dx           ! [m] grid spacing

        integer :: k

        call zremap_par_load(zr%par, filename, group)

        if (dx <= 0.0_wp) call zremap_error("zremap_init", "dx must be > 0.")

        zr%nx = nx
        zr%ny = ny
        zr%dx = dx

        ! Horizontal kernel (separable Gaussian, cut at 3 sigma)
        zr%nk = ceiling(3.0_wp*zr%par%sigma/dx)
        allocate(zr%kern(-zr%nk:zr%nk))
        do k = -zr%nk, zr%nk
            zr%kern(k) = exp(-0.5_wp*(real(k,wp)*dx/zr%par%sigma)**2)
        end do
        zr%kern = zr%kern / sum(zr%kern)

        ! Table points and the bilinear weights of the grid points
        call table_axis(zr%io, zr%ki, zr%wi, nx, zr%par%stride)
        call table_axis(zr%jo, zr%kj, zr%wj, ny, zr%par%stride)
        zr%nxo = size(zr%io)
        zr%nyo = size(zr%jo)

        ! Elevation bins
        zr%nb = nint((zr%par%z_max-zr%par%z_min)/zr%par%dz) + 1
        allocate(zr%zb(zr%nb))
        do k = 1, zr%nb
            zr%zb(k) = zr%par%z_min + real(k-1,wp)*zr%par%dz
        end do

        ! Smoothing kernel along z (cut at 3 sigma_z; a single 1 without smoothing)
        zr%nkz = ceiling(3.0_wp*zr%par%sigma_z/zr%par%dz)
        allocate(zr%kernz(-zr%nkz:zr%nkz))
        zr%kernz = 1.0_wp
        if (zr%par%sigma_z > 0.0_wp) then
            do k = -zr%nkz, zr%nkz
                zr%kernz(k) = exp(-0.5_wp*(real(k,wp)*zr%par%dz/zr%par%sigma_z)**2)
            end do
        end if

        allocate(zr%tab(zr%nb,zr%nxo,zr%nyo))
        allocate(zr%wt(zr%nb,zr%nxo,zr%nyo))
        allocate(zr%wsup(zr%nxo,zr%nyo))
        allocate(zr%z_lo(zr%nxo,zr%nyo))
        allocate(zr%z_hi(zr%nxo,zr%nyo))
        allocate(zr%f_ref(nx,ny))
        allocate(zr%z_ref(nx,ny))
        allocate(zr%w(nx,ny))
        allocate(zr%dfdz(nx,ny))

        zr%tab   = mv
        zr%wt    = 0.0_wp
        zr%wsup  = 0.0_wp
        zr%z_lo  = mv
        zr%z_hi  = mv
        zr%f_ref = mv
        zr%z_ref = mv
        zr%w     = 0.0_wp
        zr%dfdz  = zr%par%dfdz_lo

    end subroutine zremap_init

    subroutine zremap_update(zr, f_ref, z_ref, mask_ref)
        ! Build the profiles from a reference state. Points with mask_ref and
        ! valid f_ref, z_ref are the data. Can be called again whenever the
        ! reference changes.
        type(zremap_class), intent(inout) :: zr
        real(wp),           intent(in)    :: f_ref(:,:)
        real(wp),           intent(in)    :: z_ref(:,:)    ! [m]
        logical,            intent(in)    :: mask_ref(:,:)

        real(wp), allocatable :: aw(:,:,:), af(:,:,:), az(:,:,:)
        real(wp), allocatable :: sw(:), sf(:), sz(:)
        real(wp) :: w, h1, wc(4)
        integer  :: i, j, ko, jj, b0, b1, nk, kk(4), ll(4)

        call check_shape(zr, f_ref, "zremap_update", "f_ref")
        call check_shape(zr, z_ref, "zremap_update", "z_ref")
        if (size(mask_ref,1) /= zr%nx .or. size(mask_ref,2) /= zr%ny) &
            call zremap_error("zremap_update", "mask_ref is not on the zremap grid.")

        zr%f_ref = f_ref
        zr%z_ref = z_ref
        nk = zr%nk

        ! Pass 1 (x): kernel sums along each row, at the table columns.
        ! Each data point falls in two bins (hat weights).
        allocate(aw(zr%nb,zr%nxo,zr%ny), af(zr%nb,zr%nxo,zr%ny), az(zr%nb,zr%nxo,zr%ny))
        aw = 0.0_wp
        af = 0.0_wp
        az = 0.0_wp

        do j = 1, zr%ny
        do i = 1, zr%nx
            if (.not. mask_ref(i,j) .or. f_ref(i,j) == mv .or. z_ref(i,j) == mv) cycle

            call bin_weights(zr, z_ref(i,j), b0, b1, h1)

            do ko = zr%ki(max(i-nk,1)), zr%nxo
                if (zr%io(ko) - i >  nk) exit
                if (i - zr%io(ko) >  nk) cycle
                w = zr%kern(i-zr%io(ko))
                aw(b0,ko,j) = aw(b0,ko,j) + w*(1.0_wp-h1)
                af(b0,ko,j) = af(b0,ko,j) + w*(1.0_wp-h1)*f_ref(i,j)
                az(b0,ko,j) = az(b0,ko,j) + w*(1.0_wp-h1)*z_ref(i,j)
                aw(b1,ko,j) = aw(b1,ko,j) + w*h1
                af(b1,ko,j) = af(b1,ko,j) + w*h1*f_ref(i,j)
                az(b1,ko,j) = az(b1,ko,j) + w*h1*z_ref(i,j)
            end do
        end do
        end do

        ! Pass 2 (y): kernel sums along the columns at the table points, then
        ! the profile from the bin nodes.
        allocate(sw(zr%nb), sf(zr%nb), sz(zr%nb))

        do j = 1, zr%nyo
        do i = 1, zr%nxo
            sw = 0.0_wp
            sf = 0.0_wp
            sz = 0.0_wp
            do jj = max(zr%jo(j)-nk,1), min(zr%jo(j)+nk,zr%ny)
                w  = zr%kern(jj-zr%jo(j))
                sw = sw + w*aw(:,i,jj)
                sf = sf + w*af(:,i,jj)
                sz = sz + w*az(:,i,jj)
            end do

            zr%wt(:,i,j) = sw
            zr%wsup(i,j) = sum(sw)
            call calc_profile(zr%tab(:,i,j), sw > 0.0_wp, &
                              smooth_z(sw, zr%kernz), smooth_z(sf, zr%kernz), smooth_z(sz, zr%kernz), &
                              zr%zb, zr%par%dfdz_lo, zr%par%dfdz_hi)

            ! Elevation range of the data: mean elevation of the lowest and
            ! highest bin with data (before smoothing)
            zr%z_lo(i,j) = mv
            zr%z_hi(i,j) = mv
            if (zr%wsup(i,j) > 0.0_wp) then
                b0 = findloc(sw > 0.0_wp, .true., 1)
                b1 = findloc(sw > 0.0_wp, .true., 1, back=.true.)
                zr%z_lo(i,j) = sz(b0)/sw(b0)
                zr%z_hi(i,j) = sz(b1)/sw(b1)
            end if
        end do
        end do

        ! Support at the grid points
        do j = 1, zr%ny
        do i = 1, zr%nx
            call corner_weights(zr, i, j, kk, ll, wc)
            zr%w(i,j) = sum(wc)
        end do
        end do

        ! Diagnostic: gradient of the remapped field at the reference elevation
        call calc_dfdz(zr)

    end subroutine zremap_update

    subroutine zremap_map(zr, f_now, z_now, dz_ext)
        ! Remap the field to the current elevation. Optionally, dz_ext is the
        ! elevation distance beyond the data range of the profiles: z_now-z_hi
        ! above it (>0), z_now-z_lo below it (<0), 0 within it, mv without data.
        type(zremap_class), intent(in)  :: zr
        real(wp),           intent(out) :: f_now(:,:)
        real(wp),           intent(in)  :: z_now(:,:)    ! [m]
        real(wp),           intent(out), optional :: dz_ext(:,:)  ! [m]

        integer :: i, j

        call check_shape(zr, f_now, "zremap_map", "f_now")
        call check_shape(zr, z_now, "zremap_map", "z_now")
        if (present(dz_ext)) call check_shape(zr, dz_ext, "zremap_map", "dz_ext")

        do j = 1, zr%ny
        do i = 1, zr%nx
            f_now(i,j) = eval_point(zr, i, j, z_now(i,j))
        end do
        end do

        if (present(dz_ext)) then
            do j = 1, zr%ny
            do i = 1, zr%nx
                dz_ext(i,j) = calc_dz_ext(zr, i, j, z_now(i,j))
            end do
            end do
        end if

    end subroutine zremap_map

    subroutine zremap_write(zr, filename)
        ! Write the tables (on the table points) and the 2D reference fields
        ! with the support w and the gradient dfdz. x,y are distances [km]
        ! from the first grid point.
        type(zremap_class), intent(in) :: zr
        character(len=*),   intent(in) :: filename

        real(wp), allocatable :: tmp(:,:,:)
        integer :: b, i

        call nc_create(filename)
        call nc_write_dim(filename, "xc", x=[(real(i-1,wp)*zr%dx*1e-3_wp, i=1,zr%nx)], units="km")
        call nc_write_dim(filename, "yc", x=[(real(i-1,wp)*zr%dx*1e-3_wp, i=1,zr%ny)], units="km")
        call nc_write_dim(filename, "xo", x=real(zr%io-1,wp)*zr%dx*1e-3_wp, units="km")
        call nc_write_dim(filename, "yo", x=real(zr%jo-1,wp)*zr%dx*1e-3_wp, units="km")
        call nc_write_dim(filename, "zb", x=zr%zb, units="m")

        call nc_write(filename, "f_ref", zr%f_ref, dim1="xc", dim2="yc", missing_value=mv, &
                      long_name="Reference field")
        call nc_write(filename, "z_ref", zr%z_ref, dim1="xc", dim2="yc", missing_value=mv, &
                      units="m", long_name="Reference elevation")
        call nc_write(filename, "w", zr%w, dim1="xc", dim2="yc", &
                      units="1", long_name="Support (kernel weight of the data)")
        call nc_write(filename, "dfdz", zr%dfdz, dim1="xc", dim2="yc", missing_value=mv, &
                      units="m-1", long_name="Gradient of the remapped field at the reference elevation")

        call nc_write(filename, "w_sup", zr%wsup, dim1="xo", dim2="yo", &
                      units="1", long_name="Support (kernel weight of the data)")
        call nc_write(filename, "z_lo", zr%z_lo, dim1="xo", dim2="yo", missing_value=mv, &
                      units="m", long_name="Lowest elevation of the data")
        call nc_write(filename, "z_hi", zr%z_hi, dim1="xo", dim2="yo", missing_value=mv, &
                      units="m", long_name="Highest elevation of the data")

        allocate(tmp(zr%nxo,zr%nyo,zr%nb))
        do b = 1, zr%nb
            tmp(:,:,b) = zr%tab(b,:,:)
        end do
        call nc_write(filename, "f_tab", tmp, dim1="xo", dim2="yo", dim3="zb", missing_value=mv, &
                      long_name="Profile f(zb)")
        do b = 1, zr%nb
            tmp(:,:,b) = zr%wt(b,:,:)
        end do
        call nc_write(filename, "w_tab", tmp, dim1="xo", dim2="yo", dim3="zb", &
                      long_name="Kernel weight of the data in each bin")

    end subroutine zremap_write

    subroutine zremap_end(zr)
        type(zremap_class), intent(inout) :: zr

        type(zremap_class) :: zr0

        zr = zr0

    end subroutine zremap_end

    ! === Internal ===

    subroutine corner_weights(zr, i, j, kk, ll, wc)
        ! The four table points around grid point (i,j), with their bilinear
        ! weights times their support. sum(wc) is the support W at (i,j).
        type(zremap_class), intent(in)  :: zr
        integer,            intent(in)  :: i, j
        integer,            intent(out) :: kk(4), ll(4)
        real(wp),           intent(out) :: wc(4)

        integer  :: k1, k2, l1, l2

        k1 = zr%ki(i)
        k2 = min(k1+1, zr%nxo)
        l1 = zr%kj(j)
        l2 = min(l1+1, zr%nyo)
        kk = [k1, k2, k1, k2]
        ll = [l1, l1, l2, l2]
        wc = [(1.0_wp-zr%wi(i))*(1.0_wp-zr%wj(j)), zr%wi(i)*(1.0_wp-zr%wj(j)), &
              (1.0_wp-zr%wi(i))*zr%wj(j),          zr%wi(i)*zr%wj(j)]
        wc = wc * [zr%wsup(k1,l1), zr%wsup(k2,l1), zr%wsup(k1,l2), zr%wsup(k2,l2)]

    end subroutine corner_weights

    function eval_point(zr, i, j, z) result(f)
        ! Field at grid point (i,j) and elevation z: the support-weighted
        ! profiles of the surrounding table points, blended with the fallback
        ! f_ref + dfdz*(z-z_ref) (dfdz_lo below z_ref, dfdz_hi above) by
        ! alpha = min(1, W/w_full).
        type(zremap_class), intent(in) :: zr
        integer,            intent(in) :: i, j
        real(wp),           intent(in) :: z
        real(wp) :: f

        integer  :: kk(4), ll(4), c
        real(wp) :: wc(4), wsum, alpha, f_fb
        logical  :: has_fb

        f = mv
        if (z == mv) return

        call corner_weights(zr, i, j, kk, ll, wc)
        wsum = sum(wc)

        has_fb = zr%f_ref(i,j) /= mv .and. zr%z_ref(i,j) /= mv
        if (has_fb) then
            if (z < zr%z_ref(i,j)) then
                f_fb = zr%f_ref(i,j) + zr%par%dfdz_lo*(z - zr%z_ref(i,j))
            else
                f_fb = zr%f_ref(i,j) + zr%par%dfdz_hi*(z - zr%z_ref(i,j))
            end if
        end if

        if (wsum > 0.0_wp) then
            f = 0.0_wp
            do c = 1, 4
                if (wc(c) > 0.0_wp) f = f + wc(c)*eval_profile(zr, zr%tab(:,kk(c),ll(c)), z)
            end do
            f = f / wsum
            if (has_fb) then
                alpha = min(1.0_wp, wsum/zr%par%w_full)
                f = alpha*f + (1.0_wp-alpha)*f_fb
            end if
        else if (has_fb) then
            f = f_fb
        end if

    end function eval_point

    function calc_dz_ext(zr, i, j, z) result(dz)
        ! Elevation distance of z beyond the (support-weighted) data range at
        ! grid point (i,j); mv without data.
        type(zremap_class), intent(in) :: zr
        integer,            intent(in) :: i, j
        real(wp),           intent(in) :: z
        real(wp) :: dz

        integer  :: kk(4), ll(4), c
        real(wp) :: wc(4), wsum, zlo, zhi

        dz = mv
        if (z == mv) return

        call corner_weights(zr, i, j, kk, ll, wc)
        wsum = sum(wc)
        if (wsum <= 0.0_wp) return

        zlo = 0.0_wp
        zhi = 0.0_wp
        do c = 1, 4
            if (wc(c) <= 0.0_wp) cycle
            zlo = zlo + wc(c)*zr%z_lo(kk(c),ll(c))
            zhi = zhi + wc(c)*zr%z_hi(kk(c),ll(c))
        end do
        zlo = zlo / wsum
        zhi = zhi / wsum

        dz = 0.0_wp
        if (z > zhi) dz = z - zhi
        if (z < zlo) dz = z - zlo

    end function calc_dz_ext

    function eval_profile(zr, tab, z) result(f)
        ! Linear interpolation of a profile on the bins, extrapolated with the
        ! end gradients beyond them.
        type(zremap_class), intent(in) :: zr
        real(wp),           intent(in) :: tab(:)
        real(wp),           intent(in) :: z
        real(wp) :: f

        real(wp) :: s, h
        integer  :: b

        if (z <= zr%zb(1)) then
            f = tab(1) + zr%par%dfdz_lo*(z - zr%zb(1))
        else if (z >= zr%zb(zr%nb)) then
            f = tab(zr%nb) + zr%par%dfdz_hi*(z - zr%zb(zr%nb))
        else
            s = (z - zr%par%z_min)/zr%par%dz
            b = min(int(s) + 1, zr%nb-1)
            h = s - real(b-1,wp)
            f = (1.0_wp-h)*tab(b) + h*tab(b+1)
        end if

    end function eval_profile

    function smooth_z(s, kernz) result(ss)
        ! Bin sums smoothed along z (truncated at the ends: the sums of all bins
        ! within the kernel).
        real(wp), intent(in) :: s(:)
        real(wp), intent(in) :: kernz(:)        ! offsets -nkz..nkz
        real(wp) :: ss(size(s))

        integer :: nb, nkz, b, k

        nb  = size(s)
        nkz = (size(kernz)-1)/2
        do b = 1, nb
            ss(b) = 0.0_wp
            do k = max(-nkz, 1-b), min(nkz, nb-b)
                ss(b) = ss(b) + kernz(k+nkz+1)*s(b+k)
            end do
        end do

    end function smooth_z

    subroutine calc_profile(tab, has, sw, sf, sz, zb, dfdz_lo, dfdz_hi)
        ! Profile on the bins zb from the (smoothed) bin sums: nodes
        ! (z_bar,f_bar) of the bins with data `has` (kept only where z_bar
        ! increases), linear between nodes and extrapolated with dfdz_lo/hi
        ! beyond them (mv without data).
        real(wp), intent(out) :: tab(:)
        logical,  intent(in)  :: has(:)                 ! bin has data (before smoothing)
        real(wp), intent(in)  :: sw(:), sf(:), sz(:)    ! sums of w, w*f, w*z per bin
        real(wp), intent(in)  :: zb(:)
        real(wp), intent(in)  :: dfdz_lo, dfdz_hi

        real(wp) :: zn(size(zb)), fn(size(zb))
        integer  :: nb, nn, b, n

        nb = size(zb)

        nn = 0
        do b = 1, nb
            if (.not. has(b) .or. sw(b) <= 0.0_wp) cycle
            if (nn > 0) then
                if (sz(b)/sw(b) <= zn(nn)) cycle
            end if
            nn = nn + 1
            zn(nn) = sz(b)/sw(b)
            fn(nn) = sf(b)/sw(b)
        end do

        if (nn == 0) then
            tab = mv
            return
        end if

        n = 1
        do b = 1, nb
            if (zb(b) <= zn(1)) then
                tab(b) = fn(1) + dfdz_lo*(zb(b) - zn(1))
            else if (zb(b) >= zn(nn)) then
                tab(b) = fn(nn) + dfdz_hi*(zb(b) - zn(nn))
            else
                do while (zn(n+1) < zb(b))
                    n = n + 1
                end do
                tab(b) = fn(n) + (fn(n+1)-fn(n))*(zb(b)-zn(n))/(zn(n+1)-zn(n))
            end if
        end do

    end subroutine calc_profile

    subroutine calc_dfdz(zr)
        ! Gradient of the remapped field at the reference elevation (central
        ! difference over one bin).
        type(zremap_class), intent(inout) :: zr

        real(wp) :: h
        integer  :: i, j

        h = 0.5_wp*zr%par%dz
        do j = 1, zr%ny
        do i = 1, zr%nx
            if (zr%z_ref(i,j) == mv) then
                zr%dfdz(i,j) = mv
            else
                zr%dfdz(i,j) = (eval_point(zr, i, j, zr%z_ref(i,j)+h) &
                              - eval_point(zr, i, j, zr%z_ref(i,j)-h)) / (2.0_wp*h)
            end if
        end do
        end do

    end subroutine calc_dfdz
    subroutine bin_weights(zr, z, b0, b1, h1)
        ! The two bins around z and the hat weight of the upper one. Elevations
        ! outside [z_min,z_max] go to the end bin.
        type(zremap_class), intent(in)  :: zr
        real(wp),           intent(in)  :: z
        integer,            intent(out) :: b0, b1
        real(wp),           intent(out) :: h1

        real(wp) :: s

        s  = (z - zr%par%z_min)/zr%par%dz
        s  = max(0.0_wp, min(s, real(zr%nb-1,wp)))
        b0 = min(int(s) + 1, zr%nb)
        b1 = min(b0 + 1, zr%nb)
        h1 = s - real(b0-1,wp)

    end subroutine bin_weights

    subroutine table_axis(ix, k, w, n, stride)
        ! Table points along an axis of n grid points (every stride points,
        ! plus the last one), and for each grid point the table point at or
        ! before it with the bilinear weight of the next one.
        integer,  allocatable, intent(out) :: ix(:), k(:)
        real(wp), allocatable, intent(out) :: w(:)
        integer,  intent(in) :: n, stride

        integer :: no, i, m

        no = (n-1)/stride + 1
        if (1 + (no-1)*stride < n) no = no + 1

        allocate(ix(no), k(n), w(n))
        do m = 1, no
            ix(m) = min(1 + (m-1)*stride, n)
        end do

        do i = 1, n
            m = min((i-1)/stride + 1, max(no-1,1))
            k(i) = m
            if (no == 1) then
                w(i) = 0.0_wp
            else
                w(i) = real(i-ix(m),wp)/real(ix(m+1)-ix(m),wp)
            end if
        end do

    end subroutine table_axis

    subroutine check_shape(zr, var, proc, name)
        type(zremap_class), intent(in) :: zr
        real(wp),           intent(in) :: var(:,:)
        character(len=*),   intent(in) :: proc, name

        if (size(var,1) /= zr%nx .or. size(var,2) /= zr%ny) &
            call zremap_error(proc, trim(name)//" is not on the zremap grid.")

    end subroutine check_shape

    subroutine zremap_par_load(par, filename, group)
        type(zremap_par_class), intent(out) :: par
        character(len=*),       intent(in)  :: filename, group

        call nml_read(filename, group, "sigma",   par%sigma)
        call nml_read(filename, group, "z_min",   par%z_min)
        call nml_read(filename, group, "z_max",   par%z_max)
        call nml_read(filename, group, "dz",      par%dz)
        call nml_read(filename, group, "sigma_z", par%sigma_z)
        call nml_read(filename, group, "w_full",  par%w_full)
        call nml_read(filename, group, "stride",  par%stride)
        call nml_read(filename, group, "dfdz_lo", par%dfdz_lo)
        call nml_read(filename, group, "dfdz_hi", par%dfdz_hi)

        if (par%sigma <= 0.0_wp)      call zremap_error("zremap_par_load", "sigma must be > 0.")
        if (par%dz <= 0.0_wp)         call zremap_error("zremap_par_load", "dz must be > 0.")
        if (par%sigma_z < 0.0_wp)     call zremap_error("zremap_par_load", "sigma_z must be >= 0.")
        if (par%z_max <= par%z_min)   call zremap_error("zremap_par_load", "z_max must be > z_min.")
        if (par%stride < 1)           call zremap_error("zremap_par_load", "stride must be >= 1.")
        if (par%w_full <= 0.0_wp)     call zremap_error("zremap_par_load", "w_full must be > 0.")

    end subroutine zremap_par_load

    subroutine zremap_error(proc, msg)
        ! Abort with a framed message.
        character(len=*), intent(in) :: proc, msg

        write(*,*) ""
        write(*,*) "zremap:: error in "//trim(proc)
        write(*,*) "    "//trim(msg)
        error stop 1

    end subroutine zremap_error

end module zremap
