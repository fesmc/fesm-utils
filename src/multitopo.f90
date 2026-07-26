module multitopo
    ! Multi-resolution topography carrier.
    !
    ! Keeps one evolving topography (z_bed, H_ice, z_srf) consistent across two
    ! grids of different resolution, anchored to a reference state known on both
    ! grids. A driving model (ice-sheet / solid-earth) evolves the fields on one
    ! grid; multitopo_update maps that evolution onto the other grid so a
    ! consumer at either resolution sees a consistent surface.
    !
    ! The mapping is direction-agnostic (dispatched on the relative resolution of
    ! the two grids):
    !
    !   coarse -> fine : reference-anchored bedrock anomaly + a fill-level ice
    !                    reconstruction that conserves each coarse cell's mean
    !                    ice thickness. This injects the sub-grid detail the
    !                    coarse field lacks, and reproduces the physical
    !                    smooth(thick ice) <-> rough(exposed bedrock) transition.
    !   fine   -> coarse : conservative area aggregation.
    !   equal          : identity copy (single-resolution models pay nothing).
    !
    ! Band binning / lapse downscaling is deliberately NOT here -- a consumer
    ! queries the per-cell sub-grid surface (multitopo_cell_subgrid) and bins it
    ! however it likes.
    !
    ! v1: two grids, regular in x/y, sharing a projection. N-level and fields
    ! arriving on different grids are natural extensions of the pairwise API.

    use precision, only : wp, dp

    implicit none

    private

    integer, parameter, public :: MT_COARSE = 1
    integer, parameter, public :: MT_FINE   = 2

    real(dp), parameter :: FILL_TOL = 1.0e-6_dp   ! [m] fill-level bisection tol

    type multitopo_level_class
        integer               :: nx = 0, ny = 0
        real(wp)              :: dx = 0.0_wp, dy = 0.0_wp
        real(wp), allocatable :: xc(:), yc(:)          ! cell-center coords
        real(wp), allocatable :: z_bed(:,:)            ! current
        real(wp), allocatable :: H_ice(:,:)
        real(wp), allocatable :: z_srf(:,:)
        real(wp), allocatable :: z_bed_ref(:,:)        ! reference (t0)
        real(wp), allocatable :: H_ice_ref(:,:)
    end type multitopo_level_class

    type multitopo_class
        type(multitopo_level_class) :: crs   ! the coarser level
        type(multitopo_level_class) :: fin   ! the finer level
        logical :: noop = .FALSE.            ! the two grids are identical

        ! Fine -> coarse parent map + CSR grouping of fine cells per coarse cell.
        integer, allocatable :: pic(:,:), pjc(:,:)     ! (nxf,nyf) parent coarse i,j
        integer, allocatable :: cstart(:)              ! (ncrs+1) CSR offsets
        integer, allocatable :: flist(:)               ! (nfin) fine linear idx, grouped
        integer :: nsub_max = 1                        ! max fine cells in any coarse cell
    end type multitopo_class

    public :: multitopo_class
    public :: multitopo_init
    public :: multitopo_update
    public :: multitopo_cell_subgrid
    public :: multitopo_cell_bands
    public :: multitopo_end

contains

    subroutine multitopo_init(mt, xcA, ycA, zbedA, HiceA, xcB, ycB, zbedB, HiceB)
        ! Register two grids A and B with their reference topography. Which is
        ! coarse vs fine is decided from the grid spacing; identical grids -> no-op.
        implicit none
        type(multitopo_class), intent(INOUT) :: mt
        real(wp), intent(IN) :: xcA(:), ycA(:), zbedA(:,:), HiceA(:,:)
        real(wp), intent(IN) :: xcB(:), ycB(:), zbedB(:,:), HiceB(:,:)

        real(wp) :: dxA, dxB
        logical  :: same

        call multitopo_end(mt)

        dxA = grid_spacing(xcA)
        dxB = grid_spacing(xcB)

        same = (size(xcA) .eq. size(xcB)) .and. (size(ycA) .eq. size(ycB)) .and. &
               (abs(dxA - dxB) .le. epsilon(1.0_wp)*max(abs(dxA),1.0_wp))

        if (dxB .lt. dxA) then
            call level_set(mt%crs, xcA, ycA, zbedA, HiceA)
            call level_set(mt%fin, xcB, ycB, zbedB, HiceB)
        else
            call level_set(mt%crs, xcB, ycB, zbedB, HiceB)
            call level_set(mt%fin, xcA, ycA, zbedA, HiceA)
        end if

        mt%noop = same
        if (same) then
            ! Both point at the same grid; parent map is the identity.
            call level_copy(mt%fin, mt%crs)
        end if

        call build_parent_map(mt)

        return
    end subroutine multitopo_init

    subroutine multitopo_update(mt, level, z_bed, H_ice)
        ! Set the evolved (z_bed, H_ice) on the driving level and reconstruct the
        ! other level. level is MT_COARSE or MT_FINE (the grid the fields are on).
        implicit none
        type(multitopo_class), intent(INOUT) :: mt
        integer,               intent(IN)    :: level
        real(wp),              intent(IN)    :: z_bed(:,:), H_ice(:,:)

        select case (level)
            case (MT_COARSE)
                mt%crs%z_bed = z_bed
                mt%crs%H_ice = H_ice
                mt%crs%z_srf = z_bed + H_ice
                if (mt%noop) then
                    call level_copy(mt%fin, mt%crs)
                else
                    call downscale(mt)          ! coarse -> fine
                end if
            case (MT_FINE)
                mt%fin%z_bed = z_bed
                mt%fin%H_ice = H_ice
                mt%fin%z_srf = z_bed + H_ice
                if (mt%noop) then
                    call level_copy(mt%crs, mt%fin)
                else
                    call aggregate(mt)          ! fine -> coarse
                end if
            case default
                error stop "multitopo_update:: level must be MT_COARSE or MT_FINE."
        end select

        return
    end subroutine multitopo_update

    subroutine multitopo_cell_subgrid(mt, ic, jc, nsub, zsrf, fice)
        ! Return the sub-grid surface of coarse cell (ic,jc): the fine-cell
        ! surface elevations zsrf(1:nsub) and their ice fractions fice(1:nsub)
        ! (1 = ice, 0 = ice-free). A consumer bins these into elevation classes.
        !
        ! If the coarse cell has no fine cells (edge, or no-op), a single sub-cell
        ! equal to the coarse cell itself is returned, so a consumer always gets
        ! nsub >= 1.
        implicit none
        type(multitopo_class), intent(IN)  :: mt
        integer,               intent(IN)  :: ic, jc
        integer,               intent(OUT) :: nsub
        real(wp),              intent(OUT) :: zsrf(:), fice(:)

        integer :: c, k, p, fi, ifn, jfn, nxc

        nxc  = mt%crs%nx
        nsub = 0

        if (mt%noop) then
            nsub = 1
            zsrf(1) = mt%crs%z_srf(ic,jc)
            fice(1) = merge(1.0_wp, 0.0_wp, mt%crs%H_ice(ic,jc) .gt. 0.0_wp)
            return
        end if

        c = (jc-1)*nxc + ic
        do p = mt%cstart(c), mt%cstart(c+1)-1
            fi  = mt%flist(p)
            ifn = mod(fi-1, mt%fin%nx) + 1
            jfn = (fi-1)/mt%fin%nx + 1
            nsub = nsub + 1
            if (nsub .gt. size(zsrf)) error stop "multitopo_cell_subgrid:: output array too small."
            zsrf(nsub) = mt%fin%z_srf(ifn,jfn)
            fice(nsub) = merge(1.0_wp, 0.0_wp, mt%fin%H_ice(ifn,jfn) .gt. 0.0_wp)
        end do

        if (nsub .eq. 0) then       ! coarse cell with no fine coverage
            nsub = 1
            zsrf(1) = mt%crs%z_srf(ic,jc)
            fice(1) = merge(1.0_wp, 0.0_wp, mt%crs%H_ice(ic,jc) .gt. 0.0_wp)
        end if

        return
    end subroutine multitopo_cell_subgrid

    subroutine multitopo_cell_bands(mt, ic, jc, nband, zband, wband, ficeband)
        ! Equal-area elevation classes for coarse cell (ic,jc): partition the
        ! cell's sub-grid surface (sorted by elevation) into nband equal-count
        ! groups. Per band: mean elevation, area weight (group count / total),
        ! mean ice fraction. Weights sum to 1. If the cell has fewer sub-cells
        ! than nband, the surplus bands get zero weight (elevation = cell value)
        ! and a consumer skips them.
        !
        ! Band *count* and any downstream lapse downscaling are the consumer's
        ! choice; this is only the generic binning.
        implicit none
        type(multitopo_class), intent(IN)  :: mt
        integer,               intent(IN)  :: ic, jc, nband
        real(wp),              intent(OUT) :: zband(:), wband(:), ficeband(:)

        real(wp), allocatable :: z(:), f(:)
        integer  :: nsub, b, i0, i1, k, m
        real(dp) :: zs, fs

        allocate(z(mt%nsub_max), f(mt%nsub_max))
        call multitopo_cell_subgrid(mt, ic, jc, nsub, z, f)
        call sort_pair(z(1:nsub), f(1:nsub))

        zband = 0.0_wp; wband = 0.0_wp; ficeband = 0.0_wp
        do b = 1, nband
            i0 = (b-1)*nsub/nband + 1
            i1 = b*nsub/nband
            m  = i1 - i0 + 1
            if (m .le. 0) then
                zband(b)    = mt%crs%z_srf(ic,jc)   ! empty band: placeholder, w=0
                wband(b)    = 0.0_wp
                ficeband(b) = 0.0_wp
                cycle
            end if
            zs = 0.0_dp; fs = 0.0_dp
            do k = i0, i1
                zs = zs + real(z(k),dp); fs = fs + real(f(k),dp)
            end do
            zband(b)    = real(zs/real(m,dp),wp)
            wband(b)    = real(m,wp)/real(nsub,wp)
            ficeband(b) = real(fs/real(m,dp),wp)
        end do

        deallocate(z, f)

        return
    end subroutine multitopo_cell_bands

    subroutine sort_pair(z, f)
        ! Insertion sort of z ascending, carrying f. Sub-cell counts are small
        ! (tens at 8 km); revisit if very high-res sub-grids are used.
        implicit none
        real(wp), intent(INOUT) :: z(:), f(:)
        integer  :: i, j
        real(wp) :: zt, ft
        do i = 2, size(z)
            zt = z(i); ft = f(i); j = i - 1
            do while (j .ge. 1)
                if (z(j) .le. zt) exit
                z(j+1) = z(j); f(j+1) = f(j); j = j - 1
            end do
            z(j+1) = zt; f(j+1) = ft
        end do
        return
    end subroutine sort_pair

    subroutine multitopo_end(mt)
        implicit none
        type(multitopo_class), intent(INOUT) :: mt
        call level_dealloc(mt%crs)
        call level_dealloc(mt%fin)
        if (allocated(mt%pic))    deallocate(mt%pic)
        if (allocated(mt%pjc))    deallocate(mt%pjc)
        if (allocated(mt%cstart)) deallocate(mt%cstart)
        if (allocated(mt%flist))  deallocate(mt%flist)
        mt%noop = .FALSE.
        return
    end subroutine multitopo_end

    ! === internals ==========================================================

    subroutine downscale(mt)
        ! coarse -> fine: bedrock anomaly onto the fine reference, then a
        ! per-coarse-cell fill-level ice reconstruction conserving the mean.
        implicit none
        type(multitopo_class), intent(INOUT) :: mt

        integer  :: ic, jc, c, p, fi, ifn, jfn, nf, nxc
        real(wp) :: dbed, Hbar, S
        real(wp), allocatable :: bed(:)

        nxc = mt%crs%nx

        ! Bedrock: reference + coarse anomaly (piecewise constant per parent).
        do jfn = 1, mt%fin%ny
        do ifn = 1, mt%fin%nx
            ic = mt%pic(ifn,jfn); jc = mt%pjc(ifn,jfn)
            dbed = mt%crs%z_bed(ic,jc) - mt%crs%z_bed_ref(ic,jc)
            mt%fin%z_bed(ifn,jfn) = mt%fin%z_bed_ref(ifn,jfn) + dbed
        end do
        end do

        ! Ice: fill each coarse cell's fine beds to the level that matches the
        ! coarse mean thickness.
        allocate(bed(max_cell_count(mt)))
        do jc = 1, mt%crs%ny
        do ic = 1, mt%crs%nx
            c  = (jc-1)*nxc + ic
            nf = mt%cstart(c+1) - mt%cstart(c)
            if (nf .eq. 0) cycle
            Hbar = mt%crs%H_ice(ic,jc)
            do p = mt%cstart(c), mt%cstart(c+1)-1
                fi = mt%flist(p)
                ifn = mod(fi-1, mt%fin%nx) + 1
                jfn = (fi-1)/mt%fin%nx + 1
                bed(p - mt%cstart(c) + 1) = mt%fin%z_bed(ifn,jfn)
            end do
            S = fill_level(bed(1:nf), real(Hbar,dp))
            do p = mt%cstart(c), mt%cstart(c+1)-1
                fi = mt%flist(p)
                ifn = mod(fi-1, mt%fin%nx) + 1
                jfn = (fi-1)/mt%fin%nx + 1
                mt%fin%H_ice(ifn,jfn) = max(0.0_wp, S - mt%fin%z_bed(ifn,jfn))
                mt%fin%z_srf(ifn,jfn) = mt%fin%z_bed(ifn,jfn) + mt%fin%H_ice(ifn,jfn)
            end do
        end do
        end do
        deallocate(bed)

        return
    end subroutine downscale

    subroutine aggregate(mt)
        ! fine -> coarse: conservative area-average (equal-area fine cells).
        implicit none
        type(multitopo_class), intent(INOUT) :: mt

        integer  :: ic, jc, c, p, fi, ifn, jfn, nf, nxc
        real(dp) :: sbed, sice

        nxc = mt%crs%nx
        do jc = 1, mt%crs%ny
        do ic = 1, mt%crs%nx
            c  = (jc-1)*nxc + ic
            nf = mt%cstart(c+1) - mt%cstart(c)
            if (nf .eq. 0) cycle
            sbed = 0.0_dp; sice = 0.0_dp
            do p = mt%cstart(c), mt%cstart(c+1)-1
                fi = mt%flist(p)
                ifn = mod(fi-1, mt%fin%nx) + 1
                jfn = (fi-1)/mt%fin%nx + 1
                sbed = sbed + real(mt%fin%z_bed(ifn,jfn),dp)
                sice = sice + real(mt%fin%H_ice(ifn,jfn),dp)
            end do
            mt%crs%z_bed(ic,jc) = real(sbed/real(nf,dp),wp)
            mt%crs%H_ice(ic,jc) = real(sice/real(nf,dp),wp)
            mt%crs%z_srf(ic,jc) = mt%crs%z_bed(ic,jc) + mt%crs%H_ice(ic,jc)
        end do
        end do

        return
    end subroutine aggregate

    function fill_level(bed, Hbar) result(S)
        ! Smooth ice-surface level S such that mean(max(0, S-bed)) = Hbar over the
        ! sub-cell beds. Monotone in S -> bisection. Hbar<=0 -> S = min(bed)
        ! (no ice).
        implicit none
        real(wp), intent(IN) :: bed(:)
        real(dp), intent(IN) :: Hbar
        real(wp) :: S

        real(dp) :: lo, hi, mid, fmid
        integer  :: it

        lo = real(minval(bed),dp)
        if (Hbar .le. 0.0_dp) then
            S = real(lo,wp)
            return
        end if
        hi = real(maxval(bed),dp) + Hbar
        do it = 1, 100
            mid  = 0.5_dp*(lo+hi)
            fmid = mean_thickness(bed, mid) - Hbar
            if (fmid .gt. 0.0_dp) then
                hi = mid
            else
                lo = mid
            end if
            if (hi-lo .le. FILL_TOL) exit
        end do
        S = real(0.5_dp*(lo+hi),wp)

        return
    end function fill_level

    pure function mean_thickness(bed, S) result(m)
        implicit none
        real(wp), intent(IN) :: bed(:)
        real(dp), intent(IN) :: S
        real(dp) :: m
        integer  :: k
        m = 0.0_dp
        do k = 1, size(bed)
            m = m + max(0.0_dp, S - real(bed(k),dp))
        end do
        m = m/real(size(bed),dp)
        return
    end function mean_thickness

    subroutine build_parent_map(mt)
        ! Assign each fine cell to the coarse cell whose center is nearest, then
        ! build the CSR grouping (cstart/flist) of fine cells per coarse cell.
        implicit none
        type(multitopo_class), intent(INOUT) :: mt

        integer :: ifn, jfn, ic, jc, c, nxc, nyc, ncrs, fi
        integer, allocatable :: cnt(:), pos(:)
        real(wp) :: x0c, y0c, dxc, dyc

        nxc = mt%crs%nx; nyc = mt%crs%ny; ncrs = nxc*nyc
        x0c = mt%crs%xc(1); y0c = mt%crs%yc(1)
        dxc = mt%crs%dx;    dyc = mt%crs%dy

        if (allocated(mt%pic)) deallocate(mt%pic)
        if (allocated(mt%pjc)) deallocate(mt%pjc)
        allocate(mt%pic(mt%fin%nx,mt%fin%ny), mt%pjc(mt%fin%nx,mt%fin%ny))

        do jfn = 1, mt%fin%ny
        do ifn = 1, mt%fin%nx
            ic = nint((mt%fin%xc(ifn) - x0c)/dxc) + 1
            jc = nint((mt%fin%yc(jfn) - y0c)/dyc) + 1
            mt%pic(ifn,jfn) = min(max(ic,1),nxc)
            mt%pjc(ifn,jfn) = min(max(jc,1),nyc)
        end do
        end do

        ! CSR: count per coarse cell, prefix-sum, scatter.
        allocate(cnt(ncrs), pos(ncrs))
        cnt = 0
        do jfn = 1, mt%fin%ny
        do ifn = 1, mt%fin%nx
            c = (mt%pjc(ifn,jfn)-1)*nxc + mt%pic(ifn,jfn)
            cnt(c) = cnt(c) + 1
        end do
        end do

        if (allocated(mt%cstart)) deallocate(mt%cstart)
        if (allocated(mt%flist))  deallocate(mt%flist)
        allocate(mt%cstart(ncrs+1), mt%flist(mt%fin%nx*mt%fin%ny))
        mt%cstart(1) = 1
        do c = 1, ncrs
            mt%cstart(c+1) = mt%cstart(c) + cnt(c)
        end do
        pos = mt%cstart(1:ncrs)
        do jfn = 1, mt%fin%ny
        do ifn = 1, mt%fin%nx
            c  = (mt%pjc(ifn,jfn)-1)*nxc + mt%pic(ifn,jfn)
            fi = (jfn-1)*mt%fin%nx + ifn
            mt%flist(pos(c)) = fi
            pos(c) = pos(c) + 1
        end do
        end do

        mt%nsub_max = max(1, maxval(cnt))
        deallocate(cnt, pos)

        return
    end subroutine build_parent_map

    pure function max_cell_count(mt) result(n)
        implicit none
        type(multitopo_class), intent(IN) :: mt
        integer :: n, c
        n = 0
        do c = 1, size(mt%cstart)-1
            n = max(n, mt%cstart(c+1)-mt%cstart(c))
        end do
        return
    end function max_cell_count

    function grid_spacing(xc) result(dx)
        implicit none
        real(wp), intent(IN) :: xc(:)
        real(wp) :: dx
        if (size(xc) .ge. 2) then
            dx = xc(2) - xc(1)
        else
            dx = 1.0_wp
        end if
        return
    end function grid_spacing

    subroutine level_set(lv, xc, yc, z_bed, H_ice)
        implicit none
        type(multitopo_level_class), intent(INOUT) :: lv
        real(wp), intent(IN) :: xc(:), yc(:), z_bed(:,:), H_ice(:,:)
        call level_dealloc(lv)
        lv%nx = size(xc); lv%ny = size(yc)
        lv%dx = grid_spacing(xc); lv%dy = grid_spacing(yc)
        allocate(lv%xc(lv%nx), lv%yc(lv%ny))
        allocate(lv%z_bed(lv%nx,lv%ny), lv%H_ice(lv%nx,lv%ny), lv%z_srf(lv%nx,lv%ny))
        allocate(lv%z_bed_ref(lv%nx,lv%ny), lv%H_ice_ref(lv%nx,lv%ny))
        lv%xc = xc; lv%yc = yc
        lv%z_bed_ref = z_bed; lv%H_ice_ref = H_ice
        lv%z_bed = z_bed; lv%H_ice = H_ice; lv%z_srf = z_bed + H_ice
        return
    end subroutine level_set

    subroutine level_copy(dst, src)
        ! Copy the current + reference state of src into dst (used for no-op).
        implicit none
        type(multitopo_level_class), intent(INOUT) :: dst
        type(multitopo_level_class), intent(IN)    :: src
        call level_dealloc(dst)
        dst%nx = src%nx; dst%ny = src%ny; dst%dx = src%dx; dst%dy = src%dy
        allocate(dst%xc(dst%nx), dst%yc(dst%ny))
        allocate(dst%z_bed(dst%nx,dst%ny), dst%H_ice(dst%nx,dst%ny), dst%z_srf(dst%nx,dst%ny))
        allocate(dst%z_bed_ref(dst%nx,dst%ny), dst%H_ice_ref(dst%nx,dst%ny))
        dst%xc = src%xc; dst%yc = src%yc
        dst%z_bed = src%z_bed; dst%H_ice = src%H_ice; dst%z_srf = src%z_srf
        dst%z_bed_ref = src%z_bed_ref; dst%H_ice_ref = src%H_ice_ref
        return
    end subroutine level_copy

    subroutine level_dealloc(lv)
        implicit none
        type(multitopo_level_class), intent(INOUT) :: lv
        if (allocated(lv%xc))        deallocate(lv%xc)
        if (allocated(lv%yc))        deallocate(lv%yc)
        if (allocated(lv%z_bed))     deallocate(lv%z_bed)
        if (allocated(lv%H_ice))     deallocate(lv%H_ice)
        if (allocated(lv%z_srf))     deallocate(lv%z_srf)
        if (allocated(lv%z_bed_ref)) deallocate(lv%z_bed_ref)
        if (allocated(lv%H_ice_ref)) deallocate(lv%H_ice_ref)
        lv%nx = 0; lv%ny = 0
        return
    end subroutine level_dealloc

end module multitopo
