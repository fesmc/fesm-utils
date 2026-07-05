program bench_alloc
    ! Isolate the cost of per-target allocation in the map apply loop.
    !
    ! A/B on IDENTICAL weight data (a cached conservative global->NH-16KM map),
    ! same arithmetic (missing-aware area-weighted mean), differing only in how
    ! the per-target scratch is handled:
    !   A = allocate(vloc(nlink)),allocate(wloc(nlink)) inside the target loop
    !   B = index the weight arrays directly, no per-target scratch at all
    ! B is the strategy weight_map_apply uses for the MAP_WEIGHT+mean fast path;
    ! A is the pattern that made the refactored apply ~10x slower than the legacy
    ! SCRIP apply (which used a single hoisted scratch). This test guards against
    ! that regression returning.

    use coords
    implicit none

    integer, parameter :: NAPPLY = 50
    real(dp), parameter :: d2r = acos(-1.0_dp)/180.0_dp

    type(grid_class) :: g_glob, g_nh
    type(map_class)  :: map
    real(dp), allocatable :: fs(:,:), v1(:), v2(:)
    integer(8) :: c0, c1, cr
    integer    :: r
    real(dp)   :: tA, tB

    call grid_init(g_glob, name="global-1deg", mtype="latlon", units="degrees", &
                   x0=0.0_dp, dx=1.0_dp, nx=360, y0=-90.0_dp, dy=1.0_dp, ny=181)
    call grid_init(g_nh, name="NH-16KM", mtype="stereographic", units="kilometers", &
                   x0=-4000.0_dp, dx=16.0_dp, nx=500, &
                   y0=-4000.0_dp, dy=16.0_dp, ny=500, &
                   lambda=0.0_dp, phi=90.0_dp, alpha=71.0_dp)

    call map_init_conservative(map, g_glob, g_nh)

    allocate(fs(g_glob%G%nx, g_glob%G%ny))
    fs = 2.0_dp + cos(g_glob%lat*d2r)*cos(g_glob%lon*d2r)
    allocate(v1(g_glob%npts), v2(g_nh%npts))
    v1 = reshape(fs, [g_glob%npts])

    write(*,"(a)") ""
    write(*,"(a,i0,a,i0,a,i0)") "n_dst = ", map%wm%n_dst, "   n_links = ", &
        map%wm%n_links, "   n_apply = ", NAPPLY
    write(*,"(a)") ""

    call system_clock(c0, cr)
    do r = 1, NAPPLY
        call apply_A(map%wm, v1, v2)
    end do
    call system_clock(c1)
    tA = real(c1-c0,dp)/real(cr,dp)/real(NAPPLY,dp)

    call system_clock(c0, cr)
    do r = 1, NAPPLY
        call apply_B(map%wm, v1, v2)
    end do
    call system_clock(c1)
    tB = real(c1-c0,dp)/real(cr,dp)/real(NAPPLY,dp)

    write(*,"(a,f9.3,a)") "A  per-target alloc  : ", tA*1.0e3_dp, " ms/apply"
    write(*,"(a,f9.3,a)") "B  no per-target alloc: ", tB*1.0e3_dp, " ms/apply"
    write(*,"(a,f9.2,a)") "speedup B/A          : ", tA/tB, " x"
    write(*,"(a)") ""

contains

    subroutine apply_A(wm, var1, var2)
        use weight_map, only: weight_map_t
        type(weight_map_t), intent(in)  :: wm
        real(dp),           intent(in)  :: var1(:)
        real(dp),           intent(out) :: var2(:)
        integer :: k, j, j1, j2, nlink
        real(dp), allocatable :: vloc(:), wloc(:)
        real(dp) :: wsum, vsum
        var2 = -9999.0_dp
        do k = 1, wm%n_dst
            j1 = wm%dst_off(k); j2 = wm%dst_off(k+1) - 1
            nlink = j2 - j1 + 1
            if (nlink < 1) cycle
            allocate(vloc(nlink), wloc(nlink))
            do j = 1, nlink
                vloc(j) = var1(wm%src(j1+j-1))
            end do
            wloc = wm%w(j1:j2)
            wsum = 0.0_dp; vsum = 0.0_dp
            do j = 1, nlink
                if (wloc(j) > 0.0_dp) then
                    wsum = wsum + wloc(j); vsum = vsum + wloc(j)*vloc(j)
                end if
            end do
            if (wsum > 0.0_dp) var2(k) = vsum/wsum
            deallocate(vloc, wloc)
        end do
    end subroutine apply_A

    subroutine apply_B(wm, var1, var2)
        use weight_map, only: weight_map_t
        type(weight_map_t), intent(in)  :: wm
        real(dp),           intent(in)  :: var1(:)
        real(dp),           intent(out) :: var2(:)
        integer :: k, j, j1, j2
        real(dp) :: wsum, vsum, w
        var2 = -9999.0_dp
        do k = 1, wm%n_dst
            j1 = wm%dst_off(k); j2 = wm%dst_off(k+1) - 1
            if (j2 < j1) cycle
            wsum = 0.0_dp; vsum = 0.0_dp
            do j = j1, j2
                w = wm%w(j)
                if (w > 0.0_dp) then
                    wsum = wsum + w; vsum = vsum + w*var1(wm%src(j))
                end if
            end do
            if (wsum > 0.0_dp) var2(k) = vsum/wsum
        end do
    end subroutine apply_B

end program bench_alloc
