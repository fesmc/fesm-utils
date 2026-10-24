program test_coupler
    ! Exercise the coupler layer over the coords mapping core:
    !   * bilin downscale (coarse -> fine) is exact for a linear field,
    !   * the map cache builds each (src,dst,method) edge once and reuses it,
    !   * method defaults to "con", and is part of the cache key,
    !   * remap sizes the destination array on demand (2D and 3D),
    !   * coupler_prime pre-builds without adding duplicate maps,
    !   * the slot arrays grow past 16 grids / 64 maps (no fixed capacity),
    !   * remap with src == dst copies (identity shortcut, no map),
    !   * mask2 / mask_pack / reset are passed through to map_field,
    !   * coupler_take_grid moves a grid in; a from_disk name re-registered
    !     with an agreeing grid is a no-op,
    !   * method="auto" resolves to bilin when refining and con when coarsening.

    use coords
    use coupler

    implicit none

    type(coupler_class)   :: cpl
    type(grid_class)      :: gs, gt
    real(dp), allocatable :: fs(:,:), ft(:,:), fback(:,:)
    real(dp), allocatable :: fs3(:,:,:), ft3(:,:,:)
    integer  :: i, j, k, fails
    real(dp) :: emax
    character(len=*), parameter :: fldr = "maps_coupler_test"

    fails = 0
    call execute_command_line("rm -rf "//fldr)

    ! Coarse (100 km) and fine (50 km) grids over the same [0,1000] km extent.
    call grid_init(gs, name="coupler-coarse", mtype="cartesian", units="kilometers", &
                   x0=0.0_dp, dx=100.0_dp, nx=11, y0=0.0_dp, dy=100.0_dp, ny=11)
    call grid_init(gt, name="coupler-fine",   mtype="cartesian", units="kilometers", &
                   x0=0.0_dp, dx=50.0_dp,  nx=21, y0=0.0_dp, dy=50.0_dp,  ny=21)

    call coupler_init(cpl, map_fldr=fldr)
    call coupler_add_grid(cpl, "coupler-coarse", gs)
    call coupler_add_grid(cpl, "coupler-fine",   gt)

    ! Linear source field on the coarse grid.
    allocate(fs(11,11))
    do j = 1, 11
        do i = 1, 11
            fs(i,j) = flin(gs%x(i,j), gs%y(i,j))
        end do
    end do

    ! --- Test 1: bilin downscale is exact for a linear field (interior) ------
    call remap(cpl, fs, "coupler-coarse", ft, "coupler-fine", method="bilin")
    if (.not. allocated(ft)) then
        write(*,*) "FAIL: remap did not allocate the destination"; fails = fails + 1
    else if (size(ft,1) /= 21 .or. size(ft,2) /= 21) then
        write(*,*) "FAIL: destination sized ", shape(ft), " expected 21 x 21"; fails = fails + 1
    else
        emax = 0.0_dp
        do j = 2, 20
            do i = 2, 20
                emax = max(emax, abs(ft(i,j) - flin(gt%x(i,j), gt%y(i,j))))
            end do
        end do
        write(*,*) "bilin downscale interior max err =", emax
        if (emax > 1.0e-6_dp) then
            write(*,*) "FAIL: bilin downscale not exact for linear field"; fails = fails + 1
        end if
    end if

    ! --- Test 2: the map cache builds one edge and reuses it -----------------
    if (cpl%nmaps /= 1) then
        write(*,*) "FAIL: expected 1 cached map, got", cpl%nmaps; fails = fails + 1
    end if
    call remap(cpl, fs, "coupler-coarse", ft, "coupler-fine", method="bilin")
    if (cpl%nmaps /= 1) then
        write(*,*) "FAIL: repeated remap rebuilt the map (nmaps=", cpl%nmaps, ")"; fails = fails + 1
    end if

    ! --- Test 3: method defaults to "con", and is a distinct cache edge ------
    call remap(cpl, ft, "coupler-fine", fback, "coupler-coarse")   ! no method -> con
    if (cpl%nmaps /= 2) then
        write(*,*) "FAIL: default-method remap did not add a con edge (nmaps=", cpl%nmaps, ")"; fails = fails + 1
    end if
    if (trim(cpl%maps(2)%method) /= "con") then
        write(*,*) "FAIL: default method was '", trim(cpl%maps(2)%method), "', expected 'con'"; fails = fails + 1
    end if
    if (size(fback,1) /= 11 .or. size(fback,2) /= 11) then
        write(*,*) "FAIL: aggregated field sized ", shape(fback), " expected 11 x 11"; fails = fails + 1
    else
        ! Conservative aggregation of a linear field ~ value at coarse centres (interior).
        emax = 0.0_dp
        do j = 2, 10
            do i = 2, 10
                emax = max(emax, abs(fback(i,j) - flin(gs%x(i,j), gs%y(i,j))))
            end do
        end do
        write(*,*) "con aggregate interior max err  =", emax
        if (emax > 1.0e-3_dp) then
            write(*,*) "FAIL: con aggregation of linear field off by", emax; fails = fails + 1
        end if
    end if

    ! --- Test 4: coupler_prime does not duplicate an existing edge -----------
    call coupler_prime(cpl, "coupler-coarse", "coupler-fine", "bilin")
    if (cpl%nmaps /= 2) then
        write(*,*) "FAIL: prime duplicated an existing edge (nmaps=", cpl%nmaps, ")"; fails = fails + 1
    end if

    ! --- Test 5: 3D remap matches per-level 2D remap ------------------------
    allocate(fs3(11,11,3))
    do k = 1, 3
        fs3(:,:,k) = fs + real(k, dp)          ! three offset linear levels
    end do
    call remap(cpl, fs3, "coupler-coarse", ft3, "coupler-fine", method="bilin")
    if (size(ft3,1) /= 21 .or. size(ft3,2) /= 21 .or. size(ft3,3) /= 3) then
        write(*,*) "FAIL: 3D destination sized ", shape(ft3), " expected 21 x 21 x 3"; fails = fails + 1
    else
        emax = 0.0_dp
        do k = 1, 3
            call remap(cpl, fs3(:,:,k), "coupler-coarse", ft, "coupler-fine", method="bilin")
            emax = max(emax, maxval(abs(ft3(:,:,k) - ft)))
        end do
        write(*,*) "3D vs per-level 2D max diff      =", emax
        if (emax > 1.0e-12_dp) then
            write(*,*) "FAIL: 3D remap disagrees with per-level 2D remap"; fails = fails + 1
        end if
    end if

    ! --- Test 6: disk-driven grids resolved from grid_<name>.txt (no register) --
    block
        type(coupler_class)   :: cpl2
        type(grid_class)      :: gsd, gtd
        real(wp), allocatable :: fsd(:,:), ftd(:,:)
        character(len=*), parameter :: dfldr = "maps_coupler_disk_test"
        integer  :: i2, j2
        real(wp) :: emax2

        call execute_command_line("rm -rf "//dfldr//" ; mkdir -p "//dfldr)

        ! ANT-like polar_stereographic grids so the description carries a
        ! grid_mapping_name (the CF key the reader falls back on).
        call grid_init(gsd, name="disk-src", mtype="polar_stereographic", units="kilometers", &
                       lon180=.false., x0=-500.0_dp, dx=100.0_dp, nx=11, &
                       y0=-500.0_dp, dy=100.0_dp, ny=11, lambda=0.0_dp, phi=-71.0_dp)
        call grid_init(gtd, name="disk-tgt", mtype="polar_stereographic", units="kilometers", &
                       lon180=.false., x0=-500.0_dp, dx=50.0_dp,  nx=21, &
                       y0=-500.0_dp, dy=50.0_dp,  ny=21, lambda=0.0_dp, phi=-71.0_dp)

        ! Write descriptions, then strip the # header to emulate cdo-native files.
        call grid_cdo_write_desc_short(gsd, dfldr)
        call grid_cdo_write_desc_short(gtd, dfldr)
        call execute_command_line("for f in "//dfldr//"/grid_disk-*.txt; do " // &
                                  "grep -v '^#' $f > $f.tmp && mv $f.tmp $f; done")

        allocate(fsd(11,11))
        do j2 = 1, 11
            do i2 = 1, 11
                fsd(i2,j2) = flin(gsd%x(i2,j2), gsd%y(i2,j2))
            end do
        end do

        call coupler_init(cpl2, map_fldr=dfldr)      ! no coupler_add_grid
        call remap(cpl2, fsd, "disk-src", ftd, "disk-tgt", method="bilin")

        if (.not. allocated(ftd) .or. size(ftd,1) /= 21 .or. size(ftd,2) /= 21) then
            write(*,*) "FAIL: disk-driven remap wrong shape"; fails = fails + 1
        else
            emax2 = 0.0_dp
            do j2 = 2, 20
                do i2 = 2, 20
                    emax2 = max(emax2, abs(ftd(i2,j2) - flin(gtd%x(i2,j2), gtd%y(i2,j2))))
                end do
            end do
            write(*,*) "disk-driven bilin interior err  =", emax2
            if (emax2 > 1.0e-6_dp) then
                write(*,*) "FAIL: disk-driven remap inaccurate"; fails = fails + 1
            end if
        end if
        call execute_command_line("rm -rf "//dfldr)
    end block

    ! --- Test 7: growth past 16 grids / 64 maps ---------------------------------
    block
        type(coupler_class) :: cg
        type(grid_class)    :: gg
        real(dp), allocatable :: a(:,:), b(:,:), b0(:,:)
        character(len=32) :: nm(18)
        integer :: n, m, cnt
        character(len=*), parameter :: gfldr = "maps_coupler_grow_test"

        call execute_command_line("rm -rf "//gfldr)
        call coupler_init(cg, map_fldr=gfldr)
        do n = 1, 18
            write(nm(n), "(a,i0)") "grow-", n
            call grid_init(gg, name=trim(nm(n)), mtype="cartesian", units="kilometers", &
                           x0=0.0_dp, dx=100.0_dp/real(n,dp), nx=n+3, y0=0.0_dp, dy=100.0_dp/real(n,dp), ny=4)
            call coupler_add_grid(cg, trim(nm(n)), gg)
        end do
        allocate(a(4,4)); a = 1.0_dp
        a(2,3) = 5.0_dp
        call remap(cg, a, trim(nm(1)), b0, trim(nm(2)), method="nn")
        cnt = 0
        outer: do n = 1, 18
            do m = 1, 18
                if (m == n) cycle
                if (allocated(a)) deallocate(a)
                allocate(a(n+3,4)); a = real(n,dp)
                call remap(cg, a, trim(nm(n)), b, trim(nm(m)), method="nn")
                cnt = cnt + 1
                if (cnt >= 66) exit outer
            end do
        end do outer
        if (cg%ngrids /= 18 .or. cg%nmaps /= 66) then
            write(*,*) "FAIL: growth: ngrids, nmaps =", cg%ngrids, cg%nmaps, " expected 18, 66"
            fails = fails + 1
        end if
        ! the first map is still valid after both slot arrays grew
        deallocate(a); allocate(a(4,4)); a = 1.0_dp; a(2,3) = 5.0_dp
        call remap(cg, a, trim(nm(1)), b, trim(nm(2)), method="nn")
        if (any(b /= b0) .or. cg%nmaps /= 66) then
            write(*,*) "FAIL: growth: first map changed or rebuilt"; fails = fails + 1
        else
            write(*,*) "growth: 18 grids, 66 maps ok"
        end if
        call execute_command_line("rm -rf "//gfldr)
    end block

    ! --- Test 8: identity shortcut -----------------------------------------------
    block
        real(dp), allocatable :: fid(:,:), fk(:,:)
        logical,  allocatable :: m2(:,:), mp(:,:)
        real(dp) :: f2id(11,11)
        integer :: nm0

        nm0 = cpl%nmaps
        fs(3,4) = mv_dp
        allocate(m2(11,11), mp(11,11))
        call remap(cpl, fs, "coupler-coarse", fid, "coupler-coarse", method="con", mask2=m2, frac2=f2id)
        if (any(fid /= fs) .or. any(m2 .neqv. (fs /= mv_dp)) .or. cpl%nmaps /= nm0 .or. &
            any(f2id /= merge(1.0_dp, 0.0_dp, m2))) then
            write(*,*) "FAIL: identity remap is not a copy (or built a map)"; fails = fails + 1
        end if
        ! unregistered name: no registration needed
        call remap(cpl, fs, "not-registered", fid, "not-registered")
        if (any(fid /= fs) .or. cpl%ngrids /= 2) then
            write(*,*) "FAIL: identity remap on an unregistered name"; fails = fails + 1
        end if
        ! mask_pack + reset=.false. keep the destination outside the pack
        mp = .true.; mp(1:5,:) = .false.
        allocate(fk(11,11)); fk = -1.0_dp
        call remap(cpl, fs, "coupler-coarse", fk, "coupler-coarse", mask_pack=mp, reset=.false.)
        if (any(fk(1:5,:) /= -1.0_dp) .or. any(fk(6:,:) /= fs(6:,:))) then
            write(*,*) "FAIL: identity remap with mask_pack/reset=.false."; fails = fails + 1
        end if
        call remap(cpl, fs, "coupler-coarse", fid, "coupler-coarse", stat="stdev")
        if (any(fid /= merge(0.0_dp, mv_dp, fs /= mv_dp))) then
            write(*,*) "FAIL: identity remap stat=stdev is not 0"; fails = fails + 1
        end if
        call coupler_prime(cpl, "coupler-fine", "coupler-fine")
        if (cpl%nmaps /= nm0) then
            write(*,*) "FAIL: coupler_prime(src == dst) built a map"; fails = fails + 1
        end if
        fs(3,4) = flin(gs%x(3,4), gs%y(3,4))
        write(*,*) "identity shortcut ok"
    end block

    ! --- Test 9: mask2 / mask_pack pass-through ------------------------------------
    block
        type(map_class) :: mp_direct
        real(dp), allocatable :: fr(:,:), fd(:,:)
        logical,  allocatable :: m2r(:,:), m2d(:,:), mpk(:,:)
        real(dp) :: f2r(21,21), f2d(21,21)
        allocate(m2r(21,21), m2d(21,21), mpk(21,21), fd(21,21))
        mpk = .true.; mpk(:,1:4) = .false.
        fs(5,5) = mv_dp
        call remap(cpl, fs, "coupler-coarse", fr, "coupler-fine", method="con", mask2=m2r, mask_pack=mpk, &
                   frac2=f2r)
        call map_init(mp_direct, gs, gt, method="con", fldr=fldr)
        call map_field(mp_direct, "x", fs, fd, mask2=m2d, mask_pack=mpk, frac2=f2d)
        if (any(fr /= fd) .or. any(m2r .neqv. m2d) .or. count(m2r) == size(m2r) .or. any(f2r /= f2d)) then
            write(*,*) "FAIL: mask2/mask_pack/frac2 pass-through differs from map_field"; fails = fails + 1
        else
            write(*,*) "mask2/frac2 pass-through ok (", count(m2r), " of", size(m2r), " mapped)"
        end if
        ! norm="dst": the coupler supplies the registered (default) model areas
        call remap(cpl, fs, "coupler-coarse", fr, "coupler-fine", method="con", norm="dst")
        call map_field(mp_direct, "x", fs, fd, norm="dst", area1=gs%area*gs%cs%xy_conv**2, &
                       area2=gt%area*gt%cs%xy_conv**2)
        if (any(fr /= fd)) then
            write(*,*) "FAIL: norm pass-through (with registered areas) differs from map_field"
            fails = fails + 1
        end if
        fs(5,5) = flin(gs%x(5,5), gs%y(5,5))
    end block

    ! --- Test 10: coupler_take_grid; from_disk re-registration ---------------------
    block
        type(coupler_class) :: ct
        type(grid_class), allocatable :: gtake
        real(dp), allocatable :: ftk(:,:)
        character(len=*), parameter :: tfldr = "maps_coupler_take_test"

        call execute_command_line("rm -rf "//tfldr//" ; mkdir -p "//tfldr)
        call coupler_init(ct, map_fldr=tfldr)
        allocate(gtake, source=gt)
        call coupler_take_grid(ct, "coupler-fine", gtake)
        if (allocated(gtake)) then
            write(*,*) "FAIL: coupler_take_grid left the grid allocated"; fails = fails + 1
        end if
        ! coarse grid resolved from disk, then registered with the same grid
        call grid_cdo_write_desc_short(gs, tfldr)
        call remap(ct, fs, "coupler-coarse", ftk, "coupler-fine", method="bilin")
        call coupler_add_grid(ct, "coupler-coarse", gs)          ! agrees: no-op
        if (ct%ngrids /= 2 .or. .not. ct%grids(2)%from_disk) then
            write(*,*) "FAIL: take/disk registry state"; fails = fails + 1
        end if
        call remap(cpl, fs, "coupler-coarse", ft, "coupler-fine", method="bilin")
        if (any(ftk /= ft)) then
            write(*,*) "FAIL: remap through a taken / disk grid differs"; fails = fails + 1
        else
            write(*,*) "take_grid + from_disk ok"
        end if
        call execute_command_line("rm -rf "//tfldr)
    end block

    ! --- Test 11: method="auto" ------------------------------------------------------
    block
        type(coupler_class) :: ca
        real(dp), allocatable :: fa(:,:), fb(:,:)
        character(len=*), parameter :: afldr = "maps_coupler_auto_test"

        call execute_command_line("rm -rf "//afldr)
        call coupler_init(ca, map_fldr=afldr)
        call coupler_add_grid(ca, "coupler-coarse", gs)
        call coupler_add_grid(ca, "coupler-fine",   gt)
        call remap(ca, fs, "coupler-coarse", fa, "coupler-fine", method="auto")     ! refining
        call remap(ca, fa, "coupler-fine", fb, "coupler-coarse", method="auto")     ! coarsening
        if (ca%nmaps /= 2 .or. trim(ca%maps(1)%method) /= "bilin" .or. trim(ca%maps(2)%method) /= "con") then
            write(*,*) "FAIL: method=auto resolved to ", trim(ca%maps(1)%method), " / ", &
                       trim(ca%maps(2)%method), " (expected bilin / con)"
            fails = fails + 1
        else
            write(*,*) "method=auto ok (refine: bilin, coarsen: con)"
        end if
        call execute_command_line("rm -rf "//afldr)
    end block

    call execute_command_line("rm -rf "//fldr)

    if (fails > 0) stop 1
    write(*,*) "PASS: test_coupler"

contains

    pure function flin(x, y) result(f)
        real(dp), intent(in) :: x, y
        real(dp) :: f
        f = 2.0_dp + 0.003_dp*x + 0.005_dp*y
    end function flin

end program test_coupler
