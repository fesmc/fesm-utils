program test_multitopo
    ! Acceptance tests for the multitopo module: no-op, mass conservation,
    ! fill-level limiting cases (smooth thick ice / rough exposed bedrock),
    ! bedrock-anomaly mapping, and the sub-grid query.

    use precision, only : wp, dp
    use multitopo

    implicit none

    integer :: nfail

    nfail = 0

    write(*,*) "==================================================="
    write(*,*) " multitopo test suite"
    write(*,*) "==================================================="

    call test_noop()
    call test_conservation()
    call test_fill_limits()
    call test_bedrock_anomaly()
    call test_subgrid_query()

    write(*,*) "---------------------------------------------------"
    if (nfail .eq. 0) then
        write(*,*) " ALL CHECKS PASSED"
    else
        write(*,"(a,i0,a)") "  ", nfail, " CHECK(S) FAILED"
        error stop 1
    end if
    write(*,*) "==================================================="

contains

    ! --- helpers --------------------------------------------------------------

    subroutine check(name, ok)
        character(len=*), intent(IN) :: name
        logical,          intent(IN) :: ok
        write(*,"(a,a4,2x,a)") "  ", merge("ok  ","FAIL",ok), name
        if (.not. ok) nfail = nfail + 1
    end subroutine check

    ! Build a coarse (2x2, dx=4) grid and a fine (4x4, dx=2) grid sharing the
    ! same domain [0,8)x[0,8). Coarse centers {2,6}, fine centers {1,3,5,7} ->
    ! exactly 2x2 fine cells per coarse cell.
    subroutine make_grids(xcc, ycc, xcf, ycf)
        real(wp), allocatable, intent(OUT) :: xcc(:), ycc(:), xcf(:), ycf(:)
        xcc = [2.0_wp, 6.0_wp];               ycc = xcc
        xcf = [1.0_wp, 3.0_wp, 5.0_wp, 7.0_wp]; ycf = xcf
    end subroutine make_grids

    ! --- tests ----------------------------------------------------------------

    subroutine test_noop()
        type(multitopo_class) :: mt
        real(wp) :: xc(3), yc(3)
        real(wp) :: bed(3,3), ice(3,3), bed2(3,3), ice2(3,3)
        real(wp) :: zs(9), fc(9)
        integer  :: i, j, n
        logical  :: ok

        xc = [1.0_wp,2.0_wp,3.0_wp]; yc = xc
        do j=1,3; do i=1,3; bed(i,j)=real(10*i+j,wp); ice(i,j)=real(i,wp); end do; end do

        ! Identical grids -> no-op.
        call multitopo_init(mt, xc,yc,bed,ice, xc,yc,bed,ice)
        bed2 = bed + 5.0_wp; ice2 = ice + 1.0_wp
        call multitopo_update(mt, MT_COARSE, bed2, ice2)
        call multitopo_cell_subgrid(mt, 2,2, n, zs, fc)
        ok = (n .eq. 1) .and. (abs(zs(1) - (bed2(2,2)+ice2(2,2))) .lt. 1e-5_wp)
        call check("no-op: fine==coarse, single sub-cell", ok)
        call multitopo_end(mt)
    end subroutine test_noop

    subroutine test_conservation()
        type(multitopo_class) :: mt
        real(wp), allocatable :: xcc(:),ycc(:),xcf(:),ycf(:)
        real(wp) :: bedc(2,2), icec(2,2), bedf(4,4), icef(4,4)
        real(wp) :: zs(16), fc(16), meanH, bedmean
        integer  :: i,j,n
        logical  :: ok

        call make_grids(xcc,ycc,xcf,ycf)
        ! Rough fine bedrock; coarse ref bed = fine-mean; ice zero at ref.
        do j=1,4; do i=1,4; bedf(i,j)=real(100*mod(i+j,3),wp); end do; end do
        do j=1,2; do i=1,2; bedc(i,j)= sum(bedf(2*i-1:2*i,2*j-1:2*j))/4.0_wp; end do; end do
        icec = 0.0_wp; icef = 0.0_wp

        call multitopo_init(mt, xcc,ycc,bedc,icec, xcf,ycf,bedf,icef)

        ! Give the coarse cells 150 m mean ice and downscale. With no bedrock
        ! anomaly (bedc unchanged), the fine bed = reference, so per coarse cell
        ! mean(z_srf) - mean(bed) = mean(H_ice) must recover 150 m exactly.
        icec = 150.0_wp
        call multitopo_update(mt, MT_COARSE, bedc, icec)

        ok = .TRUE.
        do j=1,2; do i=1,2
            call multitopo_cell_subgrid(mt, i,j, n, zs, fc)
            if (n .ne. 4) ok = .FALSE.
            bedmean = sum(bedf(2*i-1:2*i,2*j-1:2*j))/4.0_wp
            meanH   = sum(zs(1:n))/real(n,wp) - bedmean
            if (abs(meanH - 150.0_wp) .gt. 1e-2_wp) ok = .FALSE.
        end do; end do
        call check("downscale conserves each coarse cell's mean ice (150 m)", ok)

        call multitopo_end(mt)
    end subroutine test_conservation

    subroutine test_fill_limits()
        type(multitopo_class) :: mt
        real(wp), allocatable :: xcc(:),ycc(:),xcf(:),ycf(:)
        real(wp) :: bedc(2,2), icec(2,2), bedf(4,4), icef(4,4)
        real(wp) :: zs(16), fc(16), spread, meanH
        integer  :: i,j,n
        logical  :: ok

        call make_grids(xcc,ycc,xcf,ycf)
        do j=1,4; do i=1,4; bedf(i,j)=real(200*mod(i+2*j,4),wp); end do; end do  ! 0..600 rough
        do j=1,2; do i=1,2; bedc(i,j)= sum(bedf(2*i-1:2*i,2*j-1:2*j))/4.0_wp; end do; end do
        icec=0.0_wp; icef=0.0_wp
        call multitopo_init(mt, xcc,ycc,bedc,icec, xcf,ycf,bedf,icef)

        ! Hbar = 0 -> surface = bedrock (rough): sub-grid spread = bed spread.
        icec = 0.0_wp
        call multitopo_update(mt, MT_COARSE, bedc, icec)
        call multitopo_cell_subgrid(mt, 1,1, n, zs, fc)
        spread = maxval(zs(1:n)) - minval(zs(1:n))
        ok = spread .gt. 100.0_wp          ! rough (beds span hundreds of m)
        call check("Hbar=0 -> surface follows rough bedrock", ok)

        ! Thick ice (Hbar = 5000) -> flat surface: sub-grid spread ~ 0.
        icec = 5000.0_wp
        call multitopo_update(mt, MT_COARSE, bedc, icec)
        call multitopo_cell_subgrid(mt, 1,1, n, zs, fc)
        spread = maxval(zs(1:n)) - minval(zs(1:n))
        ok = spread .lt. 1.0_wp            ! smooth
        call check("thick ice -> flat (smooth) surface", ok)

        ! And that thick-ice fill conserved the mean thickness (=5000).
        meanH = sum(zs(1:n))/real(n,wp) - sum([(bedf(1:2,1:2))])/4.0_wp
        ok = abs(meanH - 5000.0_wp) .lt. 1.0_wp
        call check("fill-level conserves mean ice thickness", ok)

        call multitopo_end(mt)
    end subroutine test_fill_limits

    subroutine test_bedrock_anomaly()
        type(multitopo_class) :: mt
        real(wp), allocatable :: xcc(:),ycc(:),xcf(:),ycf(:)
        real(wp) :: bedc(2,2), icec(2,2), bedf(4,4), icef(4,4)
        real(wp) :: zs(16), fc(16), zs0(16), delta
        integer  :: i,j,n,n0
        logical  :: ok

        call make_grids(xcc,ycc,xcf,ycf)
        do j=1,4; do i=1,4; bedf(i,j)=real(50*(i+j),wp); end do; end do
        do j=1,2; do i=1,2; bedc(i,j)= sum(bedf(2*i-1:2*i,2*j-1:2*j))/4.0_wp; end do; end do
        icec=0.0_wp; icef=0.0_wp
        call multitopo_init(mt, xcc,ycc,bedc,icec, xcf,ycf,bedf,icef)

        call multitopo_update(mt, MT_COARSE, bedc, icec)
        call multitopo_cell_subgrid(mt, 1,1, n0, zs0, fc)

        ! Uniform 300 m bedrock uplift on the coarse grid: fine surface should
        ! rise by exactly 300 m everywhere, preserving the sub-grid spread.
        delta = 300.0_wp
        call multitopo_update(mt, MT_COARSE, bedc + delta, icec)
        call multitopo_cell_subgrid(mt, 1,1, n, zs, fc)
        ok = (n .eq. n0)
        do i = 1, n
            if (abs((zs(i) - zs0(i)) - delta) .gt. 1e-3_wp) ok = .FALSE.
        end do
        call check("coarse bedrock anomaly maps 1:1 onto fine surface", ok)

        call multitopo_end(mt)
    end subroutine test_bedrock_anomaly

    subroutine test_subgrid_query()
        type(multitopo_class) :: mt
        real(wp), allocatable :: xcc(:),ycc(:),xcf(:),ycf(:)
        real(wp) :: bedc(2,2), icec(2,2), bedf(4,4), icef(4,4)
        real(wp) :: zs(16), fc(16)
        integer  :: i,j,n
        logical  :: ok

        call make_grids(xcc,ycc,xcf,ycf)
        do j=1,4; do i=1,4; bedf(i,j)=real(10*i+j,wp); end do; end do
        do j=1,2; do i=1,2; bedc(i,j)= sum(bedf(2*i-1:2*i,2*j-1:2*j))/4.0_wp; end do; end do
        icec=0.0_wp; icef=0.0_wp
        call multitopo_init(mt, xcc,ycc,bedc,icec, xcf,ycf,bedf,icef)
        call multitopo_update(mt, MT_COARSE, bedc, icec)

        ! Coarse cell (1,1) owns fine cells (1..2,1..2): surfaces = their beds.
        call multitopo_cell_subgrid(mt, 1,1, n, zs, fc)
        ok = (n .eq. 4)
        ok = ok .and. any(abs(zs(1:n) - bedf(1,1)) .lt. 1e-4_wp)
        ok = ok .and. any(abs(zs(1:n) - bedf(2,2)) .lt. 1e-4_wp)
        ok = ok .and. all(fc(1:n) .lt. 0.5_wp)      ! no ice
        call check("cell_subgrid returns the right 4 sub-cells", ok)

        call multitopo_end(mt)
    end subroutine test_subgrid_query

end program test_multitopo
