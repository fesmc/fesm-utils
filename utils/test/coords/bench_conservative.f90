program bench_conservative
    ! Performance + accuracy of coords conservative remapping vs `cdo gencon`.
    !
    ! For each regime and size:
    !   - coords: time in-process map_init_conservative (analytic polygon clip).
    !   - cdo   : time the full subprocess cost you would actually pay to get the
    !             same weights -- serialize grids (grid desc + source netcdf),
    !             run `cdo gencon`, read the SCRIP file back -- then bridge it to a
    !             weight_map via map_scrip_to_weight_map.
    !   - accuracy: map a constant field (conservation: result should be 1 where
    !             covered), the area integral where coverage is full, and the
    !             max |coords - cdo| on a smooth field where both apply.
    !
    ! Regimes: cartesian->cartesian (planar, no cdo equivalent), latlon->stereo-
    ! graphic (planar cross-system), latlon->latlon (spherical, NSUB great-circle).
    ! A short IDW section compares coords kdtree shepard vs `cdo gendis`.

    use coords
    use mapping_scrip, only : map_scrip_class, map_scrip_init_from_griddesc, &
                              map_scrip_to_weight_map, map_scrip_end
    use grid_to_cdo,   only : grid_cdo_write_desc_explicit_latlon, &
                              grid_cdo_write_desc_explicit_proj

    implicit none

    integer, parameter :: NREP = 3
    character(len=*), parameter :: fldr = "bench_tmp"
    real(dp), parameter :: d2r = acos(-1.0_dp)/180.0_dp

    call execute_command_line("mkdir -p "//fldr)

    write(*,"(a)") ""
    write(*,"(a)") "# coords conservative remapping: performance vs cdo gencon"
    write(*,"(a)") ""
    write(*,"(a)") "regime                    n_src     n_tgt   t_coords[s]  t_cdo[s]   speedup   cons_co   cons_cdo  max|co-cdo|"
    write(*,"(a)") "-----------------------------------------------------------------------------------------------------------"

    call regime_cartesian( 60)
    call regime_cartesian(120)
    call regime_cartesian(180)

    call regime_latlon_to_stereo(40.0_dp)
    call regime_latlon_to_stereo(20.0_dp)
    call regime_latlon_to_stereo(10.0_dp)

    call regime_latlon_to_latlon(1.5_dp)
    call regime_latlon_to_latlon(1.0_dp)
    call regime_latlon_to_latlon(0.75_dp)

    write(*,"(a)") ""
    write(*,"(a)") "# IDW: coords kdtree shepard vs cdo gendis (latlon -> stereographic)"
    write(*,"(a)") ""
    write(*,"(a)") "regime                    n_src     n_tgt   t_coords[s]  t_cdo[s]   speedup   max|co-cdo|"
    write(*,"(a)") "-------------------------------------------------------------------------------------------"
    call regime_idw(20.0_dp)
    call regime_idw(10.0_dp)

    write(*,"(a)") ""
    write(*,"(a)") "done."

contains

    ! ---- field helpers -------------------------------------------------------

    subroutine smooth_field(grid, f)
        type(grid_class), intent(in)  :: grid
        real(dp), allocatable, intent(out) :: f(:,:)
        allocate(f(grid%G%nx, grid%G%ny))
        if (grid%cs%is_cartesian) then
            f = 1.0_dp + 1.0e-3_dp*(grid%x + grid%y)
        else
            f = 2.0_dp + cos(grid%lat*d2r)*cos(grid%lon*d2r)
        end if
    end subroutine smooth_field

    ! ---- timing of one conservative comparison -------------------------------

    subroutine bench_row(label, tag, gs, gt, use_cdo, full_cover)
        character(len=*), intent(in) :: label        ! display label
        character(len=*), intent(in) :: tag          ! shell/filename-safe name (no spaces/special chars)
        type(grid_class), intent(in) :: gs, gt
        logical,          intent(in) :: use_cdo      ! also run + compare cdo
        logical,          intent(in) :: full_cover   ! target fully tiles source (integral check)

        type(map_class)       :: map_co, map_cdo
        type(map_scrip_class) :: mps
        real(dp), allocatable :: fs(:,:), ft_co(:,:), ft_cdo(:,:), ones_s(:,:), c_co(:,:), c_cdo(:,:)
        logical,  allocatable :: m_co(:,:), m_cdo(:,:)
        integer(8) :: c0, c1, cr
        integer    :: r
        real(dp)   :: t_co, t_cdo, cons_co, cons_cdo, maxdiff, sp
        character(len=512) :: src_nc, xnm, ynm
        character(len=64)  :: snm, tnm

        call smooth_field(gs, fs)
        allocate(ones_s(gs%G%nx,gs%G%ny)); ones_s = 1.0_dp
        allocate(ft_co(gt%G%nx,gt%G%ny), c_co(gt%G%nx,gt%G%ny), m_co(gt%G%nx,gt%G%ny))

        ! --- coords: time weight generation ---
        call system_clock(c0, cr)
        do r = 1, NREP
            call map_init_conservative(map_co, gs, gt)
        end do
        call system_clock(c1)
        t_co = real(c1-c0,dp)/real(cr,dp)/real(NREP,dp)

        ! coords accuracy
        call map_field(map_co, "c", ones_s, c_co, stat="mean", mask2=m_co)
        cons_co = 0.0_dp
        if (any(m_co)) cons_co = maxval(abs(c_co - 1.0_dp), mask=m_co)
        call map_field(map_co, "f", fs, ft_co, stat="mean", mask2=m_co)

        t_cdo    = -1.0_dp
        cons_cdo = -1.0_dp
        maxdiff  = -1.0_dp

        if (use_cdo) then
            snm = trim(tag)//"_src"
            tnm = trim(tag)//"_tgt"
            src_nc = trim(fldr)//"/src_"//trim(tag)//".nc"

            ! --- cdo: time the full serialize + gencon + read-back cost ---
            call system_clock(c0, cr)
            do r = 1, NREP
                call write_descriptions(gs, gt, snm, tnm)
                if (gs%cs%is_projection .or. gs%cs%is_cartesian) then
                    xnm = "xc"; ynm = "yc"
                else
                    xnm = "lon"; ynm = "lat"
                end if
                call grid_write(gs, fnm=trim(src_nc), xnm=trim(xnm), ynm=trim(ynm), create=.true.)
                call map_scrip_init_from_griddesc(mps, trim(snm), trim(tnm), fldr, trim(src_nc), "con", load=.false.)
            end do
            call system_clock(c1)
            t_cdo = real(c1-c0,dp)/real(cr,dp)/real(NREP,dp)

            ! bridge cdo weights -> weight_map and apply
            call map_scrip_to_weight_map(mps, map_cdo%wm)
            allocate(ft_cdo(gt%G%nx,gt%G%ny), c_cdo(gt%G%nx,gt%G%ny), m_cdo(gt%G%nx,gt%G%ny))
            call map_field(map_cdo, "c", ones_s, c_cdo, stat="mean", mask2=m_cdo)
            cons_cdo = 0.0_dp
            if (any(m_cdo)) cons_cdo = maxval(abs(c_cdo - 1.0_dp), mask=m_cdo)
            call map_field(map_cdo, "f", fs, ft_cdo, stat="mean", mask2=m_cdo)

            if (any(m_co .and. m_cdo)) then
                maxdiff = maxval(abs(ft_co - ft_cdo), mask=(m_co .and. m_cdo))
            end if
            call map_scrip_end(mps)
        end if

        if (full_cover) then
            ! integral conservation (coords), reported in the cons_co slot context note
            cons_co = max(cons_co, abs(sum(ft_co*gt%area, mask=m_co) - sum(fs*gs%area)) / abs(sum(fs*gs%area)))
        end if

        sp = -1.0_dp
        if (t_cdo > 0.0_dp .and. t_co > 0.0_dp) sp = t_cdo/t_co

        call print_row(label, gs%npts, gt%npts, t_co, t_cdo, sp, cons_co, cons_cdo, maxdiff)

    end subroutine bench_row

    ! ---- regimes -------------------------------------------------------------

    subroutine regime_cartesian(nsrc)
        integer, intent(in) :: nsrc
        type(grid_class) :: gs, gt
        character(len=64) :: lab, tag
        write(lab,"(a,i0)") "cart->cart n", nsrc
        write(tag,"(a,i0)") "cart_", nsrc
        ! source nsrc x nsrc cells over [0,L]^2; target 2x refinement (full cover)
        call grid_init(gs, name="cs", mtype="cartesian", units="kilometers", &
                       x0=0.5_dp, dx=1.0_dp, nx=nsrc, y0=0.5_dp, dy=1.0_dp, ny=nsrc)
        call grid_init(gt, name="ct", mtype="cartesian", units="kilometers", &
                       x0=0.25_dp, dx=0.5_dp, nx=2*nsrc, y0=0.25_dp, dy=0.5_dp, ny=2*nsrc)
        call bench_row(trim(lab), trim(tag), gs, gt, use_cdo=.false., full_cover=.true.)
    end subroutine regime_cartesian

    subroutine regime_latlon_to_stereo(dxkm)
        real(dp), intent(in) :: dxkm
        type(grid_class) :: gs, gt
        integer :: nx, ny
        character(len=64) :: lab, tag
        write(lab,"(a,i0)") "ll->stereo dx", nint(dxkm)
        write(tag,"(a,i0)") "ll2st_", nint(dxkm)
        ! global 2-degree source
        call grid_init(gs, name="lls", mtype="latlon", units="degrees", &
                       x0=0.0_dp, dx=2.0_dp, nx=180, y0=-90.0_dp, dy=2.0_dp, ny=91)
        ! Greenland oblique stereographic at the requested resolution
        nx = nint(1520.0_dp/dxkm); ny = nint(3020.0_dp/dxkm)
        call grid_init(gt, name="grl", mtype="stereographic", units="kilometers", &
                       dx=dxkm, nx=nx, dy=dxkm, ny=ny, &
                       lambda=-40.0_dp, phi=72.0_dp, alpha=7.5_dp)
        call bench_row(trim(lab), trim(tag), gs, gt, use_cdo=.true., full_cover=.false.)
    end subroutine regime_latlon_to_stereo

    subroutine regime_latlon_to_latlon(dxdeg)
        real(dp), intent(in) :: dxdeg
        type(grid_class) :: gs, gt
        integer :: nx, ny
        character(len=64) :: lab, tag
        write(lab,"(a,i0)") "ll->ll dx0p", nint(dxdeg*100)
        write(tag,"(a,i0)") "ll2ll_", nint(dxdeg*100)
        ! global 2-degree source -> finer global target
        call grid_init(gs, name="gs2", mtype="latlon", units="degrees", &
                       x0=0.0_dp, dx=2.0_dp, nx=180, y0=-90.0_dp, dy=2.0_dp, ny=91)
        nx = nint(360.0_dp/dxdeg); ny = nint(180.0_dp/dxdeg) + 1
        call grid_init(gt, name="gtf", mtype="latlon", units="degrees", &
                       x0=0.0_dp, dx=dxdeg, nx=nx, y0=-90.0_dp, dy=dxdeg, ny=ny)
        call bench_row(trim(lab), trim(tag), gs, gt, use_cdo=.true., full_cover=.false.)
    end subroutine regime_latlon_to_latlon

    subroutine regime_idw(dxkm)
        real(dp), intent(in) :: dxkm
        type(grid_class) :: gs, gt
        type(map_class)  :: map_co, map_cdo
        type(map_scrip_class) :: mps
        real(dp), allocatable :: fs(:,:), ft_co(:,:), ft_cdo(:,:)
        logical,  allocatable :: m_co(:,:), m_cdo(:,:)
        integer(8) :: c0, c1, cr
        integer    :: r, nx, ny
        real(dp)   :: t_co, t_cdo, maxdiff, sp
        character(len=512) :: src_nc
        character(len=64)  :: lab, tag, snm, tnm

        write(lab,"(a,i0)") "idw ll->stereo dx", nint(dxkm)
        write(tag,"(a,i0)") "idw_", nint(dxkm)
        call grid_init(gs, name="ils", mtype="latlon", units="degrees", &
                       x0=0.0_dp, dx=2.0_dp, nx=180, y0=-90.0_dp, dy=2.0_dp, ny=91)
        nx = nint(1520.0_dp/dxkm); ny = nint(3020.0_dp/dxkm)
        call grid_init(gt, name="igr", mtype="stereographic", units="kilometers", &
                       dx=dxkm, nx=nx, dy=dxkm, ny=ny, &
                       lambda=-40.0_dp, phi=72.0_dp, alpha=7.5_dp)

        call smooth_field(gs, fs)
        allocate(ft_co(gt%G%nx,gt%G%ny), m_co(gt%G%nx,gt%G%ny))

        call system_clock(c0, cr)
        do r = 1, NREP
            call map_init(map_co, gs, gt, max_neighbors=10, load=.false.)
        end do
        call system_clock(c1)
        t_co = real(c1-c0,dp)/real(cr,dp)/real(NREP,dp)
        call map_field(map_co, "f", fs, ft_co, method="shepard", mask2=m_co)

        snm = trim(tag)//"_src"; tnm = trim(tag)//"_tgt"
        src_nc = trim(fldr)//"/src_idw.nc"
        call system_clock(c0, cr)
        do r = 1, NREP
            call write_descriptions(gs, gt, snm, tnm)
            call grid_write(gs, fnm=trim(src_nc), xnm="lon", ynm="lat", create=.true.)
            call map_scrip_init_from_griddesc(mps, trim(snm), trim(tnm), fldr, trim(src_nc), "dis", load=.false.)
        end do
        call system_clock(c1)
        t_cdo = real(c1-c0,dp)/real(cr,dp)/real(NREP,dp)

        call map_scrip_to_weight_map(mps, map_cdo%wm)
        allocate(ft_cdo(gt%G%nx,gt%G%ny), m_cdo(gt%G%nx,gt%G%ny))
        call map_field(map_cdo, "f", fs, ft_cdo, stat="mean", mask2=m_cdo)
        maxdiff = -1.0_dp
        if (any(m_co .and. m_cdo)) maxdiff = maxval(abs(ft_co - ft_cdo), mask=(m_co .and. m_cdo))
        call map_scrip_end(mps)

        sp = -1.0_dp
        if (t_cdo > 0.0_dp .and. t_co > 0.0_dp) sp = t_cdo/t_co
        write(*,"(a24,1x,i9,1x,i9,2x,es10.3,1x,es10.3,2x,f7.1,3x,es10.3)") &
              adjustl(lab), gs%npts, gt%npts, t_co, t_cdo, sp, maxdiff
    end subroutine regime_idw

    ! ---- cdo helpers ---------------------------------------------------------

    subroutine write_descriptions(gs, gt, snm, tnm)
        type(grid_class), intent(in) :: gs, gt
        character(len=*), intent(in) :: snm, tnm
        if (gs%cs%is_projection) then
            call grid_cdo_write_desc_explicit_proj(gs%lon, gs%lat, trim(snm), fldr)
        else
            call grid_cdo_write_desc_explicit_latlon(gs%G%x, gs%G%y, trim(snm), fldr, .false.)
        end if
        if (gt%cs%is_projection) then
            call grid_cdo_write_desc_explicit_proj(gt%lon, gt%lat, trim(tnm), fldr)
        else
            call grid_cdo_write_desc_explicit_latlon(gt%G%x, gt%G%y, trim(tnm), fldr, .false.)
        end if
    end subroutine write_descriptions

    subroutine print_row(label, nsrc, ntgt, t_co, t_cdo, sp, cons_co, cons_cdo, maxdiff)
        character(len=*), intent(in) :: label
        integer,  intent(in) :: nsrc, ntgt
        real(dp), intent(in) :: t_co, t_cdo, sp, cons_co, cons_cdo, maxdiff
        if (t_cdo < 0.0_dp) then
            write(*,"(a24,1x,i9,1x,i9,2x,es11.3,2x,a11,1x,a8,2x,es10.2,2x,a10,2x,a11)") &
                  adjustl(label), nsrc, ntgt, t_co, "        --", "      --", cons_co, "        --", "         --"
        else
            write(*,"(a24,1x,i9,1x,i9,2x,es11.3,2x,es11.3,1x,f8.2,2x,es10.2,2x,es10.2,2x,es11.3)") &
                  adjustl(label), nsrc, ntgt, t_co, t_cdo, sp, cons_co, cons_cdo, maxdiff
        end if
    end subroutine print_row

end program bench_conservative
