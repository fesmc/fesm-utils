program test_coupler_bitident
    ! Bit-identity harness for the coupler / map_field apply surface.
    !
    ! A fixed call set (con / bilin / nn; dp / sp / int; stat variants;
    ! reset=.false.; mask_pack; fill_method / filt_method; 3D; grid <-> points)
    ! over grids that exercise every conservative generator path (same plane,
    ! lat-lon -> lat-lon incl. a lon180 seam grid and a descending Gauss grid,
    ! lat-lon -> projected, projected -> lat-lon) is run, and every result is
    ! written to its own stream file in bitident_out/<case>.bin, together with
    ! the grid geometry (lon, lat, area) and the map weights.
    !
    ! test/coords/bitident.sh <ref-commit> builds this same file at <ref-commit>
    ! and compares the stream files byte for byte. The program uses only API
    ! that already exists at the reference commit (b4c19a3).
    !
    ! Run on its own (make check), it performs the in-binary checks: remap and
    ! a direct map_field on the same map give identical results.

    use coords
    use coupler

    implicit none

    character(len=*), parameter :: fldr = "bitident_maps"
    character(len=*), parameter :: outd = "bitident_out"

    type(coupler_class) :: cpl
    type(grid_class)    :: gcc, gcf, gco, gll, gls, gga, ggd, gps, gpf
    real(dp), allocatable :: glat(:), glat_desc(:)
    integer :: fails

    fails = 0
    call execute_command_line("mkdir -p "//outd)

    ! --- grids ---------------------------------------------------------------
    call grid_init(gcc, name="bi-cart-c", mtype="cartesian", units="kilometers", &
                   x0=0.0_dp, dx=100.0_dp, nx=11, y0=0.0_dp, dy=100.0_dp, ny=11)
    call grid_init(gcf, name="bi-cart-f", mtype="cartesian", units="kilometers", &
                   x0=0.0_dp, dx=50.0_dp, nx=21, y0=0.0_dp, dy=50.0_dp, ny=21)
    call grid_init(gco, name="bi-cart-o", mtype="cartesian", units="kilometers", &
                   x0=35.0_dp, dx=70.0_dp, nx=13, y0=-20.0_dp, dy=90.0_dp, ny=9)
    call grid_init(gll, name="bi-ll2", mtype="latlon", units="degrees", planet="WGS84", &
                   x0=1.0_dp, dx=2.0_dp, nx=180, y0=-89.0_dp, dy=2.0_dp, ny=90)
    ! 0..360 centres stored as -180..180 (lon180 seam)
    call grid_init(gls, name="bi-ls", mtype="latlon", units="degrees", planet="WGS84", lon180=.true., &
                   x0=0.75_dp, dx=1.5_dp, nx=240, y0=-89.25_dp, dy=1.5_dp, ny=120)
    allocate(glat(48), glat_desc(48))
    call gaussian_latitudes_calc(48, glat)
    glat_desc = glat(48:1:-1)
    call grid_init(gga, name="bi-gauss", mtype="gaussian", units="degrees", planet="Spherical Earth", &
                   x0=0.0_dp, dx=3.75_dp, nx=96, y=glat)
    call grid_init(ggd, name="bi-gauss-desc", mtype="gaussian", units="degrees", planet="Spherical Earth", &
                   x0=0.0_dp, dx=3.75_dp, nx=96, y=glat_desc)
    call grid_init(gps, name="bi-ps", mtype="polar_stereographic", units="kilometers", planet="WGS84", &
                   lon180=.false., x0=-2000.0_dp, dx=100.0_dp, nx=41, y0=-2000.0_dp, dy=100.0_dp, ny=41, &
                   lambda=0.0_dp, phi=-71.0_dp)
    call grid_init(gpf, name="bi-ps-f", mtype="polar_stereographic", units="kilometers", planet="WGS84", &
                   lon180=.false., x0=-1000.0_dp, dx=50.0_dp, nx=41, y0=-1000.0_dp, dy=50.0_dp, ny=41, &
                   lambda=0.0_dp, phi=-71.0_dp)

    call dump_grid(gcc); call dump_grid(gcf); call dump_grid(gco)
    call dump_grid(gll); call dump_grid(gls); call dump_grid(gga); call dump_grid(ggd)
    call dump_grid(gps); call dump_grid(gpf)

    ! --- coupler -------------------------------------------------------------
    call coupler_init(cpl, map_fldr=fldr)
    call coupler_add_grid(cpl, "bi-cart-c", gcc)
    call coupler_add_grid(cpl, "bi-cart-f", gcf)
    call coupler_add_grid(cpl, "bi-cart-o", gco)
    call coupler_add_grid(cpl, "bi-ll2", gll)
    call coupler_add_grid(cpl, "bi-ls", gls)
    call coupler_add_grid(cpl, "bi-gauss", gga)
    call coupler_add_grid(cpl, "bi-gauss-desc", ggd)
    call coupler_add_grid(cpl, "bi-ps", gps)
    call coupler_add_grid(cpl, "bi-ps-f", gpf)

    call run_pair(gcc, gcf, "con");   call run_pair(gcc, gcf, "bilin"); call run_pair(gcc, gcf, "nn")
    call run_pair(gcf, gcc, "con")
    call run_pair(gcc, gco, "con");   call run_pair(gcc, gco, "bilin")
    call run_pair(gll, gls, "con");   call run_pair(gll, gls, "bilin"); call run_pair(gll, gls, "nn")
    call run_pair(gls, gll, "con")
    call run_pair(gll, gga, "con");   call run_pair(gll, gga, "bilin")
    call run_pair(ggd, gll, "con")
    call run_pair(gll, gps, "con");   call run_pair(gll, gps, "bilin"); call run_pair(gll, gps, "nn")
    call run_pair(gls, gps, "con")
    call run_pair(gps, gll, "con");   call run_pair(gps, gll, "bilin")
    call run_pair(gps, gls, "con")
    call run_pair(gps, gga, "con")
    call run_pair(gps, gpf, "con");   call run_pair(gpf, gps, "con")

    call run_points(gll)
    call run_disk(gps, gpf)

    if (fails > 0) then
        write(*,*) "FAIL: test_coupler_bitident (", fails, " in-binary mismatch(es))"
        stop 1
    end if
    write(*,*) "PASS: test_coupler_bitident (streams in "//outd//"/)"

contains

    subroutine run_pair(g1, g2, mtd)
        ! The fixed call set for one (src, dst, method) edge.
        type(grid_class), intent(in) :: g1, g2
        character(len=*), intent(in) :: mtd

        character(len=256) :: tag
        type(map_class) :: map
        real(dp), allocatable :: f1(:,:), f2(:,:), f2b(:,:), f13(:,:,:), f23(:,:,:)
        real(sp), allocatable :: s1(:,:), s2(:,:), s13(:,:,:), s23(:,:,:)
        integer,  allocatable :: i1(:,:), i2(:,:), i13(:,:,:), i23(:,:,:)
        logical,  allocatable :: m2(:,:), mp(:,:)
        integer :: nx2, ny2, i, j, k

        tag = trim(g1%name)//"__"//trim(g2%name)//"__"//trim(mtd)
        call field(g1, f1, .true.)
        nx2 = g2%G%nx; ny2 = g2%G%ny

        ! remap: dp / sp / int, stat variants
        call remap(cpl, f1, trim(g1%name), f2, trim(g2%name), method=mtd)
        call wr_dp2(trim(tag)//"__remap_dp", f2)
        call remap(cpl, f1, trim(g1%name), f2, trim(g2%name), method=mtd, stat="stdev")
        call wr_dp2(trim(tag)//"__remap_dp_stdev", f2)
        call remap(cpl, f1, trim(g1%name), f2, trim(g2%name), method=mtd, stat="count")
        call wr_dp2(trim(tag)//"__remap_dp_count", f2)

        allocate(s1(size(f1,1), size(f1,2)))
        s1 = real(f1, sp)
        call remap(cpl, s1, trim(g1%name), s2, trim(g2%name), method=mtd)
        call wr_sp2(trim(tag)//"__remap_sp", s2)

        allocate(i1(size(f1,1), size(f1,2)))
        i1 = categorical(f1)
        call remap(cpl, i1, trim(g1%name), i2, trim(g2%name), method=mtd, stat="count")
        call wr_int2(trim(tag)//"__remap_int", i2)

        ! remap: 3D
        allocate(f13(size(f1,1), size(f1,2), 3), s13(size(f1,1), size(f1,2), 3), i13(size(f1,1), size(f1,2), 3))
        do k = 1, 3
            f13(:,:,k) = f1
            where (f1 /= mv_dp) f13(:,:,k) = f1*real(k,dp) + 0.5_dp*real(k,dp)
            s13(:,:,k) = real(f13(:,:,k), sp)
            i13(:,:,k) = categorical(f13(:,:,k))
        end do
        call remap(cpl, f13, trim(g1%name), f23, trim(g2%name), method=mtd)
        call wr_dp3(trim(tag)//"__remap_dp3", f23)
        call remap(cpl, s13, trim(g1%name), s23, trim(g2%name), method=mtd)
        call wr_sp3(trim(tag)//"__remap_sp3", s23)
        call remap(cpl, i13, trim(g1%name), i23, trim(g2%name), method=mtd, stat="count")
        call wr_int3(trim(tag)//"__remap_int3", i23)

        ! direct map_field on the same edge (loaded from the coupler's cache)
        call map_init(map, g1, g2, method=mtd, fldr=fldr)
        call wr_map(trim(tag)//"__weights", map)

        allocate(f2b(nx2, ny2), m2(nx2, ny2), mp(nx2, ny2))
        call remap(cpl, f1, trim(g1%name), f2, trim(g2%name), method=mtd)
        call map_field(map, "bi", f1, f2b)
        if (any(f2b /= f2)) then
            write(*,*) "FAIL: remap /= map_field for ", trim(tag)
            fails = fails + 1
        end if

        call map_field(map, "bi", f1, f2b, mask2=m2)
        call wr_log2(trim(tag)//"__mf_mask2", m2)

        ! custom missing value
        f2 = f1
        where (f1 == mv_dp) f2 = -1.0e30_dp
        call map_field(map, "bi", f2, f2b, missing_value=-1.0e30_dp, mask2=m2)
        call wr_dp2(trim(tag)//"__mf_missing", f2b)

        ! reset=.false. keeps prefilled values where not mapped; mask_pack
        do j = 1, ny2
            do i = 1, nx2
                mp(i,j) = mod(i+j, 3) /= 0
            end do
        end do
        f2b = 7.0_dp
        call map_field(map, "bi", f1, f2b, reset=.false., mask_pack=mp, mask2=m2)
        call wr_dp2(trim(tag)//"__mf_noreset_pack", f2b)
        call wr_log2(trim(tag)//"__mf_noreset_pack_mask2", m2)

        ! fill + filter post-processing
        call map_field(map, "bi", f1, f2b, fill_method="nn")
        call wr_dp2(trim(tag)//"__mf_fill_nn", f2b)
        call map_field(map, "bi", f1, f2b, mask_pack=mp, fill_method="weighted", &
                       filt_method="gaussian", filt_par=[1.5_dp, 1.0_dp])
        call wr_dp2(trim(tag)//"__mf_fill_w_filt", f2b)

        ! sp map_field with missing value, reset=.false. and fill
        s2 = 3.0_sp
        call map_field(map, "bi", s1, s2, missing_value=-9999.0_sp, reset=.false., fill_method="nn")
        call wr_sp2(trim(tag)//"__mf_sp_fill", s2)

        ! int map_field, mask_pack
        call map_field(map, "bi", i1, i2, stat="count", mask_pack=mp, missing_value=-9999)
        call wr_int2(trim(tag)//"__mf_int_pack", i2)
    end subroutine run_pair

    subroutine run_points(g)
        ! grid -> points, points -> grid, points -> points (dp and sp).
        type(grid_class), intent(in) :: g
        type(points_class) :: p1, p2
        type(map_class)    :: map
        real(dp) :: px(12), py(12)
        real(dp), allocatable :: f1(:,:), v2(:), f2(:,:), v1(:)
        real(sp), allocatable :: s2(:), ss2(:,:)
        logical  :: m2(12)
        integer  :: k

        do k = 1, 12
            px(k) = 10.0_dp + 27.0_dp*real(k,dp)
            py(k) = -70.0_dp + 11.5_dp*real(k,dp)
        end do
        call points_init(p1, "bi-pts", mtype="latlon", units="degrees", planet="WGS84", &
                         lon180=.false., x=px, y=py)
        call points_init(p2, "bi-pts2", mtype="latlon", units="degrees", planet="WGS84", &
                         lon180=.false., x=px+3.0_dp, y=py-2.0_dp)

        call field(g, f1, .true.)
        allocate(v2(12), s2(12))

        call map_init(map, g, p1, max_neighbors=8, method="shepard")
        call map_field(map, "bi", f1, v2, mask2=m2)
        call wr_dp1("pts__grid_to_points_dp", v2)
        call map_field(map, "bi", real(f1,sp), s2)
        call wr_sp1("pts__grid_to_points_sp", s2)
        call map_field(map, "bi", f1, v2, method="nn")
        call wr_dp1("pts__grid_to_points_nn", v2)

        allocate(v1(12))
        do k = 1, 12
            v1(k) = 1.0_dp + 0.25_dp*real(k,dp)
        end do
        call map_init(map, p1, g, max_neighbors=4, method="shepard")
        allocate(f2(g%G%nx, g%G%ny), ss2(g%G%nx, g%G%ny))
        call map_field(map, "bi", v1, f2)
        call wr_dp2("pts__points_to_grid_dp", f2)
        call map_field(map, "bi", real(v1,sp), ss2, fill_method="nn")
        call wr_sp2("pts__points_to_grid_sp_fill", ss2)

        call map_init(map, p1, p2, max_neighbors=4, method="shepard")
        call map_field(map, "bi", v1, v2)
        call wr_dp1("pts__points_to_points_dp", v2)
        call map_field(map, "bi", real(v1,sp), s2)
        call wr_sp1("pts__points_to_points_sp", s2)
    end subroutine run_points

    subroutine run_disk(g1, g2)
        ! Coupler resolving both grids from grid_<name>.txt (no registration).
        type(grid_class), intent(in) :: g1, g2
        type(coupler_class) :: cpl2
        real(dp), allocatable :: f1(:,:), f2(:,:)
        character(len=*), parameter :: dfldr = "bitident_disk"

        call execute_command_line("mkdir -p "//dfldr)
        call grid_cdo_write_desc_short(g1, dfldr)
        call grid_cdo_write_desc_short(g2, dfldr)
        call coupler_init(cpl2, map_fldr=dfldr)
        call field(g1, f1, .true.)
        call remap(cpl2, f1, trim(g1%name), f2, trim(g2%name), method="con")
        call wr_dp2("disk__con", f2)
        call remap(cpl2, f1, trim(g1%name), f2, trim(g2%name), method="bilin")
        call wr_dp2("disk__bilin", f2)
    end subroutine run_disk

    subroutine field(g, f, holes)
        ! Smooth analytic field in axis units, with missing holes and a block.
        type(grid_class), intent(in) :: g
        real(dp), allocatable, intent(out) :: f(:,:)
        logical, intent(in) :: holes
        integer :: i, j, nx, ny
        nx = g%G%nx; ny = g%G%ny
        allocate(f(nx, ny))
        do j = 1, ny
            do i = 1, nx
                f(i,j) = 2.0_dp + sin(0.05_dp*g%x(i,j))*cos(0.03_dp*g%y(i,j)) + 0.001_dp*g%x(i,j)
                if (holes) then
                    if (mod(i + 2*j, 17) == 0) f(i,j) = mv_dp
                    if (i > nx/2 .and. i <= nx/2 + max(1,nx/8) .and. j > ny/3 .and. j <= ny/3 + max(1,ny/8)) &
                        f(i,j) = mv_dp
                end if
            end do
        end do
    end subroutine field

    elemental integer function categorical(v)
        real(dp), intent(in) :: v
        if (v == mv_dp) then
            categorical = -9999
        else
            categorical = mod(nint(v*3.0_dp), 5)
        end if
    end function categorical

    ! ----- stream writers --------------------------------------------------------

    integer function open_case(name) result(u)
        character(len=*), intent(in) :: name
        open(newunit=u, file=outd//"/"//trim(name)//".bin", access="stream", &
             form="unformatted", status="replace", action="write")
    end function open_case

    subroutine dump_grid(g)
        type(grid_class), intent(in) :: g
        integer :: u
        u = open_case("grid__"//trim(g%name))
        write(u) g%G%nx, g%G%ny, g%G%x, g%G%y, g%x, g%y, g%area
        if (allocated(g%lon)) write(u) g%lon, g%lat
        close(u)
    end subroutine dump_grid

    subroutine wr_map(name, map)
        character(len=*), intent(in) :: name
        type(map_class),  intent(in) :: map
        integer :: u
        u = open_case(name)
        write(u) map%wm%kind, map%wm%n_src, map%wm%n_dst, map%wm%n_links
        if (map%wm%n_links > 0) write(u) map%wm%src, map%wm%dst, map%wm%dst_off
        if (allocated(map%wm%w))    write(u) map%wm%w
        if (allocated(map%wm%dist)) write(u) map%wm%dist, map%wm%quadrant
        close(u)
    end subroutine wr_map

    subroutine wr_dp1(name, v)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: v(:)
        integer :: u
        u = open_case(name); write(u) shape(v), v; close(u)
    end subroutine wr_dp1

    subroutine wr_sp1(name, v)
        character(len=*), intent(in) :: name
        real(sp), intent(in) :: v(:)
        integer :: u
        u = open_case(name); write(u) shape(v), v; close(u)
    end subroutine wr_sp1

    subroutine wr_dp2(name, v)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: v(:,:)
        integer :: u
        u = open_case(name); write(u) shape(v), v; close(u)
    end subroutine wr_dp2

    subroutine wr_sp2(name, v)
        character(len=*), intent(in) :: name
        real(sp), intent(in) :: v(:,:)
        integer :: u
        u = open_case(name); write(u) shape(v), v; close(u)
    end subroutine wr_sp2

    subroutine wr_int2(name, v)
        character(len=*), intent(in) :: name
        integer, intent(in) :: v(:,:)
        integer :: u
        u = open_case(name); write(u) shape(v), v; close(u)
    end subroutine wr_int2

    subroutine wr_log2(name, v)
        character(len=*), intent(in) :: name
        logical, intent(in) :: v(:,:)
        integer :: u
        u = open_case(name); write(u) shape(v), merge(1, 0, v); close(u)
    end subroutine wr_log2

    subroutine wr_dp3(name, v)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: v(:,:,:)
        integer :: u
        u = open_case(name); write(u) shape(v), v; close(u)
    end subroutine wr_dp3

    subroutine wr_sp3(name, v)
        character(len=*), intent(in) :: name
        real(sp), intent(in) :: v(:,:,:)
        integer :: u
        u = open_case(name); write(u) shape(v), v; close(u)
    end subroutine wr_sp3

    subroutine wr_int3(name, v)
        character(len=*), intent(in) :: name
        integer, intent(in) :: v(:,:,:)
        integer :: u
        u = open_case(name); write(u) shape(v), v; close(u)
    end subroutine wr_int3

end program test_coupler_bitident
