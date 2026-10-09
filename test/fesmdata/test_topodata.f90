program test_topodata
    ! Exercise the topodata module against a FesmData v2 test file
    ! (NorthTest: NHT-16KM_TOPO-GEBCO2025.nc, on the grid NH-16KM):
    !   - loading the default and a full variable list, flag tables
    !   - conservative remapping NH-16KM -> NHT-32KM (fractions sum to one,
    !     the area-weighted mean bed elevation is kept within the interior of
    !     the 32 km grid; its outer ring reaches 8 km beyond the 16 km grid)
    ! The data folder is the first command-line argument, by default
    ! ../ice_data/v2/NorthTest (fesm-utils next to ice_data).

    use precision
    use constants,   only : mv
    use coordinates, only : grid_class
    use grid_cdo,    only : grid_cdo_read_desc
    use ncio,        only : nc_read
    use fesmdata,    only : flag_name
    use topodata

    implicit none

    type(topodata_class) :: td, tdr
    type(grid_class)     :: grid16, grid32
    character(len=1024)  :: fldr, path
    real(wp), allocatable :: z_bed(:,:), fsum(:,:)
    real(dp), allocatable :: w16(:,:), w32(:,:)
    logical,  allocatable :: valid(:,:)
    real(dp) :: mean16, mean32
    integer  :: nfail

    nfail = 0

    fldr = "../ice_data/v2/NorthTest"
    if (command_argument_count() .ge. 1) call get_command_argument(1, fldr)
    path = trim(fldr)//"/NHT-16KM/NHT-16KM_TOPO-GEBCO2025.nc"

    ! ========================================================================
    ! Default variables, on the grid of the file
    ! ========================================================================
    call topodata_init_arg(td, path)

    call check("size",          td%nx .eq. 551 .and. td%ny .eq. 551, nfail)
    call check("grid_src",      trim(td%par%grid_src) .eq. "NH-16KM", nfail)
    call check("product",       trim(td%par%product) .eq. "GEBCO2025", nfail)
    call check("no remap",      .not. td%remap, nfail)
    call check("default vars",  allocated(td%z_bed) .and. allocated(td%z_srf) .and. &
                                allocated(td%H_ice) .and. allocated(td%z_bed_sd), nfail)
    call check("not loaded",    .not. allocated(td%f_ocn) .and. .not. allocated(td%mask), nfail)

    allocate(z_bed(551,551))
    call nc_read(path, "z_bed", z_bed, missing_value=mv)
    call check("z_bed as file", all(td%z_bed .eq. z_bed), nfail)
    call check("H_ice >= 0",    all(td%H_ice .ge. 0.0_wp .or. td%H_ice .eq. mv), nfail)

    call topodata_end(td)

    ! ========================================================================
    ! All variables
    ! ========================================================================
    call topodata_init_arg(td, path, vars=["z_bed ", "H_ice ", "f_ocn ", "f_land", &
                                           "f_grnd", "f_flt ", "mask  ", "src_id"])

    call check("mask names",    trim(flag_name(td%tab_mask, topo_mask_grnd)) .eq. "grounded_ice" &
                                .and. size(td%tab_mask%codes) .eq. 4, nfail)
    call check("src names",     size(td%tab_src%codes) .eq. 3, nfail)
    call check("mask classes",  all(td%mask .ge. 0 .and. td%mask .le. 3), nfail)

    valid = (td%f_ocn .ne. mv)
    fsum  = td%f_ocn + td%f_land + td%f_grnd + td%f_flt
    call check("fractions sum", all(abs(fsum - 1.0_wp) .lt. 1e-3_wp .or. .not. valid), nfail)
    call check("mask ice",      all((td%mask .ge. topo_mask_grnd) .eqv. &
                                    (td%f_grnd + td%f_flt .ge. 0.5_wp)), nfail)

    call write_griddes("maps", "NH-16KM",  551, 16.0_dp)
    call grid_cdo_read_desc(grid16, "NH-16KM", "maps")
    ! Footprint of the interior of the 32 km grid (cells 2..275): 16 km cells
    ! 3..549 fully, 2 and 550 half
    w16 = edge_weights(551)*grid16%area
    mean16 = sum(td%z_bed*w16, mask=valid) / sum(w16, mask=valid)

    ! ========================================================================
    ! Conservative remapping NH-16KM -> NHT-32KM
    ! ========================================================================
    call write_griddes("maps", "NHT-32KM", 276, 32.0_dp)
    call grid_cdo_read_desc(grid32, "NHT-32KM", "maps")

    call topodata_init_arg(tdr, path, vars=["z_bed ", "H_ice ", "f_ocn ", "f_land", &
                                            "f_grnd", "f_flt ", "mask  "], grid=grid32)

    call check("remap",         tdr%remap, nfail)
    call check("remap size",    tdr%nx .eq. 276 .and. tdr%ny .eq. 276, nfail)
    call check("remap mask",    all(tdr%mask .ge. 0 .and. tdr%mask .le. 3), nfail)

    valid = (tdr%f_ocn .ne. mv)
    fsum  = tdr%f_ocn + tdr%f_land + tdr%f_grnd + tdr%f_flt
    write(*,"(a,f8.5)") "  max |sum(f) - 1|:  ", maxval(abs(fsum - 1.0_wp), mask=valid)
    call check("remap fractions", all(abs(fsum - 1.0_wp) .lt. 1e-3_wp .or. .not. valid), nfail)
    call check("remap H_ice",   all(tdr%H_ice .ge. -1e-3_wp .or. tdr%H_ice .eq. mv), nfail)

    allocate(w32(276,276))
    w32 = 0.0_dp
    w32(2:275,2:275) = grid32%area(2:275,2:275)
    mean32 = sum(tdr%z_bed*w32, mask=valid) / sum(w32, mask=valid)
    write(*,"(a,2f10.2)") "  area mean z_bed 16/32 km: ", mean16, mean32
    call check("remap mean z_bed", abs(mean32 - mean16) .lt. 1e-3_dp*abs(mean16), nfail)

    call topodata_end(td)
    call topodata_end(tdr)

    write(*,*)
    if (nfail .gt. 0) then
        write(*,"(a,i0,a)") "test_topodata: ", nfail, " check(s) FAILED"
        stop 1
    end if
    write(*,"(a)") "test_topodata: all checks passed"

contains

    function edge_weights(n) result(w)
        ! Weights of the cells of an n x n 16 km grid within the interior of
        ! the 32 km grid: 0 for the outer cells, 1/2 for the next ones
        integer, intent(in) :: n
        real(dp) :: w(n,n)
        real(dp) :: w1(n)
        integer  :: j
        w1 = 1.0_dp
        w1([1,n])   = 0.0_dp
        w1([2,n-1]) = 0.5_dp
        do j = 1, n
            w(:,j) = w1*w1(j)
        end do
    end function edge_weights

    subroutine check(label, ok, nfail)
        character(len=*), intent(in)    :: label
        logical,          intent(in)    :: ok
        integer,          intent(inout) :: nfail
        if (ok) then
            write(*,"(a,a)") "  PASS  ", label
        else
            write(*,"(a,a)") "  FAIL  ", label
            nfail = nfail + 1
        end if
    end subroutine check

    subroutine write_griddes(fldr, name, n, dx)
        ! cdo grid description of a NorthTest grid, as FesmData writes it
        character(len=*), intent(in) :: fldr, name
        integer,          intent(in) :: n
        real(dp),         intent(in) :: dx
        integer :: u
        call execute_command_line("mkdir -p "//trim(fldr))
        open(newunit=u, file=trim(fldr)//"/grid_"//trim(name)//".txt", status="replace", action="write")
        write(u,"(a)")      "gridtype = projection"
        write(u,"(a,i0)")   "gridsize = ", n*n
        write(u,"(a,i0)")   "xsize    = ", n
        write(u,"(a,i0)")   "ysize    = ", n
        write(u,"(a)")      "xname    = xc"
        write(u,"(a)")      "xunits   = km"
        write(u,"(a)")      "yname    = yc"
        write(u,"(a)")      "yunits   = km"
        write(u,"(a)")      "xfirst   = -4900.0"
        write(u,"(a,f0.1)") "xinc     = ", dx
        write(u,"(a)")      "yfirst   = -5400.0"
        write(u,"(a,f0.1)") "yinc     = ", dx
        write(u,"(a)")      "grid_mapping = crs"
        write(u,"(a)")      "grid_mapping_name = polar_stereographic"
        write(u,"(a)")      "straight_vertical_longitude_from_pole = -45.0"
        write(u,"(a)")      "latitude_of_projection_origin = 90.0"
        write(u,"(a)")      "standard_parallel = 70.0"
        write(u,"(a)")      "false_easting = 0.0"
        write(u,"(a)")      "false_northing = 0.0"
        write(u,"(a)")      "semi_major_axis = 6378137.0"
        write(u,"(a)")      "inverse_flattening = 298.25722356"
        close(u)
    end subroutine write_griddes

end program test_topodata
