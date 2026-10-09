program test_regions
    ! Exercise the regions module against FesmData v2 test files
    ! (NorthTest: NHT-16KM and NHT-32KM, regions and Zwally2012 basins):
    !   - region code arithmetic
    !   - loading, flag tables, selection expressions and named masks
    !   - online remapping NHT-16KM -> NHT-32KM, compared with the 32 km files
    ! The data folder is the first command-line argument, by default
    ! ../ice_data/v2/NorthTest (fesm-utils next to ice_data).

    use precision
    use coordinates, only : grid_class
    use grid_cdo,    only : grid_cdo_read_desc
    use ncio,        only : nc_read, nc_create, nc_write_dim, nc_write
    use regions

    implicit none

    type(regions_class) :: reg, reg32, regr, regc
    integer, allocatable :: ids(:,:)
    logical :: lmask(551,551)
    type(grid_class)    :: grid32
    character(len=1024) :: fldr, path16, path32, basins16, basins32
    integer, allocatable :: codes(:)
    integer :: nfail
    real(wp) :: agree

    nfail = 0

    fldr = "../ice_data/v2/NorthTest"
    if (command_argument_count() .ge. 1) call get_command_argument(1, fldr)

    path16   = trim(fldr)//"/NHT-16KM/NHT-16KM_REGIONS.nc"
    path32   = trim(fldr)//"/NHT-32KM/NHT-32KM_REGIONS.nc"
    basins16 = trim(fldr)//"/NHT-16KM/NHT-16KM_BASINS-{set}.nc"
    basins32 = trim(fldr)//"/NHT-32KM/NHT-32KM_BASINS-{set}.nc"

    ! ========================================================================
    ! Region code arithmetic
    ! ========================================================================
    call check("region_code 1.3.1",      region_code("1.3.1") .eq. 10301, nfail)
    call check("region_code 1.51",       region_code("1.51")  .eq. 151,   nfail)
    call check("region_code invalid",    region_code("1.0.2") .eq. 0,     nfail)
    call check("region_level",           all(region_level([0,1,99,101,9999,10301]) .eq. [0,1,1,2,2,3]), nfail)
    call check("region_ancestor",        all(region_ancestor(10301,[1,2,3,4]) .eq. [1,103,10301,0]), nfail)
    call check("region_path",            trim(region_path(10301)) .eq. "1.3.1", nfail)
    call check("in_region",              all(in_region([0,1,103,10301,10401,151],103) .eqv. &
                                             [.false.,.false.,.true.,.true.,.false.,.false.]), nfail)

    ! ========================================================================
    ! Loading on the grid of the files
    ! ========================================================================
    call regions_init_arg(reg, path16, basins16, ["Zwally2012"], &
                          mask_names=["grl_land ", "grl_shelf"], &
                          mask_exprs=["region:Greenland & zone:land             ", &
                                      "region:Greenland & zone:continental_shelf"])

    call check("size",           reg%nx .eq. 551 .and. reg%ny .eq. 551, nfail)
    call check("grid_src",       trim(reg%par%grid_src) .eq. "NHT-16KM", nfail)
    call check("no remap",       .not. reg%remap, nfail)

    codes = flag_codes(reg%tab_region, "Greenland")
    call check("name Greenland", size(codes) .eq. 1, nfail)
    if (size(codes) .eq. 1) call check("code Greenland", codes(1) .eq. 103, nfail)
    call check("level 1 name",   trim(flag_name(reg%tab_region, 1)) .eq. "Northern_Hemisphere", nfail)
    call check("level 3 name",   trim(flag_name(reg%tab_region, 10203)) .eq. "Svalbard", nfail)
    call check("zone names",     trim(flag_name(reg%tab_zone, 2)) .eq. "continental_shelf", nfail)
    call check("basin group",    reg%basins(1)%with_group, nfail)

    ! Selection expressions against the fields themselves
    call check("sel region",     all(regions_select(reg, "region:Greenland") .eqv. &
                                     reg%region_3/100**max(region_level(reg%region_3)-2,0) .eq. 103), nfail)
    call check("sel path/case",  all(regions_select(reg, " REGION : 1.3 ") .eqv. &
                                     regions_select(reg, "region:greenland")), nfail)
    call check("sel level 3",    all(regions_select(reg, "region:Svalbard,Fennoscandia") .eqv. &
                                     (reg%region_3 .eq. 10203 .or. reg%region_3 .eq. 10204)), nfail)
    call check("sel spaces",     all(regions_select(reg, "region:North America") .eqv. &
                                     in_region(reg%region_3, 101)), nfail)
    call check("sel level 1",    all(regions_select(reg, "region:Northern_Hemisphere") .eqv. &
                                     (reg%region_3 .gt. 0)), nfail)
    call check("sel and",        all(regions_select(reg, "region:Greenland & zone:land") .eqv. &
                                     (in_region(reg%region_3, 103) .and. reg%zone .eq. 3)), nfail)
    call check("sel not",        all(regions_select(reg, "~region:Greenland") .eqv. &
                                     .not. in_region(reg%region_3, 103)), nfail)
    call check("sel or",         all(regions_select(reg, "zone:0 | zone:land") .eqv. &
                                     (reg%zone .eq. 0 .or. reg%zone .eq. 3)), nfail)
    call check("sel precedence", all(regions_select(reg, "zone:land & region:Greenland | zone:open_ocean") .eqv. &
                                     ((reg%zone .eq. 3 .and. in_region(reg%region_3,103)) .or. reg%zone .eq. 0)), nfail)
    call check("sel basins",     all(regions_select(reg, "Zwally2012:11,12") .eqv. &
                                     (reg%basins(1)%basin .eq. 11 .or. reg%basins(1)%basin .eq. 12)), nfail)
    call check("sel group",      all(regions_select(reg, "zwally2012.group:1") .eqv. &
                                     (reg%basins(1)%basin_group .eq. 1)), nfail)
    call check("sel mask",       all(regions_select(reg, "Zwally2012.mask:1") .eqv. &
                                     (reg%basins(1)%basin_mask .eq. 1)), nfail)
    call check("sel all/none",   all(regions_select(reg, "all")) .and. .not. any(regions_select(reg, "none")), nfail)
    call check("basins in grl",  .not. any(regions_select(reg, "Zwally2012:11 & ~region:Greenland")), nfail)
    call check("named mask",     all(regions_mask(reg, "grl_shelf") .eqv. &
                                     (in_region(reg%region_3, 103) .and. reg%zone .eq. 2)), nfail)
    call check("named non-empty", count(regions_mask(reg, "grl_land")) .gt. 0, nfail)

    ! ========================================================================
    ! Custom basin set: a field of ids without flag attributes or grid_name
    ! (the Zwally2012 basins, with -1 for no basin)
    ! ========================================================================
    ids = merge(reg%basins(1)%basin, -1, reg%basins(1)%basin .gt. 0)
    call nc_create("test_regions_custom.nc")
    call nc_write_dim("test_regions_custom.nc", "xc", x=1.0_dp, dx=1.0_dp, nx=551)
    call nc_write_dim("test_regions_custom.nc", "yc", x=1.0_dp, dx=1.0_dp, nx=551)
    call nc_write("test_regions_custom.nc", "my_basins", ids, dim1="xc", dim2="yc")

    call regions_init_arg(regc, path16, basins16, ["Zwally2012", "custom    "], &
                          basin_paths=["                      ", "test_regions_custom.nc"], &
                          basin_vars=["         ", "my_basins"])
    call check("custom ids",      all(regc%basins(2)%basin .eq. reg%basins(1)%basin), nfail)
    call check("custom names",    all(regc%basins(2)%tab_basin%codes .eq. reg%basins(1)%tab_basin%codes) .and. &
                                  trim(regc%basins(2)%tab_basin%names(1)) .eq. "11", nfail)
    call check("basin_ids extent", all(regions_basin_ids(regc, "custom", extent=lmask) .eq. reg%basins(1)%basin) &
                                  .and. all(lmask .eqv. (reg%basins(1)%basin .gt. 0)), nfail)
    call check("basin_ids group", all(regions_basin_ids(regc, "Zwally2012.group", extent=lmask) .eq. &
                                  reg%basins(1)%basin_group) .and. all(lmask .eqv. (reg%basins(1)%basin_mask .eq. 1)), nfail)
    call check("basin_ids None",  all(regions_basin_ids(regc, "None") .eq. 0), nfail)
    call check("basin_ids domain", all(regions_basin_ids(regc, "domain", extent=lmask) .eq. 1) .and. all(lmask), nfail)
    call check("custom no group", .not. regc%basins(2)%with_group .and. .not. regc%basins(2)%with_mask, nfail)
    call check("custom select",   all(regions_select(regc, "custom:11,12") .eqv. &
                                      regions_select(regc, "Zwally2012:11,12")), nfail)
    call regions_end(regc)
    call delete_file("test_regions_custom.nc")

    ! ========================================================================
    ! Remapping NHT-16KM -> NHT-32KM, against the 32 km files (made by a
    ! dominant-class remap, so nearest neighbour agrees on most cells)
    ! ========================================================================
    call write_griddes("maps", "NHT-16KM", 551, 16.0_dp)
    call write_griddes("maps", "NHT-32KM", 276, 32.0_dp)
    call grid_cdo_read_desc(grid32, "NHT-32KM", "maps")

    call regions_init_arg(reg32, path32, basins32, ["Zwally2012"])
    call regions_init_arg(regr,  path16, basins16, ["Zwally2012"], grid=grid32)

    call check("remap",          regr%remap, nfail)
    call check("remap size",     regr%nx .eq. 276 .and. regr%ny .eq. 276, nfail)
    call check("remap filled",   all(regr%region_3 .gt. 0) .and. all(regr%zone .ge. 0), nfail)
    call check("remap nested",   all(in_region(regr%region_3, regr%region_2)), nfail)

    agree = real(count(regr%region_3 .eq. reg32%region_3),wp) / real(size(regr%region_3),wp)
    write(*,"(a,f6.3)") "  agreement region_3: ", agree
    call check("remap region_3", agree .gt. 0.95_wp, nfail)

    agree = real(count(regr%zone .eq. reg32%zone),wp) / real(size(regr%zone),wp)
    write(*,"(a,f6.3)") "  agreement zone:     ", agree
    call check("remap zone",     agree .gt. 0.95_wp, nfail)

    agree = real(count(regr%basins(1)%basin .eq. reg32%basins(1)%basin),wp) / real(size(regr%zone),wp)
    write(*,"(a,f6.3)") "  agreement basin:    ", agree
    call check("remap basin",    agree .gt. 0.95_wp, nfail)

    call check("remap select",   all(regions_select(regr, "region:Greenland") .eqv. &
                                     in_region(regr%region_3, 103)), nfail)

    call regions_end(reg)
    call regions_end(reg32)
    call regions_end(regr)

    write(*,*)
    if (nfail .gt. 0) then
        write(*,"(a,i0,a)") "test_regions: ", nfail, " check(s) FAILED"
        stop 1
    end if
    write(*,"(a)") "test_regions: all checks passed"

contains

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

    subroutine delete_file(filename)
        character(len=*), intent(in) :: filename
        integer :: u
        open(newunit=u, file=filename, status="old")
        close(u, status="delete")
    end subroutine delete_file

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

end program test_regions
