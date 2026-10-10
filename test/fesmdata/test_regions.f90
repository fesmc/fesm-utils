program test_regions
    ! Exercise the regions module against FesmData v2 test files
    ! (NorthTest: NHT-16KM and NHT-32KM, regions and Zwally2012 basins):
    !   - region code arithmetic
    !   - FesmData layers, flag tables, selection expressions and named masks
    !   - layers of other files (no flag attributes, no grid_name, given
    !     names), of the program, and an empty set
    !   - writing layers and loading them again
    !   - online remapping NHT-16KM -> NHT-32KM, compared with the 32 km files
    ! The data folder is the first command-line argument, by default
    ! ../ice_data/v2/NorthTest (fesm-utils next to ice_data).

    use precision
    use coordinates, only : grid_class
    use grid_cdo,    only : grid_cdo_read_desc
    use ncio,        only : nc_create, nc_write_dim, nc_write
    use regions

    implicit none

    type(regions_class) :: reg, reg32, regr, rege, regw
    integer, allocatable :: region(:,:), zone(:,:), basin(:,:), group(:,:), bmask(:,:), ids(:,:)
    logical :: lmask(551,551)
    type(grid_class)    :: grid16, grid32
    character(len=1024) :: fldr, path16, path32, basins16, basins32
    integer, allocatable :: codes(:)
    integer :: nfail
    real(wp) :: agree

    nfail = 0

    fldr = "../ice_data/v2/NorthTest"
    if (command_argument_count() .ge. 1) call get_command_argument(1, fldr)

    path16   = trim(fldr)//"/NHT-16KM/NHT-16KM_REGIONS.nc"
    path32   = trim(fldr)//"/NHT-32KM/NHT-32KM_REGIONS.nc"
    basins16 = trim(fldr)//"/NHT-16KM/NHT-16KM_BASINS-Zwally2012.nc"
    basins32 = trim(fldr)//"/NHT-32KM/NHT-32KM_BASINS-Zwally2012.nc"

    call write_griddes("maps", "NHT-16KM", 551, 16.0_dp)
    call write_griddes("maps", "NHT-32KM", 276, 32.0_dp)
    call grid_cdo_read_desc(grid16, "NHT-16KM", "maps")
    call grid_cdo_read_desc(grid32, "NHT-32KM", "maps")

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
    ! FesmData layers on the grid of the files
    ! ========================================================================
    call regions_init(reg)
    call regions_load_fesmdata(reg, path16)
    call regions_load_basins(reg, "Zwally2012", basins16)
    call regions_add_mask(reg, "grl_land",  "region:Greenland & zone:land")
    call regions_add_mask(reg, "grl_shelf", "region:Greenland & zone:continental_shelf")

    call check("size",           reg%nx .eq. 551 .and. reg%ny .eq. 551, nfail)
    call check("layers",         size(reg%layers) .eq. 5 .and. regions_find(reg, "ZWALLY2012.group") .eq. 4 &
                                 .and. regions_find(reg, "Zwally2012.mask") .eq. 5, nfail)
    call check("region hier",    reg%layers(1)%hier .and. .not. reg%layers(2)%hier, nfail)

    region = reg%layers(regions_find(reg, "region"))%values
    zone   = reg%layers(regions_find(reg, "zone"))%values
    basin  = reg%layers(regions_find(reg, "Zwally2012"))%values
    group  = reg%layers(regions_find(reg, "Zwally2012.group"))%values
    bmask  = reg%layers(regions_find(reg, "Zwally2012.mask"))%values

    codes = flag_codes(reg%layers(1)%tab, "Greenland")
    call check("name Greenland", size(codes) .eq. 1, nfail)
    if (size(codes) .eq. 1) call check("code Greenland", codes(1) .eq. 103, nfail)
    call check("level 1 name",   trim(flag_name(reg%layers(1)%tab, 1)) .eq. "Northern_Hemisphere", nfail)
    call check("level 3 name",   trim(flag_name(reg%layers(1)%tab, 10203)) .eq. "Svalbard", nfail)
    call check("zone names",     trim(flag_name(reg%layers(2)%tab, 2)) .eq. "continental_shelf", nfail)

    ! Selection expressions against the fields themselves
    call check("sel region",     all(regions_select(reg, "region:Greenland") .eqv. in_region(region, 103)), nfail)
    call check("sel path/case",  all(regions_select(reg, " REGION : 1.3 ") .eqv. &
                                     regions_select(reg, "region:greenland")), nfail)
    call check("sel level 3",    all(regions_select(reg, "region:Svalbard,Fennoscandia") .eqv. &
                                     (region .eq. 10203 .or. region .eq. 10204)), nfail)
    call check("sel spaces",     all(regions_select(reg, "region:North America") .eqv. in_region(region, 101)), nfail)
    call check("sel level 1",    all(regions_select(reg, "region:Northern_Hemisphere") .eqv. (region .gt. 0)), nfail)
    call check("sel and",        all(regions_select(reg, "region:Greenland & zone:land") .eqv. &
                                     (in_region(region, 103) .and. zone .eq. 3)), nfail)
    call check("sel not",        all(regions_select(reg, "~region:Greenland") .eqv. .not. in_region(region, 103)), nfail)
    call check("sel or",         all(regions_select(reg, "zone:0 | zone:land") .eqv. (zone .eq. 0 .or. zone .eq. 3)), nfail)
    call check("sel precedence", all(regions_select(reg, "zone:land & region:Greenland | zone:open_ocean") .eqv. &
                                     ((zone .eq. 3 .and. in_region(region,103)) .or. zone .eq. 0)), nfail)
    call check("sel basins",     all(regions_select(reg, "Zwally2012:11,12") .eqv. (basin .eq. 11 .or. basin .eq. 12)), nfail)
    call check("sel group",      all(regions_select(reg, "zwally2012.group:1") .eqv. (group .eq. 1)), nfail)
    call check("sel mask",       all(regions_select(reg, "Zwally2012.mask:1") .eqv. (bmask .eq. 1)), nfail)
    call check("sel all/none",   all(regions_select(reg, "all")) .and. .not. any(regions_select(reg, "none")), nfail)
    call check("basins in grl",  .not. any(regions_select(reg, "Zwally2012:11 & ~region:Greenland")), nfail)
    call check("named mask",     all(regions_mask(reg, "grl_shelf") .eqv. (in_region(region, 103) .and. zone .eq. 2)), nfail)
    call check("named non-empty", count(regions_mask(reg, "grl_land")) .gt. 0, nfail)

    call check("basin_ids",       all(regions_basin_ids(reg, "Zwally2012", extent=lmask) .eq. basin) &
                                  .and. all(lmask .eqv. (bmask .eq. 1)), nfail)
    call check("basin_ids group", all(regions_basin_ids(reg, "Zwally2012.group", extent=lmask) .eq. group) &
                                  .and. all(lmask .eqv. (bmask .eq. 1)), nfail)
    call check("basin_ids None",  all(regions_basin_ids(reg, "None") .eq. 0), nfail)
    call check("basin_ids domain", all(regions_basin_ids(reg, "domain", extent=lmask) .eq. 1) .and. all(lmask), nfail)

    ! ========================================================================
    ! Layers of another file: ids without flag attributes or grid_name (the
    ! Zwally2012 basins, with -1 for no basin), with and without given names
    ! ========================================================================
    ids = merge(basin, -9999, basin .gt. 0)
    call nc_create("test_regions_custom.nc")
    call nc_write_dim("test_regions_custom.nc", "x", x=1.0_dp, dx=1.0_dp, nx=551)
    call nc_write_dim("test_regions_custom.nc", "y", x=1.0_dp, dx=1.0_dp, nx=551)
    call nc_write("test_regions_custom.nc", "my_basins", ids, dim1="x", dim2="y")

    call regions_load_layer(reg, "custom", "test_regions_custom.nc", "my_basins")
    call regions_load_layer(reg, "named",  "test_regions_custom.nc", "my_basins", &
                            codes=[11, 12], names=["west ", "north"])

    associate(c => reg%layers(regions_find(reg, "custom")))
        call check("custom values",   all(c%values .eq. merge(basin, regions_undefined, basin .gt. 0)), nfail)
        call check("custom names",    all(c%tab%codes .eq. reg%layers(3)%tab%codes) .and. &
                                      trim(c%tab%names(1)) .eq. "11", nfail)
    end associate
    call check("custom basin_ids", all(regions_basin_ids(reg, "custom", extent=lmask) .eq. basin) &
                                  .and. all(lmask .eqv. (basin .gt. 0)), nfail)
    call check("custom select",   all(regions_select(reg, "custom:11,12") .eqv. &
                                      regions_select(reg, "Zwally2012:11,12")), nfail)
    call check("given names",     all(regions_select(reg, "named:West | NAMED:north") .eqv. &
                                      regions_select(reg, "Zwally2012:11,12")), nfail)

    ! A layer of the program, hierarchical
    call regions_add_layer(reg, "myregions", merge(10301, 201, in_region(region, 103)), &
                           codes=[1, 103, 10301, 2, 201], &
                           names=["north    ", "green    ", "greenland", "rest     ", "elsewhere"], hier=.true.)
    call check("program layer",   all(regions_select(reg, "myregions:north") .eqv. in_region(region, 103)) .and. &
                                  all(regions_select(reg, "myregions:2.1 | myregions:green") .eqv. .true.), nfail)

    ! An empty set
    call regions_init(rege, nx=5, ny=4)
    call check("empty all/none",  all(regions_select(rege, "all")) .and. .not. any(regions_select(rege, "none")) &
                                  .and. size(regions_select(rege, "all")) .eq. 20, nfail)
    call check("empty basin_ids", all(regions_basin_ids(rege, "domain") .eq. 1) .and. &
                                  all(regions_basin_ids(rege, "None") .eq. 0), nfail)
    call regions_end(rege)

    ! ========================================================================
    ! Writing layers, and loading them again (tables and hierarchy kept)
    ! ========================================================================
    call regions_write(reg, "test_regions_write.nc", grid=grid16, layers=["region   ", "myregions"])

    call regions_init(regw, grid=grid16)
    call regions_load_layer(regw, "region", "test_regions_write.nc")
    call regions_load_layer(regw, "mine",   "test_regions_write.nc", "myregions")
    call check("write region",    regw%layers(1)%hier .and. all(regw%layers(1)%values .eq. region) .and. &
                                  all(regw%layers(1)%tab%codes .eq. reg%layers(1)%tab%codes), nfail)
    call check("write select",    all(regions_select(regw, "region:Northern_Hemisphere & ~region:Greenland") .eqv. &
                                      regions_select(reg,  "region:Northern_Hemisphere & ~region:Greenland")), nfail)
    call check("write names",     all(regions_select(regw, "mine:greenland") .eqv. in_region(region, 103)), nfail)
    call regions_end(regw)

    ! ========================================================================
    ! Remapping NHT-16KM -> NHT-32KM, against the 32 km files (made by a
    ! dominant-class remap, so nearest neighbour agrees on most cells)
    ! ========================================================================
    call regions_init(reg32)
    call regions_load_fesmdata(reg32, path32)
    call regions_load_basins(reg32, "Zwally2012", basins32)

    call regions_init(regr, grid=grid32)
    call regions_load_fesmdata(regr, path16)
    call regions_load_basins(regr, "Zwally2012", basins16)
    call regions_load_layer(regr, "custom", "test_regions_custom.nc", "my_basins", grid_name="NHT-16KM")

    call check("remap size",     regr%nx .eq. 276 .and. regr%ny .eq. 276, nfail)
    call check("remap filled",   all(regr%layers(1)%values .gt. 0) .and. all(regr%layers(2)%values .ge. 0), nfail)

    agree = real(count(regr%layers(1)%values .eq. reg32%layers(1)%values),wp) / real(276*276,wp)
    write(*,"(a,f6.3)") "  agreement region:   ", agree
    call check("remap region",   agree .gt. 0.95_wp, nfail)

    agree = real(count(regr%layers(2)%values .eq. reg32%layers(2)%values),wp) / real(276*276,wp)
    write(*,"(a,f6.3)") "  agreement zone:     ", agree
    call check("remap zone",     agree .gt. 0.95_wp, nfail)

    agree = real(count(regr%layers(3)%values .eq. reg32%layers(3)%values),wp) / real(276*276,wp)
    write(*,"(a,f6.3)") "  agreement basin:    ", agree
    call check("remap basin",    agree .gt. 0.95_wp, nfail)

    call check("remap grid_name", all(regions_basin_ids(regr, "custom") .eq. regions_basin_ids(regr, "Zwally2012")), nfail)
    call check("remap select",   all(regions_select(regr, "region:Greenland") .eqv. &
                                     in_region(regr%layers(1)%values, 103)), nfail)

    call regions_end(reg)
    call regions_end(reg32)
    call regions_end(regr)
    call delete_file("test_regions_custom.nc")
    call delete_file("test_regions_write.nc")

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
