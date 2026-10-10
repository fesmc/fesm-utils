program test_htopo
    ! Exercise the htopo module against FesmData v2 test files (NorthTest).
    ! The hub grid is NHT-16KM; the topography file is on NH-16KM (the same
    ! geometry under another name), so it is remapped conservatively onto the
    ! hub. Run from the fesm-utils folder, next to ice_data.

    use precision
    use constants,      only : mv
    use ncio
    use coordinates,    only : grid_class
    use phys_constants, only : phys_const_class, phys_const_load
    use regions,        only : regions_class, regions_mask, regions_find
    use htopo

    implicit none

    character(len=*), parameter :: nml_file = "test/fesmdata/test_htopo.nml"
    integer, parameter :: n_band = 10           ! width of the missing border band [cells]

    type(htopo_class) :: ht, ht_def, ht_gaps
    type(phys_const_class) :: cnst
    real(wp), allocatable :: z_bed(:,:), zero(:,:), z_srf_exp(:,:)
    logical,  allocatable :: gap(:,:)
    integer :: nfail, j

    nfail = 0

    call phys_const_load(cnst, "par/phys_const_earth.nml", group="Earth")

    call write_griddes("maps", "NHT-16KM")
    call write_griddes("maps", "NH-16KM")
    call write_custom_basins("test_htopo_custom.nc")

    ! ========================================================================
    ! Full configuration
    ! ========================================================================
    call htopo_init(ht, nml_file, "domain", "topo", "regions", "NorthTest", "NHT-16KM", cnst)

    call check("size",          ht%nx .eq. 551 .and. ht%ny .eq. 551, nfail)
    call check("topo remapped", ht%ref%remap, nfail)

    allocate(z_bed(551,551))
    call nc_read("../ice_data/v2/NorthTest/NHT-16KM/NHT-16KM_TOPO-GEBCO2025.nc", "z_bed", z_bed, missing_value=mv)
    call check("z_bed (identity remap)", maxval(abs(ht%ref%z_bed - z_bed)) .lt. 1e-2_wp, nfail)
    call check("z_bed_sd",      maxval(ht%ref%z_bed_sd) .gt. 0.0_wp, nfail)

    call check("basins",           all(htopo_basins(ht) .eq. layer(ht%reg, "Zwally2012")) &
                                   .and. maxval(htopo_basins(ht)) .gt. 0, nfail)
    ht%par%basins = "Zwally2012.group"
    call check("basins group",     all(htopo_basins(ht) .eq. layer(ht%reg, "Zwally2012.group")), nfail)
    ht%par%basins = "Zwally2012"
    call check("custom layer",     all(layer(ht%reg, "custom") .eq. 3), nfail)
    call check("named mask",       count(regions_mask(ht%reg, "grl")) .gt. 0, nfail)

    ! Current geometry: starts from the reference; no anomaly keeps it
    call check("current = ref",    all(ht%z_bed .eq. ht%ref%z_bed) .and. all(ht%H_ice .eq. ht%ref%H_ice), nfail)
    allocate(zero(ht%nx,ht%ny)); zero = 0.0_wp
    call htopo_update(ht, zero, zero, zero)
    call check("update no anomaly", all(ht%z_bed .eq. ht%ref%z_bed) .and. &
                                    all(ht%H_ice .eq. max(ht%ref%H_ice, 0.0_wp)), nfail)
    call check("update f_grnd",     all((ht%f_grnd .eq. 1.0_wp) .eqv. (ht%z_bed .ge. 0.0_wp .or. &
                ht%H_ice - (1028.0_wp/910.0_wp)*(0.0_wp - ht%z_bed) .ge. 0.0_wp)), nfail)

    call htopo_write_init(ht, "test_htopo_out.nc", 0.0_wp)
    call htopo_write_step(ht, "test_htopo_out.nc", 0.0_wp)

    ! ========================================================================
    ! Defaults (no keys): no basins, z_bed_sd = 0
    ! ========================================================================
    call htopo_init(ht_def, nml_file, "domain_default", "topo_default", "regions_default", &
                    "NorthTest", "NHT-16KM", cnst)

    call check("default basins",    all(htopo_basins(ht_def) .eq. 0), nfail)
    call check("default z_bed_sd",  all(ht_def%ref%z_bed_sd .eq. 0.0_wp), nfail)

    ! ========================================================================
    ! Data gaps: the topography missing in a band along the x = min border
    ! ========================================================================
    allocate(gap(ht%nx,ht%ny))
    gap = .false.
    gap(1:n_band,:) = .true.
    call write_gaps_file("test_htopo_gaps.nc", ht_def, gap)
    call htopo_init(ht_gaps, nml_file, "domain_default", "topo_gaps", "regions_default", &
                    "NorthTest", "NHT-16KM", cnst)

    call check("gaps: no ice",      .not. any(ht_gaps%ref%H_ice .ne. 0.0_wp .and. gap), nfail)
    do j = 5, ht%ny-4
        if (any(ht_gaps%ref%z_bed(1:n_band,j) .ne. ht_def%ref%z_bed(n_band+1,j))) exit
    end do
    call check("gaps: nearest bed", j .gt. ht%ny-4, nfail)
    allocate(z_srf_exp(ht%nx,ht%ny))
    z_srf_exp = max(ht_gaps%ref%z_bed + ht_gaps%ref%H_ice, (1.0_wp - 910.0_wp/1028.0_wp)*ht_gaps%ref%H_ice)
    call check("gaps: z_srf",       .not. any(gap .and. abs(ht_gaps%ref%z_srf - z_srf_exp) .gt. 1e-3_wp), nfail)
    call check("gaps: valid kept",  .not. any(.not. gap .and. ht_gaps%ref%z_bed .ne. ht_def%ref%z_bed), nfail)

    call delete_file("test_htopo_gaps.nc")
    call delete_file("test_htopo_custom.nc")
    call delete_file("test_htopo_out.nc")

    write(*,*)
    if (nfail .gt. 0) then
        write(*,"(a,i0,a)") "test_htopo: ", nfail, " check(s) FAILED"
        stop 1
    end if
    write(*,"(a)") "test_htopo: all checks passed"

contains

    function layer(reg, name) result(values)
        type(regions_class), intent(in) :: reg
        character(len=*),    intent(in) :: name
        integer, allocatable :: values(:,:)
        values = reg%layers(regions_find(reg, name))%values
    end function layer

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

    subroutine write_custom_basins(filename)
        ! A custom basin field (id 3 everywhere) on the 16 km NorthTest grid,
        ! without flag attributes or grid_name
        character(len=*), intent(in) :: filename
        integer :: ids(551,551)
        ids = 3
        call nc_create(filename)
        call nc_write_dim(filename, "xc", x=1.0_dp, dx=1.0_dp, nx=551)
        call nc_write_dim(filename, "yc", x=1.0_dp, dx=1.0_dp, nx=551)
        call nc_write(filename, "my_basins", ids, dim1="xc", dim2="yc")
    end subroutine write_custom_basins

    subroutine write_gaps_file(filename, ht, gap)
        ! A FesmData-like topography file on the hub grid with all fields
        ! missing where gap.
        character(len=*),  intent(in) :: filename
        type(htopo_class), intent(in) :: ht
        logical,           intent(in) :: gap(:,:)

        call nc_create(filename)
        call nc_write_attr(filename, "grid_name", "NHT-16KM")
        call nc_write_dim(filename, "xc", x=ht%grid%G%x, units="km")
        call nc_write_dim(filename, "yc", x=ht%grid%G%y, units="km")
        call nc_write(filename, "z_bed", merge(mv, ht%ref%z_bed, gap), dim1="xc", dim2="yc", missing_value=mv)
        call nc_write(filename, "H_ice", merge(mv, ht%ref%H_ice, gap), dim1="xc", dim2="yc", missing_value=mv)
        call nc_write(filename, "z_srf", merge(mv, ht%ref%z_srf, gap), dim1="xc", dim2="yc", missing_value=mv)

    end subroutine write_gaps_file

    subroutine write_griddes(fldr, name)
        ! cdo grid description of the 16 km NorthTest grid, as FesmData writes it
        character(len=*), intent(in) :: fldr, name
        integer :: u
        call execute_command_line("mkdir -p "//trim(fldr))
        open(newunit=u, file=trim(fldr)//"/grid_"//trim(name)//".txt", status="replace", action="write")
        write(u,"(a)") "gridtype = projection"
        write(u,"(a)") "gridsize = 303601"
        write(u,"(a)") "xsize    = 551"
        write(u,"(a)") "ysize    = 551"
        write(u,"(a)") "xname    = xc"
        write(u,"(a)") "xunits   = km"
        write(u,"(a)") "yname    = yc"
        write(u,"(a)") "yunits   = km"
        write(u,"(a)") "xfirst   = -4900.0"
        write(u,"(a)") "xinc     = 16.0"
        write(u,"(a)") "yfirst   = -5400.0"
        write(u,"(a)") "yinc     = 16.0"
        write(u,"(a)") "grid_mapping = crs"
        write(u,"(a)") "grid_mapping_name = polar_stereographic"
        write(u,"(a)") "straight_vertical_longitude_from_pole = -45.0"
        write(u,"(a)") "latitude_of_projection_origin = 90.0"
        write(u,"(a)") "standard_parallel = 70.0"
        write(u,"(a)") "false_easting = 0.0"
        write(u,"(a)") "false_northing = 0.0"
        write(u,"(a)") "semi_major_axis = 6378137.0"
        write(u,"(a)") "inverse_flattening = 298.25722356"
        close(u)
    end subroutine write_griddes

end program test_htopo
