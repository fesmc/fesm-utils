module htopo
    ! Hi-resolution geometry hub of a domain (multigrid models such as yelmox).
    !
    ! The hub sits *above* the physics modules: its grid (grid_name) is the
    ! finest resolution in the setup, and it is the reference geometry that a
    ! coupler remaps *from* when a coarser module needs z_bed/H_ice/z_srf/masks.
    !
    ! Three kinds of field live here:
    !   * ref (topodata_class) -- the reference geometry z_bed, H_ice, z_srf
    !     and the bed roughness z_bed_sd, from a FesmData topography product;
    !   * reg (regions_class) -- regions, zones and basins from FesmData, and
    !     the named masks of its namelist group;
    !   * z_bed, H_ice, z_srf, f_grnd, z_sl -- the current geometry, refreshed
    !     each step from the models. On a hub finer than the ice sheet it is the
    !     reference plus the models' anomalies (htopo_update).
    ! The masks a physics module needs (e.g. where ice is allowed) are made by
    ! that module from a regions_class on its own grid; the hub's reg serves
    ! output and the basins handed to the coupler (htopo_basins).
    !
    ! Namelist groups (a domain of a multi-domain setup passes its own):
    !   group          basins (the basin set of htopo_basins)
    !   group_topo     topodata (path, vars, remap)
    !   group_regions  regions (path_regions, path_basins, basin_sets, masks, mask_*)
    ! {domain}/{grid_name} in the paths resolve to the domain and the hub grid.

    use precision
    use constants,      only : mv
    use nml
    use ncio
    use coordinates,    only : grid_class
    use grid_cdo,       only : grid_cdo_read_desc
    use interp2D,       only : fill_nearest
    use phys_constants, only : phys_const_class, phys_const_require, phys_const_get
    use topodata,       only : topodata_class, topodata_init_nml
    use regions,        only : regions_class, regions_init_nml, regions_basin_ids

    implicit none
    private

    type htopo_par_class
        character(len=256)      :: domain
        character(len=256)      :: grid_name        ! hub grid, e.g. "ANT-16KM"
        character(len=56)       :: basins           ! "<set>" or "<set>.group" of htopo_basins ("" = none)
        real(wp)                :: rho_ice          ! [kg m-3] ice density (from the domain's constants)
        real(wp)                :: rho_sw           ! [kg m-3] seawater density
    end type

    type htopo_class
        type(htopo_par_class) :: par
        type(grid_class)      :: grid           ! hub grid, from grid_<name>.txt
        integer               :: nx, ny
        ! Reference geometry and masks (static, from file).
        type(topodata_class)  :: ref
        type(regions_class)   :: reg
        ! Current geometry, refreshed from the models each step.
        real(wp), allocatable :: z_bed(:,:)     ! [m] bedrock elevation
        real(wp), allocatable :: H_ice(:,:)     ! [m] ice thickness
        real(wp), allocatable :: z_srf(:,:)     ! [m] surface elevation
        real(wp), allocatable :: f_grnd(:,:)    ! [1] grounded-ice fraction
        real(wp), allocatable :: z_sl(:,:)      ! [m] sea-surface / sea-level height
    end type

    public :: htopo_par_class, htopo_class
    public :: htopo_init, htopo_update, htopo_basins
    public :: htopo_write_init, htopo_write_step

contains

    subroutine htopo_init(htopo, filename, group, group_topo, group_regions, domain, grid_name, cnst, map_fldr)
        ! Load the hub: its grid from the disk grid table, the reference
        ! topography and the regions on that grid, and the masks derived from
        ! them.
        type(htopo_class), intent(out) :: htopo
        character(len=*),  intent(in)  :: filename       ! parameter file
        character(len=*),  intent(in)  :: group          ! hub parameters, e.g. "domain"
        character(len=*),  intent(in)  :: group_topo     ! topography, e.g. "topo"
        character(len=*),  intent(in)  :: group_regions  ! regions, e.g. "regions"
        character(len=*),  intent(in)  :: domain         ! domain name
        character(len=*),  intent(in)  :: grid_name      ! hub grid
        type(phys_const_class), intent(in) :: cnst       ! physical constants of the domain
        character(len=*),  intent(in), optional :: map_fldr

        character(len=256) :: mfldr
        integer :: k

        mfldr = "maps"
        if (present(map_fldr)) mfldr = trim(map_fldr)

        call htopo_par_load(htopo%par, filename, group, domain, grid_name)

        call phys_const_require(cnst, "htopo_init")
        call phys_const_get(cnst, "rho_ice", htopo%par%rho_ice)
        call phys_const_get(cnst, "rho_sw",  htopo%par%rho_sw)

        ! Hub grid definition (nx,ny + coordinates) from grid_<name>.txt.
        call grid_cdo_read_desc(htopo%grid, trim(htopo%par%grid_name), trim(mfldr))
        htopo%nx = htopo%grid%G%nx
        htopo%ny = htopo%grid%G%ny

        ! Reference topography and regions on the hub grid (remapped if their
        ! files are on another grid).
        call topodata_init_nml(htopo%ref, filename, group_topo, domain=htopo%par%domain, &
                               grid_name=htopo%par%grid_name, grid=htopo%grid)
        call regions_init_nml(htopo%reg, filename, group_regions, domain=htopo%par%domain, &
                              grid_name=htopo%par%grid_name, grid=htopo%grid)

        if (.not. (allocated(htopo%ref%z_bed) .and. allocated(htopo%ref%H_ice) &
                   .and. allocated(htopo%ref%z_srf))) then
            call htopo_error("htopo_init", "the topography needs z_bed, H_ice and z_srf (vars of ["// &
                             trim(group_topo)//"]).")
        end if
        if (.not. allocated(htopo%ref%z_bed_sd)) then
            allocate(htopo%ref%z_bed_sd(htopo%nx,htopo%ny))
            htopo%ref%z_bed_sd = 0.0_wp
        end if
        where (htopo%ref%z_bed_sd == mv) htopo%ref%z_bed_sd = 0.0_wp

        call htopo_fill_missing(htopo)

        ! Check par%basins against the loaded basin sets
        if (len_trim(htopo%par%basins) > 0) k = maxval(regions_basin_ids(htopo%reg, htopo%par%basins))

        ! The current geometry starts from the reference.
        htopo%z_bed = htopo%ref%z_bed
        htopo%H_ice = htopo%ref%H_ice
        htopo%z_srf = htopo%ref%z_srf
        allocate(htopo%f_grnd(htopo%nx,htopo%ny))
        allocate(htopo%z_sl(htopo%nx,htopo%ny))
        htopo%z_sl = 0.0_wp
        htopo%f_grnd = 0.0_wp
        where (calc_H_grnd(htopo%H_ice, htopo%z_bed, htopo%z_sl, &
                           htopo%par%rho_ice, htopo%par%rho_sw) >= 0.0_wp) htopo%f_grnd = 1.0_wp

    end subroutine htopo_init

    subroutine htopo_update(htopo, dz_bed, dH_ice, z_sl)
        ! Current geometry on a hub finer than the ice sheet: the hi-res reference
        ! plus the models' anomalies (on the hub grid), with ice thickness clipped
        ! at 0. Each hub cell is either fully ice-covered or ice-free, so the
        ! grounded fraction is 0 or 1 from flotation, and the surface follows from
        ! the bed, the ice and sea level.
        type(htopo_class), intent(inout) :: htopo
        real(wp),          intent(in)    :: dz_bed(:,:)   ! [m] bed displacement
        real(wp),          intent(in)    :: dH_ice(:,:)   ! [m] change in ice thickness
        real(wp),          intent(in)    :: z_sl(:,:)     ! [m] sea level

        htopo%z_bed = htopo%ref%z_bed + dz_bed
        htopo%H_ice = max(htopo%ref%H_ice + dH_ice, 0.0_wp)
        htopo%z_sl  = z_sl

        htopo%f_grnd = 0.0_wp
        where (calc_H_grnd(htopo%H_ice, htopo%z_bed, htopo%z_sl, &
                           htopo%par%rho_ice, htopo%par%rho_sw) >= 0.0_wp) htopo%f_grnd = 1.0_wp
        htopo%z_srf = calc_z_srf(htopo%H_ice, htopo%z_bed, htopo%z_sl, &
                                 htopo%par%rho_ice, htopo%par%rho_sw)

    end subroutine htopo_update

    function htopo_basins(htopo) result(basins)
        ! Basin ids of par%basins ("<set>" or "<set>.group") on the hub grid
        ! (0 = no basin; 0 everywhere without par%basins).
        type(htopo_class), intent(in) :: htopo
        integer :: basins(htopo%nx,htopo%ny)

        basins = 0
        if (len_trim(htopo%par%basins) > 0) basins = regions_basin_ids(htopo%reg, htopo%par%basins)

    end function htopo_basins

    elemental function calc_H_grnd(H_ice, z_bed, z_sl, rho_ice, rho_sw) result(H_grnd)
        ! Ice overburden relative to flotation: >= 0 grounded, < 0 floating. Above
        ! sea level, the bed's height counts too, so ice-free land is grounded.
        real(wp), intent(in) :: H_ice, z_bed, z_sl, rho_ice, rho_sw
        real(wp) :: H_grnd

        if (z_sl > z_bed) then
            H_grnd = H_ice - (rho_sw/rho_ice)*(z_sl - z_bed)
        else
            H_grnd = H_ice + (z_bed - z_sl)
        end if

    end function calc_H_grnd

    elemental function calc_z_srf(H_ice, z_bed, z_sl, rho_ice, rho_sw) result(z_srf)
        ! Surface elevation: the top of grounded ice or of floating ice in
        ! hydrostatic equilibrium, whichever is higher (sea level if ice-free ocean).
        real(wp), intent(in) :: H_ice, z_bed, z_sl, rho_ice, rho_sw
        real(wp) :: z_srf

        z_srf = max(z_bed + H_ice, z_sl + (1.0_wp - rho_ice/rho_sw)*H_ice)

    end function calc_z_srf

    subroutine htopo_fill_missing(htopo)
        ! Fill the gaps of the reference geometry (e.g. outside the coverage of
        ! the source dataset, or of the hub grid by the files): no ice, the bed
        ! from the nearest valid cell, and the surface from the bed and the ice
        ! thickness, with sea level at 0.
        type(htopo_class), intent(inout) :: htopo

        integer :: n_bed, n_ice, n_srf

        associate(z_bed => htopo%ref%z_bed, H_ice => htopo%ref%H_ice, z_srf => htopo%ref%z_srf)

        n_bed = count(z_bed == mv)
        n_ice = count(H_ice == mv)
        n_srf = count(z_srf == mv)
        if (n_bed + n_ice + n_srf == 0) return

        where (H_ice == mv) H_ice = 0.0_wp

        if (n_bed > 0) then
            if (n_bed < size(z_bed)) call fill_nearest(z_bed, mv)
            if (any(z_bed == mv)) then
                call htopo_error("htopo_fill_missing", "missing bedrock elevations could not be filled.", &
                                 "path = "//trim(htopo%ref%par%path))
            end if
        end if

        where (z_srf == mv) &
            z_srf = calc_z_srf(H_ice, z_bed, 0.0_wp, htopo%par%rho_ice, htopo%par%rho_sw)

        write(*,*) "htopo_init:: filled missing values: z_bed ", n_bed, ", H_ice ", n_ice, &
                   ", z_srf ", n_srf, " of ", size(z_bed)

        end associate

    end subroutine htopo_fill_missing

    subroutine htopo_par_load(par, filename, group, domain, grid_name)
        type(htopo_par_class), intent(out) :: par
        character(len=*),      intent(in)  :: filename, group
        character(len=*),      intent(in)  :: domain, grid_name

        par%domain    = trim(domain)
        par%grid_name = trim(grid_name)

        ! Optional: the basin set of htopo_basins (default none)
        par%basins = ""
        if (nml_has_param(filename, group, "basins")) &
            call nml_read(filename, group, "basins", par%basins)

    end subroutine htopo_par_load

    subroutine htopo_write_init(htopo, filename, time_init)
        ! Create a 2D output file on the hub grid, with the static masks.
        type(htopo_class), intent(in) :: htopo
        character(len=*),  intent(in) :: filename
        real(wp),          intent(in) :: time_init

        call nc_create(filename)
        call nc_write_dim(filename, "xc", x=htopo%grid%G%x, units="km")
        call nc_write_dim(filename, "yc", x=htopo%grid%G%y, units="km")
        call nc_write_dim(filename, "time", x=time_init, dx=1.0_wp, nx=1, &
                          units="year", unlimited=.TRUE.)

        call nc_write(filename, "region", htopo%reg%region_3, dim1="xc", dim2="yc", &
                      start=[1,1], long_name="Region code (deepest level)", units="1")
        call nc_write(filename, "zone", htopo%reg%zone, dim1="xc", dim2="yc", &
                      start=[1,1], long_name="Zone", units="1")
        call nc_write(filename, "basin", htopo_basins(htopo), dim1="xc", dim2="yc", &
                      start=[1,1], long_name="Basin ("//trim(htopo%par%basins)//")", units="1")
    end subroutine htopo_write_init

    subroutine htopo_write_step(htopo, filename, time)
        ! Append the dynamic hi-res geometry at `time`.
        type(htopo_class), intent(in) :: htopo
        character(len=*),  intent(in) :: filename
        real(wp),          intent(in) :: time

        integer :: ncid, n

        call nc_open(filename, ncid, writable=.TRUE.)
        n = nc_time_index(filename, "time", time, ncid)
        call nc_write(filename, "time", time, dim1="time", start=[n], count=[1], ncid=ncid)

        call nc_write(filename, "z_bed", htopo%z_bed, dim1="xc", dim2="yc", dim3="time", &
                      start=[1,1,n], count=[htopo%nx,htopo%ny,1], ncid=ncid, units="m", &
                      long_name="Bedrock elevation")
        call nc_write(filename, "H_ice", htopo%H_ice, dim1="xc", dim2="yc", dim3="time", &
                      start=[1,1,n], count=[htopo%nx,htopo%ny,1], ncid=ncid, units="m", &
                      long_name="Ice thickness")
        call nc_write(filename, "z_srf", htopo%z_srf, dim1="xc", dim2="yc", dim3="time", &
                      start=[1,1,n], count=[htopo%nx,htopo%ny,1], ncid=ncid, units="m", &
                      long_name="Surface elevation")
        call nc_write(filename, "f_grnd", htopo%f_grnd, dim1="xc", dim2="yc", dim3="time", &
                      start=[1,1,n], count=[htopo%nx,htopo%ny,1], ncid=ncid, units="1", &
                      long_name="Grounded-ice fraction")
        call nc_write(filename, "z_sl", htopo%z_sl, dim1="xc", dim2="yc", dim3="time", &
                      start=[1,1,n], count=[htopo%nx,htopo%ny,1], ncid=ncid, units="m", &
                      long_name="Sea-surface height")
        call nc_close(ncid)
    end subroutine htopo_write_step

    subroutine htopo_error(proc, msg, detail)
        ! Abort with a framed message.
        character(len=*), intent(in)           :: proc, msg
        character(len=*), intent(in), optional :: detail

        write(*,*) ""
        write(*,*) "htopo:: error in "//trim(proc)
        write(*,*) "    "//trim(msg)
        if (present(detail)) write(*,*) "    "//trim(detail)
        error stop 1

    end subroutine htopo_error

end module htopo
