program test_varslice_schema
    ! varslice_init_nml(..., defaults_file, defaults_group): keys not in the
    ! group come from the schema group; the arguments mode, remap, grid_src
    ! and time_ref override the group; remap = "auto" remaps only from a
    ! known, different grid.

    use precision
    use ncio
    use constants,   only : mv
    use coordinates, only : grid_class, grid_init
    use varslice,    only : varslice_class, varslice_init_nml, varslice_update, varslice_end

    implicit none

    character(len=*), parameter :: fnml = "test_ds_schema.nml", fdef = "test_ds_schema_def.nml"
    integer, parameter :: nx = 4, ny = 3
    type(varslice_class) :: vs
    type(grid_class) :: g_same, g_other
    real(wp) :: v(nx,ny)
    integer :: i, j, nfail

    nfail = 0

    do j = 1, ny; do i = 1, nx
        v(i,j) = real(10*i + j, wp)
    end do; end do
    call nc_create("test_ds_schema.nc")
    call nc_write_dim("test_ds_schema.nc", "lon", x=-10.0_dp, dx=5.0_dp, nx=nx)
    call nc_write_dim("test_ds_schema.nc", "lat", x=60.0_dp, dx=5.0_dp, nx=ny)
    call nc_write("test_ds_schema.nc", "v", v, dim1="lon", dim2="lat", units="m")
    call nc_write_attr("test_ds_schema.nc", "grid_name", "TEST-SCHEMA-SRC")
    call nc_create("test_ds_schema_tr.nc")
    call nc_write_dim("test_ds_schema_tr.nc", "xc", x=1.0_dp, dx=1.0_dp, nx=nx)
    call nc_write_dim("test_ds_schema_tr.nc", "yc", x=1.0_dp, dx=1.0_dp, nx=ny)
    call nc_write_dim("test_ds_schema_tr.nc", "time", x=0.0_dp, dx=1.0_dp, nx=2)
    call nc_write("test_ds_schema_tr.nc", "v", reshape([v, v], [nx,ny,2]), dim1="xc", dim2="yc", dim3="time")

    call grid_init(g_same, name="TEST-SCHEMA-SRC", mtype="latlon", units="degrees", &
                   x0=-10.0_dp, dx=5.0_dp, nx=nx, y0=60.0_dp, dy=5.0_dp, ny=ny)
    call grid_init(g_other, name="TEST-SCHEMA-TGT", mtype="latlon", units="degrees", &
                   x0=-7.5_dp, dx=5.0_dp, nx=nx-1, y0=62.5_dp, dy=5.0_dp, ny=ny-1)

    call write_nml()

    ! --- schema defaults: only filename, name and time in the group --------
    call varslice_init_nml(vs, fnml, "g_min", verbose=.false., defaults_file=fdef, defaults_group="tmpl")
    call check("units from the schema", trim(vs%par%units_in) .eq. "auto" .and. trim(vs%par%units_out) .eq. "auto")
    call check("scaling from the schema", vs%par%unit_scale .eq. 1.0_wp .and. vs%par%unit_offset .eq. 0.0_wp)
    call check("fill from the schema", trim(vs%par%fill) .eq. "none")
    call check("mode static from time", trim(vs%par%mode) .eq. "static" .and. .not. vs%par%with_time)
    call check("remap auto without grid: none", trim(vs%par%remap) .eq. "none")
    call varslice_update(vs)
    call check("values", vs%var(2,3,1,1) .eq. 23.0_wp)
    call varslice_end(vs)

    ! --- group keys over the schema --------------------------------------------
    call varslice_init_nml(vs, fnml, "g_keys", verbose=.false., defaults_file=fdef, defaults_group="tmpl")
    call check("scaling from the group", vs%par%unit_scale .eq. 2.0_wp .and. vs%par%unit_offset .eq. 1.0_wp)
    call check("units from the group", trim(vs%par%units_in) .eq. "m" .and. trim(vs%par%units_out) .eq. "cm")
    call check("fill const from the group", trim(vs%par%fill) .eq. "const" .and. vs%par%fill_value .eq. 7.0_wp)
    call check("explicit mode agrees", trim(vs%par%mode) .eq. "static")
    call varslice_update(vs)
    call check("scaled values", vs%var(2,3,1,1) .eq. 47.0_wp)
    call varslice_end(vs)

    ! --- remap = "auto" ------------------------------------------------------------
    call varslice_init_nml(vs, fnml, "g_min", verbose=.false., grid=g_same, &
                           defaults_file=fdef, defaults_group="tmpl")
    call check("remap auto, same grid: none", trim(vs%par%remap) .eq. "none")
    call varslice_end(vs)

    call varslice_init_nml(vs, fnml, "g_min", verbose=.false., grid=g_other, &
                           defaults_file=fdef, defaults_group="tmpl")
    call check("remap auto, other grid: bilin from grid_name", trim(vs%par%remap) .eq. "bilin" .and. &
               trim(vs%par%grid_src) .eq. "TEST-SCHEMA-SRC")
    call varslice_update(vs)
    call check("remapped shape", size(vs%var,1) .eq. nx-1 .and. size(vs%var,2) .eq. ny-1)
    call check("remapped value (bilinear midpoint)", abs(vs%var(1,1,1,1) - 0.25_wp*(11+12+21+22)) .lt. 1e-4_wp)
    call varslice_end(vs)

    ! --- argument overrides --------------------------------------------------------
    call varslice_init_nml(vs, fnml, "g_min", verbose=.false., grid=g_other, &
                           defaults_file=fdef, defaults_group="tmpl", &
                           mode="static", remap="nn", grid_src="TEST-SCHEMA-SRC", time_ref=1950.0_dp)
    call check("override remap, grid_src", trim(vs%par%remap) .eq. "nn" .and. &
               trim(vs%par%grid_src) .eq. "TEST-SCHEMA-SRC")
    call check("override time_ref", vs%par%time_ref .eq. 1950.0_dp)
    call varslice_end(vs)

    ! --- legacy path: no mode key, mode from time ---------------------------------
    call varslice_init_nml(vs, fnml, "g_legacy", verbose=.false.)
    call check("legacy: transient", trim(vs%par%mode) .eq. "transient" .and. vs%par%with_time)
    call check("legacy: remap none", trim(vs%par%remap) .eq. "none")
    call varslice_end(vs)

    if (nfail .gt. 0) then
        write(*,*) "test_varslice_schema: ", nfail, " FAILED"
        error stop 1
    end if
    write(*,*) "test_varslice_schema: all passed"

contains

    subroutine write_nml()
        integer :: u
        open(newunit=u, file=fdef, status="replace", action="write")
        write(u,"(a)") "&tmpl"
        write(u,"(a)") '    filename   = ""'
        write(u,"(a)") '    name       = ""'
        write(u,"(a)") '    mode       = "auto"'
        write(u,"(a)") '    units      = "auto" "auto"'
        write(u,"(a)") '    scaling    = 1.0 0.0'
        write(u,"(a)") '    remap      = "auto"'
        write(u,"(a)") '    grid_src   = ""'
        write(u,"(a)") '    fill       = "none"'
        write(u,"(a)") '    fill_value = 0.0'
        write(u,"(a)") '    time       = 0 0 0 0 0'
        write(u,"(a)") "/"
        close(u)

        open(newunit=u, file=fnml, status="replace", action="write")
        write(u,"(a)") "&g_min"
        write(u,"(a)") '    filename = "test_ds_schema.nc"'
        write(u,"(a)") '    name     = "v"'
        write(u,"(a)") "/"
        write(u,"(a)") "&g_keys"
        write(u,"(a)") '    filename = "test_ds_schema.nc"'
        write(u,"(a)") '    name     = "v"'
        write(u,"(a)") '    mode     = "static"'
        write(u,"(a)") '    units    = "m" "cm"'
        write(u,"(a)") '    scaling  = 2.0 1.0'
        write(u,"(a)") '    remap    = "none"'
        write(u,"(a)") '    fill     = "const"'
        write(u,"(a)") '    fill_value = 7.0'
        write(u,"(a)") "/"
        write(u,"(a)") "&g_legacy"
        write(u,"(a)") '    filename = "test_ds_schema_tr.nc"'
        write(u,"(a)") '    name     = "v"'
        write(u,"(a)") '    units    = "m" "m"'
        write(u,"(a)") '    scaling  = 1.0 0.0'
        write(u,"(a)") '    time     = 1 0 1 1 1'
        write(u,"(a)") "/"
        close(u)
    end subroutine write_nml

    subroutine check(label, ok)
        character(len=*), intent(in) :: label
        logical,          intent(in) :: ok
        write(*,"(a,a)") merge("  ok   ", "  FAIL ", ok), label
        if (.not. ok) nfail = nfail + 1
    end subroutine check

end program test_varslice_schema
