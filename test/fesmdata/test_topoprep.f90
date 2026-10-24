program test_topoprep
    ! The preparation of the reference topography: the routines of topoprep on
    ! small hand-checked cases, and htopo_init with each key on a synthetic
    ! topography (a 1 km grid written by the test; no data files needed).
    ! Run from the fesm-utils folder.

    use precision
    use ncio
    use phys_constants, only : phys_const_class, phys_const_load
    use topoprep
    use htopo

    implicit none

    character(len=*), parameter :: nml_file = "test/fesmdata/test_topoprep.nml"
    character(len=*), parameter :: grid_name = "PREP-1KM"
    integer,  parameter :: nx = 40, ny = 40
    real(wp), parameter :: rho_ice = 910.0_wp, rho_sw = 1028.0_wp
    real(wp), parameter :: dx = 1000.0_wp       ! [m] the grid spacing

    type(htopo_class) :: ht, ht_off
    type(phys_const_class) :: cnst
    real(wp) :: z_bed0(nx,ny), H_ice0(nx,ny), z_srf0(nx,ny), z_sd0(nx,ny)
    real(wp) :: z_bed(nx,ny), H_ice(nx,ny)
    integer :: nfail, n_iter, i, j

    nfail = 0

    ! ========================================================================
    ! topoprep_remove_lakes: one case per branch (sea level 0)
    ! ========================================================================
    block
        real(wp) :: Hc(5,1), zb(5,1), zs(5,1)
        !            lake      too thick  consistent  floating   floating, too thick
        zb(:,1) = [  100.0_wp,  100.0_wp,  100.0_wp, -1000.0_wp, -100.0_wp]
        Hc(:,1) = [  500.0_wp, 1000.0_wp,  800.0_wp,   300.0_wp,  300.0_wp]
        zs(:,1) = [  900.0_wp,  900.0_wp,  900.0_wp,    30.0_wp,   30.0_wp]
        call topoprep_remove_lakes(Hc, zb, zs, 0.0_wp, rho_ice, rho_sw)
        call check("lakes: grounded gap filled",    Hc(1,1) .eq. 800.0_wp, nfail)
        call check("lakes: too thick reduced",      Hc(2,1) .eq. 800.0_wp, nfail)
        call check("lakes: consistent unchanged",   Hc(3,1) .eq. 800.0_wp, nfail)
        call check("lakes: floating gap kept",      Hc(4,1) .eq. 300.0_wp, nfail)
        call check("lakes: floating, too thick",    Hc(5,1) .eq. 130.0_wp, nfail)
    end block

    ! ========================================================================
    ! topoprep_limit_gradients
    ! ========================================================================
    ! A ramp of 0.4 m/m (dx = 1 km) is below the limit of 0.5: unchanged
    do j = 1, ny
        z_bed(:,j) = [(400.0_wp*i, i=1,nx)]
    end do
    z_bed0 = z_bed
    H_ice = 0.5_wp*z_bed
    call topoprep_limit_gradients(z_bed, H_ice, 0.5_wp, dx, n_iter)
    call check("gradients: below the limit unchanged", n_iter .eq. 0 .and. all(z_bed .eq. z_bed0), nfail)

    ! A step of 1000 m over one cell (1 m/m), halfway along x: smoothed until
    ! the gradient is below the limit in the interior (i, j = 3..n-3, where it
    ! is taken); the field stays antisymmetric about the step, the far field
    ! and the border rows are untouched, H_ice (half of z_bed here) follows
    ! the bed.
    z_bed = 0.0_wp
    z_bed(nx/2+1:,:) = 1000.0_wp
    z_bed0 = z_bed
    H_ice = 0.5_wp*z_bed
    call topoprep_limit_gradients(z_bed, H_ice, 0.5_wp, dx, n_iter)
    call check("gradients: smoothed", n_iter .ge. 1 .and. n_iter .lt. 50, nfail)
    call check("gradients: below the limit", &
               maxval(abs(z_bed(4:nx-2,3:ny-3) - z_bed(3:nx-3,3:ny-3)))/dx .lt. 0.5_wp .and. &
               any(z_bed .ne. z_bed0), nfail)
    call check("gradients: antisymmetric about the step", &
               maxval(abs(z_bed(nx/2,:) + z_bed(nx/2+1,:) - 1000.0_wp)) .lt. 1e-2_wp .and. &
               maxval(abs(z_bed(nx/2-3,:) + z_bed(nx/2+4,:) - 1000.0_wp)) .lt. 1e-2_wp, nfail)
    call check("gradients: far field and rows unchanged", &
               all(z_bed(1:5,:) .eq. 0.0_wp) .and. all(z_bed(nx-4:,:) .eq. 1000.0_wp) .and. &
               all(z_bed(:,1:2) .eq. z_bed0(:,1:2)) .and. all(z_bed(:,ny-2:) .eq. z_bed0(:,ny-2:)), nfail)
    call check("gradients: H_ice follows", maxval(abs(H_ice - 0.5_wp*z_bed)) .lt. 1e-2_wp, nfail)

    ! ========================================================================
    ! htopo_init: each key on the synthetic topography
    ! ========================================================================
    call phys_const_load(cnst, "par/phys_const_earth.nml", group="Earth")
    call write_griddes("maps", grid_name)

    ! The step field: a bed step of 1000 m over one cell, ice-free except
    !   (10,10): a lake (ice 500 m, surface 200 m above the ice),
    !   (30,10): too thick ice (300 m, surface 200 m above the bed),
    !   (10,30): a 5 m and (12,30) a 20 m patch of ice.
    z_bed0 = 0.0_wp
    z_bed0(nx/2+1:,:) = 1000.0_wp
    H_ice0 = 0.0_wp
    H_ice0(10,10) = 500.0_wp
    H_ice0(30,10) = 300.0_wp
    H_ice0(10,30) = 5.0_wp
    H_ice0(12,30) = 20.0_wp
    z_srf0 = z_bed0 + H_ice0
    z_srf0(10,10) = 700.0_wp
    z_srf0(30,10) = 1100.0_wp
    z_sd0  = 50.0_wp
    call write_topo("test_topoprep_step.nc", z_bed0, H_ice0, z_srf0, z_sd0)

    ! The ramp field of 0.4 m/m, with the same ice
    do j = 1, ny
        z_bed(:,j) = [(400.0_wp*i, i=1,nx)]
    end do
    call write_topo("test_topoprep_ramp.nc", z_bed, H_ice0, z_bed + H_ice0, z_sd0)

    ! Defaults: the reference is as in the file
    call htopo_init(ht_off, nml_file, "domain", "topo_off", "regions", "Test", grid_name, cnst)
    call check("off: reference as in the file", all(ht_off%ref%z_bed .eq. z_bed0) .and. &
               all(ht_off%ref%H_ice .eq. H_ice0) .and. all(ht_off%ref%z_srf .eq. z_srf0), nfail)

    call htopo_init(ht, nml_file, "domain", "topo_sd", "regions", "Test", grid_name, cnst)
    call check("z_bed_f_sd: z_bed + f*sd", all(ht%ref%z_bed .eq. z_bed0 - 50.0_wp) .and. &
               all(ht%ref%H_ice .eq. H_ice0) .and. all(ht%ref%z_srf .eq. z_srf0), nfail)

    call htopo_init(ht, nml_file, "domain", "topo_lakes", "regions", "Test", grid_name, cnst)
    H_ice = H_ice0
    H_ice(10,10) = 700.0_wp     ! gap filled
    H_ice(30,10) = 100.0_wp     ! reduced
    call check("remove_lakes: lake and thick ice", all(ht%ref%H_ice .eq. H_ice) .and. &
               all(ht%ref%z_bed .eq. z_bed0) .and. all(ht%ref%z_srf .eq. z_srf0), nfail)

    call htopo_init(ht, nml_file, "domain", "topo_hmin", "regions", "Test", grid_name, cnst)
    H_ice = H_ice0
    H_ice(10,30) = 0.0_wp       ! thinner than 10 m
    call check("H_ice_min: thin ice removed", all(ht%ref%H_ice .eq. H_ice) .and. &
               all(ht%ref%z_bed .eq. z_bed0), nfail)

    call htopo_init(ht, nml_file, "domain", "topo_grad", "regions", "Test", grid_name, cnst)
    call check("grad_lim_zb: bed step smoothed", &
               maxval(abs(ht%ref%z_bed(4:nx-2,3:ny-3) - ht%ref%z_bed(3:nx-3,3:ny-3)))/dx .lt. 0.5_wp .and. &
               any(ht%ref%z_bed .ne. z_bed0), nfail)
    call check("grad_lim_zb: z_srf as in the file", all(ht%ref%z_srf .eq. z_srf0), nfail)
    call check("current geometry from the prepared reference", all(ht%z_bed .eq. ht%ref%z_bed) .and. &
               all(ht%H_ice .eq. ht%ref%H_ice) .and. all(ht%z_srf .eq. z_srf0), nfail)

    call htopo_init(ht, nml_file, "domain", "topo_ramp", "regions", "Test", grid_name, cnst)
    call check("grad_lim_zb: grid spacing in m (0.4 m/m kept)", all(ht%ref%z_bed .eq. z_bed), nfail)

    call delete_file("test_topoprep_step.nc")
    call delete_file("test_topoprep_ramp.nc")
    call delete_file("maps/grid_"//grid_name//".txt")

    write(*,*)
    if (nfail .gt. 0) then
        write(*,"(a,i0,a)") "test_topoprep: ", nfail, " check(s) FAILED"
        stop 1
    end if
    write(*,"(a)") "test_topoprep: all checks passed"

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

    subroutine write_topo(filename, z_bed, H_ice, z_srf, z_bed_sd)
        ! A topography file on the hub grid
        character(len=*), intent(in) :: filename
        real(wp),         intent(in) :: z_bed(:,:), H_ice(:,:), z_srf(:,:), z_bed_sd(:,:)
        real(dp) :: x(nx), y(ny)

        x = [(-20.0_dp + (i-1)*1.0_dp, i=1,nx)]
        y = [(-20.0_dp + (j-1)*1.0_dp, j=1,ny)]
        call nc_create(filename)
        call nc_write_attr(filename, "grid_name", grid_name)
        call nc_write_dim(filename, "xc", x=x, units="km")
        call nc_write_dim(filename, "yc", x=y, units="km")
        call nc_write(filename, "z_bed",    z_bed,    dim1="xc", dim2="yc")
        call nc_write(filename, "H_ice",    H_ice,    dim1="xc", dim2="yc")
        call nc_write(filename, "z_srf",    z_srf,    dim1="xc", dim2="yc")
        call nc_write(filename, "z_bed_sd", z_bed_sd, dim1="xc", dim2="yc")
    end subroutine write_topo

    subroutine write_griddes(fldr, name)
        ! cdo grid description of the 1 km test grid, as FesmData writes it
        character(len=*), intent(in) :: fldr, name
        integer :: u
        call execute_command_line("mkdir -p "//trim(fldr))
        open(newunit=u, file=trim(fldr)//"/grid_"//trim(name)//".txt", status="replace", action="write")
        write(u,"(a)") "gridtype = projection"
        write(u,"(a,i0)") "gridsize = ", nx*ny
        write(u,"(a,i0)") "xsize    = ", nx
        write(u,"(a,i0)") "ysize    = ", ny
        write(u,"(a)") "xname    = xc"
        write(u,"(a)") "xunits   = km"
        write(u,"(a)") "yname    = yc"
        write(u,"(a)") "yunits   = km"
        write(u,"(a)") "xfirst   = -20.0"
        write(u,"(a)") "xinc     = 1.0"
        write(u,"(a)") "yfirst   = -20.0"
        write(u,"(a)") "yinc     = 1.0"
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

end program test_topoprep
