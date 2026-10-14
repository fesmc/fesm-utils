program demo_zremap_grl
    ! zremap demo for the documentation: present-day Greenland SMB (MAR v3.11,
    ! 1961-1990) remapped to the surfaces of the ICE-6G_C deglaciation, i.e. the
    ! MAR reference surface plus the ICE-6G_C surface anomaly (each time slice
    ! minus present day).
    !
    ! The profiles are built from the present-day ice points only (ice fraction
    ! >= f_ice_min); f_ref is the MAR SMB on land (z_srf > 0), so that the
    ! fallback is defined on ice-free land. Output: the zremap tables and fields
    ! (zremap_write), the present-day ice mask, and per time slice z_now, the
    ! remapped f_now, dz_ext, the constant-gradient remap
    ! f_grad = f_ref + dfdz_lo*(z_now-z_ref) and the ice mask.
    !
    ! Run from the fesm-utils folder, next to ice_data:
    !   make demo-zremap && ./bin/demo_zremap_grl.x docs/figures/zremap/zremap_grl.nml

    use precision, only: wp, dp
    use constants, only: mv
    use nml
    use ncio
    use zremap

    implicit none

    character(len=512) :: file_par, file_smb, file_topo, file_out
    real(wp) :: f_ice_min, dx
    integer  :: nx, ny, nt, k, k_pd

    real(dp), allocatable :: time(:)
    real(wp), allocatable :: smb(:,:), z_srf(:,:), f_ice(:,:)
    real(wp), allocatable :: zs(:,:), zs_pd(:,:), f_ice_now(:,:)
    real(wp), allocatable :: f_ref(:,:), z_ref(:,:), z_now(:,:), f_now(:,:), dz_ext(:,:), f_grad(:,:)
    logical,  allocatable :: mask_ref(:,:), mask_now(:,:)
    integer,  allocatable :: imask(:,:)

    type(zremap_class) :: zr

    if (command_argument_count() /= 1) then
        write(*,*) "usage: demo_zremap_grl.x <parameter file>"
        error stop 1
    end if
    call get_command_argument(1, file_par)

    call nml_read(file_par, "demo", "file_smb",  file_smb)
    call nml_read(file_par, "demo", "file_topo", file_topo)
    call nml_read(file_par, "demo", "file_out",  file_out)
    call nml_read(file_par, "demo", "f_ice_min", f_ice_min)
    call nml_read(file_par, "demo", "dx",        dx)

    nx = nc_size(file_smb, "xc")
    ny = nc_size(file_smb, "yc")
    nt = nc_size(file_topo, "time")

    allocate(smb(nx,ny), z_srf(nx,ny), f_ice(nx,ny))
    allocate(zs(nx,ny), zs_pd(nx,ny), f_ice_now(nx,ny))
    allocate(f_ref(nx,ny), z_ref(nx,ny), z_now(nx,ny), f_now(nx,ny), dz_ext(nx,ny), f_grad(nx,ny))
    allocate(mask_ref(nx,ny), mask_now(nx,ny), imask(nx,ny))
    allocate(time(nt))

    ! Present-day reference (MAR)
    call nc_read(file_smb, "smb",   smb)
    call nc_read(file_smb, "z_srf", z_srf)
    call nc_read(file_smb, "mask",  f_ice)

    z_ref    = z_srf
    f_ref    = mv
    where (z_srf > 0.0_wp) f_ref = smb
    mask_ref = f_ice >= f_ice_min .and. z_srf > 0.0_wp

    ! Profiles from the present-day reference
    call zremap_init(zr, file_par, "zremap", nx, ny, dx)
    call zremap_update(zr, f_ref, z_ref, mask_ref)

    call zremap_write(zr, file_out)
    imask = merge(1, 0, mask_ref)
    call nc_write(file_out, "mask_ref", imask, dim1="xc", dim2="yc", long_name="Present-day ice")

    ! Remap to each ICE-6G_C surface: present-day reference + ICE-6G_C anomaly
    call nc_read(file_topo, "time", time)
    k_pd = minloc(abs(time), 1)
    call nc_read(file_topo, "zs", zs_pd, start=[1,1,k_pd], count=[nx,ny,1])
    call nc_write_dim(file_out, "time", x=time, units="kyr", unlimited=.true.)

    do k = 1, nt
        call nc_read(file_topo, "zs",     zs,        start=[1,1,k], count=[nx,ny,1])
        call nc_read(file_topo, "sftgif", f_ice_now, start=[1,1,k], count=[nx,ny,1])

        z_now    = z_ref + (zs - zs_pd)
        mask_now = f_ice_now >= f_ice_min

        call zremap_map(zr, f_now, z_now, dz_ext)

        f_grad = mv
        where (f_ref /= mv) f_grad = f_ref + zr%par%dfdz_lo*(z_now - z_ref)

        call nc_write(file_out, "z_now", z_now, dim1="xc", dim2="yc", dim3="time", &
                      start=[1,1,k], count=[nx,ny,1], missing_value=mv, &
                      units="m", long_name="Surface elevation")
        call nc_write(file_out, "f_now", f_now, dim1="xc", dim2="yc", dim3="time", &
                      start=[1,1,k], count=[nx,ny,1], missing_value=mv, &
                      units="mm a-1", long_name="SMB remapped to the surface")
        call nc_write(file_out, "f_grad", f_grad, dim1="xc", dim2="yc", dim3="time", &
                      start=[1,1,k], count=[nx,ny,1], missing_value=mv, &
                      units="mm a-1", long_name="SMB with a constant gradient dfdz_lo")
        call nc_write(file_out, "dz_ext", dz_ext, dim1="xc", dim2="yc", dim3="time", &
                      start=[1,1,k], count=[nx,ny,1], missing_value=mv, &
                      units="m", long_name="Elevation distance beyond the data range")
        imask = merge(1, 0, mask_now)
        call nc_write(file_out, "mask_now", imask, dim1="xc", dim2="yc", dim3="time", &
                      start=[1,1,k], count=[nx,ny,1], long_name="Ice")
    end do
    write(*,"(a,i0,a)") "remapped ", nt, " ICE-6G_C time slices"

    write(*,"(a)") "written: "//trim(file_out)

end program demo_zremap_grl
