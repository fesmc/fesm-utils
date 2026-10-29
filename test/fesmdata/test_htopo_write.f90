program test_htopo_write
    ! The precisions of the fields written by htopo_write_step, from the caller:
    !   * without `prec`, the defaults (1e-3 absolute for z_bed, H_ice, z_srf, z_sl,
    !     1e-4 for f_grnd): attributes as before, and a file identical in content to
    !     one written with the default htopo_write_prec_t();
    !   * with `prec`, each field takes its own prec_abs/prec_rel (0 = lossless).
    ! Synthetic fields (no files, no hub init); the files are written by the test and removed.
    ! Run from the fesm-utils folder.

    use precision
    use ncio
    use htopo

    implicit none

    integer, parameter :: nx = 12, ny = 9
    character(len=*), parameter :: f_def = "test_htopo_write_def.nc"
    character(len=*), parameter :: f_dft = "test_htopo_write_dft.nc"
    character(len=*), parameter :: f_cus = "test_htopo_write_cus.nc"
    character(len=8), parameter :: vars(5) = [character(len=8) :: "z_bed", "H_ice", "z_srf", "f_grnd", "z_sl"]
    type(htopo_class) :: ht
    type(htopo_write_prec_t) :: prec
    real(wp) :: a(nx,ny), b(nx,ny)
    integer :: nfail, k, i, j

    nfail = 0

    ht%nx = nx; ht%ny = ny
    allocate(ht%z_bed(nx,ny), ht%H_ice(nx,ny), ht%z_srf(nx,ny), ht%f_grnd(nx,ny), ht%z_sl(nx,ny))
    do j = 1, ny
        do i = 1, nx
            ht%z_bed(i,j)  = -1234.5678901_wp + 37.123457_wp*i - 11.98765_wp*j
            ht%H_ice(i,j)  = 2345.6789012_wp + 3.14159265_wp*i + 2.718281_wp*j
            ht%z_srf(i,j)  = ht%z_bed(i,j) + ht%H_ice(i,j)
            ht%f_grnd(i,j) = 0.123456789_wp*i/nx + 0.0123456_wp*j/ny
            ht%z_sl(i,j)   = -0.123456789_wp*i + 0.5432109_wp*j
        end do
    end do

    ! --- defaults: no prec / the default type -------------------------------------
    call new_file(f_def)
    call htopo_write_step(ht, f_def, 0.0_dp)
    call new_file(f_dft)
    call htopo_write_step(ht, f_dft, 0.0_dp, prec=htopo_write_prec_t())

    do k = 1, 5
        call expect_attr(f_def, trim(vars(k)), "precision_abs", merge(1.0e-4_dp, 1.0e-3_dp, vars(k) == "f_grnd"))
        call expect_attr(f_def, trim(vars(k)), "precision_rel", 0.0_dp)
    end do
    call check("defaults: rounded within the default precisions", within(f_def, 1.0e-3_dp, 1.0e-4_dp, 0.0_dp), nfail)
    call check("defaults: z_bed is rounded (not stored as given)", &
               nc_differs(f_def, "z_bed", ht%z_bed), nfail)
    ! the same values with and without the (default) argument
    do k = 1, 5
        call nc_read(f_def, trim(vars(k)), a, start=[1,1,1], count=[nx,ny,1])
        call nc_read(f_dft, trim(vars(k)), b, start=[1,1,1], count=[nx,ny,1])
        call check("defaults: identical with default prec, "//trim(vars(k)), all(a .eq. b), nfail)
    end do

    ! --- from the caller ------------------------------------------------------------
    prec%z_bed  = htopo_prec_t(prec_abs=1.0e-1_wp, prec_rel=0.0_wp)       ! coarser
    prec%H_ice  = htopo_prec_t(prec_abs=0.0_wp,    prec_rel=1.0e-3_wp)    ! relative
    prec%f_grnd = htopo_prec_t(prec_abs=0.0_wp,    prec_rel=0.0_wp)       ! lossless
    ! z_srf, z_sl keep the defaults
    call new_file(f_cus)
    call htopo_write_step(ht, f_cus, 0.0_dp, prec=prec)

    call expect_attr(f_cus, "z_bed",  "precision_abs", 1.0e-1_dp)
    call expect_attr(f_cus, "z_bed",  "precision_rel", 0.0_dp)
    call expect_attr(f_cus, "H_ice",  "precision_abs", 0.0_dp)
    call expect_attr(f_cus, "H_ice",  "precision_rel", 1.0e-3_dp)
    call expect_attr(f_cus, "f_grnd", "precision_abs", 0.0_dp)
    call expect_attr(f_cus, "f_grnd", "precision_rel", 0.0_dp)
    call expect_attr(f_cus, "z_srf",  "precision_abs", 1.0e-3_dp)
    call expect_attr(f_cus, "z_sl",   "precision_abs", 1.0e-3_dp)

    call nc_read(f_cus, "z_bed", a, start=[1,1,1], count=[nx,ny,1])
    call check("caller: z_bed within 1e-1", maxval(abs(a - ht%z_bed)) .le. 1.0e-1_wp, nfail)
    call check("caller: z_bed coarser than the default", &
               maxval(abs(a - ht%z_bed)) .gt. 1.0e-3_wp, nfail)
    call nc_read(f_cus, "H_ice", a, start=[1,1,1], count=[nx,ny,1])
    call check("caller: H_ice within 1e-3 relative", all(abs(a - ht%H_ice) .le. 1.0e-3_wp*abs(ht%H_ice)), nfail)
    call nc_read(f_cus, "f_grnd", a, start=[1,1,1], count=[nx,ny,1])
    call check("caller: f_grnd lossless", all(a .eq. ht%f_grnd), nfail)
    call nc_read(f_cus, "z_srf", a, start=[1,1,1], count=[nx,ny,1])
    call check("caller: z_srf default", maxval(abs(a - ht%z_srf)) .le. 1.0e-3_wp, nfail)

    call delete_file(f_def); call delete_file(f_dft); call delete_file(f_cus)

    write(*,"(a)") ""
    if (nfail .eq. 0) then
        write(*,*) "All htopo_write_step precision tests passed."
    else
        write(*,*) "htopo_write_step precision tests FAILED: ", nfail
        stop 1
    end if

contains

    subroutine new_file(fnm)
        character(len=*), intent(in) :: fnm
        call nc_create(fnm)
        call nc_write_dim(fnm, "xc", x=1.0_dp, dx=1.0_dp, nx=nx, units="km")
        call nc_write_dim(fnm, "yc", x=1.0_dp, dx=1.0_dp, nx=ny, units="km")
        call nc_write_dim(fnm, "time", x=0.0_dp, dx=1.0_dp, nx=1, units="year", unlimited=.TRUE.)
    end subroutine new_file

    ! Every field within its default precision of the data given
    logical function within(fnm, pabs, pabs_f, prel) result(ok)
        character(len=*), intent(in) :: fnm
        real(dp),         intent(in) :: pabs, pabs_f, prel
        ok = .true.
        call nc_read(fnm, "z_bed", a, start=[1,1,1], count=[nx,ny,1]);  ok = ok .and. all(abs(a - ht%z_bed)  .le. real(pabs,wp))
        call nc_read(fnm, "H_ice", a, start=[1,1,1], count=[nx,ny,1]);  ok = ok .and. all(abs(a - ht%H_ice)  .le. real(pabs,wp))
        call nc_read(fnm, "z_srf", a, start=[1,1,1], count=[nx,ny,1]);  ok = ok .and. all(abs(a - ht%z_srf)  .le. real(pabs,wp))
        call nc_read(fnm, "z_sl",  a, start=[1,1,1], count=[nx,ny,1]);  ok = ok .and. all(abs(a - ht%z_sl)   .le. real(pabs,wp))
        call nc_read(fnm, "f_grnd",a, start=[1,1,1], count=[nx,ny,1]);  ok = ok .and. all(abs(a - ht%f_grnd) .le. real(pabs_f,wp))
    end function within

    logical function nc_differs(fnm, varname, x) result(d)
        character(len=*), intent(in) :: fnm, varname
        real(wp),         intent(in) :: x(:,:)
        real(wp) :: r(nx,ny)
        call nc_read(fnm, varname, r, start=[1,1,1], count=[nx,ny,1])
        d = any(r .ne. x)
    end function nc_differs

    ! precision_abs / precision_rel present with the value given (written as real(4)), absent if 0
    subroutine expect_attr(fnm, varname, attname, val)
        character(len=*), intent(in) :: fnm, varname, attname
        real(dp),         intent(in) :: val
        double precision :: v
        logical :: ok
        if (val .eq. 0.0_dp) then
            ok = .not. nc_exists_attr(fnm, varname, attname)
        else if (.not. nc_exists_attr(fnm, varname, attname)) then
            ok = .false.
        else
            call nc_read_attr(fnm, varname, attname, v)
            ok = abs(v - val) .le. 1.0e-6_dp*val
        end if
        call check("attribute "//trim(varname)//":"//trim(attname), ok, nfail)
    end subroutine expect_attr

    subroutine check(label, ok, nfail)
        character(len=*), intent(in)    :: label
        logical,          intent(in)    :: ok
        integer,          intent(inout) :: nfail
        if (.not. ok) then
            write(*,*) "FAIL: ", trim(label)
            nfail = nfail + 1
        end if
    end subroutine check

    subroutine delete_file(fnm)
        character(len=*), intent(in) :: fnm
        integer :: u
        open(newunit=u, file=fnm, status="old")
        close(u, status="delete")
    end subroutine delete_file

end program test_htopo_write
