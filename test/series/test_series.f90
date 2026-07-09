program test_series
    ! Exercise the series (tabulated time-series reader) module:
    !   * ASCII 1-D and 2-D (multi-channel) tables;
    !   * netCDF 1-D and 2-D variables, with an optional _sd sigma companion;
    !   * clamped-linear interpolation at, between, and outside the nodes.
    ! All input files are written by the test itself, then removed.

    use precision, only: wp
    use ncio
    use series

    implicit none

    real(wp), parameter :: tol = real(1e-4,wp)
    integer :: nfail

    character(len=*), parameter :: f_asc1 = "test_series_1d.dat"
    character(len=*), parameter :: f_asc2 = "test_series_2d.dat"
    character(len=*), parameter :: f_nc   = "test_series.nc"

    nfail = 0

    write(*,*) "==================================================="
    write(*,*) " series test suite"
    write(*,*) "==================================================="
    write(*,"(a)") ""

    call test_ascii_1d()
    call test_ascii_2d()
    call test_netcdf_1d()
    call test_netcdf_2d()
    call test_sigma()

    write(*,"(a)") ""
    if (nfail == 0) then
        write(*,"(a)") " All series tests PASSED."
    else
        write(*,"(a,i0,a)") " series tests FAILED: ", nfail, " case(s)."
    end if

    call remove_file(f_asc1)
    call remove_file(f_asc2)
    call remove_file(f_nc)

    if (nfail > 0) stop 1

contains

    subroutine test_ascii_1d()
        ! time=[0,1000,2000], var=[0,10,20]
        type(series_class) :: s
        integer :: u

        open(newunit=u,file=f_asc1,status="replace",action="write")
        write(u,"(a)") "# comment line"
        write(u,"(a)") "   0.0    0.0"
        write(u,"(a)") ""
        write(u,"(a)") "1000.0   10.0"
        write(u,"(a)") "2000.0   20.0"
        close(u)

        call series_load(s,f_asc1)

        call check_int("ascii-1d  nc == 1", s%nc, 1)
        call check_int("ascii-1d  nt == 3", s%nt, 3)
        call check("ascii-1d  interp(500)  == 5 ", series_interp1(s, 500.0_wp),   5.0_wp)
        call check("ascii-1d  interp(1500) == 15", series_interp1(s,1500.0_wp),  15.0_wp)
        call check("ascii-1d  clamp(-100)  == 0 ", series_interp1(s,-100.0_wp),   0.0_wp)
        call check("ascii-1d  clamp(3000)  == 20", series_interp1(s,3000.0_wp),  20.0_wp)
    end subroutine test_ascii_1d

    subroutine test_ascii_2d()
        ! time=[0,1000], v1=[0,10], v2=[100,200]
        type(series_class) :: s
        real(wp) :: v(2)
        integer :: u

        open(newunit=u,file=f_asc2,status="replace",action="write")
        write(u,"(a)") "   0.0     0.0   100.0"
        write(u,"(a)") "1000.0    10.0   200.0"
        close(u)

        call series_load(s,f_asc2)

        call check_int("ascii-2d  nc == 2", s%nc, 2)
        v = series_interp(s, 500.0_wp)
        call check("ascii-2d  interp(500) ch1 == 5  ", v(1),   5.0_wp)
        call check("ascii-2d  interp(500) ch2 == 150", v(2), 150.0_wp)
    end subroutine test_ascii_2d

    subroutine test_netcdf_1d()
        type(series_class) :: s
        real(wp) :: time(3), co2(3)

        time = [0.0_wp, 1000.0_wp, 2000.0_wp]
        co2  = [280.0_wp, 380.0_wp, 480.0_wp]

        call nc_create(f_nc)
        call nc_write_dim(f_nc,"time",x=time,units="yr")
        call nc_write(f_nc,"co2",co2,dim1="time")

        call series_load(s,f_nc,varname="co2",time_name="time")

        call check_int("netcdf-1d nc == 1", s%nc, 1)
        call check("netcdf-1d interp(500)  == 330", series_interp1(s, 500.0_wp), 330.0_wp)
        call check("netcdf-1d interp(1500) == 430", series_interp1(s,1500.0_wp), 430.0_wp)
    end subroutine test_netcdf_1d

    subroutine test_netcdf_2d()
        ! 2-D variable (time,month): month m holds value time-slope + m
        type(series_class) :: s
        real(wp) :: time(2), month(3), dTmon(2,3)
        real(wp) :: v(3)
        integer  :: m

        time  = [0.0_wp, 1000.0_wp]
        month = [1.0_wp, 2.0_wp, 3.0_wp]
        do m = 1, 3
            dTmon(1,m) = real(m,wp)              ! t=0
            dTmon(2,m) = real(m,wp) + 10.0_wp    ! t=1000
        end do

        call nc_create(f_nc)
        call nc_write_dim(f_nc,"time", x=time, units="yr")
        call nc_write_dim(f_nc,"month",x=month)
        call nc_write(f_nc,"dTmon",dTmon,dim1="time",dim2="month")

        call series_load(s,f_nc,varname="dTmon",time_name="time")

        call check_int("netcdf-2d nc == 3", s%nc, 3)
        v = series_interp(s, 500.0_wp)   ! halfway: m + 5
        call check("netcdf-2d interp(500) m1 == 6", v(1), 6.0_wp)
        call check("netcdf-2d interp(500) m3 == 8", v(3), 8.0_wp)
    end subroutine test_netcdf_2d

    subroutine test_sigma()
        ! co2 with a companion co2_sd; sigma must interpolate, and be zero when absent.
        type(series_class) :: s
        real(wp) :: time(2), co2(2), co2_sd(2)
        real(wp) :: sg(1)

        time   = [0.0_wp, 1000.0_wp]
        co2    = [280.0_wp, 380.0_wp]
        co2_sd = [2.0_wp, 6.0_wp]

        call nc_create(f_nc)
        call nc_write_dim(f_nc,"time",x=time,units="yr")
        call nc_write(f_nc,"co2",   co2,   dim1="time")
        call nc_write(f_nc,"co2_sd",co2_sd,dim1="time")

        ! With sigma companion present
        call series_load(s,f_nc,varname="co2",time_name="time",sigma_name="co2_sd")
        call check_logical("sigma     allocated when _sd present", allocated(s%sigma))
        sg = series_interp_sig(s, 500.0_wp)
        call check("sigma     interp(500) == 4", sg(1), 4.0_wp)

        ! Without requesting sigma -> zeros
        call series_load(s,f_nc,varname="co2",time_name="time")
        call check_logical("sigma     absent when not requested ", .not. allocated(s%sigma))
        sg = series_interp_sig(s, 500.0_wp)
        call check("sigma     zero when absent    ", sg(1), 0.0_wp)
    end subroutine test_sigma

    ! ---- helpers ----

    subroutine check(name, got, expected)
        character(len=*), intent(in) :: name
        real(wp),         intent(in) :: got, expected
        logical :: ok
        character(len=4) :: tag
        ok = abs(got-expected) .le. tol
        tag = "PASS"
        if (.not. ok) then
            tag = "FAIL"
            nfail = nfail + 1
        end if
        write(*,"(a,a,a,f10.4,a,f10.4)") "  ["//tag//"] ", name, "  got=", got, " exp=", expected
    end subroutine check

    subroutine check_int(name, got, expected)
        character(len=*), intent(in) :: name
        integer,          intent(in) :: got, expected
        character(len=4) :: tag
        tag = "PASS"
        if (got /= expected) then
            tag = "FAIL"
            nfail = nfail + 1
        end if
        write(*,"(a,a,a,i0,a,i0)") "  ["//tag//"] ", name, "  got=", got, " exp=", expected
    end subroutine check_int

    subroutine check_logical(name, ok)
        character(len=*), intent(in) :: name
        logical,          intent(in) :: ok
        character(len=4) :: tag
        tag = "PASS"
        if (.not. ok) then
            tag = "FAIL"
            nfail = nfail + 1
        end if
        write(*,"(a,a)") "  ["//tag//"] ", name
    end subroutine check_logical

    subroutine remove_file(fname)
        character(len=*), intent(in) :: fname
        integer :: u, ios
        open(newunit=u, file=fname, status="old", iostat=ios)
        if (ios == 0) close(u, status="delete")
    end subroutine remove_file

end program test_series
