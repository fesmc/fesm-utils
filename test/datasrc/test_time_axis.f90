program test_time_axis
    ! The time axis of the file (varslice V2): varslice_time_decode for every
    ! calendar and `since` unit, the BP/CE forms and the errors (stat); the
    ! axis read by varslice_init_nml when time(1) < 0 or the time key is
    ! absent: a regular file gives the same vs%time as the legacy time key,
    ! multifile concatenation with different units per file, an irregular
    ! ICE-6G-like axis in ka BP served with time_ref, a monthly axis (nsub
    ! inferred), time_units = "clock", and the mode from the file.

    use precision
    use ncio
    use varslice, only : varslice_class, varslice_init_nml, varslice_update, varslice_end, &
                         varslice_time_decode, varslice_nsub

    implicit none

    character(len=*), parameter :: fnml = "test_ds_taxis.nml"
    integer, parameter :: nx = 3, ny = 2
    real(dp), parameter :: tol = 1e-9_dp
    type(varslice_class) :: vs, vl
    real(dp), allocatable :: y(:)
    real(dp) :: t9(9), tm(24)
    integer :: nfail, stat, k
    character(len=*), parameter :: ago_grp(3)   = ["ago_a", "ago_b", "ago_c"]
    character(len=*), parameter :: ago_units(3) = ["years ago           ", "years BP            ", &
                                                   "years before present"]

    nfail = 0

    ! === Decoding ===============================================================
    call dec("days since 1850-01-01", "standard",            [0.0_dp, 31.0_dp, 365.0_dp], &
             [1850.0_dp, 1850.0_dp + 31.0_dp/365.0_dp, 1851.0_dp])
    call dec("days since 1850-1-1 12:00:00", "gregorian",    [0.0_dp], [1850.0_dp + 0.5_dp/365.0_dp])
    call dec("days since 2000-01-01", "proleptic_gregorian", [366.0_dp, 366.0_dp+59.0_dp], &
             [2001.0_dp, 2001.0_dp + 59.0_dp/365.0_dp])
    call dec("days since 1500-01-01", "proleptic_gregorian", [0.0_dp], [1500.0_dp])
    call dec("hours since 1900-01-01 00:00:00", "gregorian", [24.0_dp*59.0_dp], [1900.0_dp + 59.0_dp/365.0_dp])
    call dec("minutes since 2004-03-01", "standard",         [1440.0_dp*306.0_dp], [2005.0_dp])
    call dec("seconds since 1970-01-01T00:00:00Z", "standard", [86400.0_dp*365.0_dp], [1971.0_dp])
    call dec("days since 0001-01-01", "noleap",              [3650.0_dp, 3650.0_dp+73.0_dp], [11.0_dp, 11.2_dp])
    call dec("days since 2000-03-01", "365_day",             [306.0_dp], [2001.0_dp])
    call dec("days since 2001-01-01", "all_leap",            [366.0_dp, 183.0_dp], [2002.0_dp, 2001.5_dp])
    call dec("days since 2001-01-01", "366_day",             [732.0_dp], [2003.0_dp])
    call dec("days since 1850-01-01", "360_day",             [1980.0_dp], [1855.5_dp])
    call dec("months since 1850-01-01", "360_day",           [6.0_dp, 30.0_dp], [1850.5_dp, 1852.5_dp])
    call dec("years since 1950-07-01", "",                   [-21000.0_dp, 0.5_dp], [-19050.0_dp, 1950.5_dp])
    call dec("common_years since 1-1-1", "noleap",           [10.0_dp], [11.0_dp])
    call dec("yr BP", "",                                    [21000.0_dp], [-19050.0_dp])
    call dec("a BP", "",                                     [0.0_dp], [1950.0_dp])
    call dec("years BP", "",                                 [100.0_dp], [1850.0_dp])
    call dec("ka BP", "",                                    [21.0_dp], [-19050.0_dp])
    call dec("years ago", "",                                [80000.0_dp, 0.0_dp], [-78050.0_dp, 1950.0_dp])
    call dec("Years Ago", "",                                [100.0_dp], [1850.0_dp])
    call dec("yr ago", "",                                   [21000.0_dp], [-19050.0_dp])
    call dec("a ago", "",                                    [0.0_dp], [1950.0_dp])
    call dec("years before present", "",                     [100.0_dp], [1850.0_dp])
    call dec("yr before present", "",                        [21000.0_dp], [-19050.0_dp])
    call dec("ka ago", "",                                   [21.0_dp], [-19050.0_dp])
    call dec("kyr ago", "",                                  [0.5_dp], [1450.0_dp])
    call dec("kyr BP", "",                                   [0.5_dp], [1450.0_dp])
    call dec("yr CE", "",                                    [1850.0_dp], [1850.0_dp])
    call dec("years CE", "",                                 [-500.0_dp], [-500.0_dp])
    call dec("AD", "",                                       [2000.0_dp], [2000.0_dp])

    call dec_fail("yr", "")
    call dec_fail("years", "")
    call dec_fail("ka", "")
    call dec_fail("months since 1850-01-01", "noleap")
    call dec_fail("days since 1500-01-01", "standard")
    call dec_fail("days since 1850-01-01", "julian")
    call dec_fail("fortnights since 1850-01-01", "standard")
    call dec_fail("days since 1850-02-30", "standard")

    ! === Files ==================================================================
    call write_files()
    call write_nml()

    ! A regular annual axis: legacy key and file axis give the same vs%time
    call varslice_init_nml(vl, fnml, "reg_legacy", verbose=.false.)
    call varslice_init_nml(vs, fnml, "reg_file", verbose=.false.)
    call check("regular: transient from the file", trim(vs%par%mode) .eq. "transient" .and. vs%par%with_time)
    call check("regular: vs%time equal to the legacy axis", size(vs%time) .eq. size(vl%time) .and. &
               all(vs%time .eq. vl%time))
    call varslice_update(vs, [3.5_dp], method="interp")
    call varslice_update(vl, [3.5_dp], method="interp")
    call check("regular: same interpolated slice", all(vs%var .eq. vl%var))
    call varslice_end(vs)
    call varslice_end(vl)

    ! The time key absent (legacy reader): the axis of the file
    call varslice_init_nml(vs, fnml, "reg_nokey", verbose=.false.)
    call check("no time key: axis of the file", vs%par%time_from_file .and. size(vs%time) .eq. 6 .and. &
               abs(vs%time(6) - 6.0_dp) .lt. tol)
    call varslice_end(vs)

    ! Multifile: days since 1850 (noleap) then days since 1900 (360_day), with time_ref
    call varslice_init_nml(vs, fnml, "multi", verbose=.false., time_ref=1850.0_dp)
    call check("multifile: 4 + 3 steps", size(vs%time) .eq. 7 .and. all(vs%nt_files .eq. [4,3]))
    call check("multifile: decoded and shifted", abs(vs%time(1) - 0.5_dp) .lt. tol .and. &
               abs(vs%time(4) - 3.5_dp) .lt. tol .and. abs(vs%time(5) - 50.5_dp) .lt. tol .and. &
               abs(vs%time(7) - 52.5_dp) .lt. tol)
    call varslice_update(vs, [51.0_dp], method="interp")
    call check("multifile: interp across the second file", abs(vs%var(1,1,1,1) - 51.0_wp) .lt. 1e-4_wp)
    call varslice_end(vs)

    ! Irregular ICE-6G-like axis (1 kyr, then 500 yr) in ka BP, time_ref = 1950
    call varslice_init_nml(vs, fnml, "ice6g", verbose=.false., time_ref=1950.0_dp)
    t9 = -1000.0_dp*[26.0_dp, 25.0_dp, 24.0_dp, 23.0_dp, 22.0_dp, 21.0_dp, 20.5_dp, 20.0_dp, 19.5_dp]
    call check("irregular: axis in clock years", all(abs(vs%time - t9) .lt. tol))
    call check("irregular: nsub = 1", varslice_nsub(vs) .eq. 1 .and. .not. vs%par%with_time_sub)
    call varslice_update(vs, [-20250.0_dp], method="interp")
    call check("irregular: interp in a 500 yr step", abs(vs%var(2,1,1,1) - real(lin(-20250.0_dp), wp)) .lt. 1e-3_wp)
    call varslice_update(vs, [-25500.0_dp], method="interp")
    call check("irregular: interp in a 1 kyr step", abs(vs%var(2,1,1,1) - real(lin(-25500.0_dp), wp)) .lt. 1e-3_wp)
    call varslice_end(vs)

    ! "years ago", "years BP", "years before present": the same axis, model time = -value
    do k = 1, size(ago_grp)
        call varslice_init_nml(vs, fnml, trim(ago_grp(k)), verbose=.false., time_ref=1950.0_dp)
        call check("'"//trim(ago_units(k))//"': model time = -value", size(vs%time) .eq. 5 .and. &
                   all(abs(vs%time - [-80000.0_dp, -60000.0_dp, -40000.0_dp, -20000.0_dp, 0.0_dp]) .lt. tol))
        call varslice_update(vs, [-30000.0_dp], method="interp")
        call check("'"//trim(ago_units(k))//"': interp at 30000 years ago", &
                   abs(vs%var(1,1,1,1) - 25.0_wp) .lt. 1e-3_wp)
        call varslice_end(vs)
    end do

    ! Monthly axis (days since 2000-01-01, noleap, mid-month): nsub = 12
    call varslice_init_nml(vs, fnml, "monthly", verbose=.false.)
    call check("monthly: nsub inferred = 12", varslice_nsub(vs) .eq. 12 .and. vs%par%with_time_sub)
    call check("monthly: time_sub", vs%time_sub(1) .eq. 1.0_dp .and. vs%time_sub(24) .eq. 12.0_dp)
    call varslice_update(vs, [2001.0_dp], method="interp", rep=12)
    call check("monthly: rep=12 gives the year's cycle", size(vs%var,3) .eq. 12 .and. &
               abs(vs%var(1,1,3,1) - 115.0_wp) .lt. 1e-4_wp)
    call varslice_end(vs)

    ! time_units = "clock": the values as they stand (no time_ref)
    call varslice_init_nml(vs, fnml, "clock", verbose=.false., time_ref=1950.0_dp)
    call check("clock: no decoding, no shift", all(abs(vs%time - [1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp]) .lt. tol))
    call varslice_end(vs)

    ! Mode from the file: a 2-D field is static, the climatology attribute gives clim
    call varslice_init_nml(vs, fnml, "static2d", verbose=.false.)
    call check("mode auto: static", trim(vs%par%mode) .eq. "static" .and. .not. vs%par%with_time)
    call varslice_end(vs)
    call varslice_init_nml(vs, fnml, "clim", verbose=.false.)
    call check("mode auto: clim (climatology attribute)", trim(vs%par%mode) .eq. "clim" .and. &
               varslice_nsub(vs) .eq. 12)
    call varslice_update(vs)
    call check("clim: the cycle", size(vs%var,3) .eq. 12)
    call varslice_end(vs)

    if (nfail .gt. 0) then
        write(*,*) "test_time_axis: ", nfail, " FAILED"
        error stop 1
    end if
    write(*,*) "test_time_axis: all passed"

contains

    real(dp) function lin(t)
        ! The value of cell (2,1) at clock time t (linear in t)
        real(dp), intent(in) :: t
        lin = 2.0_dp + 0.001_dp*(t + 26000.0_dp)
    end function lin

    subroutine dec(units, cal, t, yref)
        character(len=*), intent(in) :: units, cal
        real(dp),         intent(in) :: t(:), yref(:)
        logical :: ok
        call varslice_time_decode(t, units, cal, y, stat=stat)
        ok = (stat .eq. 0)
        if (ok) ok = all(abs(y - yref) .lt. tol)
        call check("decode '"//units//"' ("//cal//")", ok)
        if (.not. ok .and. stat .eq. 0) write(*,*) "      ", y
    end subroutine dec

    subroutine dec_fail(units, cal)
        character(len=*), intent(in) :: units, cal
        call varslice_time_decode([0.0_dp], units, cal, y, stat=stat)
        call check("decode '"//units//"' ("//cal//") fails", stat .ne. 0)
    end subroutine dec_fail

    subroutine write_files()
        real(wp) :: v(nx,ny,24)
        real(dp) :: tt(24)
        integer  :: i, j

        ! Regular annual axis 1..6, "years since 0-1-1"
        do k = 1, 6
            v(:,:,k) = real(k, wp)
        end do
        call ncfile("test_ds_taxis_reg.nc", v(:,:,1:6), [(real(k,dp), k=1,6)], "years since 0-1-1", "")

        ! Multifile: mid-year samples, file 1 in days since 1850 (noleap),
        ! file 2 in days since 1900 (360_day); value = year - 1850 + 0.5
        tt(1:4) = [(182.5_dp + 365.0_dp*(k-1), k=1,4)]
        do k = 1, 4
            v(:,:,k) = real(k - 1, wp) + 0.5_wp
        end do
        call ncfile("test_ds_taxis_mf_1.nc", v(:,:,1:4), tt(1:4), "days since 1850-01-01", "noleap")
        tt(1:3) = [(180.0_dp + 360.0_dp*(k-1), k=1,3)]
        do k = 1, 3
            v(:,:,k) = real(50 + k - 1, wp) + 0.5_wp
        end do
        call ncfile("test_ds_taxis_mf_2.nc", v(:,:,1:3), tt(1:3), "days since 1900-01-01", "360_day")

        ! ICE-6G-like: ka BP, 26 ... 21 every 1 kyr, then every 500 yr
        tt(1:9) = [26.0_dp, 25.0_dp, 24.0_dp, 23.0_dp, 22.0_dp, 21.0_dp, 20.5_dp, 20.0_dp, 19.5_dp]
        do k = 1, 9
            do j = 1, ny; do i = 1, nx
                v(i,j,k) = real(i + 0.001_dp*(1950.0_dp - 1000.0_dp*tt(k) - 1950.0_dp + 26000.0_dp), wp)
            end do; end do
        end do
        call ncfile("test_ds_taxis_ice6g.nc", v(:,:,1:9), tt(1:9), "ka BP", "")

        ! "Years ago" axes, 80000 ... 0 (as GLOBAL-0.5DEG_IceHistory-PaleoMIST.nc): value = 40 - t/2000
        tt(1:5) = [80000.0_dp, 60000.0_dp, 40000.0_dp, 20000.0_dp, 0.0_dp]
        do k = 1, 5
            v(:,:,k) = real(40.0_dp - tt(k)/2000.0_dp, wp)          ! 0, 10, ..., 40
        end do
        do k = 1, 3
            call ncfile("test_ds_taxis_"//ago_grp(k)//".nc", v(:,:,1:5), tt(1:5), trim(ago_units(k)), "")
        end do

        ! Monthly, two years, mid-month days since 2000-01-01 (noleap);
        ! value = 100 + 12*(year-2000) + month
        do k = 1, 24
            tm(k) = 365.0_dp*((k-1)/12) + mid_month(mod(k-1,12)+1)
            v(:,:,k) = real(100 + k, wp)
        end do
        call ncfile("test_ds_taxis_mon.nc", v, tm, "days since 2000-01-01", "noleap")

        ! Clock axis: time values 1..6 with units of no meaning
        do k = 1, 6
            v(:,:,k) = real(k, wp)
        end do
        call ncfile("test_ds_taxis_clock.nc", v(:,:,1:6), [(real(k,dp), k=1,6)], "kyr", "")

        ! Static 2-D, and a climatology with the CF attribute
        call nc_create("test_ds_taxis_2d.nc")
        call nc_write_dim("test_ds_taxis_2d.nc", "xc", x=1.0_dp, dx=1.0_dp, nx=nx)
        call nc_write_dim("test_ds_taxis_2d.nc", "yc", x=1.0_dp, dx=1.0_dp, nx=ny)
        call nc_write("test_ds_taxis_2d.nc", "v", v(:,:,1), dim1="xc", dim2="yc")
        call ncfile("test_ds_taxis_clim.nc", v(:,:,1:12), [(real(k,dp)*30.0_dp-15.0_dp, k=1,12)], &
                    "days since 0001-01-01", "360_day")
        call nc_write_attr("test_ds_taxis_clim.nc", "time", "climatology", "climatology_bounds")

    end subroutine write_files

    real(dp) function mid_month(m)
        integer, intent(in) :: m
        integer, parameter :: ML(12) = [31,28,31,30,31,30,31,31,30,31,30,31]
        mid_month = real(sum(ML(1:m-1)), dp) + 0.5_dp*ML(m)
    end function mid_month

    subroutine ncfile(fname, v, t, units, cal)
        character(len=*), intent(in) :: fname, units, cal
        real(wp),         intent(in) :: v(:,:,:)
        real(dp),         intent(in) :: t(:)
        call nc_create(fname)
        call nc_write_dim(fname, "xc", x=1.0_dp, dx=1.0_dp, nx=size(v,1))
        call nc_write_dim(fname, "yc", x=1.0_dp, dx=1.0_dp, nx=size(v,2))
        if (len_trim(cal) .gt. 0) then
            call nc_write_dim(fname, "time", x=t, units=units, calendar=cal, unlimited=.false.)
        else
            call nc_write_dim(fname, "time", x=t, units=units, unlimited=.false.)
        end if
        call nc_write(fname, "v", v, dim1="xc", dim2="yc", dim3="time")
    end subroutine ncfile

    subroutine write_nml()
        integer :: u
        open(newunit=u, file=fnml, status="replace", action="write")
        call grp(u, "reg_legacy", "test_ds_taxis_reg.nc", "time = 1 1 6 1 1")
        call grp(u, "reg_file",   "test_ds_taxis_reg.nc", "time = -1 0 0 0 0")
        call grp(u, "reg_nokey",  "test_ds_taxis_reg.nc", "")
        call grp(u, "multi",      "test_ds_taxis_mf_*.nc", "time = -1 0 0 0 0")
        call grp(u, "ice6g",      "test_ds_taxis_ice6g.nc", "time = -1 0 0 0 0")
        do k = 1, 3
            call grp(u, ago_grp(k), "test_ds_taxis_"//ago_grp(k)//".nc", "time = -1 0 0 0 0")
        end do
        call grp(u, "monthly",    "test_ds_taxis_mon.nc", "time = -1 0 0 0 0")
        call grp(u, "clock",      "test_ds_taxis_clock.nc", 'time = -1 0 0 0 0'//new_line("a")// &
                                  '    time_units = "clock"')
        call grp(u, "static2d",   "test_ds_taxis_2d.nc", "time = -1 0 0 0 0")
        call grp(u, "clim",       "test_ds_taxis_clim.nc", "time = -1 0 0 0 0")
        close(u)
    end subroutine write_nml

    subroutine grp(u, name, file, extra)
        integer,          intent(in) :: u
        character(len=*), intent(in) :: name, file, extra
        write(u,"(a)") "&"//name
        write(u,"(a)") '    filename = "'//file//'"'
        write(u,"(a)") '    name     = "v"'
        write(u,"(a)") '    units    = "1" "1"'
        write(u,"(a)") '    scaling  = 1.0 0.0'
        if (len_trim(extra) .gt. 0) write(u,"(a)") "    "//extra
        write(u,"(a)") "/"
    end subroutine grp

    subroutine check(label, ok)
        character(len=*), intent(in) :: label
        logical,          intent(in) :: ok
        write(*,"(a,a)") merge("  ok   ", "  FAIL ", ok), label
        if (.not. ok) nfail = nfail + 1
    end subroutine check

end program test_time_axis
