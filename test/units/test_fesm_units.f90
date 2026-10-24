program test_fesm_units
    ! units_conv (fesm_units): every row of the table and its inverse,
    ! identity bit-exact, failures with stat, the sp round trip of
    ! mW m-2 -> W m-2 -> mW m-2, and fesmdata_flux_conv / fesmdata_temp_ref
    ! against the factors of their implementation before fesm_units
    ! (bit patterns recorded with fesm-utils 1725385).

    use, intrinsic :: iso_fortran_env, only : int32
    use precision
    use ncio
    use fesm_units, only : units_conv
    use fesmdata,   only : fesmdata_flux_conv, fesmdata_temp_ref

    implicit none

    real(dp), parameter :: sy = 365.25_dp*86400.0_dp, rho = 917.0_dp
    character(len=*), parameter :: fname = "test_ds_units.nc"
    integer :: nfail, stat, i, k, nbad
    real(dp) :: s, o, s1, o1
    real(sp) :: x, y

    nfail = 0

    ! --- table rows and inverses ----------------------------------------------
    call row("kg m-2 s-1",   "kg m-2 yr-1", sy,             0.0_dp)
    call row("kg m-2 d-1",   "kg m-2 yr-1", sy/86400.0_dp,  0.0_dp)
    call row("mm s-1",       "kg m-2 yr-1", sy,             0.0_dp)
    call row("mm d-1",       "kg m-2 yr-1", 365.25_dp,      0.0_dp)
    call row("mm w.e. d-1",  "kg m-2 yr-1", 365.25_dp,      0.0_dp)
    call row("mm w.e. s-1",  "kg m-2 yr-1", sy,             0.0_dp)
    call row("m w.e. yr-1",  "kg m-2 yr-1", 1000.0_dp,      0.0_dp)
    call row("m i.e. yr-1",  "kg m-2 yr-1", rho,            0.0_dp, rho_ice=rho)
    call row("m ice yr-1",   "kg m-2 yr-1", rho,            0.0_dp, rho_ice=rho)
    call row("degC",         "K",           1.0_dp,         273.15_dp)
    call row("mW m-2",       "W m-2",       1.0e-3_dp,      0.0_dp)
    call row("km",           "m",           1000.0_dp,      0.0_dp)
    call row("m s-1",        "m yr-1",      sy,             0.0_dp)
    call row("m d-1",        "m yr-1",      365.25_dp,      0.0_dp)
    call row("%",            "1",           0.01_dp,        0.0_dp)
    call row("g kg-1",       "psu",         1.0_dp,         0.0_dp)

    ! Spellings, the run's year length, conversions within a class
    call row("kg m**-2 s**-1", "kg m-2 yr-1", 360.0_dp*86400.0_dp, 0.0_dp, sec_year=360.0_dp*86400.0_dp)
    call row("mm/d",         "mm yr-1",     365.25_dp,      0.0_dp)
    call row("mm w.e. a-1",  "mm d-1",      1.0_dp/365.25_dp, 0.0_dp)
    call row("degrees_Celsius", "kelvin",   1.0_dp,         273.15_dp)
    call row("m/s",          "m d-1",       86400.0_dp,     0.0_dp)
    call row("mW/m2",        "W/m2",        1.0e-3_dp,      0.0_dp)

    ! --- identity: bit-exact ----------------------------------------------------
    call units_conv("kg m-2 yr-1", "mm w.e./yr", s, o)
    call check("identity (alias)", s .eq. 1.0_dp .and. o .eq. 0.0_dp)
    call units_conv("furlong fortnight-1", "Furlong Fortnight-1", s, o)
    call check("identity (unknown units, same key)", s .eq. 1.0_dp .and. o .eq. 0.0_dp)
    call units_conv("mm d-1", "kg m-2 d-1", s, o)
    call check("identity (mm = kg m-2)", s .eq. 1.0_dp .and. o .eq. 0.0_dp)

    ! --- failures ----------------------------------------------------------------
    call units_conv("furlong", "m", s, o, stat=stat)
    call check("unknown units: stat /= 0", stat .ne. 0)
    call units_conv("K", "m", s, o, stat=stat)
    call check("different quantities: stat /= 0", stat .ne. 0)
    call units_conv("m ice yr-1", "kg m-2 yr-1", s, o, stat=stat)
    call check("m ice yr-1 without rho_ice: stat /= 0", stat .ne. 0)
    call units_conv("mm d-1", "m yr-1", s, o, stat=stat)
    call check("mass flux vs velocity: stat /= 0", stat .ne. 0)
    call units_conv("mm d-1", "kg m-2 yr-1", s, o, stat=stat)
    call check("valid: stat = 0", stat .eq. 0)

    ! --- sp round trip mW m-2 -> W m-2 (dp) -> mW m-2 ------------------------------
    call units_conv("mW m-2", "W m-2", s,  o)
    call units_conv("W m-2", "mW m-2", s1, o1)
    nbad = 0
    x = 1.0e-3_sp
    do k = 1, 200000
        y = real((real(x,dp)*s + o)*s1 + o1, sp)
        if (y .ne. x) nbad = nbad + 1
        x = x*1.0001_sp + 1.0e-4_sp
        if (x .gt. 1.0e4_sp) x = 1.0e-3_sp + real(k,sp)*1.0e-7_sp
    end do
    call check("sp round trip mW m-2 -> W m-2 -> mW m-2 (2e5 values)", nbad .eq. 0)

    ! --- fesmdata factors as before fesm_units (bit patterns) --------------------
    call nc_create(fname)
    call nc_write_dim(fname, "xc", x=1.0_wp, dx=1.0_wp, nx=1)
    i = 0
    call flux("kg m-2 s-1",   "kg m-2 yr-1", int(z"4BF0C3F0", int32))
    call flux("mm d-1",       "kg m-2 yr-1", int(z"43B6A000", int32))
    call flux("kg m-2 yr-1",  "mm d-1",      int(z"3B336D84", int32))
    call flux("kg m-2 s-1",   "mm d-1",      int(z"47A8C000", int32))
    call flux("mm w.e. a-1",  "mm d-1",      int(z"3B336D84", int32))
    call flux("mm/d",         "mm d-1",      int(z"3F800000", int32))
    call flux("mm*d**-1",     "mm d-1",      int(z"3F800000", int32))
    call flux("kg m^-2 s^-1", "mm yr-1",     int(z"4BF0C3F0", int32))
    call flux("mm s-1",       "kg m-2 d-1",  int(z"47A8C000", int32))
    call flux("mm d**-1",     "mm yr-1",     int(z"43B6A000", int32))
    call tref("K",               int(z"00000000", int32))
    call tref("degC",            int(z"43889333", int32))
    call tref("degrees C",       int(z"43889333", int32))
    call tref("degrees Celcius", int(z"43889333", int32))
    call tref("degree_celsius",  int(z"43889333", int32))
    call tref("kelvin",          int(z"00000000", int32))
    call tref("degK",            int(z"00000000", int32))
    call tref("Celsius",         int(z"43889333", int32))
    call tref("deg_C",           int(z"43889333", int32))
    call tref("degreesK",        int(z"00000000", int32))

    if (nfail .gt. 0) then
        write(*,*) "test_fesm_units: ", nfail, " FAILED"
        error stop 1
    end if
    write(*,*) "test_fesm_units: all passed"

contains

    subroutine row(from, to, scale, offset, rho_ice, sec_year)
        ! from -> to gives (scale, offset) to a relative 1e-15; to -> from
        ! its inverse.
        character(len=*), intent(in) :: from, to
        real(dp),         intent(in) :: scale, offset
        real(dp),         intent(in), optional :: rho_ice, sec_year
        real(dp) :: s, o, si, oi
        logical  :: ok
        call units_conv(from, to, s, o, rho_ice=rho_ice, sec_year=sec_year)
        call units_conv(to, from, si, oi, rho_ice=rho_ice, sec_year=sec_year)
        ok = abs(s - scale) .le. 1.0e-15_dp*abs(scale) .and. abs(o - offset) .le. 1.0e-12_dp .and. &
             abs(si - 1.0_dp/scale) .le. 1.0e-15_dp/abs(scale) .and. abs(oi + offset/scale) .le. 1.0e-12_dp
        call check(from//" <-> "//to, ok)
        if (.not. ok) write(*,*) "     ", s, o, si, oi
    end subroutine row

    subroutine flux(units, units_out, bits)
        character(len=*), intent(in) :: units, units_out
        integer(int32),   intent(in) :: bits
        character(len=8) :: vname
        real(wp) :: f
        i = i + 1
        write(vname,"(a,i0)") "p", i
        call nc_write(fname, trim(vname), [0.0_wp], dim1="xc", units=units)
        f = fesmdata_flux_conv(fname, trim(vname), units_out)
        call check("fesmdata_flux_conv "//units//" -> "//units_out//" bit-identical", &
                   transfer(f, 0_int32) .eq. bits)
    end subroutine flux

    subroutine tref(units, bits)
        character(len=*), intent(in) :: units
        integer(int32),   intent(in) :: bits
        character(len=8) :: vname
        real(wp) :: t
        i = i + 1
        write(vname,"(a,i0)") "t", i
        call nc_write(fname, trim(vname), [0.0_wp], dim1="xc", units=units)
        t = fesmdata_temp_ref(fname, trim(vname))
        call check("fesmdata_temp_ref "//units//" bit-identical", transfer(t, 0_int32) .eq. bits)
    end subroutine tref

    subroutine check(label, ok)
        character(len=*), intent(in) :: label
        logical,          intent(in) :: ok
        write(*,"(a,a)") merge("  ok   ", "  FAIL ", ok), label
        if (.not. ok) nfail = nfail + 1
    end subroutine check

end program test_fesm_units
