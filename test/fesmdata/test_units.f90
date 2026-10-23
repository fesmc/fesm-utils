program test_units
    ! fesmdata_flux_conv and fesmdata_temp_ref against the units attributes
    ! of a small file written here.

    use precision
    use ncio
    use fesmdata, only : fesmdata_flux_conv, fesmdata_temp_ref

    implicit none

    character(len=*), parameter :: fname = "test_units.nc"
    real(wp), parameter :: tol = 1e-6_wp
    integer :: nfail

    nfail = 0

    call nc_create(fname)
    call nc_write_dim(fname, "xc", x=1.0_wp, dx=1.0_wp, nx=1)
    call put(fname, "pr_yr",  "kg m-2 yr-1")
    call put(fname, "pr_d",   "mm d**-1")
    call put(fname, "pr_s",   "kg m-2 s-1")
    call put(fname, "pr_a",   "mm w.e. a-1")
    call put(fname, "pr_sl",  "mm/d")
    call put(fname, "pr_st",  "mm*d**-1")
    call put(fname, "t_k",    "K")
    call put(fname, "t_c",    "degC")
    call put(fname, "t_c2",   "degrees C")
    call put(fname, "t_c3",   "degrees Celcius")
    call put(fname, "t_c4",   "degree_celsius")

    call check("yr -> d",  fesmdata_flux_conv(fname, "pr_yr", "mm d-1"),  1.0_wp/365.25_wp, nfail)
    call check("yr -> yr", fesmdata_flux_conv(fname, "pr_yr", "mm yr-1"), 1.0_wp,           nfail)
    call check("d -> yr",  fesmdata_flux_conv(fname, "pr_d",  "kg m-2 yr-1"), 365.25_wp,    nfail)
    call check("s -> d",   fesmdata_flux_conv(fname, "pr_s",  "mm d-1"),  86400.0_wp,       nfail)
    call check("a -> d",   fesmdata_flux_conv(fname, "pr_a",  "mm d-1"),  1.0_wp/365.25_wp, nfail)
    call check("/d -> d",  fesmdata_flux_conv(fname, "pr_sl", "mm d-1"),  1.0_wp,           nfail)
    call check("*d -> d",  fesmdata_flux_conv(fname, "pr_st", "mm d-1"),  1.0_wp,           nfail)
    call check("K",        fesmdata_temp_ref(fname, "t_k"),  0.0_wp,   nfail)
    call check("degC",     fesmdata_temp_ref(fname, "t_c"),  273.15_wp, nfail)
    call check("degrees C",fesmdata_temp_ref(fname, "t_c2"), 273.15_wp, nfail)
    call check("degrees Celcius", fesmdata_temp_ref(fname, "t_c3"), 273.15_wp, nfail)
    call check("degree_celsius",  fesmdata_temp_ref(fname, "t_c4"), 273.15_wp, nfail)

    if (nfail .gt. 0) then
        write(*,*) "test_units: ", nfail, " FAILED"
        error stop 1
    end if
    write(*,*) "test_units: all passed"

contains

    subroutine put(fname, name, units)
        character(len=*), intent(IN) :: fname, name, units
        call nc_write(fname, name, [0.0_wp], dim1="xc", units=units)
    end subroutine put

    subroutine check(label, val, ref, nfail)
        character(len=*), intent(IN)    :: label
        real(wp),         intent(IN)    :: val, ref
        integer,          intent(INOUT) :: nfail
        logical :: ok
        ok = abs(val-ref) .le. tol*max(1.0_wp, abs(ref))
        write(*,"(a,a,a,2g16.8)") merge("  ok   ", "  FAIL ", ok), label, ": ", val, ref
        if (.not. ok) nfail = nfail + 1
    end subroutine check

end program test_units
