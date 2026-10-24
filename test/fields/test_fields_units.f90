program test_fields_units
    ! fesm_fields units: units_key normalisation, every alias of the design
    ! note's table equivalent to its canonical spelling, and units that need a
    ! factor or an offset kept apart.

    use fesm_fields

    implicit none

    integer :: nfail, i
    character(len=16), parameter :: rows(2, 21) = reshape([ character(len=16) :: &
        "kg m-2 yr-1", "kg m-2 a-1",   "kg m-2 yr-1", "kg/m2/yr",     &
        "kg m-2 yr-1", "mm w.e. yr-1", "kg m-2 yr-1", "mm we yr-1",   &
        "kg m-2 yr-1", "mm w.e./yr",   "kg m-2 yr-1", "mm/yr w.e.",   &
        "kg m-2",      "kg/m2",        "kg m-2",      "mm w.e.",      &
        "kg m-2",      "mm we",                                       &
        "m",           "meter",        "m",           "meters",       &
        "m",           "metre",        "m",           "metres",       &
        "m yr-1",      "m a-1",        "m yr-1",      "m/yr",         &
        "m yr-1",      "m/a",                                         &
        "K",           "kelvin",                                      &
        "W m-2",       "W/m2",         "W m-2",       "W m**-2",      &
        "1",           "-",            "1",           "fraction"      ], [2, 21])

    nfail = 0

    ! ----- units_key ---------------------------------------------------------------
    call check("units_key('kg m**-2 yr-1') = kgm-2yr-1", units_key("kg m**-2 yr-1") == "kgm-2yr-1")
    call check("units_key('kg m^-2 d^-1') = kgm-2d-1",   units_key("kg m^-2 d^-1") == "kgm-2d-1")
    call check("units_key('  K ') = k",                  units_key("  K ") == "k")
    call check("units_key('W M-2') = wm-2",              units_key("W M-2") == "wm-2")

    ! ----- aliases ----------------------------------------------------------------
    do i = 1, size(rows, 2)
        call check("'"//trim(rows(2,i))//"' == '"//trim(rows(1,i))//"'", &
                   units_equivalent(rows(2,i), rows(1,i)))
    end do
    call check("units_canonical('mm w.e./yr') = kgm-2yr-1", units_canonical("mm w.e./yr") == "kgm-2yr-1")
    call check("units_canonical('degC') = degc (own key)",   units_canonical("degC") == "degc")
    call check("spelling variants via units_key", units_equivalent("kg m**-2 yr**-1", "KG M-2 YR-1"))

    ! ----- not equivalent ---------------------------------------------------------
    call check("m /= m yr-1",           .not. units_equivalent("m", "m yr-1"))
    call check("kg m-2 yr-1 /= m yr-1", .not. units_equivalent("kg m-2 yr-1", "m yr-1"))
    call check("K /= degC",             .not. units_equivalent("K", "degC"))
    call check("kg m-2 /= kg m-2 yr-1", .not. units_equivalent("kg m-2", "kg m-2 yr-1"))
    call check("m i.e. yr-1 /= mm w.e. yr-1 (needs a factor)", &
               .not. units_equivalent("m i.e. yr-1", "mm w.e. yr-1"))

    if (nfail > 0) then
        write(*,*) "test_fields_units: ", nfail, " check(s) failed."
        stop 1
    end if
    write(*,*) "PASS: test_fields_units"

contains

    subroutine check(label, ok)
        character(len=*), intent(in) :: label
        logical,          intent(in) :: ok
        if (ok) then
            write(*,"(a)") "  ok    "//label
        else
            write(*,"(a)") "  FAIL  "//label
            nfail = nfail + 1
        end if
    end subroutine check

end program test_fields_units
