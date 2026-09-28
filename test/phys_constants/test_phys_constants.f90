program test_phys_constants
    ! Exercise the phys_constants API: load from the shipped reference set,
    ! override via a program group, build a record from parameters, and compare
    ! two records that disagree.

    use precision,      only: sp, dp
    use phys_constants

    implicit none

    type(phys_const_class) :: c_ref, c_bench, c_prog
    character(len=512), parameter :: ref = "par/phys_const_earth.nml"
    integer :: ndiff
    logical :: ok

    ok = .true.

    ! --- 1. Load the reference set on its own -------------------------------
    call phys_const_load(c_ref, ref, group="Earth")
    call phys_const_log(c_ref)

    call expect(c_ref%rho_ice, 910.0_dp,  "rho_ice")
    call expect(c_ref%rho_sw,  1028.0_dp, "rho_sw")

    ! Derived quantities follow from the primitives.
    call expect(c_ref%rho_ice_sw, 910.0_dp/1028.0_dp,  "rho_ice_sw")
    call expect(c_ref%conv_we_ie, 1000.0_dp/910.0_dp,  "conv_we_ie")
    call expect(c_ref%omega_melt, (1028.0_dp*3974.0_dp)/(910.0_dp*333500.0_dp), "omega_melt")
    call expect(c_ref%conv_km3_sle, (910.0_dp/1000.0_dp)*1.0e9_dp/(3.618e8_dp*1.0e6_dp), &
                "conv_km3_sle")

    ! --- 2. A benchmark group overriding the reference ----------------------
    ! test_phys_const_bench.nml sets only rho_ice; everything else must fall
    ! back to the Earth reference.
    call phys_const_load(c_bench, "test_phys_const_bench.nml", group="EISMINT", &
                         defaults_file=ref, defaults_group="Earth")

    call expect(c_bench%rho_ice, 917.0_dp,  "bench rho_ice (overridden)")
    call expect(c_bench%rho_sw,  1028.0_dp, "bench rho_sw (from reference)")
    call expect(c_bench%rho_ice_sw, 917.0_dp/1028.0_dp, "bench rho_ice_sw (rederived)")

    ! --- 3. Build from a program's own parameters --------------------------
    call phys_const_set(c_prog, label="test-program", source="test:constants", &
                        g=9.81_dp, T0=273.15_dp, rho_ice=910.0_dp, rho_w=1000.0_dp, &
                        rho_sw=1028.0_dp, L_ice=334000.0_dp, cp_ocn=4187.0_dp)

    call expect(c_prog%rho_ice, 910.0_dp, "prog rho_ice")

    ! Quantities not supplied stay unset, and so do the derived values needing them.
    if (phys_const_isset(c_prog%area_seasurf)) then
        write(*,*) "FAIL: area_seasurf should be unset"; ok = .false.
    end if
    if (phys_const_isset(c_prog%conv_km3_sle)) then
        write(*,*) "FAIL: conv_km3_sle should be unset (needs area_seasurf)"; ok = .false.
    end if

    ! --- 3b. phys_const_get: lookup + check + narrowing --------------------
    block
        real(sp) :: rho_ice_sp, cp_ocn_sp, omega_sp
        real(dp) :: rho_ice_dp

        ! A consumer asserts the record once, then names each quantity it needs.
        call phys_const_require(c_ref, "test:consumer")

        ! Resolved on the kind of the destination.
        call phys_const_get(c_ref, "rho_ice", rho_ice_sp)
        call phys_const_get(c_ref, "RHO_ICE", rho_ice_dp)   ! case-insensitive

        if (abs(real(rho_ice_sp, dp) - 910.0_dp) > 1.0e-3_dp) then
            write(*,*) "FAIL: get(sp) rho_ice = ", rho_ice_sp; ok = .false.
        end if
        call expect(rho_ice_dp, 910.0_dp, "get(dp) rho_ice")

        ! Derived quantities are fetched the same way as primitives.
        call phys_const_get(c_ref, "cp_ocn", cp_ocn_sp)
        call phys_const_get(c_ref, "omega_melt", omega_sp)
        if (abs(real(cp_ocn_sp, dp) - 3974.0_dp) > 1.0e-3_dp) then
            write(*,*) "FAIL: get(sp) cp_ocn = ", cp_ocn_sp; ok = .false.
        end if
        if (real(omega_sp, dp) <= 0.0_dp) then
            write(*,*) "FAIL: get(sp) omega_melt = ", omega_sp; ok = .false.
        end if
    end block

    ! --- 3b(ii). Every advertised name resolves, on a fully supplied record -
    ! This is what keeps phys_const_lookup honest: a name in PHYS_CONST_NAMES
    ! with no entry in the lookup stops the program here, and so does a quantity
    ! that the shipped reference set forgets to supply. Both are caught by the
    ! test rather than by a run.
    block
        integer  :: i
        real(dp) :: v

        do i = 1, size(PHYS_CONST_NAMES)
            ! Stops the program if the name is unrecognised or the value unset.
            call phys_const_get(c_ref, PHYS_CONST_NAMES(i), v)
            if (.not. phys_const_isset(v)) then
                write(*,*) "FAIL: ", trim(PHYS_CONST_NAMES(i)), " resolved to an unset value"
                ok = .false.
            end if
        end do
        write(*,*) "all ", size(PHYS_CONST_NAMES), " names resolve on the reference set"
    end block

    ! --- 3c. Calendar conventions are distinct, named, and exact -----------
    call expect(sec_year_365d,      31536000.0_dp,    "sec_year_365d")
    call expect(sec_year_360d,      31104000.0_dp,    "sec_year_360d")
    call expect(sec_year_366d,      31622400.0_dp,    "sec_year_366d")
    call expect(sec_year_julian,    31557600.0_dp,    "sec_year_julian")
    call expect(sec_year_gregorian, 31556952.0_dp,    "sec_year_gregorian")
    call expect(sec_year_tropical,  31556926.08_dp,   "sec_year_tropical")

    ! A 360-day calendar spanning a full tropical year has longer days.
    call expect(sec_day_360d_tropical, 31556926.08_dp/360.0_dp, "sec_day_360d_tropical")
    if (abs(sec_day_360d_tropical - sec_day) < 1.0_dp) then
        write(*,*) "FAIL: sec_day_360d_tropical should not equal the nominal day"
        ok = .false.
    end if

    ! --- 4. Compare two disagreeing records --------------------------------
    ! The program uses the rounded L_ice and cp_w in place of cp_ocn, so both
    ! should be reported; unset fields must not count as differences.
    call phys_const_compare(c_ref, c_prog, label1="reference", label2="program", n_diff=ndiff)
    if (ndiff /= 2) then
        write(*,*) "FAIL: expected 2 differences (L_ice, cp_ocn), got ", ndiff; ok = .false.
    end if

    ! --- 5. Round-trip through phys_const_write ----------------------------
    call phys_const_write(c_ref, "test_phys_const_used.nml", group="Earth")
    block
        type(phys_const_class) :: c_rt
        call phys_const_load(c_rt, "test_phys_const_used.nml", group="Earth")
        call phys_const_compare(c_ref, c_rt, label1="original", label2="round-trip", &
                                n_diff=ndiff)
        if (ndiff /= 0) then
            write(*,*) "FAIL: write/load round-trip lost precision, ", ndiff, " diffs"
            ok = .false.
        end if
    end block

    ! --- 6. An unfilled record must be rejected, not silently used ---------
    block
        type(phys_const_class) :: c_empty
        if (c_empty%initialized) then
            write(*,*) "FAIL: fresh record claims to be initialized"; ok = .false.
        end if
        if (phys_const_isset(c_empty%rho_ice)) then
            write(*,*) "FAIL: fresh record has a usable rho_ice"; ok = .false.
        end if
    end block

    write(*,*) ""
    if (ok) then
        write(*,*) "test_phys_constants: ALL CHECKS PASSED"
    else
        write(*,*) "test_phys_constants: FAILURES ABOVE"
        error stop 1
    end if

contains

    subroutine expect(got, want, name)
        real(dp), intent(IN) :: got, want
        character(len=*), intent(IN) :: name

        if (abs(got - want) > 1.0e-9_dp*max(1.0_dp, abs(want))) then
            write(*,*) "FAIL: ", trim(name), " got ", got, " want ", want
            ok = .false.
        end if

        return
    end subroutine expect

end program test_phys_constants
