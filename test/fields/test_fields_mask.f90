program test_fields_mask
    ! fesm_fields valid masks:
    !   * the published MEAN mask is the AND of the put masks over the window,
    !     and the mean at cells valid throughout is exact;
    !   * fields_get returns the mask;
    !   * a non-finite value at an invalid cell is accepted (INSTANT: stored as
    !     given; MEAN: not accumulated). At a valid cell it stops the program:
    !     see abort/fields_abort_nonfinite.f90.

    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_is_nan
    use precision, only : dp
    use fesm_fields

    implicit none

    type(field_state_t) :: st
    integer  :: hm, hi, nfail
    real(dp) :: x(2,2), y(2,2), nan, inf
    logical  :: v1(2,2), v2(2,2), vout(2,2)

    nfail = 0
    nan = ieee_value(1.0_dp, ieee_quiet_nan)
    inf = ieee_value(1.0_dp, ieee_positive_inf)

    call fields_init(st, "test:mask")
    call fields_add(st, "smb", "kg m-2 yr-1", FIELD_MEAN,    FIELD_EXTENSIVE, hm, grid="g", dims=[2,2], &
                    masked=.true.)
    call fields_add(st, "T",   "K",           FIELD_INSTANT, FIELD_INTENSIVE, hi, grid="g", dims=[2,2], &
                    masked=.true.)
    call fields_lock(st)

    ! ----- MEAN: AND over the window -------------------------------------------
    v1 = reshape([.true., .true., .false., .true.], [2,2])
    v2 = reshape([.true., .false., .false., .true.], [2,2])

    x = reshape([1.0_dp, 2.0_dp, nan, 4.0_dp], [2,2])       ! NaN at an invalid cell
    call fields_put(st, hm, x, t0=0.0_dp, t1=1.0_dp, valid=v1)
    x = reshape([3.0_dp, inf, nan, 8.0_dp], [2,2])          ! Inf, NaN at invalid cells
    call fields_put(st, hm, x, t0=1.0_dp, t1=2.0_dp, valid=v2)
    call fields_finalize_interval(st, 2.0_dp)

    call fields_get(st, hm, y, valid=vout)
    call check("MEAN mask = AND over the window", all(vout .eqv. (v1 .and. v2)))
    call check("mean exact at cells valid throughout", y(1,1) == 2.0_dp .and. y(2,2) == 6.0_dp)
    call check("non-finite at invalid cells accepted (MEAN)", .true.)

    ! The next window starts from a fully valid mask.
    call fields_put(st, hm, reshape([1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], [2,2]), t0=2.0_dp, t1=3.0_dp, &
                    valid=spread(spread(.true., 1, 2), 2, 2))
    call fields_finalize_interval(st, 3.0_dp)
    call fields_get(st, hm, y, valid=vout)
    call check("mask reset with the window", all(vout))

    ! ----- INSTANT: mask copied, invalid values stored as given ----------------
    x = reshape([270.0_dp, nan, 271.0_dp, 272.0_dp], [2,2])
    v1 = reshape([.true., .false., .true., .true.], [2,2])
    call fields_put(st, hi, x, t1=3.0_dp, valid=v1)
    call fields_get(st, hi, y, valid=vout)
    call check("INSTANT mask returned", all(vout .eqv. v1))
    call check("INSTANT valid values", y(1,1) == 270.0_dp .and. y(1,2) == 271.0_dp .and. y(2,2) == 272.0_dp)
    call check("INSTANT invalid value stored as given", ieee_is_nan(y(2,1)))

    if (nfail > 0) then
        write(*,*) "test_fields_mask: ", nfail, " check(s) failed."
        stop 1
    end if
    write(*,*) "PASS: test_fields_mask"

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

end program test_fields_mask
