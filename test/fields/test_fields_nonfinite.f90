program test_fields_nonfinite
    ! The finiteness check of fields_put (a copy of the value, or its
    ! accumulation, at valid cells; the aborts at valid cells are the
    ! abort/fields_abort_nonfinite_* programs):
    !   * finite values pass and are stored bit for bit: huge, -huge, 0, -0,
    !     tiny, the smallest subnormals, in dp and sp, INSTANT and MEAN
    !     (MEAN: acc = x*dt, with dt = 1);
    !   * NaN, +Inf and -Inf at cells marked invalid pass: an INSTANT field
    !     copies them, a MEAN field skips them.

    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf, &
                                              ieee_negative_inf, ieee_is_nan
    use precision, only : sp, dp
    use fesm_fields

    implicit none

    integer, parameter :: n = 8
    type(field_state_t) :: st
    integer  :: nfail, h_di, h_si, h_dm, h_sm, h_mi, h_mm
    real(dp) :: xd(n,1), vd(n,1), nonfin(n,1), mixd(n,1)
    real(sp) :: xs(n,1), mixs(n,1)
    logical  :: valid(n,1)

    nfail = 0

    ! huge, -huge, 0, -0, tiny, -tiny, smallest subnormal, -(smallest subnormal)
    xd(:,1) = [huge(1.0_dp), -huge(1.0_dp), 0.0_dp, sign(0.0_dp, -1.0_dp), &
               tiny(1.0_dp), -tiny(1.0_dp), tiny(1.0_dp)*epsilon(1.0_dp), -tiny(1.0_dp)*epsilon(1.0_dp)]
    xs(:,1) = [huge(1.0_sp), -huge(1.0_sp), 0.0_sp, sign(0.0_sp, -1.0_sp), &
               tiny(1.0_sp), -tiny(1.0_sp), tiny(1.0_sp)*epsilon(1.0_sp), -tiny(1.0_sp)*epsilon(1.0_sp)]

    call fields_init(st, "test:nonfinite")
    call fields_add(st, "di", "m", FIELD_INSTANT, FIELD_INTENSIVE, h_di, grid="g", dims=[n,1])
    call fields_add(st, "si", "m", FIELD_INSTANT, FIELD_INTENSIVE, h_si, grid="g", dims=[n,1])
    call fields_add(st, "dm", "m", FIELD_MEAN,    FIELD_INTENSIVE, h_dm, grid="g", dims=[n,1])
    call fields_add(st, "sm", "m", FIELD_MEAN,    FIELD_INTENSIVE, h_sm, grid="g", dims=[n,1])
    call fields_add(st, "mi", "m", FIELD_INSTANT, FIELD_INTENSIVE, h_mi, grid="g", dims=[n,1], masked=.true.)
    call fields_add(st, "mm", "m", FIELD_MEAN,    FIELD_INTENSIVE, h_mm, grid="g", dims=[n,1], masked=.true.)
    call fields_lock(st)

    ! --- Finite values: stored bit for bit ----------------------------------
    call fields_put(st, h_di, xd, t1=0.0_dp)
    call fields_get(st, h_di, vd)
    call check("dp INSTANT: huge, -huge, 0, -0, tiny, subnormals stored bit for bit", same_bits(vd, xd))

    call fields_put(st, h_si, xs, t1=0.0_dp)
    call fields_get(st, h_si, vd)
    call check("sp INSTANT: widened to dp bit for bit", same_bits(vd, real(xs, dp)))

    ! MEAN: acc = 0 + x*1 (-0 becomes +0, 0 + (-0) = +0; compared by value)
    call fields_put(st, h_dm, xd, t0=0.0_dp, t1=1.0_dp)
    call check("dp MEAN: finite values accumulated", all(st%f(h_dm)%acc(:,1,1) == xd(:,1)))
    call fields_put(st, h_sm, xs, t0=0.0_dp, t1=1.0_dp)
    call check("sp MEAN: finite values accumulated", all(st%f(h_sm)%acc(:,1,1) == real(xs(:,1), dp)))

    ! --- Non-finite values at invalid cells ---------------------------------
    nonfin(:,1) = [ieee_value(1.0_dp, ieee_quiet_nan), ieee_value(1.0_dp, ieee_positive_inf), &
                   ieee_value(1.0_dp, ieee_negative_inf), 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
    valid(:,1) = [.false., .false., .false., .true., .true., .true., .true., .true.]

    mixd = nonfin
    call fields_put(st, h_mi, mixd, t1=0.0_dp, valid=valid)
    call fields_get(st, h_mi, vd)
    call check("dp INSTANT masked: NaN, +Inf, -Inf at invalid cells copied, no abort", &
               ieee_is_nan(vd(1,1)) .and. vd(2,1) > huge(1.0_dp) .and. vd(3,1) < -huge(1.0_dp) .and. &
               all(vd(4:n,1) == nonfin(4:n,1)))

    mixs = real(nonfin, sp)
    call fields_put(st, h_mi, mixs, t1=1.0_dp, valid=valid)
    call fields_get(st, h_mi, vd)
    call check("sp INSTANT masked: the same, widened", &
               ieee_is_nan(vd(1,1)) .and. vd(2,1) > huge(1.0_dp) .and. vd(3,1) < -huge(1.0_dp) .and. &
               all(vd(4:n,1) == real(mixs(4:n,1), dp)))

    call fields_put(st, h_mm, nonfin, t0=0.0_dp, t1=1.0_dp, valid=valid)
    call check("dp MEAN masked: invalid cells skipped (acc stays 0), valid ones accumulated", &
               all(st%f(h_mm)%acc(1:3,1,1) == 0.0_dp) .and. all(st%f(h_mm)%acc(4:n,1,1) == nonfin(4:n,1)))

    call fields_put(st, h_mm, real(nonfin, sp), t0=1.0_dp, t1=2.0_dp, valid=valid)
    call check("sp MEAN masked: the same", &
               all(st%f(h_mm)%acc(1:3,1,1) == 0.0_dp) .and. all(st%f(h_mm)%acc(4:n,1,1) == 2.0_dp*nonfin(4:n,1)))

    if (nfail > 0) then
        write(*,*) "test_fields_nonfinite: ", nfail, " check(s) failed."
        stop 1
    end if
    write(*,*) "PASS: test_fields_nonfinite"

contains

    logical function same_bits(a, b)
        real(dp), intent(in) :: a(:,:), b(:,:)
        integer :: ii
        same_bits = .true.
        do ii = 1, size(a, 1)
            if (transfer(a(ii,1), 0_int64) /= transfer(b(ii,1), 0_int64)) same_bits = .false.
        end do
    end function same_bits

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

end program test_fields_nonfinite
