program test_fields_accum
    ! fesm_fields ACCUM accumulation: a field declared in integral units
    ! ("kg m-2") fed rates per year publishes sum(var*dt) exactly, with
    ! validity [ta0,ta1]; the declared units are kept.

    use precision, only : sp, dp
    use fesm_fields

    implicit none

    type(field_state_t) :: st
    integer  :: h, nfail, k
    real(dp) :: rate(3,2), y(3,2), expect(3,2), tk(4), t0, t1
    real(sp) :: rate_sp(3,2)

    nfail = 0

    call fields_init(st, "test:accum")
    call fields_add(st, "melt", "kg m-2", FIELD_ACCUM, FIELD_EXTENSIVE, h, grid="g", dims=[3,2])
    call fields_lock(st)

    tk = [0.0_dp, 0.75_dp, 1.0_dp, 3.0_dp]
    expect = 0.0_dp
    do k = 1, 3
        rate = 0.1_dp*real(k, dp) + reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp], [3,2])
        if (k == 2) then
            rate_sp = real(rate, sp)                 ! sp input widens exactly
            call fields_put(st, h, rate_sp, t0=tk(k), t1=tk(k+1))
            expect = expect + real(rate_sp, dp)*(tk(k+1) - tk(k))
        else
            call fields_put(st, h, rate, t0=tk(k), t1=tk(k+1))
            expect = expect + rate*(tk(k+1) - tk(k))
        end if
    end do

    call fields_finalize_interval(st, 3.0_dp)
    call fields_get(st, h, y, t0=t0, t1=t1)

    call check("published = sum(var*dt), exactly", all(y == expect))
    call check("validity = [0,3]", t0 == 0.0_dp .and. t1 == 3.0_dp)
    call check("declared units kept", st%f(h)%units == "kg m-2" .and. st%f(h)%ukey == "kgm-2")

    ! A second interval starts from zero.
    call fields_put(st, h, spread(spread(2.0_dp, 1, 3), 2, 2), t0=3.0_dp, t1=3.5_dp)
    call fields_finalize_interval(st, 3.5_dp)
    call fields_get(st, h, y)
    call check("second interval starts from zero", all(y == 1.0_dp))

    if (nfail > 0) then
        write(*,*) "test_fields_accum: ", nfail, " check(s) failed."
        stop 1
    end if
    write(*,*) "PASS: test_fields_accum"

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

end program test_fields_accum
