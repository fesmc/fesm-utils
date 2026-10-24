program test_fields_put_get
    ! fesm_fields INSTANT put/get round trips, 0-D, 2-D and 3-D:
    !   * dp -> dp is bit-exact; sp -> dp is the exact widening real(x, dp);
    !     dp -> sp equals real(x, sp);
    !   * validity is [t,t] without t0, [t0,t1] with it;
    !   * a rank-3 array with nz = 1 is accepted for a rank-2 field (put and get);
    !   * a put replaces the previous value.

    use precision, only : sp, dp
    use fesm_fields

    implicit none

    type(field_state_t) :: st
    integer  :: hs, h2, h3, nfail, i, j, k
    real(dp) :: x2(3,2), y2(3,2), x3(3,2,4), y3(3,2,4), z3(3,2,1), xs, ys, t0, t1
    real(sp) :: s2(3,2), r2(3,2), s3(3,2,4), r3(3,2,4), ss, rs

    nfail = 0

    call fields_init(st, "test:put_get")
    call fields_add(st, "bsl",  "m", FIELD_INSTANT, FIELD_INTENSIVE, hs)
    call fields_add(st, "z",    "m", FIELD_INSTANT, FIELD_INTENSIVE, h2, grid="g", dims=[3,2])
    call fields_add(st, "prof", "K", FIELD_INSTANT, FIELD_INTENSIVE, h3, grid="g", dims=[3,2,4], zdim="depth")
    call fields_lock(st)

    ! Non-dyadic values, so that a lossy path would show.
    do j = 1, 2
        do i = 1, 3
            x2(i,j) = 1.0_dp/real(i+3*j, dp) + 1.0e-9_dp*real(i, dp)
            do k = 1, 4
                x3(i,j,k) = sin(real(i + 10*j + 100*k, dp))
            end do
        end do
    end do
    xs = acos(-1.0_dp)/7.0_dp

    ! ----- dp -> dp ---------------------------------------------------------------
    call fields_put(st, hs, xs, t1=5.0_dp)
    call fields_get(st, hs, ys, t0=t0, t1=t1)
    call check("0-D dp -> dp bit-exact", ys == xs)
    call check("INSTANT validity [t,t] without t0", t0 == 5.0_dp .and. t1 == 5.0_dp)

    call fields_put(st, h2, x2, t1=5.0_dp)
    call fields_get(st, h2, y2)
    call check("2-D dp -> dp bit-exact", all(y2 == x2))

    call fields_put(st, h3, x3, t0=4.0_dp, t1=5.0_dp)
    call fields_get(st, h3, y3, t0=t0, t1=t1)
    call check("3-D dp -> dp bit-exact", all(y3 == x3))
    call check("INSTANT validity [t0,t1] as given", t0 == 4.0_dp .and. t1 == 5.0_dp)
    call check("published flag", st%f(h3)%published)

    ! ----- dp -> sp ---------------------------------------------------------------
    call fields_get(st, hs, rs)
    call check("0-D dp -> sp = real(x, sp)", rs == real(xs, sp))
    call fields_get(st, h2, r2)
    call check("2-D dp -> sp = real(x, sp)", all(r2 == real(x2, sp)))
    call fields_get(st, h3, r3)
    call check("3-D dp -> sp = real(x, sp)", all(r3 == real(x3, sp)))

    ! ----- sp -> dp ---------------------------------------------------------------
    ss = real(xs, sp)
    s2 = real(x2, sp)
    s3 = real(x3, sp)
    call fields_put(st, hs, ss, t1=6.0_dp)
    call fields_get(st, hs, ys)
    call check("0-D sp -> dp exact widening", ys == real(ss, dp))
    call fields_put(st, h2, s2, t1=6.0_dp)
    call fields_get(st, h2, y2, t1=t1)
    call check("2-D sp -> dp exact widening", all(y2 == real(s2, dp)))
    call check("a put replaces the previous value and validity", t1 == 6.0_dp)
    call fields_put(st, h3, s3, t1=6.0_dp)
    call fields_get(st, h3, y3)
    call check("3-D sp -> dp exact widening", all(y3 == real(s3, dp)))

    ! ----- rank-3 nz = 1 into a rank-2 field -----------------------------------------
    z3(:,:,1) = 2.0_dp*x2
    call fields_put(st, h2, z3, t1=7.0_dp)
    call fields_get(st, h2, y2)
    call check("rank-3 nz=1 put into a rank-2 field", all(y2 == z3(:,:,1)))
    z3 = 0.0_dp
    call fields_get(st, h2, z3)
    call check("rank-3 nz=1 get from a rank-2 field", all(z3(:,:,1) == 2.0_dp*x2))

    if (nfail > 0) then
        write(*,*) "test_fields_put_get: ", nfail, " check(s) failed."
        stop 1
    end if
    write(*,*) "PASS: test_fields_put_get"

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

end program test_fields_put_get
