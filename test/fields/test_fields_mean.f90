program test_fields_mean
    ! fesm_fields MEAN accumulation:
    !   * four puts of unequal length with dyadic values: the published value
    !     is exactly the weighted mean, with validity [ta0,ta1];
    !   * the open window is [ta0,ta1] and invisible until finalize;
    !   * get twice gives the same value (no reset on get);
    !   * finalize clears the window; the second interval is independent;
    !   * INSTANT fields of the same state are untouched by finalize.

    use precision, only : dp
    use fesm_fields

    implicit none

    type(field_state_t) :: st
    integer  :: hm, hs, hi, nfail, k, i, j
    real(dp) :: x(2,2), y(2,2), y_again(2,2), acc(2,2), expect(2,2)
    real(dp) :: tk(5), dt, s, s_acc, ys, t0, t1, xi

    nfail = 0

    call fields_init(st, "test:mean")
    call fields_add(st, "smb", "kg m-2 yr-1", FIELD_MEAN, FIELD_EXTENSIVE, hm, grid="g", dims=[2,2])
    call fields_add(st, "co2", "1",           FIELD_MEAN, FIELD_INTENSIVE, hs)
    call fields_add(st, "z",   "m",           FIELD_INSTANT, FIELD_INTENSIVE, hi)
    call fields_lock(st)

    call fields_put(st, hi, 3.25_dp, t1=10.0_dp)

    ! Window [10, 14] in four pieces of length 0.5, 1.25, 0.25, 2.
    tk = [10.0_dp, 10.5_dp, 11.75_dp, 12.0_dp, 14.0_dp]
    acc   = 0.0_dp
    s_acc = 0.0_dp
    do k = 1, 4
        dt = tk(k+1) - tk(k)
        do j = 1, 2
            do i = 1, 2
                x(i,j) = 0.125_dp*real(k, dp) + 0.5_dp*real(i + 2*j, dp)
            end do
        end do
        s = 0.25_dp*real(k*k, dp)
        call fields_put(st, hm, x, t0=tk(k), t1=tk(k+1))
        call fields_put(st, hs, s, t0=tk(k), t1=tk(k+1))
        acc   = acc + x*dt
        s_acc = s_acc + s*dt
    end do
    expect = acc / 4.0_dp

    call check("window [ta0,ta1] = [10,14]", st%f(hm)%ta0 == 10.0_dp .and. st%f(hm)%ta1 == 14.0_dp)
    call check("wsum = 4", st%f(hm)%wsum == 4.0_dp)
    call check("window invisible before finalize", .not. st%f(hm)%published)

    call fields_finalize_interval(st, 14.0_dp)

    call fields_get(st, hm, y, t0=t0, t1=t1)
    call check("published mean = weighted mean, exactly", all(y == expect))
    call check("validity = [10,14]", t0 == 10.0_dp .and. t1 == 14.0_dp)
    call fields_get(st, hm, y_again)
    call check("get twice gives the same", all(y_again == y))
    call fields_get(st, hs, ys)
    call check("scalar mean exact", ys == s_acc/4.0_dp)

    call check("window cleared: wsum = 0", st%f(hm)%wsum == 0.0_dp)
    call check("window cleared: acc = 0", all(st%f(hm)%acc == 0.0_dp))

    call fields_get(st, hi, xi, t1=t1)
    call check("INSTANT untouched by finalize", xi == 3.25_dp .and. t1 == 10.0_dp)

    ! Second interval [14, 15]: independent of the first.
    call fields_put(st, hm, spread(spread(3.0_dp, 1, 2), 2, 2), t0=14.0_dp, t1=14.5_dp)
    call fields_put(st, hm, spread(spread(5.0_dp, 1, 2), 2, 2), t0=14.5_dp, t1=15.0_dp)
    call fields_put(st, hs, 1.0_dp, t0=14.0_dp, t1=15.0_dp)
    call check("first value still published while the window is open", all(st%f(hm)%v(:,:,1) == expect))
    call fields_finalize_interval(st, 15.0_dp)
    call fields_get(st, hm, y, t0=t0, t1=t1)
    call check("second interval independent", all(y == 4.0_dp))
    call check("second validity = [14,15]", t0 == 14.0_dp .and. t1 == 15.0_dp)

    if (nfail > 0) then
        write(*,*) "test_fields_mean: ", nfail, " check(s) failed."
        stop 1
    end if
    write(*,*) "PASS: test_fields_mean"

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

end program test_fields_mean
