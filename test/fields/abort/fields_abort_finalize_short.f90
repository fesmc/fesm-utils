program fields_abort_finalize_short
    ! Expected abort (make check-fields): fields_finalize_interval short of t.
    ! Must stop through fesm_error with the message announced on the EXPECT line.

    use, intrinsic :: iso_fortran_env, only : output_unit
    use precision, only : dp
    use fesm_fields

    implicit none

    type(field_state_t) :: st
    integer :: h

    write(output_unit, "(a)") "EXPECT: accumulation does not reach t"
    flush(output_unit)

    call fields_init(st, "abort_test")
    call fields_add(st, "probe", "m", FIELD_MEAN, FIELD_INTENSIVE, h)
    call fields_lock(st)
    call fields_put(st, h, 1.0_dp, t0=0.0_dp, t1=1.0_dp)
    call fields_finalize_interval(st, 2.0_dp)

    write(output_unit, "(a)") "FAIL: fields_abort_finalize_short did not abort"

end program fields_abort_finalize_short
