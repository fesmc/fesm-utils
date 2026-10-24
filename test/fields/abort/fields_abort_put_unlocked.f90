program fields_abort_put_unlocked
    ! Expected abort (make check-fields): fields_put before fields_lock.
    ! Must stop through fesm_error with the message announced on the EXPECT line.

    use, intrinsic :: iso_fortran_env, only : output_unit
    use precision, only : dp
    use fesm_fields

    implicit none

    type(field_state_t) :: st
    integer :: h

    write(output_unit, "(a)") "EXPECT: state is not locked"
    flush(output_unit)

    call fields_init(st, "abort_test")
    call fields_add(st, "probe", "m", FIELD_INSTANT, FIELD_INTENSIVE, h)
    call fields_put(st, h, 1.0_dp, t1=0.0_dp)

    write(output_unit, "(a)") "FAIL: fields_abort_put_unlocked did not abort"

end program fields_abort_put_unlocked
