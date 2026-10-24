program fields_abort_invalid_handle
    ! Expected abort (make check-fields): a handle outside 1..n.
    ! Must stop through fesm_error with the message announced on the EXPECT line.

    use, intrinsic :: iso_fortran_env, only : output_unit
    use precision, only : dp
    use fesm_fields

    implicit none

    type(field_state_t) :: st
    integer :: h

    write(output_unit, "(a)") "EXPECT: invalid handle"
    flush(output_unit)

    call fields_init(st, "abort_test")
    call fields_add(st, "probe", "m", FIELD_INSTANT, FIELD_INTENSIVE, h)
    call fields_lock(st)
    call fields_put(st, h + 1, 1.0_dp, t1=0.0_dp)

    write(output_unit, "(a)") "FAIL: fields_abort_invalid_handle did not abort"

end program fields_abort_invalid_handle
