program fields_abort_read_unpublished
    ! Expected abort (make check-fields): fields_get before any value was published.
    ! Must stop through fesm_error with the message announced on the EXPECT line.

    use, intrinsic :: iso_fortran_env, only : output_unit
    use precision, only : dp
    use fesm_fields

    implicit none

    type(field_state_t) :: st
    integer :: h
    real(dp) :: x

    write(output_unit, "(a)") "EXPECT: field read before any value was published"
    flush(output_unit)

    call fields_init(st, "abort_test")
    call fields_add(st, "probe", "m", FIELD_MEAN, FIELD_INTENSIVE, h)
    call fields_lock(st)
    call fields_put(st, h, 1.0_dp, t0=0.0_dp, t1=1.0_dp)
    call fields_get(st, h, x)

    write(output_unit, "(a)") "FAIL: fields_abort_read_unpublished did not abort"

end program fields_abort_read_unpublished
