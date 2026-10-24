program fields_abort_shape
    ! Expected abort (make check-fields): a put of the wrong shape.
    ! Must stop through fesm_error with the message announced on the EXPECT line.

    use, intrinsic :: iso_fortran_env, only : output_unit
    use precision, only : dp
    use fesm_fields

    implicit none

    type(field_state_t) :: st
    integer :: h

    write(output_unit, "(a)") "EXPECT: shape mismatch"
    flush(output_unit)

    call fields_init(st, "abort_test")
    call fields_add(st, "probe", "m", FIELD_INSTANT, FIELD_INTENSIVE, h, grid="g", dims=[3,2])
    call fields_lock(st)
    call fields_put(st, h, reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp], [2,3]), t1=0.0_dp)

    write(output_unit, "(a)") "FAIL: fields_abort_shape did not abort"

end program fields_abort_shape
