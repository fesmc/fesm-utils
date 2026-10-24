program fields_abort_overlap
    ! Expected abort (make check-fields): a MEAN put overlapping the window (an exchange run twice).
    ! Must stop through fesm_error with the message announced on the EXPECT line.

    use, intrinsic :: iso_fortran_env, only : output_unit
    use precision, only : dp
    use fesm_fields

    implicit none

    type(field_state_t) :: st
    integer :: h

    write(output_unit, "(a)") "EXPECT: gap or overlap in accumulation"
    flush(output_unit)

    call fields_init(st, "abort_test")
    call fields_add(st, "probe", "m", FIELD_MEAN, FIELD_INTENSIVE, h)
    call fields_lock(st)
    call fields_put(st, h, 1.0_dp, t0=0.0_dp, t1=1.0_dp)
    call fields_put(st, h, 1.0_dp, t0=0.0_dp, t1=1.0_dp)

    write(output_unit, "(a)") "FAIL: fields_abort_overlap did not abort"

end program fields_abort_overlap
