program fields_abort_gap
    ! Expected abort (make check-fields): a MEAN put leaving a gap after the window.
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
    call fields_put(st, h, 1.0_dp, t0=1.5_dp, t1=2.0_dp)

    write(output_unit, "(a)") "FAIL: fields_abort_gap did not abort"

end program fields_abort_gap
