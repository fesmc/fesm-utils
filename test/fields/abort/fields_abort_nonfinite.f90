program fields_abort_nonfinite
    ! Expected abort (make check-fields): a NaN at a valid cell.
    ! Must stop through fesm_error with the message announced on the EXPECT line.

    use, intrinsic :: iso_fortran_env, only : output_unit
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use precision, only : dp
    use fesm_fields

    implicit none

    type(field_state_t) :: st
    integer :: h

    write(output_unit, "(a)") "EXPECT: non-finite value"
    flush(output_unit)

    call fields_init(st, "abort_test")
    call fields_add(st, "probe", "m", FIELD_INSTANT, FIELD_INTENSIVE, h, grid="g", dims=[2,2], masked=.true.)
    call fields_lock(st)
    call fields_put(st, h, reshape([1.0_dp, ieee_value(1.0_dp, ieee_quiet_nan), 3.0_dp, 4.0_dp], [2,2]), &
                    t1=0.0_dp, valid=reshape([.true., .true., .false., .true.], [2,2]))

    write(output_unit, "(a)") "FAIL: fields_abort_nonfinite did not abort"

end program fields_abort_nonfinite
