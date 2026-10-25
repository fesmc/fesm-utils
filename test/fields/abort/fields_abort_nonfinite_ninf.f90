program fields_abort_nonfinite_ninf
    ! Expected abort (make check-fields): -Inf, double precision, INSTANT at a valid cell.
    ! Must stop through fesm_error with the message announced on the EXPECT line.

    use, intrinsic :: iso_fortran_env, only : output_unit
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_negative_inf
    use precision, only : sp, dp
    use fesm_fields

    implicit none

    type(field_state_t) :: st
    integer :: h
    real(dp) :: var(2,2)

    write(output_unit, "(a)") "EXPECT: non-finite value"
    flush(output_unit)

    call fields_init(st, "abort_test")
    call fields_add(st, "probe", "m", FIELD_INSTANT, FIELD_INTENSIVE, h, grid="g", dims=[2,2])
    call fields_lock(st)
    var = reshape([1.0_dp, 2.0_dp, ieee_value(1.0_dp, ieee_negative_inf), 4.0_dp], [2,2])
    call fields_put(st, h, var, t1=0.0_dp)

    write(output_unit, "(a)") "FAIL: fields_abort_nonfinite_ninf did not abort"

end program fields_abort_nonfinite_ninf
