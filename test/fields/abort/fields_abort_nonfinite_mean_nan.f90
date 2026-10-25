program fields_abort_nonfinite_mean_nan
    ! Expected abort (make check-fields): a NaN, double precision, MEAN at a valid cell.
    ! Must stop through fesm_error with the message announced on the EXPECT line.

    use, intrinsic :: iso_fortran_env, only : output_unit
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use precision, only : sp, dp
    use fesm_fields

    implicit none

    type(field_state_t) :: st
    integer :: h
    real(dp) :: var(2,2)

    write(output_unit, "(a)") "EXPECT: non-finite value"
    flush(output_unit)

    call fields_init(st, "abort_test")
    call fields_add(st, "probe", "m", FIELD_MEAN, FIELD_INTENSIVE, h, grid="g", dims=[2,2])
    call fields_lock(st)
    var = reshape([1.0_dp, 2.0_dp, ieee_value(1.0_dp, ieee_quiet_nan), 4.0_dp], [2,2])
    call fields_put(st, h, var, t0=0.0_dp, t1=1.0_dp)

    write(output_unit, "(a)") "FAIL: fields_abort_nonfinite_mean_nan did not abort"

end program fields_abort_nonfinite_mean_nan
