program fields_abort_restart_diag_state
    ! Expected abort (make check-fields): a state field missing from the restart file
    ! (written as a diagnostic, restart = .false., read as a state field).
    ! Must stop through fesm_error with the message announced on the EXPECT line.

    use, intrinsic :: iso_fortran_env, only : output_unit
    use precision, only : dp
    use fesm_fields

    implicit none

    type(field_state_t) :: a, b
    integer :: h

    write(output_unit, "(a)") "EXPECT: restart file does not match the declaration"
    flush(output_unit)

    call fields_init(a, "abort_test")
    call fields_add(a, "H",     "m", FIELD_INSTANT, FIELD_INTENSIVE, h, grid="g", dims=[3,2])
    call fields_add(a, "probe", "m", FIELD_INSTANT, FIELD_INTENSIVE, h, restart=.false.)
    call fields_lock(a)
    call fields_put(a, 1, reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp], [3,2]), t1=1.0_dp)
    call fields_put(a, 2, 0.5_dp, t1=1.0_dp)
    call fields_restart_write(a, "test_fields_abort_restart_diag_state.nc")
    call fields_init(b, "abort_test")
    call fields_add(b, "H",     "m", FIELD_INSTANT, FIELD_INTENSIVE, h, grid="g", dims=[3,2])
    call fields_add(b, "probe", "m", FIELD_INSTANT, FIELD_INTENSIVE, h)
    call fields_lock(b)
    call fields_restart_read(b, "test_fields_abort_restart_diag_state.nc")

    write(output_unit, "(a)") "FAIL: fields_abort_restart_diag_state did not abort"

end program fields_abort_restart_diag_state
