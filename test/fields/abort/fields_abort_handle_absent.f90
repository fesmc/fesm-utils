program fields_abort_handle_absent
    ! Expected abort (make check-fields): fields_handle of an undeclared name.
    ! Must stop through fesm_error with the message announced on the EXPECT line.

    use, intrinsic :: iso_fortran_env, only : output_unit
    use precision, only : dp
    use fesm_fields

    implicit none

    type(field_state_t) :: st
    integer :: h

    write(output_unit, "(a)") "EXPECT: no field named 'probe'"
    flush(output_unit)

    call fields_init(st, "abort_test")
    call fields_add(st, "smb", "m", FIELD_INSTANT, FIELD_INTENSIVE, h)
    h = fields_handle(st, "probe")

    write(output_unit, "(a)") "FAIL: fields_abort_handle_absent did not abort"

end program fields_abort_handle_absent
