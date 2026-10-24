program fields_abort_declaration
    ! Expected abort (make check-fields): a rank-2 field without a grid.
    ! Must stop through fesm_error with the message announced on the EXPECT line.

    use, intrinsic :: iso_fortran_env, only : output_unit
    use precision, only : dp
    use fesm_fields

    implicit none

    type(field_state_t) :: st
    integer :: h

    write(output_unit, "(a)") "EXPECT: invalid declaration"
    flush(output_unit)

    call fields_init(st, "abort_test")
    call fields_add(st, "probe", "m", FIELD_INSTANT, FIELD_INTENSIVE, h, dims=[3,2])

    write(output_unit, "(a)") "FAIL: fields_abort_declaration did not abort"

end program fields_abort_declaration
