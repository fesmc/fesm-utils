program fields_abort_link_units
    ! Expected abort (make check-fields): a link between m and m yr-1.
    ! Must stop through fesm_error with the message announced on the EXPECT line.

    use, intrinsic :: iso_fortran_env, only : output_unit
    use precision, only : dp
    use fesm_fields

    implicit none

    type(field_state_t) :: s, d
    integer :: hs, hd

    write(output_unit, "(a)") "EXPECT: units differ and no equivalence is known"
    flush(output_unit)

    call fields_init(s, "abort_test:src")
    call fields_init(d, "abort_test:dst")
    call fields_add(s, "src", "m", FIELD_INSTANT, FIELD_INTENSIVE, hs)
    call fields_add(d, "probe", "m yr-1", FIELD_INSTANT, FIELD_INTENSIVE, hd)
    call fields_check_link(s, hs, d, hd)

    write(output_unit, "(a)") "FAIL: fields_abort_link_units did not abort"

end program fields_abort_link_units
