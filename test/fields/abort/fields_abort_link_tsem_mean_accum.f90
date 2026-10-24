program fields_abort_link_tsem_mean_accum
    ! Expected abort (make check-fields): a link MEAN -> ACCUM.
    ! Must stop through fesm_error with the message announced on the EXPECT line.

    use, intrinsic :: iso_fortran_env, only : output_unit
    use precision, only : dp
    use fesm_fields

    implicit none

    type(field_state_t) :: s, d
    integer :: hs, hd

    write(output_unit, "(a)") "EXPECT: time semantics incompatible"
    flush(output_unit)

    call fields_init(s, "abort_test:src")
    call fields_init(d, "abort_test:dst")
    call fields_add(s, "src", "m", FIELD_MEAN, FIELD_INTENSIVE, hs)
    call fields_add(d, "probe", "m", FIELD_ACCUM, FIELD_INTENSIVE, hd)
    call fields_check_link(s, hs, d, hd)

    write(output_unit, "(a)") "FAIL: fields_abort_link_tsem_mean_accum did not abort"

end program fields_abort_link_tsem_mean_accum
