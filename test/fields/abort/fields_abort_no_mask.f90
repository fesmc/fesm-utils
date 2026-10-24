program fields_abort_no_mask
    ! Expected abort (make check-fields): a put with a mask on an unmasked field.
    ! Must stop through fesm_error with the message announced on the EXPECT line.

    use, intrinsic :: iso_fortran_env, only : output_unit
    use precision, only : dp
    use fesm_fields

    implicit none

    type(field_state_t) :: st
    integer :: h

    write(output_unit, "(a)") "EXPECT: field has no mask"
    flush(output_unit)

    call fields_init(st, "abort_test")
    call fields_add(st, "probe", "m", FIELD_INSTANT, FIELD_INTENSIVE, h, grid="g", dims=[2,2])
    call fields_lock(st)
    call fields_put(st, h, reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], [2,2]), t1=0.0_dp, &
                    valid=reshape([.true., .true., .true., .true.], [2,2]))

    write(output_unit, "(a)") "FAIL: fields_abort_no_mask did not abort"

end program fields_abort_no_mask
