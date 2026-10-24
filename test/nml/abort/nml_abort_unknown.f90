program nml_abort_unknown
    ! Expected abort (make check-nml): nml_validate against two templates
    ! stops on a key (`probe`) that neither declares.

    use, intrinsic :: iso_fortran_env, only : output_unit
    use nml, only : nml_validate

    implicit none

    integer :: u

    write(output_unit, "(a)") "EXPECT: unknown key(s) in group &inst of test_nml_abort.nml:"
    flush(output_unit)

    open(newunit=u, file="test_nml_abort.nml", status="replace", action="write")
    write(u,"(a)") "&inst"
    write(u,"(a)") "  dt    = 1.0"
    write(u,"(a)") "  probe = 1"
    write(u,"(a)") "/"
    write(u,"(a)") "&core"
    write(u,"(a)") "  dt = 1.0"
    write(u,"(a)") "/"
    write(u,"(a)") "&back"
    write(u,"(a)") "  rf = 1.0"
    write(u,"(a)") "/"
    close(u)

    call nml_validate("test_nml_abort.nml", "inst", [character(len=32) :: "test_nml_abort.nml", &
                      "test_nml_abort.nml"], [character(len=32) :: "core", "back"])

    write(output_unit, "(a)") "FAIL: nml_abort_unknown did not abort"

end program nml_abort_unknown
