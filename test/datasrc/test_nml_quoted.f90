program test_nml_quoted
    ! nml string vectors: quoted elements keep their blanks
    ! (units = "kg m-2 yr-1" "mm d-1" is two elements); unquoted and
    ! comma-separated values as before; a single quoted string is a scalar;
    ! a comment after a group name.

    use nml, only : nml_read

    implicit none

    character(len=*), parameter :: f = "test_ds_quoted.nml"
    character(len=64) :: v2(2), v3(3), sc
    real(8) :: d(3)
    integer :: u, nfail

    nfail = 0

    open(newunit=u, file=f, status="replace", action="write")
    write(u,"(a)") "&g          ! a comment after the group name"
    write(u,"(a)") '    units  = "kg m-2 yr-1" "mm d-1"   ! two elements'
    write(u,"(a)") "    single = 'mW m-2'"
    write(u,"(a)") '    plain  = a b c'
    write(u,"(a)") '    commas = "x y", "z",  "w"'
    write(u,"(a)") '    auto   = "auto" "auto"'
    write(u,"(a)") '    nums   = 1.0, 2.5  3'
    write(u,"(a)") "/"
    close(u)

    call nml_read(f, "g", "units", v2)
    call check("quoted elements with blanks", v2(1) .eq. "kg m-2 yr-1" .and. v2(2) .eq. "mm d-1")
    call nml_read(f, "g", "single", sc)
    call check("single quoted scalar", sc .eq. "mW m-2")
    call nml_read(f, "g", "plain", v3)
    call check("unquoted elements", v3(1) .eq. "a" .and. v3(2) .eq. "b" .and. v3(3) .eq. "c")
    call nml_read(f, "g", "commas", v3)
    call check("quoted, comma-separated", v3(1) .eq. "x y" .and. v3(2) .eq. "z" .and. v3(3) .eq. "w")
    call nml_read(f, "g", "auto", v2)
    call check("quoted single words", v2(1) .eq. "auto" .and. v2(2) .eq. "auto")
    call nml_read(f, "g", "nums", d)
    call check("numbers", all(d .eq. [1.0d0, 2.5d0, 3.0d0]))

    if (nfail .gt. 0) then
        write(*,*) "test_nml_quoted: ", nfail, " FAILED"
        error stop 1
    end if
    write(*,*) "test_nml_quoted: all passed"

contains

    subroutine check(label, ok)
        character(len=*), intent(in) :: label
        logical,          intent(in) :: ok
        write(*,"(a,a)") merge("  ok   ", "  FAIL ", ok), label
        if (.not. ok) nfail = nfail + 1
    end subroutine check

end program test_nml_quoted
