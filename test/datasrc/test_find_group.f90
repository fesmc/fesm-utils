program test_find_group
    ! nml_find_group: the par file first, then the library list ("" skipped);
    ! no match gives stat = 1. nml_validate(..., nbad) counts unknown keys.

    use nml, only : nml_find_group, nml_validate

    implicit none

    character(len=*), parameter :: fpar = "test_ds_fg_par.nml", flib1 = "test_ds_fg_lib1.nml", &
                                   flib2 = "test_ds_fg_lib2.nml", fdef = "test_ds_fg_def.nml"
    character(len=256) :: files(4), f
    integer :: stat, nbad, nfail

    nfail = 0

    call write_file(fpar,  ["&g_par  ", "  a = 1 ", "/       ", "&g_both ", "  a = 2 ", "/       "])
    call write_file(flib1, ["&g_both ", "  a = 3 ", "/       ", "&g_lib  ", "  a = 4 ", "  zz = 1", "/       "])
    call write_file(flib2, ["&g_lib  ", "  a = 5 ", "/       ", "&g_two  ", "/       "])
    call write_file(fdef,  ["&tmpl   ", "  a = 0 ", "/       "])

    files = [character(len=256) :: fpar, "", flib1, flib2]

    f = nml_find_group(files, "g_par")
    call check("par only", trim(f) .eq. fpar, nfail)
    f = nml_find_group(files, "g_both")
    call check("par over library", trim(f) .eq. fpar, nfail)
    f = nml_find_group(files, "g_lib")
    call check("first library file", trim(f) .eq. flib1, nfail)
    f = nml_find_group(files, "g_two", what="test")
    call check("empty group in the second library file", trim(f) .eq. flib2, nfail)
    f = nml_find_group(files, "g_none", what="test", stat=stat)
    call check("absent: stat = 1", stat .eq. 1 .and. len_trim(f) .eq. 0, nfail)
    f = nml_find_group(files, "g_lib", stat=stat)
    call check("found: stat = 0", stat .eq. 0, nfail)

    call nml_validate(flib1, fdef, "g_lib", "tmpl", halt_on_error=.false., nbad=nbad)
    call check("nml_validate nbad = 1", nbad .eq. 1, nfail)
    call nml_validate(flib1, fdef, "g_both", "tmpl", halt_on_error=.false., nbad=nbad)
    call check("nml_validate nbad = 0", nbad .eq. 0, nfail)

    if (nfail .gt. 0) then
        write(*,*) "test_find_group: ", nfail, " FAILED"
        error stop 1
    end if
    write(*,*) "test_find_group: all passed"

contains

    subroutine write_file(fname, lines)
        character(len=*), intent(in) :: fname, lines(:)
        integer :: u, i
        open(newunit=u, file=fname, status="replace", action="write")
        do i = 1, size(lines)
            write(u,"(a)") trim(lines(i))
        end do
        close(u)
    end subroutine write_file

    subroutine check(label, ok, nfail)
        character(len=*), intent(in)    :: label
        logical,          intent(in)    :: ok
        integer,          intent(inout) :: nfail
        write(*,"(a,a)") merge("  ok   ", "  FAIL ", ok), label
        if (.not. ok) nfail = nfail + 1
    end subroutine check

end program test_find_group
