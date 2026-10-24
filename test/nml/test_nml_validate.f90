program test_nml_validate
    ! nml_group_keys / nml_group_names list a file's keys and groups (block
    ! and flat forms, each once, in file order); nml_validate against several
    ! templates accepts a key that any template declares and counts the rest
    ! (nbad); the single-template form is unchanged.

    use nml, only : nml_validate, nml_group_keys, nml_group_names, nml_key_len, nml_group_len

    implicit none

    character(len=*), parameter :: fpar = "test_nml_par.nml", fcore = "test_nml_core.nml", &
                                   fback = "test_nml_back.nml"
    character(len=nml_key_len),   allocatable :: keys(:)
    character(len=nml_group_len), allocatable :: groups(:)
    integer :: nbad, nfail

    nfail = 0

    call write_file(fpar,  [character(len=64) :: &
        "&inst   ! an instance group",            &
        "  dt       = 1.0",                        &
        "  backend  = ""yelmo""",                  &
        "  rf_const = 1e-3",                       &
        "  dt       = 2.0",                        &
        "/",                                       &
        "&empty",                                  &
        "/",                                       &
        "inst.zz = 1   flat.a = 2"])
    call write_file(fcore, [character(len=64) :: &
        "&instance",                               &
        "  dt      = 1.0",                         &
        "  backend = ""none""",                    &
        "/"])
    call write_file(fback, [character(len=64) :: &
        "&ice_yelmo",                              &
        "  rf_const = 1e-3",                       &
        "/",                                       &
        "&other",                                  &
        "  zz = 0",                                &
        "/"])

    keys = nml_group_keys(fpar, "inst")
    call check("group keys: dt backend rf_const zz (each once)", size(keys) == 4, nfail)
    if (size(keys) == 4) call check("group keys: file order", &
        keys(1) == "dt" .and. keys(2) == "backend" .and. keys(3) == "rf_const" .and. keys(4) == "zz", nfail)
    keys = nml_group_keys(fpar, "empty")
    call check("group keys: empty group", size(keys) == 0, nfail)
    keys = nml_group_keys(fpar, "absent")
    call check("group keys: absent group", size(keys) == 0, nfail)

    groups = nml_group_names(fpar)
    call check("group names: inst empty flat", size(groups) == 3, nfail)
    if (size(groups) == 3) call check("group names: comment after the group name dropped", &
        groups(1) == "inst" .and. groups(2) == "empty" .and. groups(3) == "flat", nfail)

    ! Two templates: core keys + backend keys; zz is in neither
    call nml_validate(fpar, "inst", [character(len=32) :: fcore, fback], &
                      [character(len=32) :: "instance", "ice_yelmo"], halt_on_error=.false., nbad=nbad)
    call check("multi: nbad = 1 (zz)", nbad == 1, nfail)

    ! A third template declaring zz
    call nml_validate(fpar, "inst", [character(len=32) :: fcore, fback, fback], &
                      [character(len=32) :: "instance", "ice_yelmo", "other"], halt_on_error=.false., nbad=nbad)
    call check("multi: nbad = 0 (any template declares the key)", nbad == 0, nfail)

    ! One template only
    call nml_validate(fpar, "inst", [character(len=32) :: fcore], &
                      [character(len=32) :: "instance"], halt_on_error=.false., nbad=nbad)
    call check("multi, one template: nbad = 2 (rf_const zz)", nbad == 2, nfail)

    ! defaults_groups absent: every template is searched under the user group
    call nml_validate(fpar, "flat", [character(len=32) :: fpar], halt_on_error=.false., nbad=nbad)
    call check("multi, template group = user group: nbad = 0", nbad == 0, nfail)

    ! Single-template form: duplicates counted per line, as before
    call nml_validate(fpar, fcore, "inst", "instance", halt_on_error=.false., nbad=nbad)
    call check("single: nbad = 2", nbad == 2, nfail)

    if (nfail > 0) then
        write(*,*) "test_nml_validate: ", nfail, " FAILED"
        error stop 1
    end if
    write(*,*) "test_nml_validate: all passed"

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

end program test_nml_validate
