program test_fields_declare
    ! fesm_fields declaration and lock:
    !   * handles are 1..n in declaration order (past the initial capacity);
    !   * fields_has true/false, fields_handle by name;
    !   * rank, shape, zdim (given and default) and ukey are recorded;
    !   * buffers are allocated only by fields_lock, zeroed, with acc only for
    !     MEAN/ACCUM and masks only for masked fields; size(f) == n after lock;
    !   * fields_describe writes a header and one line per field.

    use, intrinsic :: iso_fortran_env, only : output_unit
    use precision, only : dp
    use fesm_log
    use fesm_fields

    implicit none

    type(field_state_t) :: st
    integer :: h(10), i, nfail, u
    character(len=8) :: nm
    character(len=256), allocatable :: lines(:)
    logical :: none

    nfail = 0

    call fields_init(st, "test:declare")
    call fields_add(st, "bsl",   "m",          FIELD_INSTANT, FIELD_INTENSIVE, h(1))
    call fields_add(st, "H_ice", "m",          FIELD_INSTANT, FIELD_EXTENSIVE, h(2), grid="g1", dims=[4,3])
    call fields_add(st, "smb",   "mm w.e./yr", FIELD_MEAN,    FIELD_EXTENSIVE, h(3), grid="g1", dims=[4,3], &
                    masked=.true.)
    call fields_add(st, "tas",   "K",          FIELD_MEAN,    FIELD_INTENSIVE, h(4), grid="g1", dims=[4,3,12], &
                    zdim="month")
    call fields_add(st, "to",    "K",          FIELD_INSTANT, FIELD_INTENSIVE, h(5), grid="g2", dims=[2,2,5])
    do i = 6, 10
        write(nm, "(a,i0)") "x", i
        call fields_add(st, trim(nm), "1", FIELD_ACCUM, FIELD_EXTENSIVE, h(i))
    end do

    ! ----- handles -------------------------------------------------------------
    call check("handles 1..n in order", all(h == [(i, i = 1, 10)]))
    call check("n = 10 (past the initial capacity)", st%n == 10)
    call check("fields_has(smb)", fields_has(st, "smb"))
    call check("not fields_has(SMB): exact match", .not. fields_has(st, "SMB"))
    call check("not fields_has(none)", .not. fields_has(st, "none"))
    call check("fields_handle(tas) = 4", fields_handle(st, "tas") == 4)
    call check("fields_handle(x10) = 10", fields_handle(st, "x10") == 10)

    ! ----- metadata -------------------------------------------------------------
    call check("bsl rank 0, 1x1x1, no grid", st%f(1)%rank == 0 .and. st%f(1)%nx == 1 .and. &
               st%f(1)%ny == 1 .and. st%f(1)%nz == 1 .and. st%f(1)%grid == "")
    call check("H_ice rank 2, 4x3x1, grid g1", st%f(2)%rank == 2 .and. st%f(2)%nx == 4 .and. &
               st%f(2)%ny == 3 .and. st%f(2)%nz == 1 .and. st%f(2)%grid == "g1")
    call check("tas rank 3, nz 12, zdim month", st%f(4)%rank == 3 .and. st%f(4)%nz == 12 .and. &
               st%f(4)%zdim == "month")
    call check("to zdim defaults to z", st%f(5)%zdim == "z" .and. st%f(5)%nz == 5)
    call check("smb units kept as declared", st%f(3)%units == "mm w.e./yr")
    call check("smb ukey canonical", st%f(3)%ukey == "kgm-2yr-1")
    call check("tas ukey", st%f(4)%ukey == "k")
    call check("smb masked, tsem, cons", st%f(3)%masked .and. st%f(3)%tsem == FIELD_MEAN .and. &
               st%f(3)%cons == FIELD_EXTENSIVE)
    call check("H_ice not masked", .not. st%f(2)%masked)

    ! ----- no buffers before lock -------------------------------------------------
    none = .true.
    do i = 1, st%n
        none = none .and. .not. (allocated(st%f(i)%v) .or. allocated(st%f(i)%acc) .or. &
                                 allocated(st%f(i)%valid) .or. allocated(st%f(i)%acc_valid))
    end do
    call check("no buffers before lock", none)
    call check("not locked", .not. st%locked)

    ! ----- describe (before lock) -------------------------------------------------
    open(newunit=u, status="scratch", action="readwrite", form="formatted")
    call log_set_unit(u)
    call fields_describe(st)
    call read_back(u, lines)
    close(u)
    call log_set_unit(output_unit)
    call check("describe: header + 10 lines", size(lines) == 11)
    if (size(lines) >= 1) then
        call check("describe: header", index(lines(1), "INFO  fesm_fields: [test:declare] 10 fields, ") == 1)
        call check("describe: not locked", index(lines(1), "(not locked)") > 0)
    end if
    if (size(lines) >= 4) then
        call check("describe: smb line", index(lines(4), "smb") > 0 .and. index(lines(4), "MEAN") > 0 .and. &
                   index(lines(4), "unpublished") > 0 .and. index(lines(4), "masked") > 0)
    end if

    ! ----- lock -------------------------------------------------------------------
    call fields_lock(st)
    call check("locked", st%locked)
    call check("size(f) = n after lock", size(st%f) == st%n)

    call check("v allocated with the declared shape", &
               all(shape(st%f(1)%v) == [1,1,1]) .and. all(shape(st%f(2)%v) == [4,3,1]) .and. &
               all(shape(st%f(4)%v) == [4,3,12]) .and. all(shape(st%f(5)%v) == [2,2,5]))
    call check("v zeroed", all(st%f(2)%v == 0.0_dp) .and. all(st%f(4)%v == 0.0_dp))
    call check("acc only for MEAN/ACCUM", allocated(st%f(3)%acc) .and. allocated(st%f(4)%acc) .and. &
               allocated(st%f(6)%acc) .and. .not. allocated(st%f(1)%acc) .and. .not. allocated(st%f(2)%acc))
    call check("acc zeroed", all(st%f(4)%acc == 0.0_dp))
    call check("masks only for masked fields", allocated(st%f(3)%valid) .and. allocated(st%f(3)%acc_valid) .and. &
               .not. allocated(st%f(2)%valid) .and. .not. allocated(st%f(4)%acc_valid))
    call check("mask shape (nx,ny)", all(shape(st%f(3)%valid) == [4,3]))
    call check("acc_valid starts .true. (empty window)", all(st%f(3)%acc_valid))
    call check("nothing published", .not. any([(st%f(i)%published, i = 1, st%n)]))
    call check("handles survive the lock", fields_handle(st, "tas") == 4)

    if (nfail > 0) then
        write(*,*) "test_fields_declare: ", nfail, " check(s) failed."
        stop 1
    end if
    write(*,*) "PASS: test_fields_declare"

contains

    subroutine check(label, ok)
        character(len=*), intent(in) :: label
        logical,          intent(in) :: ok
        if (ok) then
            write(*,"(a)") "  ok    "//label
        else
            write(*,"(a)") "  FAIL  "//label
            nfail = nfail + 1
        end if
    end subroutine check

    subroutine read_back(unit, lines)
        ! All lines written to `unit` so far.
        integer, intent(in) :: unit
        character(len=256), allocatable, intent(out) :: lines(:)
        character(len=256) :: buf
        integer :: n, ios
        rewind(unit)
        n = 0
        do
            read(unit,"(a)",iostat=ios) buf
            if (ios /= 0) exit
            n = n + 1
        end do
        allocate(lines(n))
        rewind(unit)
        do n = 1, size(lines)
            read(unit,"(a)") lines(n)
        end do
    end subroutine read_back

end program test_fields_declare
