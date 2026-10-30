program test_fields_restart_diag
    ! fesm_fields restart with diagnostics (restart = .false.): a state with
    ! two state fields (a published INSTANT, a MEAN with an open window) and a
    ! published diagnostic is written; the file holds the state fields only.
    ! It is read into a state with one diagnostic more and into one with the
    ! diagnostic dropped: the state fields are bit-identical, a diagnostic
    ! keeps its value (unpublished after fields_lock). A file of an older
    ! declaration that holds the diagnostic as a state field reads into the
    ! new declaration: the diagnostic is ignored. A state field missing from
    ! the file, or one the declaration lacks, still stops:
    ! abort/fields_abort_restart_diag_state, abort/fields_abort_restart_extra.

    use precision, only : dp
    use ncio, only : nc_exists_var, nc_read_attr
    use fesm_fields

    implicit none

    character(len=*), parameter :: fname = "test_fields_restart_diag.nc", fold = "test_fields_restart_diag_old.nc"
    type(field_state_t) :: a, more, less, old, b
    integer  :: nfail, i, j, ival, h
    real(dp) :: x2(4,3)
    character(len=64) :: str

    nfail = 0

    do j = 1, 3
        do i = 1, 4
            x2(i,j) = sin(real(i + 7*j, dp)) * 1.0e3_dp
        end do
    end do

    ! Written: H (INSTANT), smb (MEAN, window open), resid (diagnostic).
    call declare(a, 1)
    call fields_put(a, 1, x2, t1=100.0_dp)
    call fields_put(a, 2, x2/3.0_dp, t0=99.0_dp, t1=100.0_dp)
    call fields_finalize_interval(a, 100.0_dp)
    call fields_put(a, 2, x2/7.0_dp, t0=100.0_dp, t1=100.25_dp)
    call fields_put(a, 3, 0.125_dp, t1=100.0_dp)
    call fields_restart_write(a, fname)

    call nc_read_attr(fname, "fesm_fields_n", ival)
    call check("global n: the state fields only", ival == 2)
    str = ""
    call nc_read_attr(fname, "fesm_fields_names", str)
    call check("global names: H,smb", trim(str) == "H,smb")
    call check("no variable for the diagnostic", .not. nc_exists_var(fname, "resid"))

    ! One diagnostic more (resid, extra).
    call declare(more, 2)
    call fields_restart_read(more, fname)
    call check("more: H bit-identical", all(more%f(1)%v == a%f(1)%v) .and. more%f(1)%published .and. &
               more%f(1)%t1 == 100.0_dp)
    call check("more: smb v, acc, window bit-identical", all(more%f(2)%v == a%f(2)%v) .and. &
               all(more%f(2)%acc == a%f(2)%acc) .and. more%f(2)%wsum == a%f(2)%wsum .and. &
               more%f(2)%ta0 == a%f(2)%ta0 .and. more%f(2)%ta1 == a%f(2)%ta1)
    h = fields_handle(more, "resid")
    call check("more: the diagnostics keep their value (unpublished)", .not. more%f(h)%published .and. &
               all(more%f(h)%v == 0.0_dp) .and. .not. more%f(fields_handle(more, "extra"))%published)

    ! The diagnostic dropped.
    call declare(less, 0)
    call fields_restart_read(less, fname)
    call check("less: H, smb bit-identical", all(less%f(1)%v == a%f(1)%v) .and. all(less%f(2)%v == a%f(2)%v) .and. &
               all(less%f(2)%acc == a%f(2)%acc))

    ! An older declaration: resid a state field, so in the file; read with resid a diagnostic.
    call declare(old, 1, old=.true.)
    call fields_put(old, 1, x2, t1=100.0_dp)
    call fields_put(old, 2, x2/3.0_dp, t0=99.0_dp, t1=100.0_dp)
    call fields_finalize_interval(old, 100.0_dp)
    call fields_put(old, 3, 0.125_dp, t1=100.0_dp)
    call fields_restart_write(old, fold)
    call check("old file holds resid", nc_exists_var(fold, "resid"))
    call declare(b, 1)
    call fields_restart_read(b, fold)
    call check("old file: H, smb bit-identical", all(b%f(1)%v == old%f(1)%v) .and. all(b%f(2)%v == old%f(2)%v))
    call check("old file: resid ignored (unpublished)", .not. b%f(3)%published .and. all(b%f(3)%v == 0.0_dp))

    call execute_command_line("rm -f "//fname//" "//fold)

    if (nfail > 0) then
        write(*,*) "test_fields_restart_diag: ", nfail, " check(s) failed."
        stop 1
    end if
    write(*,*) "PASS: test_fields_restart_diag"

contains

    subroutine declare(st, ndiag, old)
        ! The state fields H, smb, then ndiag diagnostics of resid, extra;
        ! old: resid is a state field (an older declaration).
        type(field_state_t), intent(out) :: st
        integer,             intent(in)  :: ndiag
        logical,             intent(in), optional :: old
        integer :: hh
        logical :: resid_state

        resid_state = .false.
        if (present(old)) resid_state = old
        call fields_init(st, "test:restart_diag")
        call fields_add(st, "H",   "m",           FIELD_INSTANT, FIELD_EXTENSIVE, hh, grid="g1", dims=[4,3])
        call fields_add(st, "smb", "kg m-2 yr-1", FIELD_MEAN,    FIELD_EXTENSIVE, hh, grid="g1", dims=[4,3])
        if (ndiag >= 1) call fields_add(st, "resid", "1", FIELD_INSTANT, FIELD_INTENSIVE, hh, restart=resid_state)
        if (ndiag >= 2) call fields_add(st, "extra", "1", FIELD_INSTANT, FIELD_INTENSIVE, hh, restart=.false.)
        call fields_lock(st)

    end subroutine declare

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

end program test_fields_restart_diag
