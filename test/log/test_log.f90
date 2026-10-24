program test_log
    ! Exercise fesm_log:
    !   * level threshold (default INFO; DEBUG hidden until enabled);
    !   * the line format `LEVEL [prefix] who: msg` and indented detail lines,
    !     read back from a scratch unit set with log_set_unit;
    !   * to_str on every kind (int, int64, sp, dp, logical, 1-D int/sp/dp);
    !   * fesm_warn (to error_unit; shown, and suppressed at LOG_ERROR);
    !   * fesm_error with `stat` present: frame on error_unit, stat = 1, no stop.
    ! Self-contained: no files or namelists are needed.

    use, intrinsic :: iso_fortran_env, only : output_unit, int64
    use precision, only : sp, dp
    use fesm_log

    implicit none

    integer :: nfail, u, stat
    character(len=256), allocatable :: lines(:)
    integer, allocatable :: iempty(:)

    nfail = 0

    ! ----- levels ------------------------------------------------------------
    call check("default level is LOG_INFO", log_get_level() == LOG_INFO)
    call check("levels ordered", LOG_DEBUG < LOG_INFO .and. LOG_INFO < LOG_WARN .and. LOG_WARN < LOG_ERROR)

    open(newunit=u, status="scratch", action="readwrite", form="formatted")
    call log_set_unit(u)

    call log_msg(LOG_DEBUG, "test_log", "hidden debug")
    call log_msg(LOG_INFO,  "test_log", "first info")
    call log_set_level(LOG_DEBUG)
    call check("log_set_level(LOG_DEBUG)", log_get_level() == LOG_DEBUG)
    call log_msg(LOG_DEBUG, "test_log", "shown debug")

    ! ----- prefix and detail -------------------------------------------------
    call log_set_prefix("grl")
    call log_msg(LOG_INFO, "varslice_update", "with prefix", &
        "a = 1"//new_line("a")//"b = 2")
    call log_set_prefix("")
    call log_msg(LOG_INFO, "test_log", "no prefix")

    call read_back(u, lines)
    call check("4 messages + 2 detail lines written", size(lines) == 6)
    if (size(lines) == 6) then
        call check("debug filtered at INFO",  index(lines(1), "hidden") == 0)
        call check("INFO line format",        lines(1) == "INFO  test_log: first info")
        call check("DEBUG shown at DEBUG",    lines(2) == "DEBUG test_log: shown debug")
        call check("prefix tag",              lines(3) == "INFO  [grl] varslice_update: with prefix")
        call check("detail line 1 indented",  lines(4) == "    a = 1")
        call check("detail line 2 indented",  lines(5) == "    b = 2")
        call check("prefix cleared",          lines(6) == "INFO  test_log: no prefix")
    end if
    close(u)
    call log_set_unit(output_unit)
    call log_set_level(LOG_INFO)

    ! ----- to_str ------------------------------------------------------------
    call check("to_str(int)",       to_str(42) == "42")
    call check("to_str(-int)",      to_str(-7) == "-7")
    call check("to_str(int64)",     to_str(123456789012_int64) == "123456789012")
    call check("to_str(sp)",        round_sp(to_str(1.5_sp)) == 1.5_sp)
    call check("to_str(dp)",        round_dp(to_str(0.1_dp)) == 0.1_dp)
    call check("to_str(sp) compact", index(to_str(2.0_sp), " ") == 0)
    call check("to_str(.true.)",    to_str(.true.)  == "T")
    call check("to_str(.false.)",   to_str(.false.) == "F")
    call check("to_str(int(:))",    to_str([1, 2, 3]) == "1 2 3")
    allocate(iempty(0))
    call check("to_str(empty)",     to_str(iempty) == "")
    call check("to_str(sp(:))",     to_str([1.0_sp, 2.5_sp]) == to_str(1.0_sp)//" "//to_str(2.5_sp))
    call check("to_str(dp(:))",     to_str([1.0_dp, -2.5_dp]) == to_str(1.0_dp)//" "//to_str(-2.5_dp))

    ! ----- warnings and errors (to error_unit; shown in the test output) -----
    write(*,"(a)") "  (expected on stderr: one WARN line, then one error frame with stat)"
    call fesm_warn("test_log", "a warning", "detail = "//to_str(3))
    call log_set_level(LOG_ERROR)
    call fesm_warn("test_log", "this warning must not appear")
    call log_set_level(LOG_INFO)

    stat = 0
    call fesm_error("test_log", "an error with stat present.", &
        "value = "//to_str(1.0e10_dp)//new_line("a")//"flag  = "//to_str(.true.), stat=stat)
    call check("fesm_error(stat=) returns stat = 1", stat == 1)

    ! ----- summary -----------------------------------------------------------
    write(*,"(a)") ""
    if (nfail == 0) then
        write(*,"(a)") " All fesm_log tests PASSED."
    else
        write(*,"(a,i0,a)") " fesm_log tests FAILED: ", nfail, " case(s)."
    end if
    if (nfail > 0) error stop 1

contains

    subroutine check(name, ok)
        character(len=*), intent(in) :: name
        logical,          intent(in) :: ok
        if (ok) then
            write(*,"(a)") "  PASS  "//name
        else
            write(*,"(a)") "  FAIL  "//name
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

    real(sp) function round_sp(s)
        character(len=*), intent(in) :: s
        read(s,*) round_sp
    end function round_sp

    real(dp) function round_dp(s)
        character(len=*), intent(in) :: s
        read(s,*) round_dp
    end function round_dp

end program test_log
