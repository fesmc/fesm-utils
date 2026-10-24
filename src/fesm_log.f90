module fesm_log
    ! Shared logging and error reporting for fesm-utils.
    !
    ! * log_msg(level, who, msg [, detail]) writes one framed line
    !       LEVEL [prefix] who: msg
    !   (plus the optional `detail`, one indented line per new_line("a")),
    !   if `level` is at or above the current threshold. DEBUG/INFO go to the
    !   log unit (default output_unit), WARN/ERROR to error_unit.
    ! * fesm_warn(proc, msg [, detail]) is log_msg at LOG_WARN.
    ! * fesm_error(proc, msg [, detail] [, stat]) writes a framed multi-line
    !   diagnostic to error_unit and aborts via `error stop 1`; if `stat` is
    !   present it sets stat = 1 and returns instead, so the caller can
    !   propagate the failure. Errors are written whatever the log level.
    ! * to_str renders scalars (integer, int64, sp, dp, logical) and 1-D arrays
    !   (integer, sp, dp) compactly (g0), for building `msg` / `detail`.
    !
    ! State (level, log unit, prefix) is module-level, i.e. process-global:
    ! logging is a property of the process, not of an object. Set it once at
    ! start-up (outside OpenMP parallel regions). Dependency-free on purpose
    ! (only iso_fortran_env and precision) so any module can use it.

    use, intrinsic :: iso_fortran_env, only : output_unit, error_unit, int64
    use precision, only : sp, dp

    implicit none

    private

    integer, parameter, public :: LOG_DEBUG = 10
    integer, parameter, public :: LOG_INFO  = 20
    integer, parameter, public :: LOG_WARN  = 30
    integer, parameter, public :: LOG_ERROR = 40

    integer,            save :: log_level  = LOG_INFO
    integer,            save :: log_unit   = output_unit
    character(len=128), save :: log_prefix = ""

    interface to_str
        module procedure to_str_int, to_str_i8, to_str_sp, to_str_dp, to_str_log
        module procedure to_str_int1, to_str_sp1, to_str_dp1
    end interface

    public :: log_set_level, log_get_level, log_set_unit, log_set_prefix
    public :: log_msg, fesm_warn, fesm_error
    public :: to_str

contains

    subroutine log_set_level(level)
        ! Threshold: messages with a level below it are not written.
        implicit none
        integer, intent(in) :: level

        if (level < LOG_DEBUG .or. level > LOG_ERROR) then
            call fesm_error("log_set_level", "invalid log level.", &
                "level = "//to_str(level)//" (use LOG_DEBUG, LOG_INFO, LOG_WARN or LOG_ERROR)")
        end if
        log_level = level

    end subroutine log_set_level

    integer function log_get_level()
        implicit none
        log_get_level = log_level
    end function log_get_level

    subroutine log_set_unit(unit)
        ! Unit for DEBUG/INFO messages (WARN/ERROR always go to error_unit).
        implicit none
        integer, intent(in) :: unit
        log_unit = unit
    end subroutine log_set_unit

    subroutine log_set_prefix(prefix)
        ! Tag written in every message, e.g. a domain or instance name;
        ! "" (the default) for none.
        implicit none
        character(len=*), intent(in) :: prefix
        log_prefix = prefix
    end subroutine log_set_prefix

    subroutine log_msg(level, who, msg, detail)
        ! One framed line `LEVEL [prefix] who: msg` (+ indented detail lines).
        implicit none
        integer,          intent(in)           :: level
        character(len=*), intent(in)           :: who     ! module, procedure or instance name
        character(len=*), intent(in)           :: msg
        character(len=*), intent(in), optional :: detail  ! extra lines (join with new_line("a"))

        integer :: unit

        if (level < log_level) return

        unit = log_unit
        if (level >= LOG_WARN) unit = error_unit

        write(unit,"(a)") level_name(level)//" "//prefix_tag()//trim(who)//": "//trim(msg)
        if (present(detail)) call write_lines(unit, detail)
        if (unit == error_unit) flush(error_unit)

    end subroutine log_msg

    subroutine fesm_warn(proc, msg, detail)
        implicit none
        character(len=*), intent(in)           :: proc
        character(len=*), intent(in)           :: msg
        character(len=*), intent(in), optional :: detail

        call log_msg(LOG_WARN, proc, msg, detail)

    end subroutine fesm_warn

    subroutine fesm_error(proc, msg, detail, stat)
        ! Framed, context-rich diagnostic on error_unit, then `error stop 1`
        ! (nonzero exit; a backtrace if built with -fbacktrace). With `stat`
        ! present, set stat = 1 and return instead.
        implicit none
        character(len=*), intent(in)            :: proc    ! failing routine, e.g. "varslice_update"
        character(len=*), intent(in)            :: msg     ! what went wrong
        character(len=*), intent(in),  optional :: detail  ! extra lines (join with new_line("a"))
        integer,          intent(out), optional :: stat

        write(error_unit,"(a)") ""
        write(error_unit,"(a)") "fesm"//trim(" "//prefix_tag())//":: error in "//trim(proc)
        write(error_unit,"(a)") "    "//trim(msg)
        if (present(detail)) call write_lines(error_unit, detail)

        if (present(stat)) then
            stat = 1
            write(error_unit,"(a)") "  returned stat = 1 to the caller."
            flush(error_unit)
            return
        end if

        write(error_unit,"(a)") "  stopped by fesm_error."
        flush(error_unit)
        error stop 1

    end subroutine fesm_error

    ! ----- internal ----------------------------------------------------------

    subroutine write_lines(unit, text)
        ! Write `text` line by line (split on new_line("a")), each indented.
        implicit none
        integer,          intent(in) :: unit
        character(len=*), intent(in) :: text

        integer :: p0, p1

        p0 = 1
        do
            p1 = index(text(p0:), new_line("a"))
            if (p1 == 0) then
                write(unit,"(a)") "    "//trim(text(p0:))
                exit
            end if
            write(unit,"(a)") "    "//trim(text(p0:p0+p1-2))
            p0 = p0 + p1
        end do

    end subroutine write_lines

    function level_name(level) result(s)
        implicit none
        integer, intent(in) :: level
        character(len=5) :: s

        if (level >= LOG_ERROR) then
            s = "ERROR"
        else if (level >= LOG_WARN) then
            s = "WARN"
        else if (level >= LOG_INFO) then
            s = "INFO"
        else
            s = "DEBUG"
        end if

    end function level_name

    function prefix_tag() result(s)
        ! "[prefix] ", or "" without a prefix.
        implicit none
        character(len=:), allocatable :: s

        if (len_trim(log_prefix) > 0) then
            s = "["//trim(log_prefix)//"] "
        else
            s = ""
        end if

    end function prefix_tag

    ! ----- to_str ------------------------------------------------------------

    function to_str_int(v) result(s)
        implicit none
        integer, intent(in) :: v
        character(len=:), allocatable :: s
        character(len=32) :: buf
        write(buf,"(i0)") v
        s = trim(buf)
    end function to_str_int

    function to_str_i8(v) result(s)
        implicit none
        integer(int64), intent(in) :: v
        character(len=:), allocatable :: s
        character(len=32) :: buf
        write(buf,"(i0)") v
        s = trim(buf)
    end function to_str_i8

    function to_str_sp(v) result(s)
        implicit none
        real(sp), intent(in) :: v
        character(len=:), allocatable :: s
        character(len=64) :: buf
        write(buf,"(g0)") v
        s = trim(adjustl(buf))
    end function to_str_sp

    function to_str_dp(v) result(s)
        implicit none
        real(dp), intent(in) :: v
        character(len=:), allocatable :: s
        character(len=64) :: buf
        write(buf,"(g0)") v
        s = trim(adjustl(buf))
    end function to_str_dp

    function to_str_log(v) result(s)
        implicit none
        logical, intent(in) :: v
        character(len=:), allocatable :: s
        if (v) then
            s = "T"
        else
            s = "F"
        end if
    end function to_str_log

    function to_str_int1(v) result(s)
        ! Space-separated values; "" for an empty array.
        implicit none
        integer, intent(in) :: v(:)
        character(len=:), allocatable :: s
        integer :: i
        s = ""
        do i = 1, size(v)
            if (i > 1) s = s//" "
            s = s//to_str_int(v(i))
        end do
    end function to_str_int1

    function to_str_sp1(v) result(s)
        implicit none
        real(sp), intent(in) :: v(:)
        character(len=:), allocatable :: s
        integer :: i
        s = ""
        do i = 1, size(v)
            if (i > 1) s = s//" "
            s = s//to_str_sp(v(i))
        end do
    end function to_str_sp1

    function to_str_dp1(v) result(s)
        implicit none
        real(dp), intent(in) :: v(:)
        character(len=:), allocatable :: s
        integer :: i
        s = ""
        do i = 1, size(v)
            if (i > 1) s = s//" "
            s = s//to_str_dp(v(i))
        end do
    end function to_str_dp1

end module fesm_log
