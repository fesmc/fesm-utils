module fesm_clock
    ! One run clock: integer ticks inside, real(dp) years at every interface,
    ! a declared calendar. Design and rationale: docs/dev/clock-design.md.
    !
    ! * clock_t     model time n [ticks], N = ticks_per_year (default 36000)
    !               ticks per year; views time, time_rel, time_cal and
    !               time_elapsed are exact affine maps of n, never sums.
    ! * cadence_t   "every dt yr" on the absolute grid phase + k*dt, or every
    !               step; due once per grid time (mark), with exact intervals.
    ! * schedule_t  output/restart times (none | const | times | file), on the
    !               step grid, with a cursor.
    !
    ! Every year value is converted once to ticks and must be a whole number
    ! of ticks (4-ulp test); every grid relation (cadence on step grid,
    ! schedule times on step grid, restart times on every cadence grid) is
    ! checked at init/start, so a run that starts can never miss a time.
    ! Errors go through fesm_error; init, start, to_ticks and restart_read take
    ! an optional `stat` (set to 1 on failure instead of stopping).
    !
    ! Free procedures throughout: no type-bound procedures.

    use, intrinsic :: iso_fortran_env, only : int64
    use precision,      only : dp
    use fesm_log,       only : LOG_INFO, log_msg, fesm_warn, fesm_error, to_str
    use nml,            only : nml_read, nml_has_param
    use ncio,           only : nc_create, nc_open, nc_close, nc_write, nc_read, nc_write_attr, nc_read_attr, &
                               nc_exists_var, nc_exists_attr
    use phys_constants, only : sec_year_360d, sec_year_365d, sec_year_julian, &
                               sec_year_gregorian, sec_year_tropical

    implicit none

    private

    integer(int64), parameter, public :: CLOCK_TICKS_PER_YEAR = 36000_int64
    integer,        parameter, public :: CADENCE_OFF   = 0
    integer,        parameter, public :: CADENCE_STEP  = 1
    integer,        parameter, public :: CADENCE_EVERY = 2

    ! |n| < 2^53: ticks convert exactly to real(dp)
    integer(int64), parameter :: N_EXACT = 9007199254740992_int64

    integer, parameter :: MONTH_END_365(0:12) = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365]

    type :: clock_entry_t                       ! registry row: describe + restart checks
        character(len=32) :: name   = ""
        character(len=8)  :: kind   = ""        ! cadence | output | restart
        character(len=8)  :: method = ""        ! schedules: none | const | times | file
        integer           :: mode   = CADENCE_OFF
        integer(int64)    :: dt = 0, phase = 0
        integer(int64), allocatable :: n(:)     ! schedule times [ticks]
    end type clock_entry_t

    type, public :: clock_t
        integer(int64)     :: tpy      = CLOCK_TICKS_PER_YEAR  ! N, ticks per year
        integer(int64)     :: n        = 0      ! model time [ticks]
        integer(int64)     :: n_init   = 0      ! this run (chunk)
        integer(int64)     :: n_end    = 0
        integer(int64)     :: dt       = 0      ! dtt [ticks]
        integer(int64)     :: k        = 0      ! steps since n_init
        integer(int64)     :: n_origin = 0      ! experiment start, kept across chunks
        character(len=8)   :: method   = ""     ! rel | cal | const
        integer(int64)     :: n_ref    = 0      ! time_ref [ticks]
        integer(int64)     :: n_const  = 0      ! frozen time_rel (const) [ticks]
        character(len=16)  :: calendar = ""     ! 360_day | 365_day | none
        integer            :: days_per_year = 0
        character(len=16)  :: year     = ""     ! phys_constants name
        real(dp)           :: sec_year = 0.0_dp ! [s]
        logical            :: started   = .false.
        logical            :: continued = .false.   ! restored at n_init from a bundle
        logical            :: branched  = .false.   ! state read from a bundle at another time
        character(len=512) :: restart_file = ""
        type(clock_entry_t), allocatable :: reg(:)
    end type clock_t

    type, public :: cadence_t
        character(len=32) :: name   = ""
        integer           :: mode   = CADENCE_OFF
        integer(int64)    :: dt     = 0         ! [ticks] (CADENCE_EVERY)
        integer(int64)    :: phase  = 0         ! [ticks], 0 <= phase < dt
        integer(int64)    :: n_last = 0         ! last mark [ticks]
        integer           :: id     = 0         ! row in clk%reg
    end type cadence_t

    type, public :: schedule_t
        character(len=32) :: name    = ""
        character(len=8)  :: method  = "none"   ! none | const | times | file
        logical           :: restart = .false.
        integer(int64), allocatable :: n(:)     ! sorted, unique, on the step grid [ticks]
        integer           :: k  = 1             ! cursor: first n(k) >= clk%n
        integer           :: id = 0
        integer(int64)    :: n_skip = -huge(1_int64)  ! last time a skipped bundle was logged
    end type schedule_t

    interface clock_init
        module procedure clock_init_nml, clock_init_args
    end interface clock_init

    interface cadence_init
        module procedure cadence_init_dt, cadence_init_step, cadence_init_per_year
    end interface cadence_init

    interface schedule_init
        module procedure schedule_init_nml, schedule_init_args
    end interface schedule_init

    interface read_key
        module procedure read_key_str, read_key_dp, read_key_int
    end interface read_key

    public :: clock_init, clock_start, clock_advance, clock_finish, clock_finished
    public :: clock_step_index, clock_time, clock_time_rel, clock_time_cal, clock_time_elapsed
    public :: clock_time_init, clock_time_end, clock_dtt
    public :: clock_to_ticks, clock_to_years, clock_is_representable, clock_restart_allowed
    public :: clock_year_fraction, clock_day_of_year, clock_month, clock_day_of_month, clock_sec_day
    public :: clock_describe, clock_restart_write, clock_restart_read
    public :: cadence_init, cadence_due, cadence_mark, cadence_interval
    public :: cadence_t_last, cadence_t_next, cadence_dt_years
    public :: cadence_restart_write, cadence_restart_read
    public :: schedule_init, schedule_due, schedule_t_next, schedule_times, schedule_size
    public :: is_multiple_of_all

contains

    ! =========================================================================
    ! Clock: init, run control
    ! =========================================================================

    subroutine clock_init_nml(clk, filename, group, defaults_file, stat)
        ! Keys: time_init, time_end, dtt, calendar, year (required; no defaults
        ! here), method (rel), time_ref (1950), ticks_per_year (36000), and for
        ! method = const exactly one of time_const, time_const_cal. A key may come
        ! from `defaults_file` (same group) instead of `filename`.
        type(clock_t),    intent(inout)         :: clk
        character(len=*), intent(in)            :: filename, group
        character(len=*), intent(in),  optional :: defaults_file
        integer,          intent(out), optional :: stat

        character(len=*), parameter :: proc = "fesm_clock::clock_init"
        character(len=*), parameter :: required(5) = [character(len=9) :: &
            "time_init", "time_end ", "dtt      ", "calendar ", "year     "]
        character(len=32) :: method, calendar, year
        real(dp) :: time_init, time_end, dtt, time_ref, tc, tcc
        integer  :: tpy, i
        logical  :: has_tc, has_tcc

        if (present(stat)) stat = 0

        do i = 1, size(required)
            if (.not. has_key(filename, group, trim(required(i)), defaults_file)) then
                call fesm_error(proc, 'key "'//trim(required(i))//'" not found in group "'//trim(group)// &
                    '" of '//trim(filename), "fesm_clock has no default for it.", stat)
                return
            end if
        end do

        call read_key(filename, group, "time_init", time_init, defaults_file)
        call read_key(filename, group, "time_end",  time_end,  defaults_file)
        call read_key(filename, group, "dtt",       dtt,       defaults_file)
        call read_key(filename, group, "calendar",  calendar,  defaults_file)
        call read_key(filename, group, "year",      year,      defaults_file)

        method = "rel"
        if (has_key(filename, group, "method", defaults_file)) &
            call read_key(filename, group, "method", method, defaults_file)
        time_ref = 1950.0_dp
        if (has_key(filename, group, "time_ref", defaults_file)) &
            call read_key(filename, group, "time_ref", time_ref, defaults_file)
        tpy = int(CLOCK_TICKS_PER_YEAR)
        if (has_key(filename, group, "ticks_per_year", defaults_file)) &
            call read_key(filename, group, "ticks_per_year", tpy, defaults_file)

        has_tc  = .false.
        has_tcc = .false.
        tc  = 0.0_dp
        tcc = 0.0_dp
        if (trim(method) == "const") then
            has_tc  = has_key(filename, group, "time_const",     defaults_file)
            has_tcc = has_key(filename, group, "time_const_cal", defaults_file)
            if (has_tc)  call read_key(filename, group, "time_const",     tc,  defaults_file)
            if (has_tcc) call read_key(filename, group, "time_const_cal", tcc, defaults_file)
        end if

        call clock_setup(clk, time_init, time_end, dtt, calendar, year, method, &
                         has_tc, tc, has_tcc, tcc, time_ref, int(tpy, int64), proc, stat)

    end subroutine clock_init_nml

    subroutine clock_init_args(clk, time_init, time_end, dtt, calendar, year, method, &
                               time_const, time_const_cal, time_ref, ticks_per_year, stat)
        ! As clock_init_nml, from arguments (tests, hosts without a namelist).
        ! time_const / time_const_cal are used with method = "const" only.
        type(clock_t),    intent(inout)         :: clk
        real(dp),         intent(in)            :: time_init, time_end, dtt
        character(len=*), intent(in)            :: calendar, year
        character(len=*), intent(in),  optional :: method
        real(dp),         intent(in),  optional :: time_const, time_const_cal, time_ref
        integer,          intent(in),  optional :: ticks_per_year
        integer,          intent(out), optional :: stat

        character(len=*), parameter :: proc = "fesm_clock::clock_init"
        character(len=32) :: meth
        real(dp)       :: tc, tcc, tref
        integer(int64) :: tpy

        meth = "rel"
        if (present(method)) meth = method
        tref = 1950.0_dp
        if (present(time_ref)) tref = time_ref
        tpy = CLOCK_TICKS_PER_YEAR
        if (present(ticks_per_year)) tpy = int(ticks_per_year, int64)
        tc  = 0.0_dp
        tcc = 0.0_dp
        if (present(time_const))     tc  = time_const
        if (present(time_const_cal)) tcc = time_const_cal

        call clock_setup(clk, time_init, time_end, dtt, calendar, year, meth, &
                         present(time_const) .and. trim(meth) == "const", tc, &
                         present(time_const_cal) .and. trim(meth) == "const", tcc, &
                         tref, tpy, proc, stat)

    end subroutine clock_init_args

    subroutine clock_setup(clk, time_init, time_end, dtt, calendar, year, method, &
                           has_tc, tc, has_tcc, tcc, time_ref, tpy, proc, stat)
        ! Shared by both clock_init forms: validate and set every field.
        type(clock_t),    intent(inout)         :: clk
        real(dp),         intent(in)            :: time_init, time_end, dtt, tc, tcc, time_ref
        character(len=*), intent(in)            :: calendar, year, method, proc
        logical,          intent(in)            :: has_tc, has_tcc
        integer(int64),   intent(in)            :: tpy
        integer,          intent(out), optional :: stat

        type(clock_t)  :: fresh
        integer(int64) :: n0

        if (present(stat)) stat = 0
        clk = fresh

        if (tpy <= 0 .or. tpy > 1000000000_int64) then
            call fesm_error(proc, "ticks_per_year = "//to_str(tpy)//" is out of range.", &
                "use a positive value below 1e9 (default "//to_str(CLOCK_TICKS_PER_YEAR)//").", stat)
            return
        end if
        clk%tpy = tpy

        select case (trim(method))
        case ("rel", "cal", "const")
            clk%method = trim(method)
        case default
            call fesm_error(proc, 'method = "'//trim(method)//'" is not recognised.', &
                "use rel, cal or const.", stat)
            return
        end select

        call set_calendar(clk, calendar, proc, stat)
        if (failed(stat)) return
        call set_year(clk, year, proc, stat)
        if (failed(stat)) return

        call conv_ticks(tpy, time_init, proc, "", "time_init", clk%n_init, stat)
        if (failed(stat)) return
        call conv_ticks(tpy, time_end, proc, "", "time_end", clk%n_end, stat)
        if (failed(stat)) return
        call conv_ticks(tpy, dtt, proc, "", "dtt", clk%dt, stat)
        if (failed(stat)) return
        call conv_ticks(tpy, time_ref, proc, "", "time_ref", clk%n_ref, stat)
        if (failed(stat)) return

        if (clk%dt <= 0) then
            call fesm_error(proc, "dtt must be positive.", "dtt = "//yr_str(dtt)//" yr", stat)
            clk%dt = 0
            return
        end if
        if (clk%n_end <= clk%n_init) then
            call fesm_error(proc, "time_end must be after time_init.", &
                "time_init = "//yr_str(time_init)//" yr, time_end = "//yr_str(time_end)//" yr", stat)
            clk%dt = 0
            return
        end if
        if (modulo(clk%n_ref, tpy) /= 0) then
            n0 = floor_div(clk%n_ref, tpy)*tpy
            call fesm_error(proc, "time_ref = "//yr_str(time_ref)//" yr is not a whole number of years.", &
                "nearest: "//pair_str(n0, n0 + tpy, tpy), stat)
            clk%dt = 0
            return
        end if
        if (modulo(clk%n_end - clk%n_init, clk%dt) /= 0) then
            n0 = clk%n_init + floor_div(clk%n_end - clk%n_init, clk%dt)*clk%dt
            call fesm_error(proc, "time_end = "//yr_str(time_end)//" yr is not on the step grid "// &
                "time_init + k*dtt.", "time_init = "//yr_str(time_init)//" yr, dtt = "//yr_str(dtt)// &
                " yr"//new_line("a")//"nearest valid time_end: "//pair_str(n0, n0 + clk%dt, tpy), stat)
            clk%dt = 0
            return
        end if

        if (clk%method == "const") then
            if (has_tc .eqv. has_tcc) then
                call fesm_error(proc, "method = const needs exactly one of time_const and time_const_cal.", &
                    "time_const freezes time_rel, time_const_cal freezes time_cal.", stat)
                clk%dt = 0
                return
            end if
            if (has_tc) then
                call conv_ticks(tpy, tc, proc, "", "time_const", clk%n_const, stat)
            else
                call conv_ticks(tpy, tcc, proc, "", "time_const_cal", clk%n_const, stat)
                clk%n_const = clk%n_const - clk%n_ref
            end if
            if (failed(stat)) then
                clk%dt = 0
                return
            end if
        end if

        clk%n        = clk%n_init
        clk%n_origin = clk%n_init
        clk%k        = 0
        allocate(clk%reg(0))

    end subroutine clock_setup

    subroutine clock_start(clk, stat)
        ! Lock the registry, check every restart time (and time_end) against
        ! every cadence grid, log the schedule table. Required before advance.
        type(clock_t), intent(inout)         :: clk
        integer,       intent(out), optional :: stat

        character(len=*), parameter :: proc = "fesm_clock::clock_start"
        integer(int64) :: p, d, nb
        integer :: r, i, j
        logical :: ok
        character(len=:), allocatable :: det, what

        if (present(stat)) stat = 0
        if (.not. is_initialised(clk, proc, stat)) return
        if (clk%started) then
            call fesm_error(proc, "the clock is already started.", stat=stat)
            return
        end if

        if (has_active_restart(clk)) then
            call restart_grid(clk, p, d, ok)
            if (.not. ok) then
                det = ""
                do i = 1, size(clk%reg)
                    if (clk%reg(i)%kind /= "cadence" .or. clk%reg(i)%mode /= CADENCE_EVERY) cycle
                    if (len(det) > 0) det = det//new_line("a")
                    det = det//'cadence "'//trim(clk%reg(i)%name)//'": every '//dt_str(clk%reg(i)%dt, clk%tpy)// &
                          " yr from "//yr_str(ticks_to_years(clk%reg(i)%phase, clk%tpy))
                end do
                call fesm_error(proc, "the cadence grids have no common time, so no restart bundle can be written.", &
                    det, stat)
                return
            end if

            do r = 1, size(clk%reg)
                if (clk%reg(r)%kind /= "restart") cycle
                do j = 1, size(clk%reg(r)%n)
                    i = offending(clk, clk%reg(r)%n(j))
                    if (i == 0) cycle
                    nb = p + floor_div(clk%reg(r)%n(j) - p, d)*d
                    if (clk%reg(r)%n(j) == clk%n_end) then
                        what = "time_end = "//t_str(clk%n_end, clk%tpy)//" yr"
                        det  = "a restart bundle is always written at time_end."//new_line("a")// &
                               "nearest valid time_end: "//pair_str(nb, nb + d, clk%tpy)
                    else
                        what = 'restart schedule "'//trim(clk%reg(r)%name)//'": time '// &
                               t_str(clk%reg(r)%n(j), clk%tpy)//" yr"
                        det  = "nearest valid: "//pair_str(nb, nb + d, clk%tpy)
                    end if
                    call fesm_error(proc, what//' is not on the grid of cadence "'//trim(clk%reg(i)%name)// &
                        '" (dt = '//dt_str(clk%reg(i)%dt, clk%tpy)//" yr, phase = "// &
                        yr_str(ticks_to_years(clk%reg(i)%phase, clk%tpy))//" yr).", det, stat)
                    return
                end do
            end do
        end if

        clk%started = .true.
        call clock_describe(clk)

    end subroutine clock_start

    subroutine clock_advance(clk)
        ! n = n + dt, k = k + 1.
        type(clock_t), intent(inout) :: clk

        if (.not. clk%started) then
            call fesm_error("fesm_clock::clock_advance", "clock_start has not been called.")
        end if
        if (clock_finished(clk)) then
            call fesm_error("fesm_clock::clock_advance", "the run is finished.", &
                "time = "//t_str(clk%n, clk%tpy)//" yr, time_end = "//t_str(clk%n_end, clk%tpy)//" yr")
        end if
        clk%n = clk%n + clk%dt
        clk%k = clk%k + 1

    end subroutine clock_advance

    subroutine clock_finish(clk)
        ! End early at the current step (kill switch): time_end = time.
        type(clock_t), intent(inout) :: clk
        clk%n_end = clk%n
    end subroutine clock_finish

    pure logical function clock_finished(clk)
        type(clock_t), intent(in) :: clk
        clock_finished = clk%n >= clk%n_end
    end function clock_finished

    pure integer(int64) function clock_step_index(clk)
        ! Steps taken in this run (resets per chunk).
        type(clock_t), intent(in) :: clk
        clock_step_index = clk%k
    end function clock_step_index

    ! =========================================================================
    ! Clock: views and run parameters [yr]
    ! =========================================================================

    pure real(dp) function clock_time(clk)
        ! Model time.
        type(clock_t), intent(in) :: clk
        clock_time = ticks_to_years(clk%n, clk%tpy)
    end function clock_time

    pure real(dp) function clock_time_rel(clk)
        ! Forcing time: model time (rel), model time - time_ref (cal), frozen (const).
        type(clock_t), intent(in) :: clk
        select case (clk%method)
        case ("cal")
            clock_time_rel = ticks_to_years(clk%n - clk%n_ref, clk%tpy)
        case ("const")
            clock_time_rel = ticks_to_years(clk%n_const, clk%tpy)
        case default
            clock_time_rel = ticks_to_years(clk%n, clk%tpy)
        end select
    end function clock_time_rel

    pure real(dp) function clock_time_cal(clk)
        ! Calendar time: always time_rel + time_ref.
        type(clock_t), intent(in) :: clk
        select case (clk%method)
        case ("cal")
            clock_time_cal = ticks_to_years(clk%n, clk%tpy)
        case ("const")
            clock_time_cal = ticks_to_years(clk%n_const + clk%n_ref, clk%tpy)
        case default
            clock_time_cal = ticks_to_years(clk%n + clk%n_ref, clk%tpy)
        end select
    end function clock_time_cal

    pure real(dp) function clock_time_elapsed(clk)
        ! Time since the experiment start (kept across chunk continuations).
        type(clock_t), intent(in) :: clk
        clock_time_elapsed = ticks_to_years(clk%n - clk%n_origin, clk%tpy)
    end function clock_time_elapsed

    pure real(dp) function clock_time_init(clk)
        type(clock_t), intent(in) :: clk
        clock_time_init = ticks_to_years(clk%n_init, clk%tpy)
    end function clock_time_init

    pure real(dp) function clock_time_end(clk)
        type(clock_t), intent(in) :: clk
        clock_time_end = ticks_to_years(clk%n_end, clk%tpy)
    end function clock_time_end

    pure real(dp) function clock_dtt(clk)
        type(clock_t), intent(in) :: clk
        clock_dtt = ticks_to_years(clk%dt, clk%tpy)
    end function clock_dtt

    ! =========================================================================
    ! Conversions
    ! =========================================================================

    integer(int64) function clock_to_ticks(clk, years, what, stat)
        ! `years` as ticks; stops (or stat = 1) if it is not a whole number of
        ! ticks. `what` names the value in the message.
        type(clock_t),    intent(in)            :: clk
        real(dp),         intent(in)            :: years
        character(len=*), intent(in)            :: what
        integer,          intent(out), optional :: stat

        call conv_ticks(clk%tpy, years, "fesm_clock::clock_to_ticks", "", what, clock_to_ticks, stat)

    end function clock_to_ticks

    pure real(dp) function clock_to_years(clk, n)
        ! real(n)/real(N): correctly rounded (both exact in dp for |n| < 2^53).
        type(clock_t),  intent(in) :: clk
        integer(int64), intent(in) :: n
        clock_to_years = ticks_to_years(n, clk%tpy)
    end function clock_to_years

    pure logical function clock_is_representable(clk, years)
        ! The test of clock_to_ticks, without stopping.
        type(clock_t), intent(in) :: clk
        real(dp),      intent(in) :: years
        real(dp) :: x
        x = years*real(clk%tpy, dp)
        clock_is_representable = .false.
        if (.not. (abs(x) < real(N_EXACT, dp))) return
        clock_is_representable = abs(years - ticks_to_years(nint(x, int64), clk%tpy)) <= 4*spacing(years)
    end function clock_is_representable

    logical function clock_restart_allowed(clk)
        ! Current time is on the grid of every registered cadence.
        type(clock_t), intent(in) :: clk
        clock_restart_allowed = offending(clk, clk%n) == 0
    end function clock_restart_allowed

    pure logical function is_multiple_of_all(n, cadences)
        ! `n` [ticks] is on the grid of every active (every-dt) cadence.
        integer(int64),  intent(in) :: n
        type(cadence_t), intent(in) :: cadences(:)
        integer :: i
        is_multiple_of_all = .true.
        do i = 1, size(cadences)
            if (cadences(i)%mode /= CADENCE_EVERY) cycle
            if (modulo(n - cadences(i)%phase, cadences(i)%dt) /= 0) then
                is_multiple_of_all = .false.
                return
            end if
        end do
    end function is_multiple_of_all

    ! =========================================================================
    ! Calendar (fixed-length years only)
    ! =========================================================================

    pure real(dp) function clock_year_fraction(clk)
        ! Position in the model-time year, [0, 1).
        type(clock_t), intent(in) :: clk
        clock_year_fraction = ticks_to_years(modulo(clk%n, clk%tpy), clk%tpy)
    end function clock_year_fraction

    integer function clock_day_of_year(clk, stat)
        ! 1..days_per_year.
        type(clock_t), intent(in)            :: clk
        integer,       intent(out), optional :: stat
        clock_day_of_year = 0
        if (.not. has_calendar(clk, "fesm_clock::clock_day_of_year", stat)) return
        clock_day_of_year = day_of_year(clk)
    end function clock_day_of_year

    integer function clock_month(clk, stat)
        ! 1..12.
        type(clock_t), intent(in)            :: clk
        integer,       intent(out), optional :: stat
        clock_month = 0
        if (.not. has_calendar(clk, "fesm_clock::clock_month", stat)) return
        clock_month = month_of(clk, day_of_year(clk))
    end function clock_month

    integer function clock_day_of_month(clk, stat)
        ! 1..30 (360_day), 1..31 (365_day).
        type(clock_t), intent(in)            :: clk
        integer,       intent(out), optional :: stat
        integer :: doy, m
        clock_day_of_month = 0
        if (.not. has_calendar(clk, "fesm_clock::clock_day_of_month", stat)) return
        doy = day_of_year(clk)
        m   = month_of(clk, doy)
        if (clk%days_per_year == 360) then
            clock_day_of_month = doy - 30*(m - 1)
        else
            clock_day_of_month = doy - MONTH_END_365(m - 1)
        end if
    end function clock_day_of_month

    real(dp) function clock_sec_day(clk, stat)
        ! Day length [s]: sec_year/days_per_year.
        type(clock_t), intent(in)            :: clk
        integer,       intent(out), optional :: stat
        clock_sec_day = 0.0_dp
        if (.not. has_calendar(clk, "fesm_clock::clock_sec_day", stat)) return
        clock_sec_day = clk%sec_year/real(clk%days_per_year, dp)
    end function clock_sec_day

    ! =========================================================================
    ! Describe
    ! =========================================================================

    subroutine clock_describe(clk)
        ! The schedule table of the run, through log_msg (LOG_INFO).
        type(clock_t), intent(in) :: clk

        character(len=1), parameter :: nl = new_line("a")
        character(len=:), allocatable :: head, det, names, sdt, sph, sfirst, slast
        integer(int64) :: p, d, first, last, cnt, tpy
        integer :: i
        logical :: ok

        tpy = clk%tpy
        head = trim(clk%method)//", model time "//t_str(clk%n_init, tpy)//" -> "//t_str(clk%n_end, tpy)// &
               " yr, dtt = "//dt_str(clk%dt, tpy)//" yr ("//to_str((clk%n_end - clk%n_init)/max(clk%dt, 1_int64))// &
               " steps)"

        select case (clk%method)
        case ("cal")
            det = "forcing  time_cal = model time; time_rel = time_cal - "//yr_str(ticks_to_years(clk%n_ref, tpy))
        case ("const")
            det = "forcing  time_rel frozen at "//t_str(clk%n_const, tpy)//"; time_cal = "// &
                  t_str(clk%n_const + clk%n_ref, tpy)
        case default
            det = "forcing  time_rel = model time; time_cal = time_rel + "//yr_str(ticks_to_years(clk%n_ref, tpy))
        end select
        det = det//nl//"calendar "//trim(clk%calendar)//"; year = "//trim(clk%year)//" ("// &
              yr_str(clk%sec_year, 2)//" s); tick = 1/"//to_str(tpy)//" yr"
        det = det//nl//"elapsed  origin "//t_str(clk%n_origin, tpy)
        if (clk%continued) then
            det = det//" (continued from "//trim(clk%restart_file)//")"
        else if (clk%branched) then
            det = det//" (this run; branched from "//trim(clk%restart_file)//")"
        else
            det = det//" (this run)"
        end if
        det = det//nl//pad("name", 10)//pad("kind", 9)//pad("dt [yr]", 10)//pad("phase", 7)// &
              pad("first", 10)//pad("last", 10)//"count"

        names = ""
        if (allocated(clk%reg)) then
            do i = 1, size(clk%reg)
                associate (e => clk%reg(i))
                    first = 0
                    last  = 0
                    sph = "-"
                    if (e%kind == "cadence") then
                        select case (e%mode)
                        case (CADENCE_STEP)
                            sdt   = "step"
                            first = clk%n_init + clk%dt
                            last  = clk%n_end
                            cnt   = (clk%n_end - clk%n_init)/clk%dt
                        case (CADENCE_EVERY)
                            sdt   = dt_str(e%dt, tpy)
                            sph   = yr_str(ticks_to_years(e%phase, tpy))
                            first = e%phase + (floor_div(clk%n_init - e%phase, e%dt) + 1)*e%dt
                            last  = e%phase + floor_div(clk%n_end - e%phase, e%dt)*e%dt
                            cnt   = max(0_int64, (last - first)/e%dt + 1)
                            if (len(names) > 0) names = names//", "
                            names = names//trim(e%name)
                        case default
                            sdt = "off"
                            cnt = 0
                        end select
                    else
                        select case (e%method)
                        case ("const")
                            sdt = dt_str(e%dt, tpy)
                            sph = "0"
                        case ("none")
                            sdt = "none"
                        case default
                            sdt = "list"
                        end select
                        cnt = size(e%n)
                        if (cnt > 0) then
                            first = e%n(1)
                            last  = e%n(cnt)
                        end if
                    end if
                    if (cnt > 0) then
                        sfirst = t_str(first, tpy)
                        slast  = t_str(last, tpy)
                    else
                        sfirst = "-"
                        slast  = "-"
                    end if
                    det = det//nl//pad(trim(e%name), 10)//pad(trim(e%kind), 9)//pad(sdt, 10)//pad(sph, 7)// &
                          pad(sfirst, 10)//pad(slast, 10)//to_str(cnt)
                end associate
            end do
        end if

        call restart_grid(clk, p, d, ok)
        if (len(names) == 0) then
            det = det//nl//"restart grid: every step"
        else if (ok) then
            det = det//nl//"restart grid: every "//dt_str(d, tpy)//" yr from "//yr_str(ticks_to_years(p, tpy))// &
                  " ("//names//")"
        else
            det = det//nl//"restart grid: empty, the cadence grids have no common time ("//names//")"
        end if

        call log_msg(LOG_INFO, "fesm_clock", head, det)

    end subroutine clock_describe

    ! =========================================================================
    ! Clock restart
    ! =========================================================================

    subroutine clock_restart_write(clk, filename)
        ! Clock state as dp scalars [ticks] (exact, |n| < 2^53) and the
        ! attribute clock_method. Creates the file if it does not exist.
        type(clock_t),    intent(in) :: clk
        character(len=*), intent(in) :: filename

        integer :: ncid

        call ensure_file(filename)
        call nc_open(filename, ncid, writable=.true.)
        call nc_write(filename, "clock_ticks_per_year", real(clk%tpy, dp), long_name="ticks per year", &
                      ncid=ncid)
        call nc_write(filename, "clock_n",        real(clk%n, dp),        units="ticks", long_name="model time", &
                      ncid=ncid)
        call nc_write(filename, "clock_n_origin", real(clk%n_origin, dp), units="ticks", &
                      long_name="experiment start", ncid=ncid)
        call nc_write(filename, "clock_dtt",      real(clk%dt, dp),       units="ticks", long_name="time step", &
                      ncid=ncid)
        call nc_write_attr(filename, ncid, "clock_method", trim(clk%method))
        call nc_close(ncid)

    end subroutine clock_restart_write

    subroutine clock_restart_read(clk, filename, stat)
        ! Bundle time == time_init: continuation (n_origin from the bundle,
        ! continued = .true.); any other time: branch (n_origin = n_init). A
        ! bundle written at another ticks_per_year is converted exactly or
        ! rejected. Call after clock_init, before clock_start.
        type(clock_t),    intent(inout)         :: clk
        character(len=*), intent(in)            :: filename
        integer,          intent(out), optional :: stat

        character(len=*), parameter :: proc = "fesm_clock::clock_restart_read"
        integer(int64) :: tpy_b, n, n_org, ndt, nb
        character(len=32) :: meth_b
        real(dp) :: x
        logical :: ok

        if (present(stat)) stat = 0
        if (.not. is_initialised(clk, proc, stat)) return
        if (clk%started) then
            call fesm_error(proc, "restart_read must come before clock_start.", trim(filename), stat)
            return
        end if
        if (.not. has_vars(filename, [character(len=20) :: "clock_ticks_per_year", "clock_n", &
                           "clock_n_origin", "clock_dtt"], proc, stat)) return

        call read_bundle_tpy(filename, tpy_b, proc, stat)
        if (failed(stat)) return

        call read_bundle_ticks(filename, "clock_n", tpy_b, clk%tpy, n, proc, stat)
        if (failed(stat)) return
        call read_bundle_ticks(filename, "clock_n_origin", tpy_b, clk%tpy, n_org, proc, stat)
        if (failed(stat)) return

        clk%restart_file = filename
        if (n == clk%n_init) then
            clk%continued = .true.
            clk%branched  = .false.
            clk%n_origin  = n_org
            call log_msg(LOG_INFO, "fesm_clock", "continuation from "//trim(filename)//" at "// &
                t_str(n, clk%tpy)//" yr; elapsed origin "//t_str(n_org, clk%tpy)//" yr")
        else
            clk%continued = .false.
            clk%branched  = .true.
            clk%n_origin  = clk%n_init
            call log_msg(LOG_INFO, "fesm_clock", "branch from "//trim(filename)//": bundle time "// &
                t_str(n, clk%tpy)//" yr, time_init "//t_str(clk%n_init, clk%tpy)//" yr; elapsed origin reset")
        end if

        if (tpy_b /= clk%tpy) call log_msg(LOG_INFO, "fesm_clock", "bundle ticks_per_year "// &
            to_str(tpy_b)//" converted to "//to_str(clk%tpy))
        call nc_read(filename, "clock_dtt", x)
        call dp_ticks(x, nb, ok)
        if (ok) call convert_ticks(nb, tpy_b, clk%tpy, ndt, ok)
        if (ok .and. ndt /= clk%dt) call log_msg(LOG_INFO, "fesm_clock", "dtt changed from "// &
            dt_str(ndt, clk%tpy)//" to "//dt_str(clk%dt, clk%tpy)//" yr")
        if (nc_exists_attr(filename, "clock_method")) then
            meth_b = ""
            call nc_read_attr(filename, "clock_method", meth_b)
            if (trim(meth_b) /= trim(clk%method)) call log_msg(LOG_INFO, "fesm_clock", &
                "method changed from "//trim(meth_b)//" to "//trim(clk%method))
        end if

    end subroutine clock_restart_read

    ! =========================================================================
    ! Cadences
    ! =========================================================================

    subroutine cadence_init_dt(c, clk, name, dt, phase, stat)
        ! Every `dt` yr on the absolute grid phase + k*dt; dt <= 0: off.
        type(cadence_t),  intent(inout)         :: c
        type(clock_t),    intent(inout)         :: clk
        character(len=*), intent(in)            :: name
        real(dp),         intent(in)            :: dt
        real(dp),         intent(in),  optional :: phase
        integer,          intent(out), optional :: stat

        character(len=*), parameter :: proc = "fesm_clock::cadence_init"
        character(len=:), allocatable :: ctx
        integer(int64) :: ndt, nph

        if (present(stat)) stat = 0
        ctx = 'cadence "'//trim(name)//'"'
        if (.not. can_register(clk, proc, ctx, stat)) return

        if (dt <= 0.0_dp) then
            call cadence_setup(c, clk, name, CADENCE_OFF, 0_int64, 0_int64, proc, ctx, stat)
            return
        end if
        call conv_ticks(clk%tpy, dt, proc, ctx, "dt", ndt, stat, per_year_hint=.true.)
        if (failed(stat)) return
        nph = 0
        if (present(phase)) call conv_ticks(clk%tpy, phase, proc, ctx, "phase", nph, stat)
        if (failed(stat)) return
        call cadence_setup(c, clk, name, CADENCE_EVERY, ndt, nph, proc, ctx, stat)

    end subroutine cadence_init_dt

    subroutine cadence_init_step(c, clk, name, every_step, stat)
        ! Due at every step, whatever dtt is in this chunk; .false.: off.
        type(cadence_t),  intent(inout)         :: c
        type(clock_t),    intent(inout)         :: clk
        character(len=*), intent(in)            :: name
        logical,          intent(in)            :: every_step
        integer,          intent(out), optional :: stat

        character(len=*), parameter :: proc = "fesm_clock::cadence_init"
        character(len=:), allocatable :: ctx

        if (present(stat)) stat = 0
        ctx = 'cadence "'//trim(name)//'"'
        if (.not. can_register(clk, proc, ctx, stat)) return
        if (every_step) then
            call cadence_setup(c, clk, name, CADENCE_STEP, 0_int64, 0_int64, proc, ctx, stat)
        else
            call cadence_setup(c, clk, name, CADENCE_OFF, 0_int64, 0_int64, proc, ctx, stat)
        end if

    end subroutine cadence_init_step

    subroutine cadence_init_per_year(c, clk, name, per_year, phase, stat)
        ! dt = N/per_year exactly (sub-annual cadences without decimal
        ! literals); per_year <= 0: off.
        type(cadence_t),  intent(inout)         :: c
        type(clock_t),    intent(inout)         :: clk
        character(len=*), intent(in)            :: name
        integer,          intent(in)            :: per_year
        real(dp),         intent(in),  optional :: phase
        integer,          intent(out), optional :: stat

        character(len=*), parameter :: proc = "fesm_clock::cadence_init"
        character(len=:), allocatable :: ctx
        integer(int64) :: nph

        if (present(stat)) stat = 0
        ctx = 'cadence "'//trim(name)//'"'
        if (.not. can_register(clk, proc, ctx, stat)) return

        if (per_year <= 0) then
            call cadence_setup(c, clk, name, CADENCE_OFF, 0_int64, 0_int64, proc, ctx, stat)
            return
        end if
        if (modulo(clk%tpy, int(per_year, int64)) /= 0) then
            call fesm_error(proc, ctx//": per_year = "//to_str(per_year)//" does not divide ticks_per_year = "// &
                to_str(clk%tpy)//".", "1/"//to_str(per_year)//" yr is not a whole number of ticks.", stat)
            return
        end if
        nph = 0
        if (present(phase)) call conv_ticks(clk%tpy, phase, proc, ctx, "phase", nph, stat)
        if (failed(stat)) return
        call cadence_setup(c, clk, name, CADENCE_EVERY, clk%tpy/int(per_year, int64), nph, proc, ctx, stat)

    end subroutine cadence_init_per_year

    pure logical function cadence_due(c, clk)
        ! n > n_last and n on the cadence grid (or every step).
        type(cadence_t), intent(in) :: c
        type(clock_t),   intent(in) :: clk
        cadence_due = .false.
        if (clk%n <= c%n_last) return
        select case (c%mode)
        case (CADENCE_STEP)
            cadence_due = .true.
        case (CADENCE_EVERY)
            cadence_due = modulo(clk%n - c%phase, c%dt) == 0
        end select
    end function cadence_due

    subroutine cadence_mark(c, clk)
        ! n_last = n; call after the component reached t1.
        type(cadence_t), intent(inout) :: c
        type(clock_t),   intent(in)    :: clk
        c%n_last = clk%n
    end subroutine cadence_mark

    subroutine cadence_interval(c, clk, t0, t1)
        ! [t0, t1] = [last mark, now] [yr], each evaluated once from ticks.
        type(cadence_t), intent(in)  :: c
        type(clock_t),   intent(in)  :: clk
        real(dp),        intent(out) :: t0, t1
        t0 = ticks_to_years(c%n_last, clk%tpy)
        t1 = ticks_to_years(clk%n, clk%tpy)
    end subroutine cadence_interval

    pure real(dp) function cadence_t_last(c, clk)
        type(cadence_t), intent(in) :: c
        type(clock_t),   intent(in) :: clk
        cadence_t_last = ticks_to_years(c%n_last, clk%tpy)
    end function cadence_t_last

    pure real(dp) function cadence_t_next(c, clk)
        ! Next grid time after now [yr]; huge() for an off cadence.
        type(cadence_t), intent(in) :: c
        type(clock_t),   intent(in) :: clk
        select case (c%mode)
        case (CADENCE_STEP)
            cadence_t_next = ticks_to_years(clk%n + clk%dt, clk%tpy)
        case (CADENCE_EVERY)
            cadence_t_next = ticks_to_years(c%phase + (floor_div(clk%n - c%phase, c%dt) + 1)*c%dt, clk%tpy)
        case default
            cadence_t_next = huge(1.0_dp)
        end select
    end function cadence_t_next

    pure real(dp) function cadence_dt_years(c, clk)
        ! dt [yr]: the cadence's, dtt for every step, 0 when off.
        type(cadence_t), intent(in) :: c
        type(clock_t),   intent(in) :: clk
        select case (c%mode)
        case (CADENCE_STEP)
            cadence_dt_years = ticks_to_years(clk%dt, clk%tpy)
        case (CADENCE_EVERY)
            cadence_dt_years = ticks_to_years(c%dt, clk%tpy)
        case default
            cadence_dt_years = 0.0_dp
        end select
    end function cadence_dt_years

    subroutine cadence_restart_write(c, filename)
        ! cadence_<name>_{n_last,dt,phase} [ticks]; dt = 0: every step, -1: off.
        ! Write to the bundle that holds the clock (clock_restart_write).
        type(cadence_t),  intent(in) :: c
        character(len=*), intent(in) :: filename

        character(len=:), allocatable :: base
        integer :: ncid

        base = "cadence_"//trim(c%name)
        call ensure_file(filename)
        call nc_open(filename, ncid, writable=.true.)
        call nc_write(filename, base//"_n_last", real(c%n_last, dp), units="ticks", ncid=ncid)
        call nc_write(filename, base//"_dt",     real(dt_code(c), dp), units="ticks", ncid=ncid)
        call nc_write(filename, base//"_phase",  real(c%phase, dp), units="ticks", ncid=ncid)
        call nc_close(ncid)

    end subroutine cadence_restart_write

    subroutine cadence_restart_read(c, clk, filename, stat)
        ! Continuation: restore n_last (checked == time_init); branch or a
        ! cadence absent from the bundle: n_last = time_init. A changed dt or
        ! phase is logged (the new grid was checked by cadence_init). Call after
        ! clock_restart_read.
        type(cadence_t),  intent(inout)         :: c
        type(clock_t),    intent(in)            :: clk
        character(len=*), intent(in)            :: filename
        integer,          intent(out), optional :: stat

        character(len=*), parameter :: proc = "fesm_clock::cadence_restart_read"
        character(len=:), allocatable :: base, ctx
        integer(int64) :: tpy_b, n_last, ndt, nph
        real(dp) :: x

        if (present(stat)) stat = 0
        base = "cadence_"//trim(c%name)
        ctx  = 'cadence "'//trim(c%name)//'"'
        c%n_last = clk%n_init

        if (.not. clk%continued) return
        if (.not. nc_exists_var(filename, base//"_n_last")) then
            call log_msg(LOG_INFO, "fesm_clock", ctx//": not in "//trim(filename)//"; starts at time_init")
            return
        end if
        if (.not. has_vars(filename, [character(len=64) :: "clock_ticks_per_year", base//"_dt", &
                           base//"_phase"], proc, stat)) return

        call read_bundle_tpy(filename, tpy_b, proc, stat)
        if (failed(stat)) return
        call read_bundle_ticks(filename, base//"_n_last", tpy_b, clk%tpy, n_last, proc, stat)
        if (failed(stat)) return
        call read_bundle_ticks(filename, base//"_phase", tpy_b, clk%tpy, nph, proc, stat)
        if (failed(stat)) return
        call nc_read(filename, base//"_dt", x)
        if (x > 0.0_dp) then
            call read_bundle_ticks(filename, base//"_dt", tpy_b, clk%tpy, ndt, proc, stat)
            if (failed(stat)) return
        else
            ndt = nint(x, int64)
        end if

        if (ndt /= -1 .and. c%mode /= CADENCE_OFF .and. n_last /= clk%n_init) then
            call fesm_error(proc, ctx//": last mark in "//trim(filename)//" ("//t_str(n_last, clk%tpy)// &
                " yr) is not time_init ("//t_str(clk%n_init, clk%tpy)//" yr).", &
                "restart bundles are written only where every cadence fired.", stat)
            return
        end if
        if (ndt /= dt_code(c) .or. (c%mode == CADENCE_EVERY .and. nph /= c%phase)) then
            call log_msg(LOG_INFO, "fesm_clock", ctx//": grid changed from "//grid_str(ndt, nph, clk%tpy)// &
                " to "//grid_str(dt_code(c), c%phase, clk%tpy))
        end if

    end subroutine cadence_restart_read

    ! =========================================================================
    ! Schedules
    ! =========================================================================

    subroutine schedule_init_nml(s, clk, name, filename, group, restart, defaults_file, stat)
        ! Keys as timeout: method = none | const | times | file, and dt
        ! (const), times (one string of values and t0:dt:t1 ranges), file
        ! (lines of the same; # comments). restart = .true.: a restart schedule.
        type(schedule_t), intent(inout)         :: s
        type(clock_t),    intent(inout)         :: clk
        character(len=*), intent(in)            :: name, filename, group
        logical,          intent(in),  optional :: restart
        character(len=*), intent(in),  optional :: defaults_file
        integer,          intent(out), optional :: stat

        character(len=*), parameter :: proc = "fesm_clock::schedule_init"
        character(len=:), allocatable :: ctx, key
        character(len=32)  :: method
        character(len=512) :: str
        real(dp) :: dt
        integer(int64), allocatable :: lst(:)
        integer(int64) :: ndt
        integer :: nl
        logical :: is_rst

        if (present(stat)) stat = 0
        ctx = 'schedule "'//trim(name)//'"'
        if (.not. can_register(clk, proc, ctx, stat)) return
        is_rst = .false.
        if (present(restart)) is_rst = restart
        ndt = 0

        if (.not. has_key(filename, group, "method", defaults_file)) then
            call fesm_error(proc, ctx//': key "method" not found in group "'//trim(group)//'" of '// &
                trim(filename), "use none, const, times or file.", stat)
            return
        end if
        call read_key(filename, group, "method", method, defaults_file)

        allocate(lst(16))
        nl = 0
        select case (trim(method))
        case ("none")
        case ("const", "times", "file")
            key = "dt"
            if (trim(method) == "times") key = "times"
            if (trim(method) == "file")  key = "file"
            if (.not. has_key(filename, group, key, defaults_file)) then
                call fesm_error(proc, ctx//': method = '//trim(method)//' needs key "'//key//'" in group "'// &
                    trim(group)//'" of '//trim(filename), stat=stat)
                return
            end if
            select case (trim(method))
            case ("const")
                call read_key(filename, group, "dt", dt, defaults_file)
                call const_times(clk, dt, proc, ctx, lst, nl, ndt, stat)
            case ("times")
                str = ""
                call read_key(filename, group, "times", str, defaults_file)
                call parse_times(clk, str, proc, ctx, "times", lst, nl, stat)
            case ("file")
                str = ""
                call read_key(filename, group, "file", str, defaults_file)
                call parse_time_file(clk, trim(str), proc, ctx, lst, nl, stat)
            end select
            if (failed(stat)) return
        case default
            call fesm_error(proc, ctx//': method = "'//trim(method)//'" is not recognised.', &
                "use none, const, times or file.", stat)
            return
        end select

        call schedule_setup(s, clk, name, trim(method), is_rst, lst(1:nl), ndt, proc, ctx, stat)

    end subroutine schedule_init_nml

    subroutine schedule_init_args(s, clk, name, dt, times, restart, stat)
        ! From arguments: dt > 0 as method const, times as method times;
        ! neither (or dt <= 0): none.
        type(schedule_t), intent(inout)         :: s
        type(clock_t),    intent(inout)         :: clk
        character(len=*), intent(in)            :: name
        real(dp),         intent(in),  optional :: dt
        real(dp),         intent(in),  optional :: times(:)
        logical,          intent(in),  optional :: restart
        integer,          intent(out), optional :: stat

        character(len=*), parameter :: proc = "fesm_clock::schedule_init"
        character(len=:), allocatable :: ctx
        character(len=8) :: method
        integer(int64), allocatable :: lst(:)
        integer :: nl, i
        integer(int64) :: n, ndt
        logical :: is_rst

        if (present(stat)) stat = 0
        ctx = 'schedule "'//trim(name)//'"'
        if (.not. can_register(clk, proc, ctx, stat)) return
        is_rst = .false.
        if (present(restart)) is_rst = restart
        if (present(dt) .and. present(times)) then
            call fesm_error(proc, ctx//": give dt or times, not both.", stat=stat)
            return
        end if

        allocate(lst(16))
        nl = 0
        ndt = 0
        method = "none"
        if (present(dt)) then
            if (dt > 0.0_dp) then
                method = "const"
                call const_times(clk, dt, proc, ctx, lst, nl, ndt, stat)
                if (failed(stat)) return
            end if
        else if (present(times)) then
            method = "times"
            do i = 1, size(times)
                call conv_ticks(clk%tpy, times(i), proc, ctx, "time", n, stat)
                if (failed(stat)) return
                call push(lst, nl, n)
            end do
        end if

        call schedule_setup(s, clk, name, trim(method), is_rst, lst(1:nl), ndt, proc, ctx, stat)

    end subroutine schedule_init_args

    logical function schedule_due(s, clk) result(due)
        ! Now is a scheduled time, or the run is finished (kill switch: an
        ! output schedule is due; a restart schedule only where every cadence
        ! fired, otherwise the skipped bundle is logged). Advances the cursor.
        type(schedule_t), intent(inout) :: s
        type(clock_t),    intent(in)    :: clk

        integer :: i

        due = .false.
        if (.not. allocated(s%n) .or. s%method == "none") return

        do while (s%k <= size(s%n))
            if (s%n(s%k) >= clk%n) exit
            s%k = s%k + 1
        end do
        if (s%k <= size(s%n)) due = s%n(s%k) == clk%n
        if (due .or. .not. clock_finished(clk)) return

        if (.not. s%restart) then
            due = .true.
            return
        end if
        i = offending(clk, clk%n)
        if (i == 0) then
            due = .true.
        else if (s%n_skip /= clk%n) then
            s%n_skip = clk%n
            call fesm_warn("fesm_clock", 'schedule "'//trim(s%name)//'": no restart bundle written at '// &
                t_str(clk%n, clk%tpy)//" yr (run ended early).", 'not on the grid of cadence "'// &
                trim(clk%reg(i)%name)//'" ('//grid_str(clk%reg(i)%dt, clk%reg(i)%phase, clk%tpy)//")")
        end if

    end function schedule_due

    real(dp) function schedule_t_next(s, clk)
        ! First scheduled time after now [yr]; huge() if none is left.
        type(schedule_t), intent(in) :: s
        type(clock_t),    intent(in) :: clk
        integer :: i
        schedule_t_next = huge(1.0_dp)
        if (.not. allocated(s%n)) return
        do i = max(1, s%k), size(s%n)
            if (s%n(i) > clk%n) then
                schedule_t_next = ticks_to_years(s%n(i), clk%tpy)
                return
            end if
        end do
    end function schedule_t_next

    function schedule_times(s, clk) result(t)
        ! All scheduled times [yr].
        type(schedule_t), intent(in) :: s
        type(clock_t),    intent(in) :: clk
        real(dp), allocatable :: t(:)
        integer :: i
        allocate(t(schedule_size(s)))
        do i = 1, size(t)
            t(i) = ticks_to_years(s%n(i), clk%tpy)
        end do
    end function schedule_times

    pure integer function schedule_size(s)
        type(schedule_t), intent(in) :: s
        schedule_size = 0
        if (allocated(s%n)) schedule_size = size(s%n)
    end function schedule_size

    ! =========================================================================
    ! Internal: setup helpers
    ! =========================================================================

    subroutine set_calendar(clk, calendar, proc, stat)
        type(clock_t),    intent(inout)         :: clk
        character(len=*), intent(in)            :: calendar, proc
        integer,          intent(out), optional :: stat

        if (present(stat)) stat = 0
        select case (trim(calendar))
        case ("360_day")
            clk%calendar = "360_day"
            clk%days_per_year = 360
        case ("365_day", "noleap")
            clk%calendar = "365_day"
            clk%days_per_year = 365
        case ("none")
            clk%calendar = "none"
            clk%days_per_year = 0
        case ("gregorian", "standard", "proleptic_gregorian", "julian")
            call fesm_error(proc, 'calendar = "'//trim(calendar)//'" is not supported.', &
                "its year length varies, so a day is not a fixed fraction of a year."//new_line("a")// &
                "use 360_day, 365_day (noleap) or none.", stat)
        case default
            call fesm_error(proc, 'calendar = "'//trim(calendar)//'" is not recognised.', &
                "use 360_day, 365_day (noleap) or none; fesm_clock has no default.", stat)
        end select

    end subroutine set_calendar

    subroutine set_year(clk, year, proc, stat)
        type(clock_t),    intent(inout)         :: clk
        character(len=*), intent(in)            :: year, proc
        integer,          intent(out), optional :: stat

        if (present(stat)) stat = 0
        select case (trim(year))
        case ("360d")
            clk%sec_year = sec_year_360d
        case ("365d")
            clk%sec_year = sec_year_365d
        case ("julian")
            clk%sec_year = sec_year_julian
        case ("gregorian")
            clk%sec_year = sec_year_gregorian
        case ("tropical")
            clk%sec_year = sec_year_tropical
        case default
            call fesm_error(proc, 'year = "'//trim(year)//'" is not recognised.', &
                "use 360d, 365d, julian, gregorian or tropical (phys_constants sec_year_*); "// &
                "fesm_clock has no default.", stat)
            return
        end select
        clk%year = trim(year)

    end subroutine set_year

    subroutine cadence_setup(c, clk, name, mode, ndt, nph, proc, ctx, stat)
        ! Validate the cadence grid against the step grid, reset c, register.
        type(cadence_t),  intent(inout)         :: c
        type(clock_t),    intent(inout)         :: clk
        character(len=*), intent(in)            :: name, proc, ctx
        integer,          intent(in)            :: mode
        integer(int64),   intent(in)            :: ndt, nph
        integer,          intent(out), optional :: stat

        type(cadence_t) :: fresh
        integer(int64)  :: ph, nb

        if (present(stat)) stat = 0
        ph = 0
        if (mode == CADENCE_EVERY) then
            if (ndt <= 0 .or. modulo(ndt, clk%dt) /= 0) then
                nb = floor_div(ndt, clk%dt)*clk%dt
                if (nb > 0) then
                    call fesm_error(proc, ctx//": dt = "//dt_str(ndt, clk%tpy)//" yr is not a multiple of dtt = "// &
                        dt_str(clk%dt, clk%tpy)//" yr.", "nearest: "//pair_str(nb, nb + clk%dt, clk%tpy), stat)
                else
                    call fesm_error(proc, ctx//": dt = "//dt_str(ndt, clk%tpy)//" yr is not a multiple of dtt = "// &
                        dt_str(clk%dt, clk%tpy)//" yr.", "nearest: "//dt_str(clk%dt, clk%tpy), stat)
                end if
                return
            end if
            ph = modulo(nph, ndt)
            if (modulo(ph - clk%n_init, clk%dt) /= 0) then
                nb = clk%n_init + floor_div(ph - clk%n_init, clk%dt)*clk%dt
                call fesm_error(proc, ctx//": the grid phase + k*dt is not on the step grid time_init + k*dtt.", &
                    "phase = "//yr_str(ticks_to_years(ph, clk%tpy))//" yr, time_init = "// &
                    t_str(clk%n_init, clk%tpy)//" yr, dtt = "//dt_str(clk%dt, clk%tpy)//" yr"//new_line("a")// &
                    "nearest valid phase: "//pair_str(nb, nb + clk%dt, clk%tpy), stat)
                return
            end if
        end if

        c = fresh
        c%name   = name
        c%mode   = mode
        if (mode == CADENCE_EVERY) c%dt = ndt
        c%phase  = ph
        c%n_last = clk%n_init
        call register(clk, "cadence", name, mode, c%dt, ph, "", c%id)

    end subroutine cadence_setup

    subroutine schedule_setup(s, clk, name, method, is_rst, raw, ndt, proc, ctx, stat)
        ! Clip to [time_init, time_end], check the step grid, add time_init
        ! (output, const) and time_end, drop time_init (restart), sort, unique.
        ! ndt: the const spacing [ticks] (0 for lists), for describe.
        type(schedule_t), intent(inout)         :: s
        type(clock_t),    intent(inout)         :: clk
        character(len=*), intent(in)            :: name, method, proc, ctx
        logical,          intent(in)            :: is_rst
        integer(int64),   intent(in)            :: raw(:), ndt
        integer,          intent(out), optional :: stat

        type(schedule_t) :: fresh
        integer(int64), allocatable :: lst(:)
        integer(int64) :: nb
        integer :: i, nl

        if (present(stat)) stat = 0
        allocate(lst(max(16, size(raw) + 2)))
        nl = 0
        do i = 1, size(raw)
            if (raw(i) < clk%n_init .or. raw(i) > clk%n_end) cycle
            if (modulo(raw(i) - clk%n_init, clk%dt) /= 0) then
                nb = clk%n_init + floor_div(raw(i) - clk%n_init, clk%dt)*clk%dt
                call fesm_error(proc, ctx//": time "//t_str(raw(i), clk%tpy)//" yr is not on the step grid.", &
                    "time_init = "//t_str(clk%n_init, clk%tpy)//" yr, dtt = "//dt_str(clk%dt, clk%tpy)//" yr"// &
                    new_line("a")//"nearest valid: "//pair_str(nb, nb + clk%dt, clk%tpy), stat)
                return
            end if
            call push(lst, nl, raw(i))
        end do
        if (method /= "none") then
            if (method == "const" .and. .not. is_rst) call push(lst, nl, clk%n_init)
            call push(lst, nl, clk%n_end)
        end if
        call sort_unique(lst, nl)
        if (is_rst .and. nl > 0) then
            if (lst(1) == clk%n_init) then
                lst(1:nl-1) = lst(2:nl)
                nl = nl - 1
            end if
        end if

        s = fresh
        s%name    = name
        s%method  = method
        s%restart = is_rst
        s%n       = lst(1:nl)
        s%k       = 1
        if (is_rst) then
            call register(clk, "restart", name, CADENCE_OFF, ndt, 0_int64, method, s%id, s%n)
        else
            call register(clk, "output", name, CADENCE_OFF, ndt, 0_int64, method, s%id, s%n)
        end if

    end subroutine schedule_setup

    subroutine const_times(clk, dt, proc, ctx, lst, nl, ndt, stat)
        ! The absolute multiples of dt in [time_init, time_end]; ndt = dt [ticks].
        type(clock_t),    intent(in)              :: clk
        real(dp),         intent(in)              :: dt
        character(len=*), intent(in)              :: proc, ctx
        integer(int64),   intent(inout), allocatable :: lst(:)
        integer,          intent(inout)           :: nl
        integer(int64),   intent(out)             :: ndt
        integer,          intent(out),   optional :: stat

        integer(int64) :: k, k0, k1, nb

        if (present(stat)) stat = 0
        ndt = 0
        if (dt <= 0.0_dp) then
            call fesm_error(proc, ctx//": method = const needs dt > 0.", "dt = "//yr_str(dt)//" yr", stat)
            return
        end if
        call conv_ticks(clk%tpy, dt, proc, ctx, "dt", ndt, stat)
        if (failed(stat)) return
        if (modulo(ndt, clk%dt) /= 0) then
            nb = floor_div(ndt, clk%dt)*clk%dt
            call fesm_error(proc, ctx//": dt = "//dt_str(ndt, clk%tpy)//" yr is not a multiple of dtt = "// &
                dt_str(clk%dt, clk%tpy)//" yr.", "nearest: "//pair_str(max(nb, clk%dt), &
                max(nb, clk%dt) + clk%dt, clk%tpy), stat)
            return
        end if
        k0 = -floor_div(-clk%n_init, ndt)       ! ceiling
        k1 = floor_div(clk%n_end, ndt)
        do k = k0, k1
            call push(lst, nl, k*ndt)
        end do

    end subroutine const_times

    subroutine parse_times(clk, str, proc, ctx, src, lst, nl, stat)
        ! Tokens separated by blanks or commas: a value, or t0:dt:t1. Ranges
        ! are generated in ticks and clipped to [time_init, time_end].
        type(clock_t),    intent(in)              :: clk
        character(len=*), intent(in)              :: str, proc, ctx, src
        integer(int64),   intent(inout), allocatable :: lst(:)
        integer,          intent(inout)           :: nl
        integer,          intent(out),   optional :: stat

        integer :: p0, p1, ic1, ic2, ios
        integer(int64) :: n, n0, dn, n1, k, k0, k1
        real(dp) :: x, y, z
        character(len=:), allocatable :: tok

        if (present(stat)) stat = 0
        p0 = 1
        do
            do while (p0 <= len(str))
                if (.not. is_sep(str(p0:p0))) exit
                p0 = p0 + 1
            end do
            if (p0 > len(str)) exit
            p1 = p0
            do while (p1 < len(str))
                if (is_sep(str(p1+1:p1+1))) exit
                p1 = p1 + 1
            end do
            tok = str(p0:p1)
            p0 = p1 + 1

            ic1 = index(tok, ":")
            ic2 = index(tok, ":", back=.true.)
            if (ic1 == 0) then
                read(tok, *, iostat=ios) x
                if (ios /= 0) then
                    call fesm_error(proc, ctx//': cannot read "'//tok//'" in '//src//".", stat=stat)
                    return
                end if
                call conv_ticks(clk%tpy, x, proc, ctx, "time", n, stat)
                if (failed(stat)) return
                call push(lst, nl, n)
            else
                ios = 1
                if (ic2 > ic1 + 1 .and. ic1 > 1 .and. ic2 < len(tok)) then
                    read(tok(1:ic1-1), *, iostat=ios) x
                    if (ios == 0) read(tok(ic1+1:ic2-1), *, iostat=ios) y
                    if (ios == 0) read(tok(ic2+1:), *, iostat=ios) z
                end if
                if (ios /= 0) then
                    call fesm_error(proc, ctx//': cannot read "'//tok//'" in '//src//" (use t0:dt:t1).", stat=stat)
                    return
                end if
                call conv_ticks(clk%tpy, x, proc, ctx, "range start", n0, stat)
                if (failed(stat)) return
                call conv_ticks(clk%tpy, y, proc, ctx, "range step", dn, stat)
                if (failed(stat)) return
                call conv_ticks(clk%tpy, z, proc, ctx, "range end", n1, stat)
                if (failed(stat)) return
                if (dn <= 0 .or. n1 < n0) then
                    call fesm_error(proc, ctx//': range "'//tok//'" needs dt > 0 and t1 >= t0.', stat=stat)
                    return
                end if
                ! Points inside [time_init, time_end] only.
                k0 = max(0_int64, -floor_div(-(clk%n_init - n0), dn))
                k1 = floor_div(min(n1, clk%n_end) - n0, dn)
                if (k1 > k0 .and. modulo(dn, clk%dt) /= 0) then
                    n = max(clk%dt, floor_div(dn, clk%dt)*clk%dt)
                    call fesm_error(proc, ctx//': range "'//tok//'" in '//src//": dt = "//dt_str(dn, clk%tpy)// &
                        " yr is not a multiple of dtt = "//dt_str(clk%dt, clk%tpy)//" yr.", &
                        "nearest: "//pair_str(n, n + clk%dt, clk%tpy), stat)
                    return
                end if
                do k = k0, k1
                    call push(lst, nl, n0 + k*dn)
                end do
            end if
        end do

    end subroutine parse_times

    subroutine parse_time_file(clk, filename, proc, ctx, lst, nl, stat)
        ! Lines of values or t0:dt:t1; blank lines and # comments skipped.
        type(clock_t),    intent(in)              :: clk
        character(len=*), intent(in)              :: filename, proc, ctx
        integer(int64),   intent(inout), allocatable :: lst(:)
        integer,          intent(inout)           :: nl
        integer,          intent(out),   optional :: stat

        character(len=1024) :: line
        integer :: io, ios, ic, iline
        logical :: ex

        if (present(stat)) stat = 0
        inquire(file=filename, exist=ex)
        if (.not. ex) then
            call fesm_error(proc, ctx//": time file not found.", "file = "//filename, stat)
            return
        end if
        open(newunit=io, file=filename, status="old", action="read")
        iline = 0
        do
            read(io, "(a)", iostat=ios) line
            if (ios /= 0) exit
            iline = iline + 1
            ic = index(line, "#")
            if (ic > 0) line(ic:) = ""
            if (len_trim(line) == 0) cycle
            call parse_times(clk, trim(line), proc, ctx, filename//" line "//to_str(iline), lst, nl, stat)
            if (failed(stat)) exit
        end do
        close(io)

    end subroutine parse_time_file

    ! =========================================================================
    ! Internal: registry
    ! =========================================================================

    subroutine register(clk, kind, name, mode, dt, phase, method, id, n)
        ! Add (or replace, same kind and name) a registry row.
        type(clock_t),    intent(inout)        :: clk
        character(len=*), intent(in)           :: kind, name, method
        integer,          intent(in)           :: mode
        integer(int64),   intent(in)           :: dt, phase
        integer,          intent(out)          :: id
        integer(int64),   intent(in), optional :: n(:)

        type(clock_entry_t), allocatable :: tmp(:)
        integer :: i, m

        if (.not. allocated(clk%reg)) allocate(clk%reg(0))
        id = 0
        do i = 1, size(clk%reg)
            if (clk%reg(i)%kind == kind .and. clk%reg(i)%name == name) id = i
        end do
        if (id == 0) then
            m = size(clk%reg)
            allocate(tmp(m + 1))
            tmp(1:m) = clk%reg
            call move_alloc(tmp, clk%reg)
            id = m + 1
        end if
        clk%reg(id)%name   = name
        clk%reg(id)%kind   = kind
        clk%reg(id)%method = method
        clk%reg(id)%mode   = mode
        clk%reg(id)%dt     = dt
        clk%reg(id)%phase  = phase
        if (allocated(clk%reg(id)%n)) deallocate(clk%reg(id)%n)
        if (present(n)) then
            clk%reg(id)%n = n
        else
            allocate(clk%reg(id)%n(0))
        end if

    end subroutine register

    integer function offending(clk, n)
        ! First registered every-dt cadence whose grid misses n, or 0.
        type(clock_t),  intent(in) :: clk
        integer(int64), intent(in) :: n
        integer :: i
        offending = 0
        if (.not. allocated(clk%reg)) return
        do i = 1, size(clk%reg)
            if (clk%reg(i)%kind /= "cadence" .or. clk%reg(i)%mode /= CADENCE_EVERY) cycle
            if (modulo(n - clk%reg(i)%phase, clk%reg(i)%dt) /= 0) then
                offending = i
                return
            end if
        end do
    end function offending

    logical function has_active_restart(clk)
        type(clock_t), intent(in) :: clk
        integer :: i
        has_active_restart = .false.
        if (.not. allocated(clk%reg)) return
        do i = 1, size(clk%reg)
            if (clk%reg(i)%kind == "restart" .and. clk%reg(i)%method /= "none") has_active_restart = .true.
        end do
    end function has_active_restart

    subroutine restart_grid(clk, p, d, ok)
        ! Intersection of the every-dt cadence grids (CRT): p + k*d; the step
        ! grid when there is none. ok = .false.: empty intersection.
        type(clock_t),  intent(in)  :: clk
        integer(int64), intent(out) :: p, d
        logical,        intent(out) :: ok
        integer :: i
        logical :: first

        d  = max(clk%dt, 1_int64)
        p  = modulo(clk%n_init, d)
        ok = .true.
        first = .true.
        if (.not. allocated(clk%reg)) return
        do i = 1, size(clk%reg)
            if (clk%reg(i)%kind /= "cadence" .or. clk%reg(i)%mode /= CADENCE_EVERY) cycle
            if (first) then
                p = clk%reg(i)%phase
                d = clk%reg(i)%dt
                first = .false.
            else
                call crt(p, d, clk%reg(i)%phase, clk%reg(i)%dt, ok)
                if (.not. ok) return
            end if
        end do
    end subroutine restart_grid

    subroutine crt(p, d, p2, d2, ok)
        ! {p + k*d} intersect {p2 + k*d2} -> {p + k*d}, 0 <= p < d.
        integer(int64), intent(inout) :: p, d
        integer(int64), intent(in)    :: p2, d2
        logical,        intent(out)   :: ok
        integer(int64) :: g, m2, t, a

        g  = gcd64(d, d2)
        ok = modulo(p2 - p, g) == 0
        if (.not. ok) return
        m2 = d2/g
        if (d/g > huge(1_int64)/d2) then
            ok = .false.                         ! period beyond int64: no common time in range
            return
        end if
        a = modulo(d/g, m2)
        t = mulmod(modulo((p2 - p)/g, m2), inv_mod(a, m2), m2)
        p = p + d*t
        d = d/g*d2
        p = modulo(p, d)
    end subroutine crt

    ! =========================================================================
    ! Internal: conversion and checks
    ! =========================================================================

    subroutine conv_ticks(tpy, years, proc, ctx, name, n, stat, per_year_hint)
        ! years -> ticks with the 4-ulp test; message with the nearest
        ! representable values on failure.
        integer(int64),   intent(in)            :: tpy
        real(dp),         intent(in)            :: years
        character(len=*), intent(in)            :: proc, ctx, name
        integer(int64),   intent(out)           :: n
        integer,          intent(out), optional :: stat
        logical,          intent(in),  optional :: per_year_hint

        character(len=:), allocatable :: what, det
        real(dp) :: x
        integer(int64) :: n0
        logical :: hint

        if (present(stat)) stat = 0
        n = 0
        what = name
        if (len(ctx) > 0) what = ctx//": "//name
        x = years*real(tpy, dp)
        if (.not. (abs(x) < real(N_EXACT, dp))) then
            call fesm_error(proc, what//" is outside the exact range of the clock.", &
                name//" = "//to_str(years)//" yr; |t| < 2^53 ticks = "// &
                yr_str(real(N_EXACT, dp)/real(tpy, dp), 0)//" yr", stat)
            return
        end if
        n = nint(x, int64)
        if (abs(years - ticks_to_years(n, tpy)) <= 4*spacing(years)) return

        n0 = floor(x, int64)
        det = name//" = "//yr_str(years, 8)//" yr, 1 tick = 1/"//to_str(tpy)//" yr ("//name//" = "// &
              yr_str(x, 3)//" ticks)"//new_line("a")//"nearest representable: "// &
              yr_str(ticks_to_years(n0, tpy), 8)//" yr ("//to_str(n0)//" ticks), "// &
              yr_str(ticks_to_years(n0 + 1, tpy), 8)//" yr ("//to_str(n0 + 1)//" ticks)"
        hint = .false.
        if (present(per_year_hint)) hint = per_year_hint
        det = det//fraction_hint(years, tpy, hint)
        n = 0
        call fesm_error(proc, what//" is not a whole number of ticks.", det, stat)

    end subroutine conv_ticks

    function fraction_hint(years, tpy, per_year) result(s)
        ! Hint lines when `years` is a truncated simple fraction k/m (m <= 100).
        real(dp),       intent(in) :: years
        integer(int64), intent(in) :: tpy
        logical,        intent(in) :: per_year
        character(len=:), allocatable :: s
        character(len=40) :: buf
        integer(int64) :: k, m

        s = ""
        do m = 2, 100
            k = nint(years*real(m, dp), int64)
            if (k == 0) cycle
            if (abs(years - real(k, dp)/real(m, dp)) > 1.0e-4_dp*max(1.0_dp, abs(years))) cycle
            if (modulo(tpy, m) /= 0) return
            write(buf, "(g0.16)") real(k, dp)/real(m, dp)
            s = new_line("a")//to_str(k)//"/"//to_str(m)//" yr is a whole number of ticks: write it with "// &
                "16 significant digits, "//trim(adjustl(buf))
            if (per_year .and. k == 1) then
                select case (m)
                case (2)
                    s = s//new_line("a")//"for halves of a year use per_year = 2"
                case (3)
                    s = s//new_line("a")//"for thirds of a year use per_year = 3"
                case (4)
                    s = s//new_line("a")//"for quarters of a year use per_year = 4"
                case default
                    s = s//new_line("a")//"for steps of 1/"//to_str(m)//" yr use per_year = "//to_str(m)
                end select
            end if
            return
        end do

    end function fraction_hint

    logical function is_initialised(clk, proc, stat)
        type(clock_t),    intent(in)            :: clk
        character(len=*), intent(in)            :: proc
        integer,          intent(out), optional :: stat
        if (present(stat)) stat = 0
        is_initialised = clk%dt > 0
        if (.not. is_initialised) call fesm_error(proc, "the clock is not initialised (clock_init).", stat=stat)
    end function is_initialised

    logical function can_register(clk, proc, ctx, stat)
        type(clock_t),    intent(in)            :: clk
        character(len=*), intent(in)            :: proc, ctx
        integer,          intent(out), optional :: stat
        can_register = is_initialised(clk, proc, stat)
        if (.not. can_register) return
        if (clk%started) then
            can_register = .false.
            call fesm_error(proc, ctx//": registered after clock_start.", &
                "initialise every cadence and schedule before clock_start.", stat)
        end if
    end function can_register

    logical function has_calendar(clk, proc, stat)
        type(clock_t),    intent(in)            :: clk
        character(len=*), intent(in)            :: proc
        integer,          intent(out), optional :: stat
        if (present(stat)) stat = 0
        has_calendar = clk%days_per_year > 0
        if (.not. has_calendar) call fesm_error(proc, "calendar = "//trim(clk%calendar)// &
            ": no day or month is defined.", "set calendar = 360_day or 365_day for calendar queries.", stat)
    end function has_calendar

    pure integer function day_of_year(clk)
        type(clock_t), intent(in) :: clk
        day_of_year = int(modulo(clk%n, clk%tpy)*int(clk%days_per_year, int64)/clk%tpy) + 1
    end function day_of_year

    pure integer function month_of(clk, doy)
        type(clock_t), intent(in) :: clk
        integer,       intent(in) :: doy
        integer :: m
        if (clk%days_per_year == 360) then
            month_of = (doy - 1)/30 + 1
        else
            month_of = 12
            do m = 1, 12
                if (doy <= MONTH_END_365(m)) then
                    month_of = m
                    return
                end if
            end do
        end if
    end function month_of

    pure integer(int64) function dt_code(c)
        ! Stored cadence dt: ticks, 0 every step, -1 off.
        type(cadence_t), intent(in) :: c
        select case (c%mode)
        case (CADENCE_EVERY)
            dt_code = c%dt
        case (CADENCE_STEP)
            dt_code = 0
        case default
            dt_code = -1
        end select
    end function dt_code

    ! =========================================================================
    ! Internal: restart file helpers
    ! =========================================================================

    subroutine ensure_file(filename)
        character(len=*), intent(in) :: filename
        logical :: ex
        inquire(file=filename, exist=ex)
        if (.not. ex) call nc_create(filename)
    end subroutine ensure_file

    logical function has_vars(filename, names, proc, stat)
        character(len=*), intent(in)            :: filename, names(:), proc
        integer,          intent(out), optional :: stat
        logical :: ex
        integer :: i
        if (present(stat)) stat = 0
        has_vars = .false.
        inquire(file=filename, exist=ex)
        if (.not. ex) then
            call fesm_error(proc, "restart file not found.", "file = "//trim(filename), stat)
            return
        end if
        do i = 1, size(names)
            if (.not. nc_exists_var(filename, trim(names(i)))) then
                call fesm_error(proc, trim(filename)//' has no variable "'//trim(names(i))//'".', stat=stat)
                return
            end if
        end do
        has_vars = .true.
    end function has_vars

    subroutine read_bundle_tpy(filename, tpy_b, proc, stat)
        character(len=*), intent(in)            :: filename, proc
        integer(int64),   intent(out)           :: tpy_b
        integer,          intent(out), optional :: stat
        real(dp) :: x
        logical :: ok
        if (present(stat)) stat = 0
        call nc_read(filename, "clock_ticks_per_year", x)
        call dp_ticks(x, tpy_b, ok)
        if (.not. ok .or. tpy_b <= 0) then
            call fesm_error(proc, trim(filename)//": clock_ticks_per_year = "//to_str(x)//" is not valid.", stat=stat)
        end if
    end subroutine read_bundle_tpy

    subroutine read_bundle_ticks(filename, name, tpy_b, tpy, n, proc, stat)
        ! Read a tick count written at tpy_b and convert it exactly to tpy.
        character(len=*), intent(in)            :: filename, name, proc
        integer(int64),   intent(in)            :: tpy_b, tpy
        integer(int64),   intent(out)           :: n
        integer,          intent(out), optional :: stat
        real(dp) :: x
        integer(int64) :: nb
        logical :: ok
        if (present(stat)) stat = 0
        n = 0
        call nc_read(filename, name, x)
        call dp_ticks(x, nb, ok)
        if (.not. ok) then
            call fesm_error(proc, trim(filename)//": "//name//" = "//to_str(x)//" is not a tick count.", stat=stat)
            return
        end if
        call convert_ticks(nb, tpy_b, tpy, n, ok)
        if (.not. ok) then
            call fesm_error(proc, trim(filename)//": "//name//" = "//to_str(nb)//" ticks at "//to_str(tpy_b)// &
                " per year cannot be converted exactly to "//to_str(tpy)//" per year.", &
                name//" = "//yr_str(ticks_to_years(nb, tpy_b), 10)//" yr", stat)
        end if
    end subroutine read_bundle_ticks

    pure subroutine dp_ticks(x, n, ok)
        ! A dp holding a whole tick count (|n| < 2^53) -> int64.
        real(dp),       intent(in)  :: x
        integer(int64), intent(out) :: n
        logical,        intent(out) :: ok
        n  = 0
        ok = abs(x) < real(N_EXACT, dp)
        if (.not. ok) return
        ok = same(x, aint(x))
        if (ok) n = int(x, int64)
    end subroutine dp_ticks

    pure subroutine convert_ticks(nb, tpy_b, tpy, n, ok)
        ! nb ticks at tpy_b per year -> n ticks at tpy per year, if exact.
        integer(int64), intent(in)  :: nb, tpy_b, tpy
        integer(int64), intent(out) :: n
        logical,        intent(out) :: ok
        integer(int64) :: g, a, b
        n = 0
        g = gcd64(tpy, tpy_b)
        a = tpy_b/g
        b = tpy/g
        ok = modulo(nb, a) == 0
        if (.not. ok) return
        ok = abs(nb/a) < N_EXACT/b
        if (ok) n = (nb/a)*b
    end subroutine convert_ticks

    ! =========================================================================
    ! Internal: namelist keys
    ! =========================================================================

    logical function has_key(filename, group, name, defaults_file)
        character(len=*), intent(in)           :: filename, group, name
        character(len=*), intent(in), optional :: defaults_file
        has_key = nml_has_param(filename, group, name)
        if (.not. has_key .and. present(defaults_file)) has_key = nml_has_param(defaults_file, group, name)
    end function has_key

    logical function in_defaults(group, name, defaults_file)
        character(len=*), intent(in)           :: group, name
        character(len=*), intent(in), optional :: defaults_file
        in_defaults = .false.
        if (present(defaults_file)) in_defaults = nml_has_param(defaults_file, group, name)
    end function in_defaults

    subroutine read_key_str(filename, group, name, value, defaults_file)
        character(len=*), intent(in)           :: filename, group, name
        character(len=*), intent(inout)        :: value
        character(len=*), intent(in), optional :: defaults_file
        if (in_defaults(group, name, defaults_file)) then
            call nml_read(filename, group, name, value, defaults_file=defaults_file)
        else
            call nml_read(filename, group, name, value)
        end if
    end subroutine read_key_str

    subroutine read_key_dp(filename, group, name, value, defaults_file)
        character(len=*), intent(in)           :: filename, group, name
        real(dp),         intent(inout)        :: value
        character(len=*), intent(in), optional :: defaults_file
        if (in_defaults(group, name, defaults_file)) then
            call nml_read(filename, group, name, value, defaults_file=defaults_file)
        else
            call nml_read(filename, group, name, value)
        end if
    end subroutine read_key_dp

    subroutine read_key_int(filename, group, name, value, defaults_file)
        character(len=*), intent(in)           :: filename, group, name
        integer,          intent(inout)        :: value
        character(len=*), intent(in), optional :: defaults_file
        if (in_defaults(group, name, defaults_file)) then
            call nml_read(filename, group, name, value, defaults_file=defaults_file)
        else
            call nml_read(filename, group, name, value)
        end if
    end subroutine read_key_int

    ! =========================================================================
    ! Internal: integer arithmetic, lists
    ! =========================================================================

    pure real(dp) function ticks_to_years(n, tpy)
        integer(int64), intent(in) :: n, tpy
        ticks_to_years = real(n, dp)/real(tpy, dp)
    end function ticks_to_years

    pure integer(int64) function floor_div(a, b)
        ! floor(a/b) for b > 0.
        integer(int64), intent(in) :: a, b
        floor_div = (a - modulo(a, b))/b
    end function floor_div

    pure recursive integer(int64) function gcd64(a, b) result(g)
        integer(int64), intent(in) :: a, b
        if (b == 0) then
            g = abs(a)
        else
            g = gcd64(b, modulo(a, b))
        end if
    end function gcd64

    pure integer(int64) function mulmod(a, b, m)
        ! a*b mod m without overflow, 0 <= a, b < m < 2^62.
        integer(int64), intent(in) :: a, b, m
        integer(int64) :: x, y
        mulmod = 0
        x = a
        y = b
        do while (y > 0)
            if (modulo(y, 2_int64) == 1) mulmod = modulo(mulmod + x, m)
            x = modulo(2*x, m)
            y = y/2
        end do
    end function mulmod

    pure integer(int64) function inv_mod(a, m)
        ! Inverse of a modulo m (gcd(a, m) = 1); 0 for m = 1.
        integer(int64), intent(in) :: a, m
        integer(int64) :: r0, r1, s0, s1, q, tmp
        r0 = m
        r1 = modulo(a, m)
        s0 = 0
        s1 = 1
        do while (r1 /= 0)
            q   = r0/r1
            tmp = r0 - q*r1
            r0  = r1
            r1  = tmp
            tmp = s0 - q*s1
            s0  = s1
            s1  = tmp
        end do
        inv_mod = modulo(s0, m)
    end function inv_mod

    subroutine push(lst, nl, v)
        integer(int64), allocatable, intent(inout) :: lst(:)
        integer,                     intent(inout) :: nl
        integer(int64),              intent(in)    :: v
        integer(int64), allocatable :: tmp(:)
        if (.not. allocated(lst)) allocate(lst(16))
        if (nl == size(lst)) then
            allocate(tmp(max(16, 2*size(lst))))
            tmp(1:nl) = lst(1:nl)
            call move_alloc(tmp, lst)
        end if
        nl = nl + 1
        lst(nl) = v
    end subroutine push

    subroutine sort_unique(lst, nl)
        ! Heap sort of lst(1:nl), then drop duplicates.
        integer(int64), intent(inout) :: lst(:)
        integer,        intent(inout) :: nl
        integer :: i, j, root, child, nk
        integer(int64) :: tmp

        do i = nl/2, 1, -1
            call sift(i, nl)
        end do
        do i = nl, 2, -1
            tmp = lst(1)
            lst(1) = lst(i)
            lst(i) = tmp
            call sift(1, i - 1)
        end do
        nk = min(nl, 1)
        do j = 2, nl
            if (lst(j) /= lst(nk)) then
                nk = nk + 1
                lst(nk) = lst(j)
            end if
        end do
        nl = nk

    contains

        subroutine sift(start, nend)
            integer, intent(in) :: start, nend
            root = start
            do
                child = 2*root
                if (child > nend) exit
                if (child < nend) then
                    if (lst(child + 1) > lst(child)) child = child + 1
                end if
                if (lst(root) >= lst(child)) exit
                tmp = lst(root)
                lst(root) = lst(child)
                lst(child) = tmp
                root = child
            end do
        end subroutine sift

    end subroutine sort_unique

    pure logical function is_sep(ch)
        character(len=1), intent(in) :: ch
        is_sep = ch == " " .or. ch == "," .or. ch == achar(9)
    end function is_sep

    pure logical function same(a, b)
        ! Exact equality of two (non-NaN) reals, for round-trip tests.
        real(dp), intent(in) :: a, b
        same = .not. (a < b .or. a > b)
    end function same

    logical function failed(stat)
        ! An optional stat was returned non-zero (without stat, fesm_error stops).
        integer, intent(in), optional :: stat
        failed = .false.
        if (present(stat)) failed = stat /= 0
    end function failed

    ! =========================================================================
    ! Internal: formatting
    ! =========================================================================

    function yr_str(x, maxd, mind) result(s)
        ! Shortest fixed-point form of x with mind..maxd decimals that reads
        ! back as x (else maxd decimals, trailing zeros dropped).
        real(dp), intent(in)           :: x
        integer,  intent(in), optional :: maxd, mind
        character(len=:), allocatable :: s
        character(len=64) :: buf, fmt
        real(dp) :: y
        integer :: nd, nd0, nd1, ios, i

        nd0 = 0
        nd1 = 10
        if (present(mind)) nd0 = mind
        if (present(maxd)) nd1 = max(maxd, nd0)
        do nd = nd0, nd1
            write(fmt, "(a,i0,a)") "(f0.", nd, ")"
            write(buf, fmt) x
            read(buf, *, iostat=ios) y
            if (ios == 0 .and. same(y, x)) exit
        end do
        s = trim(adjustl(buf))
        if (index(s, ".") > 0) then
            ! drop trailing zeros beyond nd0 decimals, and a bare "."
            i = len(s)
            do while (i > index(s, ".") + nd0 .and. s(i:i) == "0")
                i = i - 1
            end do
            s = s(1:i)
            if (s(len(s):len(s)) == ".") s = s(1:len(s)-1)
        end if
        if (len(s) == 0) s = "0"
        if (s(1:1) == ".") s = "0"//s
        if (len(s) > 1) then
            if (s(1:2) == "-.") s = "-0"//s(2:)
        end if
        if (s == "-0") s = "0"
    end function yr_str

    function t_str(n, tpy) result(s)
        ! A time [ticks] in years, at least one decimal (-21000.0).
        integer(int64), intent(in) :: n, tpy
        character(len=:), allocatable :: s
        s = yr_str(ticks_to_years(n, tpy), mind=1)
    end function t_str

    function dt_str(n, tpy) result(s)
        ! An interval [ticks] in years: a short decimal, else 1/m yr.
        integer(int64), intent(in) :: n, tpy
        character(len=:), allocatable :: s
        real(dp) :: x, y
        integer :: ios
        x = ticks_to_years(n, tpy)
        s = yr_str(x, 6)
        read(s, *, iostat=ios) y
        if ((ios /= 0 .or. .not. same(y, x)) .and. n > 0 .and. n < tpy) then
            if (modulo(tpy, n) == 0) s = "1/"//to_str(tpy/n)
        end if
    end function dt_str

    function pair_str(n1, n2, tpy) result(s)
        integer(int64), intent(in) :: n1, n2, tpy
        character(len=:), allocatable :: s
        s = yr_str(ticks_to_years(n1, tpy))//", "//yr_str(ticks_to_years(n2, tpy))
    end function pair_str

    function grid_str(ndt, nph, tpy) result(s)
        ! A stored cadence grid for messages.
        integer(int64), intent(in) :: ndt, nph, tpy
        character(len=:), allocatable :: s
        if (ndt > 0) then
            s = "every "//dt_str(ndt, tpy)//" yr from "//yr_str(ticks_to_years(nph, tpy))
        else if (ndt == 0) then
            s = "every step"
        else
            s = "off"
        end if
    end function grid_str

    function pad(s, w) result(r)
        character(len=*), intent(in) :: s
        integer,          intent(in) :: w
        character(len=:), allocatable :: r
        if (len(s) >= w) then
            r = s//" "
        else
            r = s//repeat(" ", w - len(s))
        end if
    end function pad

end module fesm_clock
