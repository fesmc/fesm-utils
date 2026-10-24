# `fesm_clock` — design note

One run clock for fesm-utils: integer ticks inside, `dp` years at every
interface, a declared calendar. Settled in kryos `docs/design.md` §7, §12 row 4
and decision 8; this note fixes the interface.

**Status:** implemented (commit pending): `src/fesm_clock.f90`, tests
`test/clock/` (`make check-clock`), page `docs/packages/clock.qmd`. Sections
2, 3 and 6 show the free-procedure interface as implemented.

It replaces three things that fail today: `cadence_due` (`mod(nint(t*100),
nint(dt*100))` on `time_elapsed`: silent with non-multiples, 0.01 yr resolution,
overflow past 21 Myr), the 1e-6 tolerance of `timeout_check`, and the chunk rule
"boundaries must be multiples of `dt_clim`, `dt_isos`, `tm_*`, counted from the
start of each run" (yelmox `docs/run-chunks.md`).


**Review outcome (2026-10-24, Alex).** Free-procedure style throughout, no
type-bound procedures or `class(...)`: `clock_init` (generic: namelist /
arguments), `clock_start`, `clock_advance`, `clock_finished`, `clock_time(clk)`,
`clock_time_rel`, `clock_time_cal`, `clock_time_elapsed`, `clock_to_ticks`,
`clock_to_years`, `clock_describe`, `clock_restart_write/read`,
`clock_day_of_year` ...; `cadence_init(c, clk, ...)`, `cadence_due(c, clk)`,
`cadence_mark`, `cadence_interval`; `schedule_init`, `schedule_due`,
`schedule_times`. Q1: N = 36000. Q2: no defaults inside `fesm_clock`
(`calendar` and `year` required; the kryos defaults file sets `365_day` and
`tropical`). Q3: stop at init. Q4: continuation detected by bundle time ==
`time_init`, logged. Q5: `calendar` takes the CF meaning in kryos par files;
yelmox is not changed. §5 superseded: `tstep_class` and `timeout` stay
untouched (no wrappers, no behaviour change for yelmox); kryos uses
`fesm_clock` only; deprecation of the old two is a later, separate step.

## 1. Concepts

**Tick.** Model time is `n` ticks, `integer(int64)`, with `N` ticks per year,
fixed for the run (`ticks_per_year`, default **36000**).

| dt [yr] | 1 | 0.5 | 0.2 | 0.1 | 0.05 | 0.01 | 1/12 | 1/360 | 0.001 | 1/32 |
|---|---|---|---|---|---|---|---|---|---|---|
| ticks at N = 36000 | 36000 | 18000 | 7200 | 3600 | 1800 | 360 | 3000 | 100 | 36 | 1125 |

The smallest N with every required value whole is lcm(1, 2, 5, 10, 20, 100,
12, 360) = 1800; 36000 = 20 × 1800 adds 0.001 yr and 1/16, 1/32 yr, and makes a
360-day day 100 ticks. Range: int64 holds ±2.6e14 yr; the binding limit is
exact conversion to `dp`, |n| < 2^53, i.e. |t| < 2.5e11 yr. A ±1e8 yr paleo run
uses |n| ≤ 3.6e12. Not representable: a 365-day day (1/365 yr; see open
question 1).

**Origin.** Tick 0 is model time 0. Cadence grids and schedule times are
absolute (`phase + k·dt` from 0), never counted from `time_init`; this is what
makes chunks with different `dtt` align.

**Three views, one counter.** The clock stores only `n`; the views are exact
affine maps, so no view is ever accumulated.

| `method` | model time `time()` | `time_rel()` (forcing time) | `time_cal()` |
|---|---|---|---|
| `rel` | n/N | = time | time_rel + time_ref |
| `cal` | n/N | time − time_ref | = time |
| `const` | n/N | frozen at `time_const` | time_rel + time_ref |

`time_elapsed() = (n − n_origin)/N`, where `n_origin` is the experiment start
and survives chunk restarts (section 4.7). `time_ref` (default 1950) is stored
in ticks and must be a whole number of years, so the year fraction is the same
in every view. In `const`, the frozen forcing time is given by exactly one of
`time_const` (a `time_rel`) or `time_const_cal` (a `time_cal`).

**Calendar (a service).** Sub-annual position comes from the run's declared
calendar applied to the model-time year fraction `modulo(n, N)/N`, never from a
tick count. Only fixed-length years are supported: `360_day` (30-day months),
`365_day` (alias `noleap`), `none` (calendar queries stop with an error).
`gregorian`/`julian` are rejected at init: their year length varies, so a day
is not a fixed fraction of a year (PISM's adapter does its own CF calendar).

**Year length (a run constant).** `year = "tropical"` selects `sec_year` from
the named lengths in `phys_constants` (`360d | 365d | julian | gregorian |
tropical`). It converts rates (`kg m-2 s-1` ↔ `yr-1`) and is independent of the
calendar: a 360-day calendar on a tropical year has `sec_day =
sec_day_360d_tropical` (CLIMBER-X, smbpal).

## 2. Types

```fortran
module fesm_clock
    use, intrinsic :: iso_fortran_env, only: int64
    use precision, only: dp
    implicit none
    private

    integer(int64), parameter :: CLOCK_TICKS_PER_YEAR = 36000_int64
    integer, parameter :: CADENCE_OFF = 0, CADENCE_STEP = 1, CADENCE_EVERY = 2

    type :: clock_entry_t                     ! registry row: describe + restart checks
        character(len=32) :: name = ""
        character(len=8)  :: kind = ""        ! cadence | output | restart
        integer           :: mode = CADENCE_OFF
        integer(int64)    :: dt = 0, phase = 0
        integer(int64), allocatable :: n(:)   ! schedule times
    end type

    type, public :: clock_t
        integer(int64)    :: tpy      = CLOCK_TICKS_PER_YEAR  ! N, ticks per year
        integer(int64)    :: n        = 0     ! model time [ticks]
        integer(int64)    :: n_init   = 0     ! this run (chunk)
        integer(int64)    :: n_end    = 0
        integer(int64)    :: dt       = 0     ! dtt [ticks]
        integer(int64)    :: k        = 0     ! steps since n_init
        integer(int64)    :: n_origin = 0     ! experiment start, kept across chunks
        character(len=8)  :: method   = ""    ! rel | cal | const
        integer(int64)    :: n_ref    = 0     ! time_ref [ticks]
        integer(int64)    :: n_const  = 0     ! frozen time_rel (const) [ticks]
        character(len=16) :: calendar = ""    ! 360_day | 365_day | none
        integer           :: days_per_year = 0
        character(len=16) :: year     = ""    ! phys_constants name
        real(dp)          :: sec_year = 0.0_dp
        logical           :: started   = .false.
        logical           :: continued = .false.   ! restored at n_init from a bundle
        logical           :: branched  = .false.   ! read from a bundle at another time
        character(len=512) :: restart_file = ""
        type(clock_entry_t), allocatable :: reg(:)
    end type

    type, public :: cadence_t
        character(len=32) :: name   = ""
        integer           :: mode   = CADENCE_OFF
        integer(int64)    :: dt     = 0       ! [ticks] (CADENCE_EVERY)
        integer(int64)    :: phase  = 0       ! [ticks], 0 <= phase < dt; grid phase + k*dt
        integer(int64)    :: n_last = 0       ! last mark [ticks]
        integer           :: id     = 0       ! row in clk%reg
    end type

    type, public :: schedule_t
        character(len=32) :: name    = ""
        character(len=8)  :: method  = "none" ! none | const | times | file
        logical           :: restart = .false.
        integer(int64), allocatable :: n(:)   ! sorted, unique, on the step grid
        integer           :: k  = 1           ! cursor: first n(k) >= clk%n
        integer           :: id = 0
        integer(int64)    :: n_skip = -huge(1_int64)  ! last skipped bundle logged
    end type

    interface clock_init     ! namelist | arguments
    interface cadence_init   ! dt [, phase] | every_step | per_year [, phase]
    interface schedule_init  ! namelist | arguments
end module
```

No lead/lag field on `cadence_t`: whether an exchanged field is valid at `t0`
or `t1` of the interval is declared in the kryos wiring table (§7), not here.
The registry holds copies of definitions only; state (`n_last`, cursor) stays
in the objects the caller owns.

## 3. Public interface

All times in and out are `real(dp)` years. Ticks appear only in
`clock_to_ticks`, `clock_to_years`, `is_multiple_of_all` and restart files.
Free procedures; the object is the first argument.

| Call | Semantics |
|---|---|
| `call clock_init(clk, filename, group [, defaults_file, stat])` | Read `time_init, time_end, dtt, calendar, year` (required, no defaults), `method` (rel), `time_ref` (1950), `ticks_per_year` (36000), and for `const` exactly one of `time_const`, `time_const_cal`; a key may come from `defaults_file` (same group); validate; `n = n_origin = n_init`. |
| `call clock_init(clk, time_init, time_end, dtt, calendar, year [, method, time_const, time_const_cal, time_ref, ticks_per_year, stat])` | Same from arguments (tests, hosts without a namelist). |
| `call clock_start(clk [, stat])` | Lock the registry; with an active restart schedule, check every restart time and `time_end` against every cadence grid (an empty intersection stops); log `clock_describe`. Required before the first `clock_advance`. |
| `call clock_advance(clk)` | `n = n + dt`, `k = k + 1`. Error if `clock_finished` or not started. |
| `clock_finished(clk)` | `n >= n_end`. |
| `call clock_finish(clk)` | End early at the current step (kill switch): `n_end = n`. |
| `clock_step_index(clk)` | `k` (`int64`), steps taken in this run (resets per chunk). |
| `clock_time(clk)`, `clock_time_rel(clk)`, `clock_time_cal(clk)`, `clock_time_elapsed(clk)` | Views of section 1, `dp` years. |
| `clock_time_init(clk)`, `clock_time_end(clk)`, `clock_dtt(clk)` | Run parameters, `dp` years. |
| `clock_to_ticks(clk, years, what [, stat])` | `int64`; stops if `years` is not a whole number of ticks (section 4.4); `what` names the value in the message. |
| `clock_to_years(clk, n)` | `real(n,dp)/real(N,dp)`, correctly rounded. |
| `clock_is_representable(clk, years)` | The test of `clock_to_ticks`, without stopping. |
| `clock_year_fraction(clk)`, `clock_day_of_year(clk [, stat])`, `clock_month(clk [, stat])`, `clock_day_of_month(clk [, stat])` | Calendar position of model time; `day_of_year` in `1..days_per_year`. With calendar `none` all but `clock_year_fraction` stop. |
| `clk%sec_year`, `clock_sec_day(clk [, stat])` | Year length; `sec_year/days_per_year`. |
| `call cadence_init(c, clk, name, dt [, phase, stat])` | Cadence every `dt` yr on the absolute grid `phase + k·dt`; `dt <= 0` = off. `n_last = n_init`. |
| `call cadence_init(c, clk, name, every_step=.true. [, stat])` | Due at every step, whatever `dtt` is in this chunk. |
| `call cadence_init(c, clk, name, per_year=12 [, phase, stat])` | `dt = N/12` exactly; for sub-annual cadences without decimal literals. |
| `cadence_due(c, clk)` | `n > n_last` and (`every_step` or `modulo(n − phase, dt) == 0`). |
| `call cadence_mark(c, clk)` | `n_last = n`. Call after the component reached `t1`. |
| `call cadence_interval(c, clk, t0, t1)` | `t0 = to_years(n_last)`, `t1 = to_years(n)`. |
| `cadence_t_last(c, clk)`, `cadence_t_next(c, clk)`, `cadence_dt_years(c, clk)` | Last mark, next grid time (`huge` when off), `dt` in years. |
| `call schedule_init(s, clk, name, filename, group [, restart, defaults_file, stat])` | Schedule from keys `method, dt, times, file` (as `timeout`). |
| `call schedule_init(s, clk, name [, dt= \| times=, restart, stat])` | Same from arguments: `dt > 0` is `const`, `times` a list, neither `none`. |
| `schedule_due(s, clk)` | `n` is in the list, or `clock_finished(clk)` (section 4.2); advances the cursor (`s` is `intent(inout)`). |
| `schedule_t_next(s, clk)`, `schedule_times(s, clk)`, `schedule_size(s)` | Next time after now (`huge` if none); all times as `dp` years; count. |
| `is_multiple_of_all(n, cadences)` | Pure: `n` is on the grid of every active cadence. |
| `clock_restart_allowed(clk)` | `n` is on the grid of every registered cadence. |
| `call clock_restart_write(clk, filename)` / `clock_restart_read(clk, filename [, stat])` | Clock state; `read` sets `continued` (section 4.7). |
| `call cadence_restart_write(c, filename)` / `cadence_restart_read(c, clk, filename [, stat])` | Cadence state, keyed by `name`, in the bundle that holds the clock; read after the clock. |
| `call clock_describe(clk)` | The schedule table of section 4.9, through `log_msg`. |

Errors go through `fesm_log`: `call fesm_error("fesm_clock::"//proc, msg,
detail [, stat])`; every `init`, `clock_start`, `clock_to_ticks`, the
restart reads and the calendar queries take an optional `stat` that they
forward, so tests can check rejections. Decision 8 spells `due(cadence)`,
`from_years`, `is_output_time`, `is_restart_time`; here they are
`cadence_due(c, clk)`, `clock_to_ticks`, `schedule_due(s, clk)` on an output
or restart schedule.

## 4. Semantics

### 4.1 Grids and validation

Every check runs at init; a run that passes can never miss a time.

| Check | Where | Rule |
|---|---|---|
| Every year value (`time_init`, `time_end`, `dtt`, `time_ref`, `time_const`, `dt`, `phase`, schedule times) | `clock_to_ticks` | whole number of ticks |
| Step grid | `clock_init` | `dtt > 0`, `time_end > time_init`, `modulo(n_end − n_init, dt) == 0`, \|n\| < 2^53 |
| Cadence grid ⊂ step grid | `cadence_init` | `modulo(dt_c, dt) == 0` and `modulo(phase − n_init, dt) == 0` |
| Schedule times on the step grid | `schedule_init` | `modulo(t − n_init, dt) == 0`; times outside `[time_init, time_end]` are dropped |
| Restart times ⊂ every cadence grid | `clock_start` | with an active restart schedule: every restart time, and `time_end`, on the grid of every active cadence |

Nothing is rounded or skipped silently: `dt_isos = 15` with `dtt = 10` stops
(today it fires only at multiples of 30). `describe` prints the restart grid
(the intersection of the cadence grids by the CRT, e.g. "every 100 yr from 0");
an empty intersection (two cadences of 10 yr with phases 0 and 5) stops.

### 4.2 When things are due

- A cadence starts with `n_last = n_init`: the state at `time_init` comes from
  the component's `init` or its restart, so **no cadence is due before the
  first `advance`**. After `advance`, a cadence is due when `n` is on its grid.
  Its first interval is `[time_init, first grid time]`, shorter than `dt` when
  `time_init` is off its grid; `describe` shows it.
- `due` is false once marked (`n == n_last`), so calling it again in the same
  step, or after a restart, never re-steps: `update` stays idempotent at
  `dt = 0`.
- An output schedule includes `time_init` (`const`; `times`/`file` if listed)
  and always `time_end`; it is due at the top of the loop, before `advance`,
  as in yelmox today. A restart schedule excludes `time_init` and always
  includes `time_end`.
- After `clk%finish()` (kill switch), every output schedule is due once. A
  restart schedule is due only if `clk%restart_allowed()`; otherwise it logs
  that no bundle is written and why (offending cadence). Forcing-increment
  restarts (`coefs_restart_due`) use the same test and defer to the next
  allowed step.

### 4.3 What a component is told

Only years: the target `t1 = clk%time()` (or `time_rel()` per the clock the
component declares in kryos) and, for advance-by-`dt` backends, `t1 − t0` from
`c%interval`. Both are `n·tick` evaluated once, never sums, so `t1` of one
interval is bit-identical to `t0` of the next.

### 4.4 Converting years from a namelist

`dt_isos = 10.0` is read as `dp` by `nml_read` and converted once:

```fortran
x = years*real(clk%tpy, dp)
n = nint(x, int64)
ok = abs(years - clk%to_years(n)) <= 4*spacing(years)
```

The 4-ulp test accepts any decimal literal that denotes a whole number of
ticks (`0.1`, `-21e3`, `0.0833333333333333` for 1/12) and rejects anything off
by more than rounding; at N = 36000 two neighbouring ticks differ by far more
than 4 ulp for |t| < 1e10 yr, so there is no ambiguity. On failure:

```
fesm_clock::cadence_init: cadence "geo": dt is not a whole number of ticks
  dt = 0.33330000 yr, 1 tick = 1/36000 yr (dt = 11998.8 ticks)
  nearest representable: 0.33327778 yr (11998 ticks), 0.33330556 yr (11999 ticks)
  for thirds of a year use per_year = 3
```

The grid errors read the same way: `dt_isos = 15 yr is not a multiple of
dtt = 10 yr; nearest: 10, 20`, and `time_end = 105 yr is not on the grid of
cadence "geo" (dt = 10 yr); nearest valid time_end: 100, 110`.

### 4.5 Precision of `to_years`

`real(n,dp)/real(N,dp)`: both operands are exact in `dp` (|n| < 2^53), so the
IEEE division returns the correctly rounded value of the rational n/N — the
best any `dp` can do, and the same number the literal `0.1` parses to. The
split `real(n/N) + real(mod(n,N))/N` is worse (two roundings) and is not used.
Years are never accumulated; Kahan summation and `round_time` go.

### 4.6 Negative times

Ticks are signed; grid tests use `modulo` (floor semantics), so
`modulo(n − phase, dt)` is correct for `n < 0` and the grid is the same on
both sides of 0. Paleo runs use `method = rel` with negative `time_rel`
(years relative to `time_ref`, i.e. −21000 = 21 ka BP); kyr is a display
choice of `describe` only. `time_bp` for insolation is `-clk%time_rel()`.

### 4.7 Chunks, `dtt` changes, restarts

A chunk starts with `time_init` = previous `time_end` (run-chunks.md). Because
grids are absolute, the cadence and output times of a chain are those of one
long run whatever `dtt` each chunk uses, provided each chunk passes 4.1; the
alignment rules of run-chunks.md become init checks.

`clk%restart_read`:

- **continuation** (bundle time == `time_init`): `continued = .true.`,
  `n_origin` from the bundle, so `time_elapsed` (relax/opt windows) runs on
  across chunks instead of restarting per chunk;
- **branch** (any other time, e.g. a spin-up bundle into a transient):
  `continued = .false.`, `n_origin = n_init`, logged as a branch.

`c%restart_read` restores `n_last` on a continuation and sets `n_init` on a
branch. Since restarts are only written where every cadence fired,
`n_last == n_init` always holds on a continuation; it is stored and checked
anyway, so the rule can be relaxed later without a format change. A changed
`dt` or `phase` between chunks is allowed and logged; the new grid is checked
as in 4.1.

File content (ncio, scalar `dp`, exact since |n| < 2^53; ncio has no
`int64`): `clock_ticks_per_year`, `clock_n`, `clock_n_origin`, `clock_dtt`,
attribute `clock_method`; per cadence `cadence_<name>_{n_last,dt,phase}`. A
bundle written at another N is converted exactly or rejected.

### 4.8 Schedules

`method = none | const | times | file`, keys as `timeout`. `const` lists the
absolute multiples of `dt` in `[time_init, time_end]` (today `time_init +
k·dt`; identical when `time_init` is a multiple of `dt`). `times` is read as
one string and parsed like a file line, so `times = 0 500 1000:1000:15000`
works; `file` lines are a value or `t0:dt:t1`, `#` comments. Parsing is in
ticks (exact ranges, no `ceiling` overshoot), the list is sorted and
de-duplicated as integers, stored allocatable (no `nmax`, no 25 MB stack
vector), and `due` advances a cursor (O(1) per step).

### 4.9 `describe`

Written by `clk%start()`; one table for the whole run:

```
fesm_clock: rel, model time -21000.0 -> 0.0 yr, dtt = 5 yr (4200 steps)
  forcing  time_rel = model time; time_cal = time_rel + 1950
  calendar 365_day; year = tropical (31556926.08 s); tick = 1/36000 yr
  elapsed  origin -21000.0 (continued from restart-...)
  name      kind     dt [yr]   phase  first     last      count
  ice       cadence  step      -      -20995.0  0.0       4200
  geo       cadence  10        0      -20990.0  0.0       2100
  out_2D    output   1000      0      -21000.0  0.0       22
  restart   restart  5000      0      -20000.0  0.0       5
  restart grid: every 10 yr from 0 (geo)
```

## 5. Relationship to `tstep_class` and `timeout`

Both become thin, deprecated wrappers during migration (no new features), and
are removed once yelmox runs on kryos. Their only consumers are yelmox
(`yelmox.f90`, `kryos_*.f90`, `yelmox_bipolar.f90`, `tests/test_timestepping.f90`).
`timer` (profiling) stays separate and unchanged.

- `tstep_class` keeps its public `real(dp)` fields, which consumers read, and
  gains `type(clock_t) :: clk`. `tstep_init_nml` maps `tstep_method →
  method`, `tstep_const → time_const` (`time_const_cal` when `cal=.true.`),
  `time_ref`. `tstep_update(ts, dt)` checks `to_ticks(dt) == clk%dt`, calls
  `clk%advance()` and refreshes `time, time_cal, time_rel, time_elapsed, n,
  is_finished` from the views. Kahan, `comp_*`, `round_time` go.
- `timeout_class` wraps a `schedule_t` (it stores N); `timeout_check(tm,
  time)` converts `time` with the 4-ulp rule and tests membership. No 1e-6.
- yelmox `cadence_due` is replaced by `cadence_t`.

Behaviour changes, to be listed in the migration commit: views always satisfy
`time_cal = time_rel + time_ref` (today `rel` with `cal=.true.` freezes
`time_cal` at `tstep_const`, and `const` without `const_cal` sets `time_cal =
-1e8`); `time_end` must be reached exactly (today overshoot ends the run);
`const` output lists are absolute; `time_elapsed` continues across chunks.

## 6. Example

```fortran
call clock_init(clk, path_par, "ctrl")
call cadence_init(c_ice, clk, "ice", every_step=.true.)
call cadence_init(c_geo, clk, "geo", dt=10.0_dp)
call schedule_init(s_2D,  clk, "out_2D",  path_par, "tm_2D")
call schedule_init(s_rst, clk, "restart", path_par, "tm_rst", restart=.true.)
if (restart /= "None") call run_restore(restart, dom, clk, c_ice, c_geo)
call clock_start(clk)                                   ! validate, log the table

do
    if (schedule_due(s_2D, clk))  call write_2D(dom, clock_time(clk))
    if (schedule_due(s_rst, clk)) call run_restart_write(dom, clk, c_ice, c_geo)
    if (clock_finished(clk)) exit

    call clock_advance(clk)
    call forcing_update(frc, clock_time_rel(clk))

    if (cadence_due(c_geo, clk)) then
        call cadence_interval(c_geo, clk, t0, t1)
        call geo_update(dom%geo, t0, t1)                ! integrates [t0, t1]
        call cadence_mark(c_geo, clk)
    end if
    if (cadence_due(c_ice, clk)) then
        call ice_update(dom%ice, clock_time(clk))
        call cadence_mark(c_ice, clk)
    end if
end do
```

## 7. Tests (`test/clock/`, `make test-clock`)

| Test | Checks |
|---|---|
| `test_clock_exact` | 1e6 steps of `dtt = 0.1` from −130000: `time()` equals `to_years(n_init + k·dt)` bit for bit at every step, ends exactly on `time_end`; same with `dtt = 1/12` via `per_year` |
| `test_clock_chunks` | one run 0–200 (`dtt = 1`) against 0–100 (`dtt = 1`) + restart + 100–200 (`dtt = 0.5`): `geo` (10) and `clim` (5) fire at identical times with identical intervals; `time_elapsed` continues |
| `test_clock_reject` | `dt = 0.3333`, `dt_isos = 15` with `dtt = 10`, schedule time 105, `time_end = 105` with `geo = 10`, disjoint phases, `gregorian`: each fails with the value and the nearest valid values (via the optional `stat`) |
| `test_clock_restart` | write at 100, read: continuation restores `n_origin`, `n_last`; branch resets them; bundle at N = 1800 read at N = 36000; a non-convertible N is rejected |
| `test_clock_forcing` | `rel`, `cal`, `const` (`time_const`, `time_const_cal`): views and `time_cal = time_rel + time_ref` |
| `test_clock_negative` | −1e8 to −1e8 + 1000, `dtt = 0.5`: grids, `modulo`, `to_years(−n) == −to_years(n)` |
| `test_clock_schedule` | file with values, `t0:dt:t1`, comments, duplicates, out-of-range times; `time_end` added; off-grid time rejected; `times` string in a namelist |
| `test_clock_calendar` | `day_of_year`, `month` at known fractions for `360_day`, `365_day`; `none` stops on query |
| `test_clock_convert` | `to_ticks` of 0.1, 0.01, −21e3, 16-digit 1/12 accepted, 15-digit rejected |

## 8. Open questions

1. **Default N**: 36000 (proposed), 1800 (the minimum), or 2628000 = 36000 × 73,
   which also makes a 365-day day whole (|t| < 3.4e9 yr still exact)?
2. **Defaults for `calendar` and `year`**: `365_day` + `tropical`, `360_day`
   (smbpal, rembo, insol, CLIMBER-X), or no default, as `phys_constants` has
   none?
3. **`time_end` off the restart grid**: stop at init (proposed, follows
   decision 8, since a bundle is always written at `time_end`), or warn and
   skip the final bundle?
4. **Continuation vs branch** detected by bundle time == `time_init`
   (proposed), or an explicit key (`restart_mode = continue | branch`)?
5. **Key names**: `calendar` (day calendar, CF meaning) collides with yelmox
   `ctrl.calendar` (logical "timeline in calendar years", here `method = cal`
   / `time_const_cal`). Take the CF meaning and rename the yelmox key?
