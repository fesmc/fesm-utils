# `fesm_fields` — design note

Status: proposal for review (2026-10-24), no code yet. Row 2 of the fesm-utils
additions in `kryos/docs/design.md` §12; the exchange of §5 there. Precondition:
`fesm_log` (errors and log lines).

## 1. Purpose and scope

A generic container for fields that cross a component boundary: `dp` buffers
with the metadata an exchange layer needs (units, grid tag, validity interval,
time semantics, conservation flag), resolved to integer handles at init.

- **Model-agnostic.** fesm-utils knows no field names. The vocabulary (`H_ice`,
  `smb`, `dz_bed`, ...), the wiring tables and the lags live in kryos
  (`kryos_fields_vocab`, `kryos_wiring`).
- **No remapping, no physics, no conversion factors.** `fesm_fields` does not
  `use coupler`; the exchange in kryos combines the two (§5). Unit handling is
  limited to recognising equal units spelled differently (§3.5).
- One module, `src/fesm_fields.f90`; depends on `precision` (`sp`, `dp`),
  `fesm_log`, `ncio`. Public interfaces never use `wp`.

## 2. Types

```fortran
integer, parameter :: FIELD_INSTANT = 1     ! value at t1 (t0 = t1)
integer, parameter :: FIELD_MEAN    = 2     ! time mean over [t0, t1]
integer, parameter :: FIELD_ACCUM   = 3     ! time integral over [t0, t1]

integer, parameter :: FIELD_INTENSIVE = 1   ! per-point quantity (T_srf, z_srf); remap default bilin
integer, parameter :: FIELD_EXTENSIVE = 2   ! integral to conserve (H_ice, smb); remap default con

integer, parameter :: FIELD_NAME_LEN = 32, FIELD_UNITS_LEN = 32, FIELD_GRID_LEN = 64

type :: field_t
    ! --- Metadata: set by add, fixed after lock ------------------------------
    character(len=FIELD_NAME_LEN)  :: name  = ""   ! e.g. "smb"; unique within a state
    character(len=FIELD_UNITS_LEN) :: units = ""   ! canonical string as declared, "kg m-2 yr-1"
    character(len=FIELD_UNITS_LEN) :: ukey  = ""   ! normalised + equivalence-resolved (3.5)
    character(len=FIELD_GRID_LEN)  :: grid  = ""   ! coupler grid tag; "" for rank 0
    character(len=FIELD_NAME_LEN)  :: zdim  = ""   ! name of axis 3, rank 3 only ("month", "depth")
    integer :: rank = 0                            ! 0 | 2 | 3
    integer :: nx = 1, ny = 1, nz = 1              ! rank 0: 1,1,1; rank 2: nz = 1
    integer :: tsem = FIELD_INSTANT
    integer :: cons = FIELD_INTENSIVE
    logical :: masked = .false.                    ! carries a valid mask

    ! --- Published value: what get returns -----------------------------------
    logical  :: published = .false.
    real(dp) :: t0 = 0.0_dp, t1 = 0.0_dp           ! validity interval [yr]
    real(dp), allocatable :: v(:,:,:)              ! (nx,ny,nz), owned
    logical,  allocatable :: valid(:,:)            ! (nx,ny) iff masked

    ! --- Open accumulation window: MEAN and ACCUM only ------------------------
    real(dp), allocatable :: acc(:,:,:)            ! sum of var*(t1-t0) over puts
    logical,  allocatable :: acc_valid(:,:)        ! AND of the put masks
    real(dp) :: wsum = 0.0_dp                      ! sum of weights (t1-t0) [yr]
    real(dp) :: ta0 = 0.0_dp, ta1 = 0.0_dp         ! window; empty while wsum == 0
end type

type :: field_state_t
    character(len=64) :: owner = ""                ! "ice_grl:exports"; in every message
    integer :: n = 0
    logical :: locked = .false.
    type(field_t), allocatable :: f(:)             ! f(h), h = 1..n in declaration order
contains
    procedure :: init, add, lock, handle, has
    generic   :: put => put_0d_dp, put_0d_sp, put_2d_dp, put_2d_sp, put_3d_dp, put_3d_sp
    generic   :: get => get_0d_dp, get_0d_sp, get_2d_dp, get_2d_sp, get_3d_dp, get_3d_sp
    procedure :: finalize_interval, describe, restart_write, restart_read
end type
```

Decisions:

1. **One buffer shape.** Every field stores `v(nx,ny,nz)`; rank is metadata.
   Scalars (`bsl`, CO2, indices) are `v(1,1,1)` with `grid = ""`. Monthly
   climate (`tas`, `pr`) and ocean profiles are rank 3 with `zdim` naming
   axis 3. No rank 1 (open question 3).
2. **Components are public** (Fortran has no `protected` components). The rule
   is "read `f(h)` freely, write only through `put`". The exchange reads
   `f(h)%v` and `f(h)%grid` directly to avoid a copy (§5).
3. **Type-bound procedures** (`call ex%put(...)`), as the kryos sketch and
   `component_t` use, rather than the `phys_const_*` free-procedure style
   (open question 4).
4. **Two phases per state**: declare (`init`, `add`) then `lock`, which
   allocates every buffer once. `add` after `lock` and `put`/`get` before it
   are errors. While declaring, `f(:)` grows by doubling, but holds no
   buffers yet, so growth copies metadata only.
5. **No integer fields.** Masks cross as fractions in [0,1] (`f_ice`,
   `f_grnd`) or as the logical `valid`; integer codes stay inside backends.

## 3. Procedures

### 3.1 Declaration

```fortran
subroutine init(st, owner)
    class(field_state_t), intent(inout) :: st
    character(len=*),     intent(in)    :: owner         ! instance + role, for messages

subroutine add(st, name, units, tsem, cons, h, grid, dims, zdim, masked)
    class(field_state_t), intent(inout) :: st
    character(len=*),     intent(in)    :: name, units
    integer,              intent(in)    :: tsem, cons    ! FIELD_INSTANT|MEAN|ACCUM, FIELD_INTENSIVE|EXTENSIVE
    integer,              intent(out)   :: h             ! handle, = st%n after the call
    character(len=*), intent(in), optional :: grid       ! required iff dims present
    integer,          intent(in), optional :: dims(:)    ! absent: scalar; [nx,ny]; [nx,ny,nz]
    character(len=*), intent(in), optional :: zdim       ! rank 3 only; default "z"
    logical,          intent(in), optional :: masked     ! default .false.; not for rank 0

subroutine lock(st)                                      ! allocate v (+ acc for MEAN/ACCUM, + masks); v = 0
integer function handle(st, name)                        ! index of name; fesm_error if absent
logical function has(st, name)                           ! never errors
```

`add` returns the handle as an argument (not a function result) so that a
call with side effects never sits in an expression. Names match exactly after `trim`.

### 3.2 `put`: copy in

```fortran
subroutine put_2d_dp(st, h, var, t1, t0, valid)          ! 6 specifics: {0d,2d,3d} x {dp,sp}
    class(field_state_t), intent(inout) :: st
    integer,  intent(in) :: h
    real(dp), intent(in) :: var(:,:)                     ! real(sp) in the _sp specifics
    real(dp), intent(in) :: t1                           ! end of validity [yr]; always dp
    real(dp), intent(in), optional :: t0                 ! start; default t1
    logical,  intent(in), optional :: valid(:,:)         ! required iff field masked
```

Call by keyword: `call ex%put(c%h_smb, smb, t0=ta, t1=tb)`.
| tsem | effect of `put` |
|---|---|
| INSTANT | `v = var`, `[t0,t1]` recorded as given, `published = .true.`; replaces the previous value |
| MEAN, ACCUM | `acc += var*(t1-t0)`, `wsum += t1-t0`, `acc_valid = acc_valid .and. valid`; window extends to `[ta0, t1]`. Requires `t0 < t1` and, unless the window is empty, `t0 == ta1` exactly. `v` is untouched until `finalize_interval` |

- The contiguity test is exact equality. Times come from `fesm_clock`
  (`n*tick`, design decision 8), so a sub-stepping adapter takes its sub-step
  ends from the clock and sets the last one to `t_target`; a gap or overlap
  (for instance an exchange run twice for one producer interval) is an error,
  not a silent double count.
- A rank-3 `var` with `size(var,3) == 1` is accepted for a rank-2 field (the
  exchange works on `v(:,:,:)` throughout); otherwise ranks must match.
- `sp` input is widened with `real(var, dp)`, exact.
- Values at valid cells must be finite (`ieee_is_finite`, in the copy loop);
  values at invalid cells are stored as given and carry no meaning.

### 3.3 `get`: copy out

```fortran
subroutine get_2d_sp(st, h, var, t0, t1, valid)         ! 6 specifics, as put
    class(field_state_t), intent(in)  :: st              ! intent(in): get cannot reset anything
    integer,  intent(in)  :: h
    real(sp), intent(out) :: var(:,:)                    ! caller's array; shape must match
    real(dp), intent(out), optional :: t0, t1            ! validity of the value returned
    logical,  intent(out), optional :: valid(:,:)
```

Returns the **published** value; the open window is invisible. `dp` to `sp`
narrows with `real(v, sp)` (no range check). Error if nothing was published.

Decision 6: **accumulators reset in `finalize_interval`, never in `get`.**
One published value is read by several consumers (the ice import, the output
writer, a diagnostic), and with lag 1 it is read one step after it was
produced. Reset-on-get (the `study/fesm-utils.md` §6.2 sketch) would give the
second reader zeros. `intent(in)` on `get` makes this a compile-time property.

### 3.4 `finalize_interval`

```fortran
subroutine finalize_interval(st, t)
    class(field_state_t), intent(inout) :: st
    real(dp),             intent(in)    :: t                 ! end of the interval being closed
```

For every MEAN and ACCUM field: requires `wsum > 0` and `ta1 == t`; then
`v = acc/wsum` (MEAN) or `v = acc` (ACCUM), `valid = acc_valid`,
`[t0,t1] = [ta0,ta1]`, `published = .true.`, and the window is cleared
(`acc = 0`, `acc_valid = .true.`, `wsum = 0`). INSTANT fields are untouched.
kryos calls `exports%finalize_interval(t_target)` after each `update`, and
`imports%finalize_interval(t_target)` before it when imports accumulate (§6).

An ACCUM field is declared in integral units (`"kg m-2"`) and fed rates in
those units per year; its published value is `sum(var*dt)`.

### 3.5 Units

```fortran
pure function units_key(units) result(key)               ! lower case; drop " ", "*", "^"
function units_canonical(units) result(key)              ! units_key, then the alias table
logical function units_equivalent(a, b)                  ! units_canonical(a) == units_canonical(b)
subroutine field_check_link(src, hs, dst, hd)            ! errors on any incompatibility below
```

`units_key` moves here from `fesmdata` (`fesmdata.f90:334`, private today);
fesmdata then uses this one. Alias table (identities only, factor 1, offset 0):

| canonical | accepted spellings (before `units_key`) |
|---|---|
| `kg m-2 yr-1` | `kg m-2 a-1`, `kg/m2/yr`, `mm w.e. yr-1`, `mm we yr-1`, `mm w.e./yr`, `mm/yr w.e.` |
| `kg m-2` | `kg/m2`, `mm w.e.`, `mm we` |
| `m` | `meter`, `meters`, `metre`, `metres` |
| `m yr-1` | `m a-1`, `m/yr`, `m/a` |
| `K` | `kelvin` |
| `W m-2` | `W/m2`, `W m**-2` |
| `1` | `-`, `fraction` |

`mm w.e. = kg m-2` is an identity by definition of the water equivalent
(`rho_w = 1000`), which is why yelmox's mm w.e./yr exchange is bit-identical
under decision 3. `m i.e.`, `degC`, `mW m-2`, `s-1` rates are **not** listed:
they need a factor or an offset, and conversions belong in the adapters
(design §5). Decision 7: the exchange never scales, so it stays bit-exact.

`field_check_link` is called by the kryos wiring at init for every row and
stops on: units not equivalent; rank or `nz` differ; `cons` differs;
time semantics not in the allowed table:

| source \ destination | INSTANT | MEAN | ACCUM |
|---|---|---|---|
| INSTANT | ok | error | error |
| MEAN | ok (last interval's mean) | ok (accumulates) | error |
| ACCUM | error | error | ok (sums) |

### 3.6 `describe`

`subroutine describe(st, level)` (`level` default `LOG_INFO`).
A header (`[ice_grl:exports] 14 fields, 6.9 MB`) and one `log_msg` line per
field: handle, name, units, grid, shape, tsem, cons, validity (or
`unpublished`), mask. kryos logs every state at init, next to the wiring table.

### 3.7 Restart

```fortran
subroutine restart_write(st, filename)                   ! creates or overwrites filename
subroutine restart_read(st, filename)                    ! into a locked state with the same declaration
```

Written per field with `ncio`: `<name>` (`v`, double), `<name>_acc` (MEAN
and ACCUM), `<name>_valid` and `<name>_acc_valid` (int 0/1, masked fields);
variable attributes `units`, `tsem`, `cons`, `grid`, `published`, `t0`, `t1`,
`ta0`, `ta1`, `wsum`. Dimensions `<grid>_x`, `<grid>_y`, `<zdim>`; scalars on a
length-1 dimension `one`. `restart_read` checks every declared field is present
with the same shape, units key, `tsem` and `cons`, and restores all of the
above, so a restored state is bit-identical to the written one. Extra
variables in the file are an error (a renamed or dropped field must not pass
silently). This is what makes the exchange accumulators part of the restart
bundle (design §7): with lags, a restart time can hold a published value not
yet consumed and an import window already open.

## 4. Memory and performance policy

Decision 8: **owned `dp` copies, not pointers into backend arrays.**

- Backends differ in kind (Yelmo, FastIsostasy, chion `sp`; CISM, VILMA2 `dp`):
  a pointer cannot alias across kinds, so a `dp` exchange needs a copy anyway.
- Backend arrays sit inside derived types without `target`, and some are
  reallocated (Yelmo regridding, chion); pointers into them dangle silently.
- A lagged consumer needs the value from the previous step after the producer
  has moved on. A snapshot is the semantics, not a cost.

Rules:

- Buffers are allocated once in `lock`; `put`, `get`, `finalize_interval`
  and `restart_read` never allocate. The exchange keeps its remap target as a
  persistent workspace in the wiring row (coupler's `alloc_*` reallocates only
  on a shape change, so after the first call it does not).
- No `class(*)`, no polymorphic components (design decision 1); the 12
  `put`/`get` specifics are written out or fypp-generated like `ncio`.
- Per-call checks (handle, shape, mask, finiteness) are O(1) plus one pass
  fused with the copy; always on.
- Cost per exchange: one remap (O(links), indirect addressing) plus two
  streaming copies of `8 nx ny nz` bytes (the `put` into the destination, the
  backend's `get`; the source is read in place). Small next to the remap.
- Memory: `8 nx ny nz` bytes per field, twice for MEAN/ACCUM. A 381×381 grid:
  1.2 MB per 2-D field, 14 MB per monthly field; a typical instance holds a
  few tens of MB in its two states.

## 5. Use by the wiring layer (kryos, not fesm-utils)

```fortran
type :: wire_t                                         ! one row of a fixed wiring table
    integer :: hs, hd                                  ! handles, resolved at init
    character(len=16) :: method = ""                   ! "" = from cons: con | bilin
    integer :: lag = 0
    real(dp), allocatable :: work(:,:,:)               ! persistent remap target
end type

subroutine exchange(src, w, dst, cpl)                  ! kryos_wiring
    type(field_state_t), intent(in)    :: src
    type(wire_t),        intent(inout) :: w
    type(field_state_t), intent(inout) :: dst
    type(coupler_class), intent(inout) :: cpl
    associate (fs => src%f(w%hs), fd => dst%f(w%hd))
    if (.not. fs%published) call fesm_error("exchange", "source not published", ...)
    if (fs%rank == 0) then                              ! scalars bypass remap
        call dst%put(w%hd, fs%v(1,1,1), t0=fs%t0, t1=fs%t1)
    else
        call remap(cpl, fs%v, fs%grid, w%work, fd%grid, method=method_of(w, fs%cons))
        call dst%put(w%hd, w%work, t0=fs%t0, t1=fs%t1)  ! valid: forwarded once coupler row 3 lands
    end if
    end associate
end subroutine
```

At init kryos resolves `w%hs = src%handle(name_s)`, `w%hd = dst%handle(name_d)`
(a missing producer is the named gap of design decision 7), calls
`field_check_link`, and checks `lag` against the cadences. Transforms such as
`z_bed = z_bed_ref(grid) + con(dz_bed)` are a wiring-row option applied to
`w%work` before `put`; they are not in `fesm_fields`. Until the coupler gains
mask pass-through (design §12 row 3), wiring a masked field between two
different grids is an init error in kryos.

## 6. Interaction with the clock

- `put` sets validity: INSTANT `[t1,t1]` (or the interval given, for a MEAN
  source landing in an INSTANT import); MEAN/ACCUM the closed window after
  `finalize_interval`.
- A consumer due at `t_target` reads the last **published** value with its
  `[t0,t1]`. kryos compares `t1` with the time the row's lag implies (lag 0:
  end of the producer step at or before `t_target`; lag 1: one producer step
  earlier) and stops on a mismatch, which catches a producer that forgot to
  `put`. `fesm_fields` itself knows nothing about lags or cadences.
- Consumer slower than producer on a MEAN edge (ice `dt = 5`, surface
  `dt = 1`): the import is declared MEAN, kryos runs the exchange after every
  producer step, each `put` extends the import window by one producer
  interval, and the consumer's `imports%finalize_interval(t_target)` yields
  the mean over its own step. Equal cadences need no import accumulation
  (import declared INSTANT, `MEAN -> INSTANT` in the table). Open question 1.
- Restarts happen only where every cadence divides the restart tick; the
  state file (§3.7) holds whatever is published or open at that time.

## 7. Error cases

All through `fesm_error(proc, msg, detail)`, which ends in `error stop 1`; no
`stat` arguments (every case is a wiring or backend bug). `detail` always names
`owner` and the field.

| Case | `proc` | `msg` | `detail` adds |
|---|---|---|---|
| `add` after `lock`; `put`/`get` before | `fesm_fields_add` / `_put` | "state is locked" / "state is not locked" | |
| duplicate name | `_add` | "duplicate field name" | existing handle |
| bad `tsem`/`cons`; `dims` size not 2 or 3, or ≤ 0; `grid` missing for rank ≥ 2 or given for rank 0; `zdim` or `masked` on the wrong rank | `_add` | "invalid declaration" | the offending argument and value |
| name absent | `_handle` | "no field named 'smb'" | the declared names |
| handle outside 1..n | `_put` / `_get` | "invalid handle" | `h`, `n` |
| shape or rank mismatch | `_put` / `_get` | "shape mismatch" | declared `[nx,ny,nz]`, given shape |
| mask missing on a masked field, or given on an unmasked one | `_put` | "valid mask required" / "field has no mask" | |
| MEAN/ACCUM without `t0`, or `t0 >= t1` | `_put` | "invalid interval" | `t0`, `t1` |
| interval does not continue the window | `_put` | "gap or overlap in accumulation" | window `[ta0,ta1]`, given `[t0,t1]`, full precision |
| non-finite value | `_put` | "non-finite value" | count, first index |
| window empty or not ending at `t` | `_finalize_interval` | "accumulation does not reach t" | `ta1`, `t`, `wsum` |
| read before publish | `_get` | "field read before any value was published" | |
| incompatible link | `field_check_link` | "units differ and no equivalence is known" / "rank differs" / "time semantics incompatible" / "conservation flag differs" | both sides: owner, name, units, key, shape, tsem |
| restart mismatch | `_restart_read` | "restart file does not match the declaration" | file, field, expected vs found |

## 8. Tests (`test/fields/`)

Run by `make check` (row 10). Abort cases are separate small programs that
`make check` runs expecting a nonzero exit and the field name on stderr.

| Program | Asserts |
|---|---|
| `test_fields_declare` | handles are 1..n in order; `has` true/false; rank, shape, `zdim`, `ukey` recorded; buffers allocated only after `lock`, zeroed |
| `test_fields_put_get` | 0-D/2-D/3-D round trips: dp→dp bit-exact, sp→dp exact widening, dp→sp equals `real(x,sp)`; INSTANT validity `[t,t]`; rank-3 `nz=1` put into a rank-2 field |
| `test_fields_mean` | four puts of unequal length with dyadic values: published mean equals the weighted mean exactly; window `[ta0,ta1]`; `get` twice gives the same; window cleared; second interval independent of the first |
| `test_fields_accum` | published value equals `sum(var*dt)` exactly; units of the declaration kept |
| `test_fields_mask` | MEAN mask is the AND over the window; `valid` returned by `get`; non-finite at an invalid cell accepted, at a valid cell rejected (abort program) |
| `test_fields_units` | `units_key` examples (`"kg m**-2 yr-1"` → `"kgm-2yr-1"`); each alias row equivalent to its canonical; `m` ≠ `m yr-1`, `kg m-2 yr-1` ≠ `m yr-1`, `K` ≠ `degC` |
| `test_fields_link` | a compatible pair passes; every row of the tsem table behaves as listed; units, rank, `nz`, `cons` mismatches abort (one abort program each) |
| `test_fields_restart` | a state with published INSTANT, published MEAN and an open window written and read into a fresh identical declaration: every component bit-identical; a changed declaration (units, shape, extra field) aborts |
| `test_fields_exchange` | two states on two small grids in a `coupler`: `con` of an extensive field conserves the area integral within the coupler's tolerance; rank 0 bypasses remap bit-exactly; MEAN→MEAN with a slow consumer gives the mean over two producer intervals; the workspace is not reallocated after the first call |
| `test_fields_errors` | one abort program per row of §7 not covered above (duplicate name, add after lock, gap, overlap, read before publish, finalize short of `t`) |

## 9. Open questions for review

1. **Slow consumers of MEAN fields.** Support them now as in §6 (import
   declared MEAN, exchange after each producer step), or forbid
   `consumer dt > producer dt` on MEAN edges at init until a run needs it?
   Does any yelmox phase 1 configuration run ice with a longer `dt` than the
   surface?
2. **Restart file.** The exchange state as its own file,
   `restart-<time>/<instance>_fields.nc` (proposed: `fesm_fields` owns the
   whole file, no name clashes with backend variables), or as prefixed
   variables inside `<instance>_restart.nc`?
3. **Rank 1.** obm box values and per-basin profiles are 1-D. Add rank 1
   (no grid, no remap, like rank 0) now, or at phase 3 when obm arrives?
4. **Call style.** Type-bound (`call st%put(...)`, proposed, matches
   `component_t`) or the fesm-utils free-procedure style
   (`call fields_put(st, ...)`, as `phys_const_*`, `coupler_*`)?
