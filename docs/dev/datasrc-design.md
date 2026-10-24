# `fesm_datasrc` — design note

Status: implemented (commit pending), 2026-10-24. Row 5 of the fesm-utils
additions in `kryos/docs/design.md` §12; the data backend of §8 there.
Preconditions: `fesm_log`, `fesm_fields` ([fields-design.md](fields-design.md):
`field_state_t`, `add`, `put(h, var, t1, t0, valid)`, public `f(h)%units,
%rank, %nx, %ny, %nz, %grid, %tsem, %cons, %masked`, `units_canonical`), and
the varslice steps of §4. Times come from `fesm_clock`
([clock-design.md](clock-design.md)).


**Review outcome (2026-10-24, Alex).** Free-procedure style throughout, no
type-bound procedures or `class(...)`: `datasrc_init`, `datasrc_update`
(generic: `t` / `t0,t1`), `datasrc_fill`, `datasrc_has`, `datasrc_describe`,
`datasrc_get` (generic), `datasrc_end`; §3 is superseded accordingly, and the
`fesm_fields` / `fesm_clock` calls follow their free-procedure names. Q1:
`series` for 0-D. Q2: explicit binding only. Q3: stop. Q4: the run's
`sec_year` when given. Q5: accessor (V4). Q6: V6 after phase 1. Q7:
`hold`/`error` only. Q8, Q9: kryos backend concerns, not datasrc. Q10:
decoding in varslice now. Q11: same-dataset slice; kryos `design.md`
decision 4 amended to say so.

**Implementation details (2026-10-24)** that refine or differ from the text
below; §3 has been updated to match.

- Steps as committed: N1, U1, V1, V2, **N2**, V3, V4, V5, datasrc. N2 (new):
  `nml` string vectors split at blanks outside quotes (`units = "kg m-2
  yr-1" "mm d-1"` was read as `kg`, `m-2`), and a comment after a group name
  is allowed. N1 also adds `nml_validate(..., nbad)`; `nml_find_group(files,
  group [, what, stat])` takes the par file as `files(1)`.
- U1: `units_conv` converts between any two units of a class (every row
  inverts), matched after `units_canonical`, rates written `<q> <t>-1` or
  `<q>/<t>`; units with the same canonical key give exactly `1, 0` (also
  unknown ones). `fesmdata_flux_conv`/`_temp_ref` give bit-identical factors
  (checked against values recorded before the change).
- V1: the schema keys read by varslice are `filename, name, units, scaling,
  time, mode, remap, grid_src, fill, fill_value, time_units, calendar`.
  varslice's own `remap = "auto"` remaps with `bilin`; datasrc resolves
  `auto` itself (`con` for `FIELD_EXTENSIVE`) and passes the result. New
  `par%mode` (`static | clim | transient`), checked against a legacy `time`.
- V2: decoding is `varslice_time_decode(t, units, calendar, y [, stat])`
  (public). `standard`/`gregorian` stop for dates before 1582-10-15;
  `proleptic_gregorian` is exact at all dates. A missing calendar attribute
  is `standard`; `julian` and others stop. With a file axis, `mode = auto`
  takes the last dimension as time when its coordinate is named `time`, has
  `axis = "T"`, `standard_name = "time"` or decodable units; the CF
  `climatology` attribute makes it `clim`. `nsub` > 1 only when every year
  of the axis holds the same number of steps; it then needs a whole-year
  `time_ref`. New `par%time_name`; bounds in `vs%time_bnds(2,nt)`.
- V4: `varslice_get` also has 1-D specifics; the caller's array is shape
  checked.
- V5: with `rep = nsub`, each year's step stands for `[year, year+1)` (cell
  weights), so a calendar year gives exactly its cycle. A single time gives
  the interpolated value (cells: the cell holding it). The slab read of
  `varslice_update` is one helper shared with `interval_mean`.
- datasrc: `cpl` is not an argument (V6 deferred). `datasrc_update` takes an
  optional `stat` too. Scalars are interpolated and averaged in dp by datasrc
  from the series samples (`series_interp` returns sp). The monthly cycle of
  a transient axis is served for `floor(t)`. Missing values left after
  `fill`: `valid = .false.` (value 0) for a masked field, else update stops.
  `ref_group` is read as a static group; `anomaly_ref` is not available for
  const fields and `group` not for scalars. Errors inside `varslice`/`nml`
  (a malformed file or group, a file not on the grid) always stop; they are
  tested as expected-abort programs (`test/datasrc/abort/`).
- Tests: §9 as listed, `test_units` as `test/units/test_fesm_units.f90`
  (`make check-units`); in `test/datasrc/` also `test_find_group` (N1),
  `test_nml_quoted` (N2), `test_varslice_schema` (V1), `test_varslice_units_attr`
  (V3), `test_varslice_get` (V4), `test_varslice_interval` (V5) and
  `test_examples` (§8 a–d). `test_varslice_legacy` compares hashes of the bit
  patterns of `vs%var`/`vs%time` for 31 cases (transient, monthly, clim,
  static with fill, multifile, remap) against `test_varslice_legacy.ref`,
  recorded before V1.

## 1. Scope

`datasrc_t` serves a list of named fields from files into a `field_state_t`.
A kryos `data` backend is `component_t + datasrc_t`: the component declares
its exports, then hands their names to `datasrc_t`.

Does: bind each name to a *field group* (file, variable, mode) found in the
par file or a shared library file; read it with `varslice` (2-D/3-D) or
`series` (0-D), with time axis and units taken from the file; remap it once
at init onto the instance grid; convert to the units declared in the state or
stop at init; evaluate at a time or over an interval, optionally as an anomaly
against a reference; `put` each field with its validity.

Does not: know component or field names (names, units, shape, tsem come from
the caller's state); do physics (no flotation, no masks from thickness, no sign
conventions: a melt file needs `scaling = -1 0`); choose a clock (§6); remap
between instances (the exchange does); hold restart state (it is a function
of time).

It replaces `ice_bnd`, `with_ice_sheet = False`, FastIsostasy `bsl method =
file|const`, the geothermal and sediments readers (`study/yelmox.md` §8 item
6) and CLIMBER-X `fake_*`, keeping their good idea (a consumer cannot tell
read from computed) without the central fake module, the per-reader time
interpolation, or data doubling as nudging targets (those are imports of the
backend that wants them, read by its own `datasrc_t`).

## 2. Namelist schema

Two levels. The **instance group** (schema owned by the kryos data backend,
`input/data_defaults.nml`, one defaults group per component) binds export
names to field groups, one key per export name, default `""` (unbound):

```fortran
&ice_bg                     ! data ice on the globe (background of design §6)
    backend = "data"
    grid    = "GLOBAL-1DEG"
    dt      = 100.0
    clock   = "forcing"
    H_ice   = "ice6g_c_H_ice"
    z_srf   = "ice6g_c_z_srf"
    f_grnd  = "ice6g_c_f_grnd"
/
```

A **field group** binds one field to a file. It is read with `defaults_group
= "datasrc_field"` of the datasrc schema (shipped as
`par/datasrc_defaults.nml`, synchronised copy in the host's `input/`, as
chion/VILMA2 do), so it holds only what differs, usually `filename` and `name`:

```fortran
&datasrc_field
    filename     = ""            ! keywords {domain} {grid_name} {key}; a glob = multifile
    name         = ""            ! variable (unused for an ascii series)
    mode         = "auto"        ! auto | const | static | clim | transient
    value        = 0.0           ! const: the value, in the declared units
    time_method  = "auto"        ! auto (from tsem) | interp | mean | exact
    out_of_range = "error"       ! error | hold
    units        = "auto" "auto" ! file units (checked), units after scaling; auto = attribute
    scaling      = 1.0 0.0       ! var*scale + offset, before unit conversion
    remap        = "auto"        ! auto | none | con | bilin | nn | shepard | quadrant
    grid_src     = ""            ! file grid; "" = its global attribute grid_name
    fill         = "none"        ! none | nearest | const
    fill_value   = 0.0
    anomaly_ref  = "none"        ! none | slice | group
    ref_time     = 0.0           ! slice: reference time [yr, instance clock]
    ref_group    = ""            ! group: a static field group of the same quantity
    time         = -1 0 0 0 0    ! legacy "active x0 x1 dx nsub"; active < 0: axis from the file
    time_units   = ""            ! "" = attribute of the time variable; "clock" = clock years
    calendar     = ""            ! "" = attribute of the time variable
/
```

- `mode = auto`: `transient` if the variable has a time dimension, else
  `static`; `clim` must be explicit (or the time variable has the CF
  `climatology` attribute). A legacy `time` with `active >= 0` fixes the mode
  (1: transient; 0 with `nsub > 1`: clim; else static); a conflicting
  explicit `mode` stops.
- `remap = auto`: remap when the file grid (`grid_src`, else `grid_name`) is
  known and differs from the instance grid, `con` for fields declared
  `FIELD_EXTENSIVE`, `bilin` otherwise; else read as is (size checked,
  `varslice_check_grid`).
- Today's library groups (`yelmox/input/clim_db.nml`: `units`, `scaling`,
  `time` all set) are valid field groups unchanged.

**Lookup.** A field group comes from the par file if it has it, else from the
first file of the run's `var_defs` list that does (yelmox `var_defs_find`,
`yelmox_nml.f90:36`, moved to `nml`, step N1); the par file never names
itself, so its groups override the library. Each field group is checked with
`nml_validate(gfile, defaults_file, group, "datasrc_field")`; the instance
group is validated by its owner.

## 3. Types and interface

```fortran
type :: datasrc_field_t                        ! private: one bound field
    character(len=32)   :: name  = ""          ! export name (the caller's)
    character(len=128)  :: group = ""          ! field group
    character(len=1024) :: gfile = "", file = ""   ! file the group was found in; data file
    integer             :: h = 0, rank = 0     ! state handle; rank 0 | 2 | 3
    integer             :: nx = 1, ny = 1, nz = 1, tsem = FIELD_INSTANT, nsub = 1
    logical             :: masked = .false., identity = .true.
    character(len=16)   :: mode = "", method = "", out_of_range = "", anomaly_ref = ""
    character(len=128)  :: units_data = "", units_decl = ""
    real(dp)            :: scale = 1.0_dp, offset = 0.0_dp   ! after scaling -> declared
    real(dp)            :: scaling(2), value = 0.0_dp      ! series scaling; const value
    real(dp)            :: trange(2) = 0.0_dp  ! coverage of the file [clock yr]
    character(len=32)   :: remap = "none"
    character(len=256)  :: grid_src = "", grid_tgt = ""
    type(varslice_class)  :: vs                ! rank 2, 3
    type(series_class)    :: ser               ! rank 0
    real(dp), allocatable :: ref(:,:,:)        ! anomaly reference, declared units
    logical,  allocatable :: ref_ok(:,:,:)
    real(dp), allocatable :: buf(:,:,:)        ! current value (nx,ny,nz); rank 0: (1,1,1)
    logical,  allocatable :: valid(:,:)        ! masked fields only
end type

type, public :: datasrc_t
    character(len=128) :: group    = ""        ! instance group
    real(dp)           :: time_ref = 0.0_dp    ! calendar year of the caller's time zero
    real(dp)           :: t0 = 0.0_dp, t1 = 0.0_dp   ! last update
    logical            :: updated  = .false.
    type(datasrc_field_t), allocatable :: f(:) ! bound fields only
end type

interface datasrc_update                       ! datasrc_update_t, datasrc_update_iv
interface datasrc_get                          ! datasrc_get_0d, _2d, _3d
```

```fortran
subroutine datasrc_init(ds, par_file, group, state, names, grid, defaults_file, &
                        required, var_defs, domain, subs, time_ref, time_span, &
                        rho_ice, sec_year, stat)
    type(datasrc_t),     intent(out) :: ds
    character(len=*),    intent(in)  :: par_file, group     ! the instance group
    type(field_state_t), intent(in)  :: state               ! declared (locked) exports
    character(len=*),    intent(in)  :: names(:)            ! fields this backend may serve
    type(grid_class),    intent(in)  :: grid                ! instance grid: {grid_name}, remap target
    character(len=*),    intent(in)  :: defaults_file       ! holds [datasrc_field]
    character(len=*),    intent(in), optional :: required(:)    ! must be bound (default: all)
    character(len=*),    intent(in), optional :: var_defs(:)    ! library files after par_file
    character(len=*),    intent(in), optional :: domain, subs(:,:)
    real(dp),            intent(in), optional :: time_ref       ! default 0: calendar years
    real(dp),            intent(in), optional :: time_span(2)   ! run interval: coverage checked here
    real(dp),            intent(in), optional :: rho_ice        ! enables "m ice yr-1"
    real(dp),            intent(in), optional :: sec_year       ! default 365.25 d (FesmData)
    integer,             intent(out), optional :: stat          ! tests: fesm_error returns

call datasrc_update(ds, t [, stat])        ! INSTANT fields at t; stops if a MEAN/ACCUM field is bound
call datasrc_update(ds, t0, t1 [, stat])   ! INSTANT at t1; MEAN/ACCUM: mean rate over [t0,t1], t0 < t1
call datasrc_fill(ds, state)               ! put every bound field with [t0,t1] of the last update
datasrc_has(ds, name)                      ! is name bound
call datasrc_describe(ds [, level])        ! per field: group, file, mode, method,
                                           ! units -> declared (scale, offset), coverage, grid_src -> grid
call datasrc_get(ds, name, var)            ! dp copy of buf (tests, backend-internal use); 0-D, 2-D, 3-D
call datasrc_end(ds)
```

**`init`, per name**: binding key in the instance group, else unbound (stop
if required: "field `z_srf` of [ice_bg] is not bound: set `z_srf =
"<group>"`"); find and validate the field group; `const` fills `buf`;
rank 0 opens a `series`, rank 2/3 a `varslice` (with the arguments of §4);
check shape and `f(h)%grid == grid%name`; derive the unit conversion (§5);
check coverage against `time_span` when `out_of_range = error`; build the
anomaly reference; read `static` and `clim` once. Errors are
`fesm_error("datasrc_init", msg, detail)` naming instance, field, group, file.

**`datasrc_update` and `datasrc_fill`**:

| mode | read | `buf` |
|---|---|---|
| const | never | `value` |
| static | at init | the field |
| clim | at init | `(nx,ny,nsub)`, or the annual mean with `time_method = mean` |
| transient | each new `(t0,t1)` | INSTANT: value at `t1`; MEAN/ACCUM: mean over `[t0,t1]` |

`fill` puts every bound field at every call with the interval of the last
update (INSTANT `t0 = t1`), static ones included, because kryos checks each
published `t1` against the time the wiring lag implies (fields-design §6).
For ACCUM fields (declared in integral units, fed rates) `buf` is the mean
rate and the units check is against `<declared> yr-1`; the state integrates.
`update` is a no-op for an unchanged `(t0,t1)` (varslice's cache), so it is
idempotent at `dt = 0` (design §7). Missing values left after `fill`: a field
declared `masked` gets `valid = .false.` there; otherwise init/update stops
with the field and time named (as `ice_bnd` does today).

**0-D.** Read with `series` (ascii or netCDF, irregular axes already; the
sea-level curves and indices are ascii today). `datasrc` replaces `ser%time`
by the decoded axis (§6), checks the range itself (`series` clamps), and
requires `units(1)` for ascii. Rank 1 is not supported (as in `fesm_fields`).

## 4. Changes to `varslice` and two neighbours

One commit each, in this order. *Legacy path* = `varslice_init_nml` without
the new optional arguments on a group with a `time` key: every yelmox call
(`ice_bnd`, `snapshots`, `climate_terms`, esm) and every `clim_db.nml` group.

| # | Change | Compatibility | yelmox bit-identity |
|---|---|---|---|
| N1 | `nml_find_group(files, group [, what, stat]) -> file`: par file first (`files(1)`), then the list; stops listing all files. yelmox `var_defs_find` becomes a call to it | new procedure | identical lookup |
| U1 | `units_conv(from, to, scale, offset, rho_ice, sec_year, stat)` in `fesm_units`, on `fesm_fields`' `units_canonical` (§5); `fesmdata_flux_conv`/`_temp_ref` call it | fesmdata API unchanged | identical factors |
| V1 | `varslice_init_nml(..., defaults_file, defaults_group)` reads every key through a schema; argument overrides `mode`, `remap`, `grid_src`, `time_ref` | without them: today's mandatory keys | unchanged path |
| V2 | Time axis from the file when `time(1) < 0` or the key is absent: read the time coordinate (the variable's last dimension) of each file, decode `units`/`calendar` (§6), concatenate over multifile, require strictly increasing. Irregular axes allowed (`get_indices` already handles them). `nsub` inferred as steps per whole year. New `par%time_units`, `par%calendar` | groups with `time(1) >= 0` unchanged | unchanged path |
| V3 | `par%units_file`: the variable's `units` attribute (`""` without); stored only, the check is `datasrc`'s | metadata | none |
| V4 | dp access: generic `varslice_get(vs, var [, k])`, dp/sp, 0-D..3-D, natural rank (as `v1/v2/v3`) | additive | none |
| V5 | Method `interval_mean`: overlap-weighted mean over `[t0,t1]` of the file's piecewise representation (cell means over CF `time_bnds` when the time variable has `bounds`, else the linear interpolant, partial segments at both ends); with `rep = nsub`; `hold` for the part outside the axis. `range_mean` (equal weights over the samples in the closed range, both ends counted, fails with none inside) stays | new method | none |
| V6 (deferred, after phase 1) | Optional `cpl`: build the map through the coupler cache (`coupler_prime`), one map per (src, dst, method) per run, grid names from the coupler registry, instead of a private `map_init` | absent: today's private map | unchanged path |

Not changed: one variable per object; `var` stays `real(wp)` (= sp), hence
V4, not dp storage: the files are float32, a km-scale bed differenced in sp
keeps ~0.1 mm, and dp storage doubles every yelmox forcing's memory (`buf`,
`ref` and the state are dp). `ncio` needs nothing (`nc_read_attr`,
`nc_exists_attr`, dp `nc_read`). The unchanged path is verified by
`test_varslice_clim`, `test_varslice_legacy` (§9) and the yelmox
Greenland-paleo case (ice_bnd + snapshots) on Levante.

## 5. Units

Declared units come from the state (`f(h)%units`; kryos: m, K, W m-2,
kg m-2 yr-1, m yr-1, 1). Per field, at init:

1. `u_file` = `units(1)` if not `auto`, else `par%units_file`; neither: stop;
   both and not equivalent: stop (a wrong attribute or group is caught).
2. `u_data` = `units(2)` if not `auto`, else `u_file`: today's meaning of
   `units = in out` with `scaling`, so `clim_db`'s `"degC" "K"` with offset
   273.15 is not converted twice.
3. `units_conv(u_data, declared, scale, offset, ...)`; failure stops ("no
   conversion from 'mm ice/yr' to 'kg m-2 yr-1' for `smb` of [surface_data],
   group ..., file ..."). Identity skips the arithmetic (bit-exact copy).

`units_canonical` (`fesm_fields`) recognises spellings; `fesm_units` adds the
factors that `fesm_fields` leaves to adapters (its own module so kryos
adapters can use it too). Small and explicit, no udunits:

| Class | From | To | Factor / offset |
|---|---|---|---|
| mass flux | kg m-2 {s,d}-1, mm [w.e.] {s,d}-1 | kg m-2 yr-1 | × `sec_year`/{1, 86400} |
| | m w.e. yr-1 | | × 1000 |
| | m i.e. yr-1, m ice yr-1 | | × `rho_ice` (absent: stop) |
| temperature | degC | K | + 273.15 |
| heat flux | mW m-2 | W m-2 | × 1e-3 |
| length | km | m | × 1000 |
| velocity | m s-1, m d-1 | m yr-1 | × `sec_year`/{1, 86400} |
| dimensionless | % | 1 | × 0.01 |
| salinity | g kg-1 | psu | 1 |

`sec_year`: the run's (clock-design §1) when passed, else 365.25 d as in
`fesmdata_flux_conv`. Gate note: `Q_geo` in mW m-2 → W m-2 in the state →
mW m-2 in the Yelmo adapter: `(x·1e-3)·1e3` in dp rounds back to the sp `x`
(dp error ≪ half an sp ulp; `test_units` sweeps sp values), so phase 1 stays
bit-identical.

## 6. Time

**Clock.** The instance key `clock = "forcing" | "model"` (kryos data
defaults, default `forcing`, as every time-dependent input in yelmox
`_forcing-design.md` §3) is read by the kryos backend. It passes
`clk%time_rel()` with `t_ref = clk%time_cal() − clk%time_rel()` (1950), or for
`model` (the esm exception: calendar periods) `clk%time_cal()` with `t_ref =
0`; `t_span` is the run interval on the same view. `datasrc` decodes each
file axis to calendar years `y` and uses `t = y − time_ref`; a legacy `time`
axis and `time_units = "clock"` are clock years as they stand.

**Decoding** (V2):

| `units` of the time variable | Calendar year |
|---|---|
| `{seconds,minutes,hours,days} since <date>` | by `calendar`: `standard`/`gregorian`/`proleptic_gregorian` (exact day count; before 1582-10-15: stop), `noleap`/`365_day`, `all_leap`/`366_day`, `360_day`; year + day-of-year/days-in-year |
| `months since <date>` | `360_day` only, else stop |
| `years since <date>`, `common_years since <date>` | year of date + value (calendar years) |
| `yr BP`, `a BP`, `years BP`, `ka BP`, `kyr BP` | 1950 − value × {1, 1000} |
| `yr CE`, `years CE`, `AD` | value |
| `yr`, `years`, `ka` with no reference | stop: set `time_units` |

File calendars are decoded independently of the run calendar (the clock
accepts only `360_day`/`365_day`; CMIP files are often `gregorian`).

**Outside the axis**, per field: `out_of_range = error` (default) stops, at
init against `time_span` and at every update; `hold` takes the first or last
sample (varslice `extrap`), for an interval mean the part outside. No linear
extrapolation: an extrapolated ice history or sea-level curve is never what
is meant.

**Method**, `time_method = auto`: INSTANT → `interp`; MEAN, ACCUM →
`interval_mean` (V5). `exact` stops when no sample matches; explicit `interp`
on a MEAN field takes the value at `(t0+t1)/2`.

**Sub-annual.** A climatology (12 steps, no years) is served as rank 3
`(nx,ny,12)` (`zdim = "month"`), step k = month k; `time_method = mean`
serves its annual mean (varslice's collapse). A 3-D climatology is rank 4:
only its annual mean. A transient monthly axis is served `(nx,ny,12)` for the
year of `t`, interpolated between years (varslice `rep = nsub`).

## 7. Grids

The instance grid is the target of the varslice remap, built once at init.
By default the kryos data backend takes its consumer's grid in the same
domain (data ice, data geothermal → the ice grid): the only remap is at init
and the exchange is an identity (design §12 row 3). An instance may set its
own `grid` (the globe's background ice on GLOBAL-1DEG, a coarse product for
several consumers); the exchange then remaps every step with the declared
`cons`. With `cpl` (V6) the maps come from the run's coupler cache and
`grid_src` may name any grid the coupler knows; without it, as today.

## 8. Worked examples

**(a) Background ice on the globe from ICE-6G_C**, instance `[ice_bg]` of §2.
Library entries (product name illustrative); the 1 kyr → 500 yr axis is read
from the file, which varslice cannot do today (`clim_db.nml:631`):

```fortran
&ice6g_c_H_ice
    filename = "fesmdata/Global/{grid_name}/{grid_name}_IceHistory-ICE6GC.nc"
    name     = "H_ice"
/
! &ice6g_c_z_srf, &ice6g_c_f_grnd alike (name = "z_srf", "f_grnd"; f_grnd: fill = "const")
```
```fortran
! ice_data_t extends ice_t with a datasrc_t `src`; ice_t has declared and locked self%exports
call datasrc_init(self%src, par_file, group, self%exports, ICE_EXPORT_NAMES, grid, run%datasrc_defaults, &
                   required=[character(len=8) :: "H_ice", "z_srf", "f_grnd"], var_defs=run%var_defs, &
                   time_ref=t_ref, time_span=t_span)
call datasrc_describe(self%src)
! update(t0, t1): geometry is INSTANT, so the values at t1
call datasrc_update(self%src, t0, t1)
call datasrc_fill(self%src, self%exports)
```

**(b) Static geothermal.** `[geothermal_grl]` binds `Q_geo = "q_geo_shapiro"`;
a constant flux would be `mode = "const"  value = 0.042` (replaces
`ghf_const`; sediments `method = 0` is `value = 0.0`):

```fortran
&q_geo_shapiro
    filename = "fesmdata/{domain}/{grid_name}/{grid_name}_GHF-Shapiro2004.nc"
    name     = "ghf"
    fill     = "nearest"    ! file units "mW m-2" -> W m-2 (x 1e-3)
/
```
```fortran
call datasrc_init(self%src, par_file, group, self%exports, ["Q_geo"], grid, run%datasrc_defaults, &
                   domain=dom%name, var_defs=run%var_defs)
call datasrc_update(self%src, t)          ! nothing read; validity [t,t]
call datasrc_fill(self%src, self%exports)
```

**(c) Scalar `bsl` curve for the globe's `forcing`.** `[forcing]` binds
`bsl = "bsl_spratt2016"`:

```fortran
&bsl_spratt2016
    filename   = "input/sealevel_spratt2016.dat"   ! ascii: time value
    units      = "m" "m"                           ! ascii has no attribute
    time_units = "clock"                           ! first column is forcing time
/
```
```fortran
call datasrc_init(self%src, par_file, "forcing", self%exports, ["bsl"], globe%grid, run%datasrc_defaults, &
                   time_ref=t_ref, time_span=t_span)
call datasrc_update(self%src, run%clk%time_rel())   ! INSTANT: series value at t
call datasrc_fill(self%src, self%exports)           ! put(h_bsl, value, t1=t)
```

**(d) Transient `geo` from a relative sea-level history, as anomalies**
(decision 4). `rsl` = sea surface − solid surface; the bed anomaly is
`−(rsl(t) − rsl(0))`, `dz_ss` and `bsl` are 0 (`rsl` holds the barystatic part):

```fortran
&geo_grl
    backend = "data"
    clock   = "forcing"
    dz_bed  = "ice6g_c_rsl_dz_bed"
    dz_ss   = "zero"
    bsl     = "zero"
/
&ice6g_c_rsl_dz_bed
    filename    = "fesmdata/{domain}/{grid_name}/{grid_name}_RSL-ICE6GC.nc"
    name        = "rsl"
    scaling     = -1.0 0.0
    anomaly_ref = "slice"
    ref_time    = 0.0
/
&zero
    mode  = "const"
    value = 0.0
/
```
```fortran
call datasrc_init(self%src, par_file, group, self%exports, [character(len=8) :: "dz_bed", "dz_ss", "bsl"], &
                   grid, run%datasrc_defaults, domain=dom%name, var_defs=run%var_defs, &
                   time_ref=t_ref, time_span=t_span)
call datasrc_update(self%src, run%clk%time_rel())
call datasrc_fill(self%src, self%exports)
```

Consumers form `z_sl − z_bed = −z_bed_ref(grid) + (rsl(t) − rsl(0))`: the
hub's high-resolution reference plus the dataset's change, never its coarse
absolute bed. One `const` group serves 0-D and 2-D fields.

## 9. Tests (`test/datasrc/`, `make test-datasrc`, part of `make check`)

Synthetic netCDF/ascii files written by the tests with `ncio`; failures via `stat`.

| Program | Checks |
|---|---|
| `test_units` | every table row and its inverse; identity bit-exact; unknown units and `m ice yr-1` without `rho_ice` give `stat /= 0`; sp round trip of mW m-2 → W m-2 → mW m-2 |
| `test_time_axis` | each calendar and `since` unit, BP/CE forms, multifile concatenation, an irregular ICE-6G-like axis, non-monotonic stops; a regular file via legacy `time` and via the file gives equal `vs%time` |
| `test_static` | `const`, `static` with `remap = auto` from `grid_name`, `fill`; shape and grid mismatch stop |
| `test_transient` | `interp` exact at samples, linear between; `interval_mean` equals the analytic mean for linear data and the overlap-weighted mean with `time_bnds`; `hold` vs `error`; coverage stop at init; repeated update reads nothing |
| `test_clim` | 12-month 2-D served `(nx,ny,12)`; `mean` collapse; rank-4 climatology only as mean |
| `test_scalar` | ascii and netCDF series; `units` required for ascii; `fill` into a `field_state_t`, validity as in §3 |
| `test_anomaly` | `slice` gives exactly 0 at `ref_time`; `group` reference; conversion applied before differencing |
| `test_bind` | lookup order (par file over library); unbound required field stops; unknown key stops; a `clim_db` group accepted unchanged; MEAN field with `update(t)` stops |
| `test_varslice_legacy` | byte-equal `vs%var`, `vs%time` for a `clim_db` excerpt before and after V1–V5 |

## 10. Open questions for review

1. **0-D through `series`** (proposed) or varslice (one reader, ascii to add)?
2. **Binding**: explicit group per field (proposed, as snapshots do), or also
   a default group `<instance>_<field>` when the key is unset?
3. **Unit mismatch** between `units(1)` and the attribute: stop (proposed), or
   allow the group to override, with a warning, for files with wrong attributes?
4. **Year length** for d/s → yr: the run's `sec_year` when given (proposed),
   or always 365.25 d to match how FesmData products were built?
5. **varslice dp storage** later, or the accessor (V4) for good?
6. **V6 now or after phase 1?** Without it each field holds its own map copy.
7. **Beyond the axis**: `hold`/`error` enough, or a `fill_value` (e.g. `bsl =
   0` after the record)?
8. **Geothermal extras** (`f_stdev` × error field, `ghf_min` clamp): a generic
   `clip = lo hi` key, and `f_stdev` in the kryos geothermal backend (it would
   read `Q_geo_err` as a backend-internal field, so `init` must accept a name
   with no state handle, with units given by the caller)?
9. **ice_bnd `scaled`** (coefficient × one slice): a generic `scale_by =
   "<scalar>"` with the scalar passed to `update`, or an experiment plugin?
   The strip (fixed/relaxation cells) belongs to the Yelmo adapter either way.
10. **Calendar decoding** in varslice now; later in `fesm_clock`?
11. **Anomaly reference**: same-dataset slice (proposed) or the domain
    reference topography, as decision 4 literally reads (which mixes datasets)?
