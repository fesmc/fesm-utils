# Changelog

All notable changes to fesm-utils are documented here.

## Unreleased

### Added
- **fesm_log**: shared logging and error reporting. Levels (`LOG_DEBUG` ... `LOG_ERROR`,
  `log_set_level`), `log_set_unit`, `log_set_prefix`, `log_msg`, `fesm_warn`,
  `fesm_error(proc, msg [, detail] [, stat])` (framed message on `error_unit`, then
  `error stop 1`, or `stat = 1` and return), and `to_str`. `make test-log`.
- **fesm_clock**: one run clock, integer ticks inside (`int64`, N = 36000 per year) and
  `dp` years at every interface. `clock_t` (views `clock_time`, `_rel`, `_cal`, `_elapsed`
  for `method = rel | cal | const`; calendar `360_day | 365_day | none` and `year`
  required, no defaults), `cadence_t` (absolute grids `phase + k*dt`, every step, or
  `per_year`; `cadence_due`/`_mark`/`_interval`), `schedule_t` (`none | const | times |
  file`, keys as `timeout`). Every year value must be a whole number of ticks and every
  grid relation is checked at init (errors give the nearest valid values); restart
  bundles detect continuation vs branch and convert other N exactly. `clock_describe`
  logs the schedule table. `timestepping`/`timeout` are unchanged. `make check-clock`.
- **fesm_fields**: containers for fields exchanged between components (`field_t`,
  `field_state_t`; `dp` buffers with units, grid, validity, INSTANT/MEAN/ACCUM semantics,
  conservation flag, optional mask). Free procedures `fields_init`, `fields_add`,
  `fields_lock`, `fields_handle`, `fields_has`, `fields_put`/`fields_get` (generic over
  `dp`/`sp` and 0-D/2-D/3-D), `fields_finalize_interval`, `fields_check_link`,
  `fields_describe`, `fields_restart_write`/`_read`, and `units_key`, `units_canonical`,
  `units_equivalent`. `make check-fields` (with expected-abort programs, `check_abort`).
- **nml**: `nml_find_group(files, group [, what, stat])`, the first file of a list that
  holds a group (par file first, then the library list; stops listing every file, or
  `stat = 1`); replaces yelmox `var_defs_find`. `nml_validate(..., nbad)` returns the
  number of unknown keys.
- **nml**: multi-template `nml_validate(filename, group, defaults_files(:)
  [, defaults_groups(:), halt_on_error, nbad])`, a key is valid if any template declares
  it; unknown keys stop through `fesm_error` naming the group and every template searched
  (the single-template form is unchanged). `nml_group_keys(filename, group)` and
  `nml_group_names(filename)` list a group's keys and a file's groups. `nml` now uses
  `fesm_log`. `make check-nml`; `make check-nml-compat ref=<commit> [files=...]`
  (`test/nml/nml_compat.sh`) diffs every key read as a string vector against a reference
  commit (not part of `make check`). Run against `v1.4` on 34 yelmox/yelmo/chion
  namelist files (7248 keys): 17 keys differ, all `units` vectors with quoted
  elements containing blanks (now `"kg m-2 yr-1" "mm d-1"` gives 2 elements, before
  5). In v1.4 varslice only stored and printed `units_in/out` (conversion is
  `scaling`), so no numerical result of yelmox changes.
- **coupler**: `coupler_method(cpl, src, dst [, method, norm])`, the method `remap` would
  use (`"auto"` resolved) without building a map.
- **fesm_units**: `units_conv(from, to, scale, offset [, rho_ice, sec_year, stat])`, the
  factor and offset between two units on `units_canonical`: mass flux (kg m-2 or mm w.e.
  per s/d/yr, m w.e., m ice with `rho_ice`), degC/K, mW m-2/W m-2, km/m, velocity per
  s/d/yr, %/1, g kg-1/psu; the run's `sec_year` or 365.25 d. `make check-units`
  (`test_fesm_units`).
- **varslice**: `varslice_init_nml(..., defaults_file, defaults_group)` reads every key
  through a schema group (keys absent from the group come from the schema); the
  arguments `mode`, `remap`, `grid_src`, `time_ref` override the group; `remap = "auto"`
  remaps (bilin) only from a known grid other than the target; `par%mode` (`static |
  clim | transient`), checked against a legacy `time` key. Without the new arguments the
  legacy path is unchanged (`test_varslice_legacy`: byte-identical `var`/`time` against a
  reference recorded before the change). `make check-datasrc`.
- **varslice**: the time axis from the file when `time(1) < 0` or the `time` key is
  absent: the coordinate of the variable's last dimension of each file, CF `units` and
  `calendar` decoded (`varslice_time_decode`: `{seconds..days} since <date>` in the
  standard/gregorian, proleptic_gregorian, noleap/365_day, all_leap/366_day and 360_day
  calendars, `months since` (360_day), `years since`, BP and CE forms), minus `time_ref`,
  concatenated over multifile sets, strictly increasing (irregular axes allowed); `n_sub`
  inferred; `mode = auto` resolved from the file (CF `climatology` attribute: clim). New
  `par%time_units`, `par%calendar` (keys; `"clock"`: values as they stand), `par%time_name`.
- **varslice**: `par%units_file`, the `units` attribute of the variable (`""` without);
  stored only.
- **varslice**: `varslice_get(vs, var [, k])`, generic over `dp`/`sp` and 0-D to 3-D: a
  copy of the slice at its natural rank (as `v1`/`v2`/`v3`) into the caller's array, shape
  checked; storage stays `real(wp)`.
- **varslice**: `method = "interval_mean"`, the overlap-weighted mean over `[t0,t1]` of the
  piecewise representation of the file (cell means over CF `time_bnds`, read into
  `vs%time_bnds`, else the linear interpolant; partial segments at both ends; the first or
  last sample held outside the axis); `rep = n_sub` averages each sub-annual step over the
  years. `range_mean` is unchanged. The slab read of `varslice_update` is now one helper
  (`read_slab`) shared by both.
- **fesm_datasrc**: fields served from files into a `field_state_t` (the data backend of
  kryos). Free procedures `datasrc_init`, `datasrc_update` (`t` | `t0, t1`),
  `datasrc_fill`, `datasrc_has`, `datasrc_describe`, `datasrc_get` (0-D/2-D/3-D),
  `datasrc_end`. Names are bound explicitly to field groups (par file over `var_defs`),
  read with the schema `par/datasrc_defaults.nml` (`[datasrc_field]`; hosts keep a copy in
  `input/`): modes const/static/clim/transient, varslice (2-D/3-D, remapped once at init)
  or series (0-D), the file's time axis and units, conversion to the declared units with
  `units_conv` (units(1) vs the attribute: stop), `interp`/`interval_mean`/`exact`,
  `out_of_range = error | hold`, anomalies against a slice of the same data or a static
  group, masked fields with validity. `make check-datasrc` (§9 tests, the worked
  examples, expected-abort programs).

- **coords**: `grid_axis_class` carries the cell edges `x_bnds(0:nx)`, `y_bnds(0:ny)` and
  `explicit_bnds`, filled at `grid_init`; `grid_area` and the conservative generators read
  them (bit-identical to the previous midpoint edges). `axis_midpoint_bounds` is public.
- **coupler**: `coupler_take_grid(cpl, name, grid)` moves a grid in without a copy;
  `remap` passes `missing_value, mask2, reset, mask_pack, fill_method, filt_method,
  filt_par` through to `map_field` (the `int` specifics: `missing_value, mask2, reset,
  mask_pack`); `method="auto"` (`bilin` when refining, `con` when coarsening, by median
  model cell area).
- **coords**: `map_postprocess` (the fill/filter pass of `map_field`) is public.
- **coords**: map cache identity. New cache files carry the grid signatures and the
  method/generator/cutoff/format as global attributes (`coords_src_sig`, ...); `map_init`
  checks them on load (`on_mismatch="error"` default, or `"rebuild"`), loads files
  without signatures unverified (logged once), and always stops on a file whose sizes do
  not match the grids. Opt-in `hash_names=.true.` (`map_init`, `coupler_init`) tags the
  file names with `grid_hash8`. New module `grid_signature` (`grid_signature_str`,
  `grid_signature_match`, `grid_hash8`). `make check` runs `test_grid_signature`.
- **coords**: normalisation and coverage for conservative maps. The coords generator
  stores the full cell areas and pre-prune link sums (`wfull_src/dst`, `wcov_src/dst` in
  `weight_map_t`, `mp_wfull_*`/`mp_wcov_*` in the cache). `weight_map_apply` and the 8 real
  `map_field` specifics take `norm` (`"frac"` default | `"dst"` | `"sum"`), `frac2`
  (covered fraction, `real(dp)`), `area1`, `area2` (model cell areas). Default paths are
  unchanged. `coupler::remap` (dp/sp) passes `norm` and `frac2` through; `method="auto"`
  resolves to `con` for `norm /= "frac"`. `make check` runs `test_coupler_global`.
- **coupler**: global <-> regional exchange. Registered grids carry a model cell area
  (default plane `dx·dy` for projected grids, `grid%area·xy_conv²` otherwise;
  `coupler_set_area`), supplied as `area1`/`area2` for `remap(..., norm="dst"|"sum")`.
  `coupler_integral` (with `describe=.true.` it logs the ellipsoidal-area integral as a
  diagnostic), `coupler_merge_begin/add/end` (`coupler_merge_t`; extensive/intensive,
  owner weights, order independent, overlap check), `coupler_owner_reduce`,
  `coupler_split` (optionally as an anomaly against references).
- **coords**: Gaussian grids with bounds from the Gauss weights: `grid_init_gaussian(grid,
  name, nlon, nlat [, planet, lon180, x0])`; `grid_init(..., x_bnds=, y_bnds=)` (validated);
  for explicit-bounds lat-lon/Gauss grids `grid%area` is the latitude-band area (= the
  quadrature weight `R²Δλw_j` on a sphere). `grid_cdo` writes/reads such bounds as cdo
  `xbounds`/`ybounds`. `make check` runs `test_gauss_bounds`.
- **test**: `test_coupler_bitident` + `test/coords/bitident.sh <ref-commit>`, a bit-identity
  gate for the coupler/`map_field` apply surface against a reference commit.

### Changed
- **varslice**: time units `years ago`, `yr ago`, `a ago`, `years before present`, `yr before present`
  (and `ka ago`, `kyr ago`) are equivalent to `years BP` etc.: positive values are before
  present, year = 1950 - t (model time -t with `time_ref = 1950`). Reads e.g. FesmData
  `GLOBAL-0.5DEG_IceHistory-PaleoMIST.nc` (`years ago`).
- **htopo**: `htopo_write_step(htopo, filename, time [, prec])` takes the precisions of the
  written fields from the caller: `prec` (`htopo_write_prec_t`, one `htopo_prec_t(prec_abs,
  prec_rel)` per field `z_bed`, `H_ice`, `z_srf`, `f_grnd`, `z_sl`; 0 = lossless). The defaults
  are the former fixed values (1e-3 absolute, 1e-4 for `f_grnd`); without `prec` the file is
  byte-identical to before. `make check-htopo-write`.
- **coords**: optional sphere radius `radius` [m] for `grid_init_gaussian` and `grid_init` (opts form), e.g.
  `radius = 6371000.0_dp` (VILMA2, CLIMBER-X) instead of the table's 6.371221E6; areas and the cache
  signature use it (unchanged for a given radius). The radius is a field of the planet
  (`planet_init(planet, name [, radius])`; the name stays `"Spherical Earth"`) that grids and points
  derived from a grid, points, the cdo grid description (`# planet_radius`) and a map namelist
  (`planet_radius = <m>`) keep. Spherical planets only; the default is unchanged.
- **coords**: `grid_init_gaussian` defaults to `planet = "Spherical Earth"` (it was WGS84, against its
  documentation), so Gauss cell areas are `R^2*dlon*w_j`. Other grids still default to WGS84;
  pass `planet = "Spherical Earth"` for spherical models (lon-lat areas differ by -0.45% at the
  equator to +0.89% at the poles). Cached maps from such grids must be regenerated.
- **nml**: string vectors split at blanks outside quotes, so quoted elements keep their
  blanks: `units = "kg m-2 yr-1" "mm d-1"` gives `kg m-2 yr-1`, `mm d-1` (before, every
  quote was removed first and the value split at each blank: `kg`, `m-2`). Unquoted and
  single-string values are read as before. A comment after a group name
  (`&ice_bg   ! data ice`) is allowed (it was read as part of the name).
- **fesmdata**: `fesmdata_flux_conv` and `fesmdata_temp_ref` are calls to `units_conv`
  (identical factors, checked bit-for-bit); `fesmdata_flux_conv` also accepts `kg/m2/s`
  and `m w.e.` rates.
- **fesmdata**: `units_key` now comes from `fesm_fields` and also drops `^`, so
  `kg m^-2 s^-1` is recognised (it was rejected before); other units are unchanged.
- **coords**: `weight_map_prune` drops links below `cutoff·max(w)` of their target (was
  `cutoff·Σw`, which drops every link of a target with more than `1/cutoff` near-equal
  links). Changes newly generated `con`/bilinear maps where slivers were dropped (in the
  bit-identity harness: the lat-lon <-> polar-stereographic `con` maps and a polar-
  stereographic -> lat-lon `bilin` map); existing cache files load unchanged.
- **coords**: `conservative_spherical` drops repeated vertices of a cell polygon (a cell with
  two corners at a pole). Changes newly generated projected -> lat-lon `con` maps whose
  target polar row reaches +-90 (in the harness: a pole-covering polar-stereographic grid
  onto 2 deg and 1.5 deg lon-lat grids).
- **coords**: `conservative_spherical` skips candidate pairs whose bounding spherical caps
  do not meet (ANT-16KM -> T64 generation 7.2 s -> 3.8 s serial). Such pairs cannot
  overlap; only clip-noise links of disjoint polygons (<= 1e-10 of a cell) disappear, so
  pruned maps are unchanged except for targets whose only links were such noise (polar
  cells next to mis-built seam source cells in the harness). `bench_conservative` has an
  ANT-16KM -> T64 regime.
- **coupler**: no fixed capacity (`GRID_MAX`/`MAP_MAX` removed); grids and maps live in
  slot arrays that grow by moving their components, and maps are built from the
  registered grids without copying them. A grid read from `grid_<name>.txt` is read once.
- **coupler**: `remap` with `src == dst` is a copy (no map); `coupler_prime` with
  `src == dst` is a no-op. Not bit-identical for a caller that remapped a grid onto itself
  through the coupler (`(w·v)/w` can differ from `v` by 1 ulp); no downstream caller does
  (yelmox short-circuits before the coupler).
- **varslice, coupler, fesmdata, regions, topodata, htopo, zremap**: errors go through
  `fesm_error`; the private `*_error`/`to_str` helpers are removed. The frame header is now
  `fesm:: error in <proc>`. `coupler`, `htopo` and `zremap` errors now go to `error_unit`
  (were stdout); `coupler` stops with `error stop 1` (was `stop 1`).
- **ncio**: attribute routines on an open file: `nc_write_attr`, `nc_read_attr`,
  `nc_exists_attr` and `nc_size_attr` accept `(filename, ncid, [varname,] name, ...)`
  (`filename` is only for messages; the file stays open, in data mode); `nc_exists_var`
  and `nc_dims` take an optional `ncid`. The calls without `ncid` are unchanged.
- **fesm_fields, fesm_clock**: `fields_restart_write`/`_read`, `clock_restart_write` and
  `cadence_restart_write` open the file once and write or read through that handle (they
  opened and closed it for every variable and attribute: 81 opens for a 6-field state, 205
  for write plus read; now 1 each). Same variables, attributes, values and dimension order
  (`h5diff` finds no difference); only the order in which HDF5 lists a variable's attributes
  differs. `fields_restart_write` checks the dimension lengths before creating the file.
- **fesm_fields**: the finiteness check of `fields_put` is `abs(x) <= huge(x)` (same
  classification: NaN and +-Inf fail) instead of `ieee_is_finite`, which ifx calls out of
  line per value (3.3 vs 0.56 ns). `make check-fields` covers NaN, +-Inf, +-huge, +-0 and
  subnormals (`test_fields_nonfinite`, `fields_abort_nonfinite_*`).

### Fixed
- **ncio**: `nc_read` returns the values as stored. The reader set every value below `1e-30` in magnitude to
  exactly zero (`NCIO_TOL_UNDERFLOW`, removed), so a restart lost the tiny values a model had written (decayed
  Maxwell memories below `1e-30`). `nc_write` never flushed. `test_ncio_tiny` (double `1e-35`, `3e-120`,
  `7e-300`, float `1e-35`, `1e-38`; 0-3-D, slices, with and without `missing_value`; bit-identical) fails on the old
  reader.
- **ncio**: `nc_read` keeps the sign of a stored `-0.0`. The scale/offset step ran as `x*1+0` for every variable
  with `missing_value`, which gives `+0.0`; it now runs only if the variable is packed (`scale_factor /= 1` or
  `add_offset /= 0`). `test_ncio_tiny` covers `-0.0` (double and float, 1-D/2-D/scalar, with `missing_value`).
- **coords**: `conservative_planar` (lat-lon source onto a projected target) stopped under ifx with
  `forrtl 151, allocatable array is already allocated` at the allocate of its OpenMP region: the pre-pass that marks
  the candidate source cells allocated and freed `cidx`, `cd2`, which the region then lists as `private`. The
  pre-pass has its own arrays (`cidx0`, `cd20`). No change in results. `test_conservative_project` (no cdo; two
  lat-lon maps onto a polar-stereographic grid, constant preserved) fails on the old code with ifx.
- **variable_io**: a bad `prec` string (and the other table errors, an unknown variable, a
  rank/dims mismatch in `var_io_write`) stopped with exit status 0 (bare `stop`); they now `error stop 1`.
  `test_var_io` checks the exit status of a child process.
- **coupler**: `coupler_merge_end` logs the double-counted amount only for coverage above
  `1 + 1e-9` (it logged at INFO for round-off such as 1.0000000000000107); the hard error above
  `overlap_tol` is unchanged.
- **coords**: `conservative_spherical` builds each source cell polygon by sampling its
  edges in the source grid's plane (NSUB = 8 per edge) and inverse-projecting every
  sample; it interpolated the 4 inverse-projected corners in lon/lat, so projected cells
  across the 0/360 meridian or around the pole were mis-built (a garbage meridian strip in
  every pole-covering ANT/GRL map). The loop is source-major (each source polygon built
  once; links of every target in ascending source order for any thread count). Changes
  newly generated projected -> lon-lat/Gauss `con` maps (in the bit-identity harness:
  `bi-ps` onto `bi-ll2`, `bi-ls`, `bi-gauss`); existing cache files load unchanged.
- **coords**: `conservative_spherical` clips against each target cell as 8 convex strips in
  longitude; the cell's poleward parallel (great-circle chords) is concave and the clip
  lost the slivers behind it (1.3e-3 of a 1 deg cell at 70S, 2.4e-3 for T64). For a target
  that covers the source, `sum_k phi_jk = 1` to ~1e-12 on a sphere and ~1e-8 on WGS84.
  Target longitude bounds come from the unwrapped axis (fixes the seam cell of a lon180
  grid built from 0..360 centres) and latitude bounds are clamped to +-90. Changes the
  same newly generated maps as the previous entry. ANT-16KM -> T64: 3.8 s -> 4.3 s
  serial, 1.6 s -> 0.8 s on 10 threads.
- **coords**: `polar_stereographic` on a spherical planet is centred at the pole (with
  `phi` the standard parallel), as on the ellipsoid and in the cdo/CF description; it was
  centred at `phi`. Changes the coordinates of spherical polar-stereographic grids with
  `|phi| < 90`, and every map that uses them; WGS84 grids are unchanged.

## [v1.4] - 2026-10-23

### Added
- **tsgen**: window modifier of the value, every method (replaces snapclim's `f_hol`):
  `f' = f + w(t)*(win_fac - 1)*(f - win_ref)` for `win_t0 < t < win_t1`, only where
  `f > win_ref` with `win_above`; `w` ramps linearly over `win_taper` inside both edges.
  Applied before the noise to `f_now` (`f_mean_now` unchanged) and in `tsgen_tabulate`.
  Inactive by default (`win_t1 <= win_t0`). With `defaults_file` the `win_*` keys must
  be in the template group; without it they are optional.
- **tsgen**: the restart holds the random number generator state (`rng_seed`,
  `random_seed(get=)`) when the object draws noise; `tsgen_restart_read` restores it, so
  the noise after a restart is the noise without one. A file without the state, or with
  one of another size (compiler), keeps the current state with a message.
- **nml**: `nml_has_group(filename, group)`: does the file hold the group (block or flat
  lines); for looking a group up in a list of files.
- **tsgen**: `tsgen_init(..., defaults_file, defaults_group)`: keys not set in the par group fall
  back to a template group in a defaults file (as `nml_read`); `norm_t0`/`norm_t1` stay optional.
- **tsgen**: normalisation of a `series` into an index, `f = (S - S0)/(S1 - S0)`, with
  `norm_t0`/`norm_t1` (one time: the value there; two times: the mean over that window);
  active when both are given.
- **series**: `series_mean(ser, tw)` (exact window mean of the interpolant) and
  `series_normalize(ser, v0, v1)`. `make test-series` covers both.
- **fesmdata**: `fesmdata_flux_conv(fname, varname, units_out)`, the factor of a water flux
  from the units attribute of the variable to `units_out` (kg m-2 or mm w.e. per s, d or yr;
  1 yr = 365.25 d), and `fesmdata_temp_ref(fname, varname)`, the reference temperature of
  the units (0 for K, 273.15 for degC); both stop on missing or unknown units. `make test-units`.
- **interp2D**: `fill_nearest_levels(z, missing_value)` fills a field on levels (e.g. ocean by
  depth): each level by `fill_nearest`, else from the level above.
- **varslice**: `fill = "nearest"` fills a level without data from the level above
  (`fill_nearest_levels`; it stayed missing).
- **varslice**: `remap` from a projected grid: `grid_src` with a cdo description
  (`grid_<grid_src>.txt` next to the file, as FesmData writes it, else in `maps/`) is read
  from it (`fesmdata_grid_read`); otherwise from the lon/lat of the file as before.
  New `fesmdata_grid_desc_exists`.
- **varslice**: climatologies (a sub-annual cycle without years, `time = 0 ... nsub`, `nsub > 1`)
  of 3D fields (`to(x,y,depth,month)`; was "4D array without time dimension is not yet
  supported"); a `range_mean` update with `rep = 1` gives the mean over the cycle (also for 2D
  fields with a cycle), `varslice_sub_mean` averages it. `make test-varslice-clim`.

### Fixed
- **ncio**: `nc_read` without `count` stops when a file dimension covered by the array
  has another length than the array (it read the leading block silently, e.g. the corner
  of a finer grid); further file dimensions (e.g. time) are still read at index 1.
- **varslice**: `varslice_init_nml` with a target `grid` and no remapping stops when the
  variable's horizontal dimensions are not the grid's.
- **varslice**: data without a time dimension is read once (`varslice_update` re-read it
  on every call); callers must not modify `vs%var`.
- **gaussian_filter**: `filter_gaussian_fast` without `mask` smooths using all points
  (it returned the field unchanged); the reflected border now also fills the corners;
  points outside `mask` stay unchanged (they were filled from valid neighbours), as in
  `filter_gaussian`.
- **gaussian_filter**: `filter_gaussian` without `dx` uses `sigma` in grid points (it used
  an undefined value); `hgrad` sets the last x boundary (`du(nx,:)`, was `du(ny,:)`).

### Added
- **zremap**: elevation remapping of a 2D field (e.g. SMB, surface temperature).
  `zremap_update` builds local f(z) profiles from a reference state (f_ref, z_ref,
  mask_ref) in a Gaussian neighbourhood of width `sigma` (cf. Goelzer et al., 2020,
  with a radius instead of basins); `zremap_map` evaluates them at z_now. Tables every
  `stride` points, smoothing of the profiles along z (`sigma_z`), end gradients
  `dfdz_lo`/`dfdz_hi`. The interpolated profile P is blended with the fallback
  f_ref + dfdz*(z_now-z_ref) (dfdz_lo below z_ref, dfdz_hi above) by alpha =
  min(1, W/w_full), W the data support (kernel weight of the data, 0-1), so the field
  fades continuously to the fallback away from the data. On the data points P enters
  as an anomaly, f_ref + alpha*[P(z_now) - P(z_ref)] + ..., so f_ref is reproduced
  exactly at z_ref; elsewhere as alpha*P(z_now) + ... (f_ref there may be gap-filled). Optional `dz_ext` from `zremap_map`: elevation distance
  beyond the data range. `zremap_write` writes the tables, the support `w` and the
  gradient field `dfdz`. Includes a `make test-zremap` self-check.
- **varslice**: filling of missing values, optional group keys `fill = "nearest"`
  (nearest valid cell, per horizontal slice) or `fill = "const"` with `fill_value`;
  the `fill` argument of `varslice_update` overrides it.
- **timer**: `timer_init(tmr, labels)` registers the components up front;
  `timer_print_totals` and `timer_write_totals` report the time accumulated by
  each component since the reset (with its share of the total), the latter as a
  table file (header from the labels with `init=.TRUE.`, one line per call).
- **regions**: named categorical layers (regions, zones, basins, ...) in a
  `regions_class`, from any file (`regions_load_layer`; namelist `layers`,
  `path_<layer>`, ...), from the program (`regions_add_layer`), or from FesmData
  v2 files (`regions_load_fesmdata`: hierarchical `region` and `zone`;
  `regions_load_basins`: `<set>`, `<set>.group`, `<set>.mask`), with names from
  `codes`/`names`, `flag_values`/`flag_meanings` or the codes; masks from
  selection expressions (`regions_select(reg, "region:Greenland & zone:land")`,
  named masks, `regions_mask`), `regions_basin_ids`, `regions_write` (layers
  with their tables, to load again), region-code arithmetic as in FesmUtils.jl,
  and nearest-neighbour remapping onto a target grid (`regions_init(reg, grid)`).
- **topodata**: FesmData v2 topography products (`<GRID>_TOPO-<product>.nc`) in a
  `topodata_class` (`topodata_init_nml`, `topodata_init_arg`): only the listed
  variables are loaded (default `z_bed z_srf H_ice z_bed_sd`); optional remapping
  onto a target grid (real fields `con` by default, `mask`/`src_id` nearest
  neighbour). Other files (e.g. ISMIP7) via `names` (variables of the file);
  a file without `grid_name` must be on the target grid.
- **htopo**: the hi-res geometry hub of a domain, moved from yelmox and built on
  `topodata` and `regions`: `htopo_init(htopo, filename, group, group_topo,
  group_regions, domain, grid_name, cnst)` (one set of groups per domain),
  `htopo_basins`. The physics masks (yelmox's `ice_codes`/`relax_codes`/
  `region_codes`) move to the modules, as selection expressions.
  Optional preparation of the reference topography on the hub grid, from keys of the
  topography group, each off by default: `z_bed_f_sd`, `remove_lakes`, `H_ice_min`,
  `grad_lim_zb` (the steps of Yelmo's present-day topography, in that order).
- **topoprep**: `topoprep_remove_lakes` and `topoprep_limit_gradients`, ported from Yelmo
  (`remove_englacial_lakes`, `adjust_topography_gradients`) for the hub. `make check-topoprep`.
- **fesmdata**: helpers shared by the readers of FesmData v2 and other files
  (flag tables, the grid of a file from `grid_<GRID>.txt`, `fesmdata_map_init`:
  how a field gets onto a target grid, reading with remapping, path
  placeholders).

### Changed
- **timestepping**, **timeout**: model time in double precision. `tstep_class`
  times (`time`, `time_init`, `time_end`, `time_elapsed`, `time_cal`,
  `time_rel`, `time_ref`, `time_const_*`), the `tstep_init` time arguments
  (`time_init`, `time_end`, `time_ref`, `const_cal`, `const_rel`) and the
  `timeout` times (`timeout_init` `time_init`/`time_end`, `timeout_check` `time`,
  `tm%times`) and the step (`tstep_update` `dt`, `tstep_init` namelist `dtt`)
  are `real(dp)`. In single precision the spacing of the time is 0.016 yr at
  1e5 yr. Callers pass `real(dp)` times and steps.
- **tsgen**, **series**: times and time intervals in double precision. The `time`
  arguments of `tsgen_init`, `tsgen_update`, `tsgen_tabulate`, `tsgen_write`,
  `tsgen_write_step` and `tsgen_restart_write`, the state (`time_init`,
  `time_prev`, `dt`, history `time(:)`), the `dt_init`/`dt_ramp`/`dt_conv`/`dt_ave`
  parameters, the series time axis (`ser%time`) and the `time` argument of
  `series_interp`, `series_interp1`, `series_interp_sig` are `real(dp)`. Forcing
  values and rates stay `wp`.
- **varslice**, **htopo**: time in double precision. The `time` argument of
  `varslice_update`, `time_par` of `varslice_init_arg`, `time` of `print_var_range`
  (the time axis and `time_par` already were) and the times of `htopo_write_init` and
  `htopo_write_step` are `real(dp)`.
- **subgrid**: the single-precision wrappers (`calc_subgrid_array`, `_mask`, `_quad`)
  promote only the 3x3 neighbourhood of the cell, not the whole field (O(nx*ny) per
  call). Results unchanged; Yelmo ANT-8KM with `bmb_gl_method = "pmpt"`: 10-yr loop
  196 -> 120 s (16 threads).
- **varslice**: optional online remapping (`remap`, `grid_src` group keys; target
  grid passed to `varslice_init_nml`): regular lon-lat files via a cached map,
  curvilinear files (2D lon/lat, e.g. tripolar oceans) as a point set.
- **varslice**: `range_*` on a variable with a sub-annual axis selects whole years
  also with `rep = 1` (before: the range missed year 0 and returned/stopped on
  missing values). A file longer than the declared time axis is accepted (its
  leading steps are used, with a note).
- **mapping**: map caches are written under a per-process name and renamed into
  place, so concurrent runs sharing `maps/` no longer read half-written files.
- **nml**: `nml_has_param` (is a parameter set?).
- **timestepping**: `tstep_update` advances the time keepers on every call,
  including the first one. Before, the first call left `time == time_init`, so a
  driver loop of the form `do while (.not. ts%is_finished); call tstep_update(...)`
  ran its first pass as a zero-length step at `time_init`. Drivers that wrote the
  initial state on that pass should write the state at the top of the loop, before
  calling `tstep_update`, and exit there once `ts%is_finished`.
- **varslice**: `varslice_update` with a `range_*` method (`range_mean`,
  `range_sd`, `range_min`, `range_max`, `range_sum`) stops with an error when the
  requested time range is not within the time axis of the file. Before, it
  silently returned a field of missing values.

### Added
- **subgrid**: `calc_subgrid_array_quad(vint,v,nxi,i,j,im1,ip1,jm1,jp1)`, the
  subgrid values of a cell from the bilinear interpolation of `v` between cell
  centres (per quadrant: centre, face midpoints, corner mean), at the centres of
  an `nxi` x `nxi` partition of the cell. Unlike `calc_subgrid_array`, which
  interpolates between the four corner means, the field passes through `v(i,j)`,
  so a positive centre surrounded by negative neighbours keeps a positive area
  (e.g. a grounded cell next to deep ocean). Same quadrants as the grounded
  fractions of Leguy et al. (2021).
- **coords**: transverse Mercator projection (`mtype = "transverse_mercator"`,
  ellipsoidal, Snyder 1987), e.g. for UTM grids. `lambda` is the central meridian,
  `phi` the latitude of origin, the new `k0` the scale factor (default 1) and
  `x_e`/`y_n` the false easting/northing [m]. CF netCDF grids and `nc_write_map`
  (new optional `k0`) read/write its CF keys, with false easting/northing in the
  units of the axes. Grid descriptions (`grid_cdo_*`) write it as a PROJ string
  (`proj_params`; cdo has no CF `transverse_mercator`), which cdo uses for its
  weights. Within 1 mm / 1e-8 deg of pyproj in a UTM zone (`test_proj`).
- **varslice**: `varslice_sub_mean(vs)` returns the mean of the current slice
  over its time-like axis (e.g. the sub-annual cycle after `rep=varslice_nsub(vs)`)
  as an `(x,y,z)` field, for any rank: the time-like axis is axis `ndim`, so a
  2D field gives `(nx,ny,1)`. Missing values are skipped.
- **gaussian_quadrature**: optional `act` mask in `gq2D_to_nodes_aa` and
  `gq3D_to_nodes_aa` (cells carrying a value; 2D horizontal mask for 3D), as for
  the `_acx`/`_acy` routines: the corner means use only the active cells, the
  plain mean if none is active. Without `act`, unchanged.
- **timestepping**: `tstep_init` is a generic interface. Besides the explicit
  form, `tstep_init(ts, filename, group, dtt [, time_ref, cal])` reads
  `tstep_method`, `tstep_const`, `time_init`, `time_end` and `dtt` from a
  namelist group.

### Fixed
- **coords/conservative**: conservative remapping (`gen="coords"`) from a
  rotated-pole source grid (e.g. RACMO) onto a projected target now places the
  source cells correctly. The cross-system planar path treated the source axis
  values as geographic lon/lat, but for a rotated grid they are *rotated*
  lon/lat; the corners are now inverse-rotated to geographic before projection,
  so a rotated source covers its target instead of producing an empty map.
  Adds `make test-coords` self-check `test_conservative_rotated`.

## [v1.3] - 2026-07-15

### Added
- **series** (tabulated time-series reader): load a scalar (`nc=1`) or
  multi-channel (`nc>1`, e.g. 12 monthly) forcing/index curve from an ASCII
  table or netCDF variable (auto-detected by extension), with clamped-linear
  interpolation and an optional per-time standard-deviation companion. The
  lightweight 0-D cousin of `varslice`. Includes a `make test-series`
  self-check.
- **tsgen** (time-series generator): transient scalar forcing driven by time
  and/or model response. Time-driven methods (`const`, `ramp-time`,
  `ramp-time-step`, `ramp-slope`, `sin`) evaluate an analytic series; the
  `series` method reads a tabulated curve from file (via the `series` module,
  scalar or multi-channel); response-driven methods (`exp`, `PI42`, `H312b`,
  `H312PID`, `H321PID`, `PID1`) run stateful feedback controllers. For a
  multi-channel series, `f_now` is the channel mean and per-channel values are
  in `vec%f_now(:)`. Includes `tsgen_tabulate` / `tsgen_write` diagnostics and
  a `make test-tsgen` self-check. Supersedes the former (unbuilt) `hyster`
  module.
- **tsgen** `tsgen_write_step` (append a per-step 1-D diagnostic record) and
  `tsgen_restart_write` / `tsgen_restart_read` (persist and restore controller
  state across a restart).
- **root_finder** promoted into the build with hardened bracketing/convergence
  and a `make test-roots` self-check.
- **varslice** natural-rank accessors `v1`/`v2`/`v3` (1-D/2-D/3-D views onto the
  rank-4 `var`) and a `varslice_nsub` accessor for the declared sub-annual
  period.
- **distances** `compute_distance_to_mask` gains optional periodic-`x`/`y`
  support.
- **error handling**: framed, context-rich abort messages across the library.
- **docs**: Quarto documentation site with per-package pages, generated figures,
  and gh-pages publishing from `dev`.

### Changed
- **tsgen** `tsgen_init` argument `label` renamed to `group` and used verbatim as
  the namelist group name (default `"tsgen"`), instead of forcing a `tsgen_`
  prefix. Callers now fully control the group name (e.g. `&snp2_idx_at`). The one
  no-argument call path is unaffected (still reads `&tsgen`).
- **tsgen** namelist: the feedback-convergence tolerance is renamed `eps` →
  `tol`; `eps` now denotes the realized noise sample (was `eta_now`
  internally). Existing tsgen namelists must rename their `eps` key to `tol`.
- **varslice** namelist format condensed (breaking): the per-variable group now
  uses bundled keys — `units = "<in>" "<out>"`, `scaling = <scale> <offset>`,
  and `time = <active> <x0> <x1> <dx> <n_sub>` — replacing `units_in`,
  `units_out`, `unit_scale`, `unit_offset`, `with_time`, and `time_par`.
  `active` (1.0/0.0) folds in `with_time` and lets a static field document the
  period it represents (`x0..x1` preserved); `n_sub` is an optional trailing
  value. Existing namelists must be updated. Internals and the `varslice_class`
  API are unchanged.
- **varslice** `varslice_init_nml`/`varslice_par_load`/`parse_path` gain an
  optional `subs(:,2)` argument of extra `{key}`→value path substitutions
  (`subs(k,1)`=key without braces, `subs(k,2)`=value), applied after
  `{domain}`/`{grid_name}`. Non-breaking. Enables per-call tokens such as
  `{snapshot}`, `{gcm}`, `{experiment}` without a bespoke path parser.
- **ncio** `nc_write` no longer writes the `actual_range` attribute.
- **ncio** wrapper families are generated from a single `ncio.fypp` template
  (fypp 3.2 vendored under `tools/`); `ncio.f90` is now a generated artifact.
- **subgrid** `calc_subgrid_array`/`_cell`/`_mask` made precision-generic
  (single/double).
- **nml** uses `error stop` on read failures so callers exit with a nonzero
  status.
- **staggering** marked a work-in-progress stub (not exported, not built).

### Fixed
- **series** `series_load_ascii` now skips a plain (non-`#`/`!`) header row such as
  `time  value`, so ascii series files need not comment their header line.
- **varslice** `varslice_update` allocates the missing-values fallback for `ndim == 4`
  (3-D field + time) instead of leaving `vs%var` unallocated — a no-index-found update
  on a 4-D field previously segfaulted on `vs%var = mv`.
- **varslice** correctness fixes in init and time-slicing; `get_matching_files`
  now uses a per-process temp file (avoids collisions between concurrent runs).
- **varslice** avoid an `ifx` compilation error from the `get_pid()` intrinsic
  usage.

### Performance
- **ncio** skip the pre-write copy when no scale/offset packing is active, and
  skip `redef`/`enddef` when overwriting an existing variable.
- **varslice** cache per-file time lengths at init, and collapse the per-rank
  reduction into a single `(space, time)` loop.

### Removed
- **interp2D** dead code removed.
- **hyster** module removed (superseded by `tsgen`).

## [v1.2] - 2026-07-08

### Added
- **coords/mapping subsystem** (major): coordinate geometry, projections
  (oblimap, rotated-pole), Gaussian-latitude grids, and an array-backed kdtree.
- **Weight-map remapping**: nearest-neighbour, bilinear, and conservative
  kernels with SCRIP-superset weight-map file I/O (read/write + round-trip).
- **CDO interop**: read/write cdo grid descriptions, online cdo map generation,
  and cached grid→grid maps (load-if-exists, save-on-generate).
- **Multigrid coupler**: coupling layer over coords mapping (grid registry,
  lazy map cache, remap).
- Package-wide constants module; ncio informative errors, vector attributes,
  and transpose helpers; nml defaults overlay + validation.

### Changed
- **Flattened layout**: `utils/` hoisted to repo root; `fftw/`, `lis/`, `SHTns/`
  nested under container dirs.
- **Build**: configme + `make`-driven build replaces `build.py`; SHTns added as
  a build component.
- netCDF4 output now defaults to deflate level 1 (lower memory use).

### Performance
- Parallelized conservative weight generation and `weight_map_apply`.
- Analytic separable weights and structured grid-locate fast paths.
- Removed per-target allocation in apply (~10× faster remap).

### Fixed
- Conservative maps for descending-lat and cross-system (projected) targets;
  non-monotonic (lon180-wrapped) source longitude axes.
- ncio stack-temp overflow on large 3-D writes.

## [v1.1] - 2026-06-18

- Config cleanup (many configs moved to `config/legacy/`).
- Unified `build.py` driver with per-machine TOML configs.
- nml/ncio robustness fixes; lis 2.1.11 with autotools modernized for icx/ifx.

## [v0.1.1] - 2025-03-20

- `varslice`: make `get_matching_files` public; subgrid allocation fixes.

## [v0.1] - 2025-03-06

- Initial import from `climber-x-exlib`.
