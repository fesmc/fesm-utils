# Changelog

All notable changes to fesm-utils are documented here.

## [Unreleased]

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
