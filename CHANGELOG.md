# Changelog

All notable changes to fesm-utils are documented here.

## [Unreleased]

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
