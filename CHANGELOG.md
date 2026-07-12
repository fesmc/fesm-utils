# Changelog

All notable changes to fesm-utils are documented here.

## [Unreleased]

### Changed
- **tsgen** `tsgen_init` argument `label` renamed to `group` and used verbatim as
  the namelist group name (default `"tsgen"`), instead of forcing a `tsgen_`
  prefix. Callers now fully control the group name (e.g. `&snp2_idx_at`). The one
  no-argument call path is unaffected (still reads `&tsgen`).
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

### Changed
- **tsgen** namelist: the feedback-convergence tolerance is renamed `eps` →
  `tol`; `eps` now denotes the realized noise sample (was `eta_now`
  internally). Existing tsgen namelists must rename their `eps` key to `tol`.

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
