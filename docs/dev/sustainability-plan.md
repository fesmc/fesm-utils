# fesm-utils — sustainability / modernization plan

Working plan for improving the `src/` library. We go step by step; each step is
its own commit (no mixing concerns). Status is updated as we land each one.

Ordering rationale: safest, highest-confidence changes first (dead code, then the
mechanical coupler fix), then the larger error-handling work.

---

## Step 1 — Delete dead code  *(done)*

Zero-caller routines in `src/coords/interp2D.f90`, verified by grep across
`src/` and `test/`:

- `filter_poisson_dble_new` (lines 1762–1873) — not in the `filter_poisson`
  interface (`_float`/`_dble` only), no references. ~112 lines.
- `interp_nearest_dble0` (lines 337–416) — not in the `interp_nearest` interface
  (`_dble`/`_int` only), no references. ~80 lines.
- `fill_nearest_dble_new` (lines 620–723) — listed in the `fill_nearest`
  interface (line 34), but the two live generic calls (`mapping.f90:1188`,
  `mapping_scrip.f90:493`) pass `var, missing_value, mask` only and bind to
  `fill_nearest_dble`; they cannot bind to `_dble_new` (mandatory
  `fill_value/fill_dist/n/dx`). Delete the body **and** interface line 34. ~104 lines.

Also sweep obviously-orphaned commented-out call lines in the same neighborhoods
only if trivially safe (e.g. `interp2D.f90:1618`).

**Verify:** `make fesmutils-static` builds; `make test-coords` compiles.

---

## Step 2 — coupler precision  *(done — comment only)*

The original survey flagged `coupler.f90`'s `real(sp)` args as a `wp`-genericity
gap. On inspection this was a **misdiagnosis**: the public `remap` generic already
overloads both concrete kinds plus int (`remap_{2D,3D}_{dp,sp,int}`), so it
accepts `dp` and `sp` caller arrays *independent of* the compile-time `wp` — more
general than "follow `wp`", not less.

Changing the `_sp` wrappers to `real(wp)` would be a **bug**: under a future
`wp = dp` build, `remap_*_sp` becomes `real(dp)` and collides with `remap_*_dp`
(ambiguous generic → compile error), and drops `sp` input support.

Action taken: no code change; added a comment on the `remap` interface
documenting why the `dp`/`sp` overloads are intentional and must not be
genericized to `wp`. The real precision-switch obstacles live in the parked
`real(8)`/`real(kind=8)` cleanup.

---

## Step 3 — Error handling with rich context  *(varslice done; pilot)*

**Decisions (from review):** keep aborting behaviour but with rich context; use
`error stop` (nonzero exit); leave `ncio` self-contained; **no new shared error
module yet** — build a self-contained, private helper inside `varslice` first and
evaluate before generalizing; `-fbacktrace` is left to downstream build flags.

**Done in `varslice.f90`:**
- Added private, dependency-free `varslice_error(proc, msg, detail)` +
  a `to_str` generic (int/sp/dp/logical + 1-D arrays).
- Converted **all ~22 bare `stop` error sites** to `varslice_error`, routing
  every error to `error_unit` with a uniform frame + `error stop 1`, preserving
  (and formatting) all the per-site context values.
- Remaining `write(*,*)` in varslice are informational/verbose/loading output
  (not error paths) — intentionally left on stdout.
- Verified: `make fesmutils-static` + `make test-coords` build clean;
  `test_varslice_map` PASS.

**Verdict on the pattern:** `proc + msg + free-form detail` (with `to_str`) fits
the varied per-site dumps far better than a fixed context struct. `varslice_error`
+ `to_str` are ~120 self-contained lines; if we want this in other modules,
promoting them to a shared module later is mechanical.

**Rolled out to the other heavy `stop` users (module-local, per decision):**
- `coords/ncio_interp.f90` — 7 sites → `ncio_interp_error` + int `to_str`.
- `gaussian_quadrature.f90` — 9 sites → `gaussian_quadrature_error`. Vendored
  CISM-style code: kept its existing `print*` numeric diagnostics (indices,
  Jacobian/basis values) and replaced only the terminal message + `stop` with a
  framed `error stop`.
- `coords/coordinates.f90` — 8 sites → `coordinates_error` + int `to_str`
  (also added the missing `error_unit` import; errors previously went to stdout).

Each module carries its own private copy of the helper (deliberate duplication to
avoid interdependency). All error `stop`s across these four modules are now
`error stop` with framed, `error_unit`-routed, context-rich messages. Build clean;
`test-coords` + `test_varslice_map`/`mapinit`/`bilinear`/`conservative`/`ccsm3_grl`
all PASS.

**Note on duplication:** the `*_error` + `to_str` helpers are now near-identical
across four modules (~30–120 lines each). If a fifth module needs it, that's the
signal to reconsider a shared `error_handling` module vs. continued copies.

**Remaining `stop` users (not yet done, lighter):** `mapping_scrip`, `grid_cdo`,
`planet`, `nml`, `timeout`, `timestepping`, `hyster`.

### Original design notes (kept for reference)

Goal: library errors should carry **as much context as possible** to help find
bugs, and stop halting the host process by default.

Current state (from the survey): ~105 bare `stop`, 0 `error stop`, no
status-return convention (except `ncio`'s optional `stat`/`iostat`); diagnostics
go to stdout 3:1 over stderr, inconsistent even within a file.

Proposed direction (to finalize interactively before coding):

- A small shared `error_handling` module providing one reporting entry point,
  modeled on `ncio`'s `nc_ctx` idea but package-wide. It captures:
  module/routine name, a message, and optional structured context (file, var,
  dim, index/bounds, offending value(s), expected vs got). Writes a formatted,
  multi-line diagnostic to `error_unit`.
- An optional `stat` out-argument convention so callers can propagate failures
  instead of aborting; default behaviour (no `stat`) aborts via `error stop`
  with a nonzero code so exit status is meaningful.
- Migrate modules in batches (start with the heaviest / most user-facing:
  `varslice`, `coords/ncio_interp`, `coordinates`), replacing bare `stop` +
  `write(*,*)` with the shared reporter.

Open design questions to settle at Step 3:
1. Exact context struct fields and the formatted message layout.
2. `error stop` vs returnable `stat` default per call site.
3. Whether to fold `ncio`'s existing `nc_ctx`/`nc_check` into the shared module
   or keep ncio self-contained and mirror the convention.

---

## Parked (later steps, from the full analysis)

- fypp / include-file templating of `ncio.f90` (~3500 lines of hand-rolled
  type×rank variants) and `mapping.f90` precision wrappers (~300 lines).
- `make check` target that *runs* the test `.x` binaries and asserts exit codes,
  plus GitHub Actions CI (gfortran + netCDF).
- `fpm.toml` (and/or CMake) for general installability.
- Slim vendored upstream trees (submodules / fetch-on-build) — 4235 tracked files.
- Split large modules via F2008 submodules; unify the two remap apply paths onto
  `weight_map` as the sole engine.
- Cosmetic `real(8)`/`real(kind=8)` → `real(dp)`/`real(sp)` in
  `gaussian_quadrature.f90`, `gaussian_filter.f90`.
- Document runtime shell/`cdo` dependencies (`execute_command_line`).
