# coords — remaining work & the precision/constants proposal

Handoff notes for a future session. The coords library itself is functionally
complete and gate-tested (`make test-coords`, 16 programs), the conservative
cross-system performance fix has landed, and performance vs `cdo` is documented
in [coords-performance.md](coords-performance.md). What follows is **not yet done**.

The headline item is the **package-wide precision/constants consolidation** (§1);
§2 lists smaller interpolation-parity gaps; §3 the held `dev` merge and a known
accuracy caveat.

---

## 1. Precision / constants consolidation (the main proposal)

### Problem — the current situation in `fesm-utils/utils`

Kind (`dp`/`sp`/`wp`) and constant (`mv`/`pi`/…) definitions are fragmented:

- `precision.f90` mixes **kinds** (`dp, sp, wp`) with **values** (`mv`, `TOL`,
  `TOL_UNDERFLOW`). "precision" is the wrong home for `mv`/`TOL`.
- `coords/coord_constants.f90` (coords-only) holds the domain values
  (`MISSING_VALUE_DEFAULT`/`mv`, `pi`, `degrees_to_radians`, `radians_to_degrees`,
  `ERR_DIST`, `ERR_IND`) and correctly sources `dp`/`sp` from `precision`.
- **~12 modules redefine their own `dp`/`sp`/`wp` locally** instead of
  `use precision`: `root_finder`, `gaussian_filter`, `timestepping`, `timeout`,
  `esm`, `hyster`, and the coords legacy set `interp1D`, `interp2D`, `grid_to_cdo`,
  `mapping_scrip` (and historically `esm_base`). Only ~6 modules actually
  `use precision` (`distances`, `gaussian_quadrature`, `derivatives`, `subgrid`,
  `varslice`, `staggering`).
- **Two `mv`s with different kinds**: `precision%mv = -9999.0_sp` vs
  `coord_constants%mv = -9999.0_dp`. Numerically equal, but two sources of truth.

### Proposed target shape

- **`precision`** → kinds **only** (`dp, sp, wp`). Nothing else.
- **`constants`** (package-wide) → all values: one `mv`, `pi`,
  `degrees_to_radians`/`radians_to_degrees`, `ERR_DIST`, `ERR_IND`, `TOL`,
  `TOL_UNDERFLOW`. Generalize/rename `coord_constants` and drop the `coord_`
  scoping (these constants are universal, not coords-specific); optionally keep a
  thin `coord_constants` that just re-exports for back-compat.
- Every module: `use precision` (kinds) + `use constants` (values); delete the
  ~12 local copies.

### Critical safety check FIRST (read-only, no edits)

Switching a module from a local `wp` to `precision%wp` is only safe if its local
`wp` already equals `precision%wp` (currently `= sp`). If any module set
`wp = dp` locally and relies on it, the swap **silently changes its working
precision** — a real correctness risk, not a compile error.

Do this audit before touching anything:

```sh
cd fesm-utils/utils/src
# every local kind definition and its value
grep -rnE ':: *(dp|sp|wp) *= *kind|parameter *:: *(dp|sp|wp)' . | sort
# every local constant definition (mv/pi/missing) and its value/kind
grep -rnE 'parameter *:: *(mv|MV|pi|MISSING_VALUE|TOL)' .
```

Tabulate each module's local `wp` value. Proceed only once all are confirmed
`= sp` (or handle the exceptions explicitly).

### Blast radius & sequencing

- Touches **~every module in `utils`** (including `conservative.f90`). The
  conservative fix is now merged, so that conflict is clear.
- Touches **production modules consumed directly by climber-x / yelmo**
  (`interp1D`, `interp2D`, `gaussian_filter`, …). It is API-compatible **only if**
  the kind values are identical everywhere (they should be: `dp = kind(1.d0)`,
  `sp = kind(1.0)`, `wp = sp`). Do it on a branch, build, and ideally smoke-test a
  downstream build before merging.
- This is a **fesm-utils-wide cleanup, broader than the coords project** — treat
  it as its own reviewed change, not folded into coords.

### Suggested steps

1. Read-only `wp`/constant audit (above); confirm consistency.
2. Add/extend `constants` as the single source of `mv`/`pi`/`TOL`/…; standardize
   on `dp` internally and expose `sp` forms only where public APIs require them
   (resolves the two-`mv` issue).
3. Strip `precision` to kinds only.
4. Per module (one or a few per commit): replace local kind defs with
   `use precision`, local constants with `use constants`; delete the duplicates.
5. After each batch: `python config.py /tmp/coords_frag.mk && make fesmutils-static`
   then `make test-coords` (16 green). (Worktree build is not configme-driven; see
   the memory note / coords-design.md for the fragment trick.)

---

## 2. Interpolation parity gaps (close only where needed)

The grid→grid path is at full parity and richer than the old `coordinates`
library (adds conservative count/stdev, weighted/poisson fill-smoothing). Gaps:

- **fill / smoothing / `mask_pack` / `reset` live only on `map_field_grid_grid`.**
  `points→grid` has a 2-D target and should support them but doesn't;
  `grid→points` / `points→points` (1-D targets) could still take `mask_pack`/`reset`.
- **sp/integer generics exist only for `grid→grid`**; the other three combos are
  dp-only (old lib had integer/float/double for all four).
- **`border` option not ported.** Old `border` relaxed the neighbor criterion so
  domain-edge cells still fill. coords' default already uses whatever valid
  neighbors exist (graceful degradation), so this is likely subsumed — confirm on a
  real case before adding.
- **Dropped by design** (not interpolation per se): `subset`/`subset2` region
  extraction, `coordinates_sigma`.

In-tree consumers (`varslice`, `ncio_interp`) are all grid→grid, so they have
everything. The gaps only bite if something maps onto a `points→grid` target and
wants fill/smoothing.

---

## 3. Other remaining items

- **Merge `coords-dev` → `dev`** when ready (currently held by request).
  `origin/coords-dev` is stale; local `coords-dev` is the source of truth.
- **ll→stereo accuracy caveat**: `max|coords−cdo|` ~2e-1 for projected targets vs
  ~3e-3 for lat-lon→lat-lon. Conservation is exact in both, so it's an
  edge/longitude-wrap interpolation-detail difference (the benchmark uses the
  deliberately simple grid descriptions with the known-broken `wraplon` disabled).
  Worth a per-point look before relying on projected-target accuracy.
