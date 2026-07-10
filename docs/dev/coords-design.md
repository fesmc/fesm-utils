# `coords` — design & implementation plan

A clean rebuild of the `coordinates` library as the `coords` module set inside
`fesm-utils/utils`. It keeps the user-facing API of the old library
(`points_class` / `grid_class` / `map_class`, `*_init`, `map_field`) while
replacing the internal storage, neighbor search, and conservative remapping with
a single clean design.

This document is the agreed plan. It was produced by stress-testing the design
decision-by-decision; the numbered decisions below are settled.

## Goals

- **One storage way.** The old library had three conservative implementations,
  two `pt_wts_class` definitions, and two on-disk formats encoding the same
  sparse source→target weight matrix. `coords` has a single in-memory weight
  type, a single apply path, and a single file format.
- **Fast.** Replace the brute-force `O(N·M)` neighbor search with a KD-tree, and
  the subgrid-sampling conservative weights with analytic polygon clipping —
  removing the runtime `cdo` dependency once validated.
- **Clean objects from the start.** Coordinate-system metadata is defined once
  and shared, not duplicated across three types.
- **Call-compatible.** Downstream code keeps `call grid_init(...)` /
  `call map_field(...)`; the redesign is internal. (Direct struct-field reads
  like `grid%mtype` move under `grid%cs%mtype` — a one-time downstream sweep.)

## Location & build

- Lives in `utils/src/coords/` (the existing `coordinates-light/` folder renamed
  and grown). Builds into `libfesmutils.a` via the same Makefile strategy as the
  rest of `utils`.
- Reuses utils modules directly (no copies): `precision`, `ncio`, `interp1D`,
  `interp2D`, `index`, `gaussian_filter`, `grid_to_cdo`.
- Umbrella module `coords` re-exports the public surface; users `use coords`.

## Object model

```
coord_system   ! mtype, units, planet, proj, xy_conv,
               !   is_cartesian, is_projection, is_lon180   — defined ONCE,
               !   embedded as %cs in points/grid/map
grid_axis      ! nx, ny, x0, dx, y0, dy, x(:), y(:)
               !   supports explicit non-uniform axes (Gaussian latitudes)
points         ! npts; x/y/lon/lat/area/dx/dy(:); cs
grid           ! 2D x/y/lon/lat/area(:,:); grid_axis; cs; cached points view
weight_map     ! the single sparse store (internal)
map_class      ! public wrapper: weight_map + target geometry + names + cs
```

### `weight_map` — the single store

```fortran
integer, parameter :: MAP_DISTANCE = 1   ! links carry raw distances + positions;
                                         !   weights computed at apply time
integer, parameter :: MAP_WEIGHT   = 2   ! links carry final weights; apply reduces

type weight_map
   integer :: kind                       ! MAP_DISTANCE | MAP_WEIGHT
   integer :: n_src, n_dst, n_links
   integer,  allocatable :: src(:), dst(:)     ! SCRIP-style link addresses
   real(dp), allocatable :: w(:)               ! final weight (MAP_WEIGHT)
   ! present only when kind == MAP_DISTANCE:
   real(dp), allocatable :: x(:), y(:)         ! source-neighbor positions
   real(dp), allocatable :: dist(:)
   integer,  allocatable :: quadrant(:), border(:)
end type
```

- `MAP_DISTANCE` holds distances **and neighbor positions**, so every
  distance-based method — `nn`, `quadrant`, `shepard`, and **bilinear** — is
  derived uniformly at apply time. No method-specific α/iquad fields.
- `MAP_WEIGHT` holds final weights (conservative overlap areas, or an imported
  SCRIP matrix); apply just reduces them (`mean`/`count`/`stdev`).
- **All geometry and all weights/distances are `dp`.** Internal accumulation is
  always `dp`; the public `map_field` stays generic over the field element type
  (`integer`/`sp`/`dp`). This removes the old `sp`/`dp` mixing and the
  conservative-sum drift.

### Apply semantics (dynamic masks)

Source masks change between calls, so `MAP_DISTANCE` weights are computed
**per call** from `(positions, current mask, method)` — never cached across
calls. Bilinear under a masked corner:

- **default:** graceful degradation to inverse-distance over the surviving
  quadrant neighbors (true bilinear when all 4 present; missing only when 0
  neighbors valid). The returned mask reports which cells fell back.
- **strict mode** (per-call argument): missing unless all 4 corners valid.

## File format

One NetCDF layout, a **SCRIP superset**:

- Standard SCRIP core (`src_address`, `dst_address`, `remap_matrix`) — fully
  interoperable with `cdo` and other tools.
- Optional `coords` extension variables (`mp_x`, `mp_y`, `mp_dist`,
  `mp_quadrant`, `mp_border`) written for `MAP_DISTANCE` maps, so a `coords`
  reload restores the exact distance map and keeps change-method-for-free.

Writer: `call map_write(wm, filename, scrip_only=.false., method=...)`

- `MAP_WEIGHT` → always a pure, standard SCRIP file (`scrip_only`/`method`
  ignored).
- `MAP_DISTANCE`, `scrip_only=.false.` (default) → SCRIP core (baked from
  `method`, default = the map's own method) **plus** extensions; the file is
  still valid SCRIP for other tools.
- `MAP_DISTANCE`, `scrip_only=.true.` → standard SCRIP core only, single method,
  `num_wgts = 1`.

Reading a plain CDO SCRIP file yields a `MAP_WEIGHT`. `mapping_scrip` is kept as
the SCRIP reader feeding `weight_map`; a new in-house SCRIP **writer** is added
(prerequisite for dropping `cdo`).

## Conservative remapping (staged)

All combinations are required: `cart→cart`, `latlon→latlon`, `cart↔latlon`.

- **Stage A — import CDO.** `map_init(method="con")` runs `cdo gencon` and
  imports the result into a `MAP_WEIGHT`. Works for every case on day one.
- **Stage B — planar analytic clip.** Sutherland–Hodgman in the target plane;
  covers any Cartesian/projected target (incl. `latlon`-source → projected).
  Area via shoelace. Drops `cdo` for projected targets.
- **Stage C — spherical clip.** Great-circle edge / point-side tests, area via
  geodesic `planet_area`; unlocks `latlon`/Gaussian targets and full `cdo`
  independence. `cdo` covers `latlon` targets until C lands.

## Projections & grids

- Port `projection_oblimap2` verbatim (stereographic, polar-stereo, LAEA) and
  add **rotated-pole** to its dispatch table.
- **Gaussian grids are a `latlon` type** with a non-uniform latitude axis
  (geodesic distances, `is_cartesian=.false.`). Cell corners/areas are derived
  from actual adjacent-axis midpoints, so uniform and Gaussian grids share one
  area path (correct for conservative). Add a **Gaussian-latitude generator**
  (Legendre roots + weights) so a Gaussian grid can be defined from a resolution.

## Port vs. implement

| Bucket | Modules | Action |
|---|---|---|
| Reuse from utils | precision, ncio, interp1D, interp2D, index, gaussian_filter, grid_to_cdo | `use`, no copy |
| Port verbatim | geodesic, projection_oblimap2 | copy, rename header only |
| Port + light refactor | planet (data-driven, dp), polygons (extend for clipping), grid_gen, loess, interp_time | |
| Redesign (keep call API) | coordinates (geometry), coordinates_mapping (`map_init`/`map_field`) | new store + KD-tree |
| Collapse to one | 3 conservative paths → one conservative module | unified `MAP_WEIGHT` |
| Adapt | mapping_scrip | reader → `weight_map`; add SCRIP writer |
| New | kdtree, polygon clip (planar+spherical), unified apply+I/O, rotated-pole, Gaussian-lat generator | implement |
| Drop | mod_toms526, coordinates_sigma, subset, subset2 | not ported |

## Milestones

Each is gated on its validation before the next; validation is self-contained
(no dependency on the old `coordinates` repo).

1. **Scaffold + geometry.** `coords/` folder + Makefile wiring; `coord_system` /
   `points` / `grid` / `grid_axis`; port geodesic / oblimap (+rotated-pole) /
   planet / polygons; Gaussian-latitude generator. Grids build, project, write.
   *Gate: ported `test_proj` / `test_geodinverse`.*
2. **Unified IDW mapping.** `weight_map` + KD-tree + `map_init` / `map_field`
   (nn/bilinear/quadrant/shepard) with per-call mask semantics; SCRIP-superset
   reader/writer. *Gate: in-harness brute-force neighbor-set equality;
   CCSM3→GRL regression.*
3. **Conservative A→B.** CDO import into `MAP_WEIGHT`; then planar analytic clip.
   *Gate: area conservation + `cdo gencon` diff (Cartesian targets).*
4. **Conservative C.** Spherical clip → `latlon`/Gaussian-target conservative;
   retire `cdo`. *Gate: area conservation + `cdo` diff (latlon/Gaussian targets).*
5. **Migrate `varslice`** to the unified `map_field`; remove the legacy SCRIP
   applier. (Separate, isolated change.)

## Status (2026-06)

All milestones are implemented and gate-tested (`make test-coords`, 16 programs):

1. Scaffold + geometry — **done**.
2. Unified IDW mapping (`weight_map`, k-d tree, `map_init`/`map_field`,
   SCRIP-superset I/O) — **done**.
3. Conservative A→B and 4. Conservative C — **done** (planar same-system +
   cross-system; spherical great-circle clip via NSUB=8 subsampling).
5. `varslice` migrated to `map_field`; the legacy SCRIP applier (`map_scrip_field`)
   is retained only as a test reference. The `nc_read_interp` family moved to a new
   `ncio_interp` module (also on `map_field`).

Stage-A bridge: `map_scrip_to_weight_map` imports a plain CDO SCRIP map into the
unified store; `map_read` loads one straight into a `map_class`.

**Performance vs `cdo`** is measured in [coords-performance.md](coords-performance.md)
(`make bench-coords`). Summary: correctness matches `cdo` and planar/cartesian
remapping (which `cdo` cannot do natively) is fast and exact; but conservative
weight generation is not yet speed-competitive with `cdo`'s C implementation
(~25× slower for lat-lon→lat-lon), and the cross-system projected path has a known,
fixable O(n_src·n_tgt) search inefficiency. The runtime `cdo` dependency is removed.

## Validation

Self-contained only:

- **KD-tree IDW:** a ~50-LOC brute-force neighbor search lives in the test
  program (not the library); on small grids the KD-tree must return identical
  neighbor sets and bit-equal distances.
- **Conservative:** area-conservation invariant (Σ weight·area matches the source
  integral) plus a diff against `cdo gencon` on the test grids.
- **Projections/geodesics:** ported `test_proj` / `test_geodinverse`,
  self-checking against known values.
- **Regression grids:** CCSM3-T42, GRL-20KM, MAR, ETOPO.
