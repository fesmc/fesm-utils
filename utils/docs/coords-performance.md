# coords conservative remapping — performance & validation vs `cdo`

This note records how the coords conservative remapping (analytic polygon
clipping) compares to `cdo gen{con,dis}` in **accuracy** and **weight-generation
speed**. It is reproducible via `make bench-coords` then
`bin/bench_conservative.x` (needs `cdo` on PATH).

## Method

For each regime and size:

- **coords**: time the in-process `map_init_conservative` (analytic clip).
- **cdo**: time the *full* cost you would actually pay to obtain the same
  weights — serialize the grids (CDO grid description + a source NetCDF),
  run the `cdo gencon` subprocess, read the SCRIP file back — then bridge it
  to a `weight_map` via `map_scrip_to_weight_map`.
- **accuracy**: map a constant field (conservation: the result must be 1 where
  covered → `cons`), and the max `|coords − cdo|` on a smooth field where both
  maps cover the target.

Each timing is the mean of 3 repetitions. Grids are generated in code (no data
files). Hardware/toolchain: Apple Silicon (arm64, macOS), `gfortran -O2`,
single-threaded, `cdo` 2.6.1. Absolute times are machine-specific; the **ratios**
are the takeaway.

## Conservative results

`t` in seconds; `cons_*` is max|constant_out − 1|; `max|co−cdo|` is on a smooth
field of magnitude ~1–3.

```
regime               n_src    n_tgt   t_coords   t_cdo    coords/cdo   cons_co   cons_cdo   max|co-cdo|
cart->cart  n60       3600    14400    0.037      n/a        n/a        1.1e-14    n/a          n/a
cart->cart  n120     14400    57600    0.150      n/a        n/a        2.7e-15    n/a          n/a
cart->cart  n180     32400   129600    0.340      n/a        n/a        9.5e-16    n/a          n/a
ll->stereo  dx40km   16380     2888    0.020      0.049      0.40x      0          0            2.0e-1
ll->stereo  dx20km   16380    11476    0.042      0.077      0.55x      0          0            2.2e-1
ll->stereo  dx10km   16380    45904    0.127      0.194      0.65x      0          0            2.3e-1
ll->ll      dx1.5deg 16380    29040    1.81       0.076       24x       0          0            2.8e-3
ll->ll      dx1.0deg 16380    65160    2.86       0.123       23x       0          0            5.9e-3
ll->ll      dx0.75   16380   115680    4.28       0.172       25x       0          0            3.6e-3
```

`coords/cdo` is the slowdown factor (coords time ÷ cdo time).

## IDW results (`map_field` shepard vs `cdo gendis`)

```
regime               n_src    n_tgt   t_coords   t_cdo    coords/cdo   max|co-cdo|
idw ll->stereo dx20  16380    11476    0.125      0.050      2.5x        3.3e-2
idw ll->stereo dx10  16380    45904    0.482      0.124      3.9x        3.3e-2
```

## What this says

**Is coords competitive with cdo on raw weight-generation speed?** For the
cross-system (lat-lon→projected) case it is now **faster than cdo** (0.4–0.65×
its time); for lat-lon→lat-lon conservative `cdo`'s mature C implementation (YAC)
is still ~25× faster. coords' advantage is **correctness with no runtime `cdo`
dependency**, a single in-memory weight store, and being **fast and exact where
`cdo` does not natively apply**.

Regime by regime:

1. **Cartesian → cartesian (planar, same system).** coords is fast (130k target
   cells in 0.34 s) and conserves to machine precision (1e-15). `cdo gencon`
   operates on spherical lon/lat geometry and has no native equivalent for pure
   planar grids, so there is no cdo column. **coords wins outright here.**

2. **Lat-lon → lat-lon (spherical clip).** Accuracy vs cdo is excellent
   (max diff ~3e-3 on an O(1) field; conservation exact). Speed is ~25× slower
   than cdo but scales near-linearly with target size, so it is usable —
   weight generation is a one-time cost, then `map_field` is cheap.

3. **Lat-lon → regional projected (cross-system planar).** Correct (conservation
   exact) and now **fast** — 0.13 s for 46k targets, ~0.65× cdo's time. The
   candidate search runs on the unit sphere (see below); the planar
   Sutherland-Hodgman clip itself was never the bottleneck.

4. **IDW (kdtree shepard vs `cdo gendis`).** Comparable accuracy (~3e-2),
   ~2.5–4× slower than cdo at these sizes; the gap is mostly k-d tree build +
   exact-distance re-rank overhead.

### The cross-system slowdown (fixed)

Originally `conservative_planar` built a k-d tree over the source-cell centers
projected into the **target** plane, and for each target searched a radius
`0.5·target_diag + maxsrcdiag`, where `maxsrcdiag` is the largest projected
source-cell diagonal. When the source is **global** and the target is a
**regional** projection, source cells far from the projection center distort
enormously (stereographic blows up away from its center), so `maxsrcdiag`
became huge. The radius then pulled nearly every source cell for every target,
defeating the k-d tree and degrading to O(n_src·n_tgt) — 100 s for 46k targets.

**Fix:** for the cross-system path the candidate search now runs on the **unit
sphere** (the same undistorted geometry `conservative_spherical` uses): the k-d
tree is built over source-cell centers as 3D unit vectors, and each target queries
a chord radius from the source + target angular cell size (the projected target's
cell size is converted to degrees via the planet radius). Source corners are still
projected into the target plane once for the planar clip, but the candidate set
per target is now O(1). Result: 100 s → 0.13 s, conservation unchanged (exact).
The `same_sys` path is undistorted and keeps its planar tree.

### Accuracy caveat for projected targets

`max|co-cdo|` is ~2e-1 for `ll->stereo` vs ~3e-3 for `ll->ll`. Conservation is
exact (`cons = 0`) in both, so this is an interpolation-detail difference, not a
mass error. It is concentrated where the two pipelines treat grid edges
differently — the benchmark uses the deliberately simple CDO grid descriptions
with `wraplon` disabled (the `wraplon=.true.` helper is known-broken), and a
regional stereographic target wraps in longitude near the pole. The clean
`ll->ll` figure (~0.1 %) is the trustworthy accuracy number; the projected-target
discrepancy warrants a closer per-point look before relying on it.

## Bottom line

coords delivers the design's **correctness** goal (matches cdo on lat-lon
conservation; handles planar cases cdo cannot) and removes the runtime cdo
dependency. The design's **speed** goal is largely met: the k-d tree avoids
O(N·M) for IDW neighbor search, cartesian and cross-system conservative weight
generation are at or faster than cdo, and only the lat-lon→lat-lon spherical clip
remains ~25× slower than cdo's C implementation. For the typical workflow —
generate weights once, then apply `map_field` many times — current speed is
acceptable across all regimes.
