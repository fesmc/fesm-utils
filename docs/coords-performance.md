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
regime               n_src    n_tgt   t_coords   t_cdo    speedup    cons_co   cons_cdo   max|co-cdo|
cart->cart  n60       3600    14400    0.00039    n/a        n/a      1.1e-14    n/a          n/a
cart->cart  n120     14400    57600    0.00080    n/a        n/a      2.7e-15    n/a          n/a
cart->cart  n180     32400   129600    0.0016     n/a        n/a      9.5e-16    n/a          n/a
ll->stereo  dx40km   16380     2888    0.018      0.060      3.3x     0          0            2.0e-1
ll->stereo  dx20km   16380    11476    0.038      0.072      1.9x     0          0            2.2e-1
ll->stereo  dx10km   16380    45904    0.117      0.180      1.5x     0          0            2.3e-1
ll->ll      dx1.5deg 16380    29040    0.00074    0.076      102x     0          0            1.0e-4
ll->ll      dx1.0deg 16380    65160    0.0012     0.114       94x     0          0            6.8e-14
ll->ll      dx0.75   16380   115680    0.0019     0.167       88x     0          0            1.3e-13
```

`speedup` is cdo time ÷ coords time (>1 means coords is faster). The cart->cart
and ll->ll rows now use the **analytic separable** path (`conservative_planar_reg`
/ `conservative_latlon`): axis-aligned cells make the overlap the product of two
1-D interval overlaps, so there is no polygon clip and no k-d tree. ll->stereo is
the cross-system planar clip (lat-lon source projected into the target plane),
which is not separable and keeps the Sutherland-Hodgman clip.

The ll->ll change is the headline: from ~24x slower than cdo to ~90-100x
**faster**, and agreement with cdo tightened from ~3e-3 to 1e-4..1e-13 (the
analytic overlap computes the same lon x sin(lat) rectangle area cdo's YAC does).

## IDW results (`map_field` shepard vs `cdo gendis`)

```
regime               n_src    n_tgt   t_coords   t_cdo    speedup    max|co-cdo|
idw ll->stereo dx20  16380    11476    0.029      0.061      2.1x       3.2e-2
idw ll->stereo dx10  16380    45904    0.092      0.115      1.3x       3.3e-2
```

Shepard still uses the k-d tree (it is a genuine neighbourhood kernel, not
separable); it is now roughly on par with cdo at these sizes.

## Structured nn / bilinear results (`map_init` grid-locate vs `cdo gennn`/`genbil`)

nn and bilinear grid->grid maps with a regular source now locate each target
directly in the source grid (bisection on the axis centres) instead of a k-d
tree — see `map_init_structured`. Bilinear bakes the exact index-space weights
(matching `cdo genbil`); nn takes the nearest cell corner. Both directions of the
lat-lon <-> stereographic pair are shown (the Rtopo <-> NH-16KM shapes).

```
regime               n_src    n_tgt   t_coords   t_cdo    speedup    max|co-cdo|
bil ll->stereo dx20  16380    11476    0.011      0.046      4.1x       2.0e-3
bil stereo->ll dx20  11476    16380    0.0068     0.053      7.8x       1.0e-4
nn  ll->stereo dx20  16380    11476    0.028      0.046      1.6x       2.7e-2 (NN tie-break)
nn  stereo->ll dx20  11476    16380    0.0086     0.085      9.9x       0.0
```

Bilinear matches cdo genbil to 1e-4..2e-3 and runs 4-8x faster. nn matches cdo
exactly where the nearest source is unambiguous; the 2.7e-2 in the ll->stereo row
is boundary tie-breaking (two sources near-equidistant), inherent to nearest-
neighbour and also present between any two independent nn implementations.

## What this says

**Is coords competitive with cdo on raw weight-generation speed?** Yes — for the
grid pairs the coupler actually uses (axis-aligned regular grids) coords is now at
or faster than cdo across every conservative and interpolation regime, having
replaced the general polygon-clip / k-d-tree work with the structured shortcuts
cdo's YAC uses internally. coords' additional advantage is **correctness with no
runtime `cdo` dependency**, a single in-memory weight store, and handling cases
`cdo` does not natively apply (pure planar grids).

Regime by regime:

1. **Cartesian → cartesian (planar, same system).** Analytic separable rectangle
   overlap: 130k target cells in ~0.0016 s, conserving to machine precision
   (1e-15). `cdo gencon` operates on spherical lon/lat geometry and has no native
   equivalent for pure planar grids, so there is no cdo column. **coords wins
   outright here.**

2. **Lat-lon → lat-lon.** Analytic separable lon × sin(lat) overlap: ~90–100×
   **faster** than cdo, agreement 1e-4..1e-13 (it computes the same rectangle area
   YAC does), conservation exact. This replaced the great-circle polygon clip that
   was ~25× *slower* than cdo.

3. **Lat-lon → regional projected (cross-system planar).** Correct (conservation
   exact) and fast — ~0.12 s for 46k targets, ~1.5× cdo's time. Not separable (the
   projected source cells are general quadrilaterals), so it keeps the planar
   Sutherland-Hodgman clip with a unit-sphere candidate search (see below).

4. **nn / bilinear (structured grid-locate vs `cdo gennn`/`genbil`).** Regular
   source grids locate the target directly (no k-d tree): bilinear is exact
   index-space (matches genbil to 1e-4, 4–8× faster), nn is the nearest cell corner
   (matches gennn exactly modulo boundary ties, up to ~10× faster). See
   `map_init_structured`.

5. **IDW / shepard (kdtree vs `cdo gendis`).** A genuine neighbourhood kernel, so
   it keeps the k-d tree; now roughly on par with cdo (~1.3–2×).

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

### The spherical-clip slowdown (fixed)

The lat-lon→lat-lon conservative case originally ran the general great-circle
polygon clip (`conservative_spherical`) for every source/target cell pair, at
~25× cdo. But lat-lon cells are axis-aligned in (lon,lat), so their overlap is
the separable product of a longitude overlap and a sin(lat)-band overlap — no
clipping needed. `conservative_latlon` computes this analytically with a
separable 1-D candidate search (no k-d tree), which is what makes it ~90–100×
faster than cdo. `conservative_planar_reg` does the same for same-Cartesian-plane
grids. A related units bug in the `conservative_spherical` candidate radius
(a projected source's cell size was used as degrees) was also fixed; that path
now only runs for the genuinely non-separable mixes (e.g. projected source →
lat-lon target).

## Bottom line

coords delivers the design's **correctness** goal (matches cdo on lat-lon
conservation and bilinear; handles planar cases cdo cannot) and removes the
runtime cdo dependency. The **speed** goal is now fully met for the regular-grid
pairs the coupler uses: analytic separable conservative weights and structured
nn/bilinear locate replace the polygon clip and k-d tree with the same shortcuts
cdo uses, so coords is at or faster than cdo everywhere except where a genuine
neighbourhood search (shepard) or a non-separable cross-system clip is required,
and even those are now on par. For the typical workflow — generate weights once,
then apply `map_field` many times — weight generation is no longer a bottleneck.
