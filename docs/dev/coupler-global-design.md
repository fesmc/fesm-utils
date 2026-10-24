# coupler: global <-> regional exchange — design note

Status: implemented through step 7 (branch phase0-coupler). Follow-up on
phase0-kryos: the two geometry problems of `conservative_spherical` found in
step 7 (seam/pole source cells, concave target cells) are fixed, and `Σ_k φ_jk
= 1` holds for a covering target (§11).

Scope: make conservative exchange between N regional projected grids and one
global grid (lon-lat or Gauss-Legendre, e.g. VILMA2 T64) correct and
first-class in `coords` + `coupler`, and close the small coupler gaps. Source of
the requirements: kryos `docs/design.md` §5, §6, §12 row 3, decision 4;
kryos `docs/study/fesm-utils.md` §2.3–2.7, §6.3, §8.6.

**Hard constraint.** Every existing call keeps its default behaviour
bit-identical (downstream refactors are gated on it). Every new behaviour is
behind a new optional argument or a new procedure. Where a change cannot be
bit-identical, it is listed explicitly below.

Notation: map links `(j,k)`, source `j`, target `k`, baked weight `w_jk`
(overlap area in the *map metric*, §1.1). `V_k` = links of `k` with `w_jk > 0`
and `v_j /= missing_value`.


**Review outcome (2026-10-24, Alex).** Q1: the default model area of a
projected grid is the plane `dx·dy` (what Yelmo and FastIsostasy integrate),
so merges and `bsl` conserve the volume the ice model conserves; the
ellipsoidal (`grid%area`) integral is logged beside it as a diagnostic. Q3:
prune fix accepted. Q2: no raw sum. Q4: accept unverified (logged). Q5:
coords-only. Q6: explicit bounds only. Q7: accumulator. Q8: 1e-2, stop. Q9:
add `method="auto"` to the coupler (extensive: `con` both ways; intensive:
`bilin` when refining, `con` when coarsening, by median model area). Q10: ok.
Q11: `dp`. Style: free procedures as already written (no type-bound
procedures anywhere in fesm-utils).

## 1. Normalisation and coverage

### 1.1 Area metrics in use today

Each conservative generator measures overlap in its own metric, and none of
them is `grid%area`:

| Path (`conservative.f90`) | When | Unit of `w_jk` |
|---|---|---|
| `conservative_planar_reg` | same Cartesian plane | axis units² in the projection plane (no scale factor) |
| `conservative_planar` | other source -> projected target | target-plane axis units² of NSUB-subsampled projected source cells |
| `conservative_latlon` | ll/Gauss -> ll/Gauss | `Δλ[deg]·Δsinφ`, spherical band, poles clamped, no `R²` |
| `conservative_spherical` | projected -> ll/Gauss | Karney area (target planet) of great-circle polygons, NSUB=8 samples per edge in each grid's plane; target cells as 8 convex strips (§11) |

`grid%area` (`grid_area`) is: plane area for Cartesian; the Karney area of the 4
inverse-projected corners for projected grids (in axis units², i.e. km² for km
grids); the Karney area of the 4 corners for lon-lat (geodesic edges, not a
latitude band; polar edges extrapolated, not clamped). For a polar
stereographic ice grid the plane metric differs from `grid%area` by the local
`k²` (several % across Antarctica). Yelmo integrates volume as `Σ H dx dy`
(plane metric).

Consequence: "conservative" today means conservative in the map metric only.
A consistent integral check must state the metric. The design below keeps the
map metric internal (for partitions and coverage) and does all mass accounting
in a *model area* chosen per grid (§2.1).

### 1.2 What the map must carry

`weight_map_t` gains four `real(dp), allocatable` vectors, filled by the
conservative generator (cheap there: the cell polygons are already built) and
written to the cache as coords extension variables (`mp_wfull_src`, ...):

```fortran
real(dp), allocatable :: wfull_src(:)   ! (n_src) full source-cell area, map metric
real(dp), allocatable :: wfull_dst(:)   ! (n_dst) full target-cell area, map metric
real(dp), allocatable :: wcov_src(:)    ! (n_src) sum_k w_jk before pruning
real(dp), allocatable :: wcov_dst(:)    ! (n_dst) sum_j w_jk before pruning
```

At load, the post-prune sums `s_j = Σ_k w_jk`, `s_k = Σ_j w_jk` are computed from
the CSR, and `q_j = wcov_src_j/s_j`, `r_k = wcov_dst_k/s_k` restore the coverage
removed by `weight_map_prune`. With them, coverage and mass are exact even on a
pruned map. A cache file without these variables (every existing file, and cdo
SCRIP files) loads as today; requesting `norm /= "frac"` or `frac2` on it stops
with "regenerate <file>". Not inferred from geometry at load, to avoid a second
copy of each generator's cell construction.

Prerequisite fix (found by inspection, not yet reproduced): `weight_map_prune`
drops links with `w < 1e-4·Σw`. A target with more than 10⁴ near-equal links
(a 1 km source onto a ~2° cell) loses **all** its links. Proposed criterion:
`w < cutoff·maxval(w over the target)`. This only affects newly generated maps;
existing caches are unchanged (see §10, Q3).

### 1.3 `norm` and `frac2`

```fortran
subroutine weight_map_apply(wm, var1, var2, kernel, stat, missing_value, radius, mask2, &
                            norm, frac2, area1, area2)
    character(len=*), optional, intent(in)  :: norm      ! "frac" (default) | "dst" | "sum"
    real(dp),         optional, intent(out) :: frac2(:)  ! (n_dst) covered fraction in [0,1]
    real(dp),         optional, intent(in)  :: area1(:)  ! (n_src) model area of source cells
    real(dp),         optional, intent(in)  :: area2(:)  ! (n_dst) model area of target cells
```

Definitions (`φ_jk = q_j·w_jk/wfull_src_j`, the fraction of source cell `j`
lying in target `k`; `Σ_k φ_jk` = covered fraction of `j`):

| Output | Value at target `k` |
|---|---|
| `frac2` | `min(1, r_k·Σ_{V_k} w_jk / wfull_dst_k)`; 0 where no valid link, outside `mask_pack`, or not mapped |
| `norm="frac"` | `Σ_{V_k} w v / Σ_{V_k} w` — today's code path, untouched |
| `norm="dst"`, no areas | `r_k·Σ_{V_k} w v / wfull_dst_k` = `frac2·mean` (map metric; preserves constants × coverage) |
| `norm="dst"`, `area1`+`area2` | `Σ_{V_k} φ_jk·area1_j·v_j / area2_k` (exact in the model metric) |
| `norm="sum"` | `Σ_{V_k} φ_jk·α_j·v_j`, `α = area1` if present else 1 (amount per target) |

Rules:
- Partially covered target: `frac` returns the covered mean (integral inflated
  by `1/frac2`, as today); `dst` returns integral / full target area;
  `frac2` says how much was covered.
- Masked source (value `== missing_value`): excluded from `V_k`. It contributes
  neither mass nor coverage. This is the same rule for every `norm`.
- `mask2 = (Σ_{V_k} w > 0)` for every `norm` (unchanged meaning); where
  `.not. mask2`, `var2` is `missing_value` (or kept with `reset=.false.`).
- Conservation (`dst` with areas): `Σ_k area2_k v2_k = Σ_j area1_j v_j Σ_k φ_jk`.
  For a global target `Σ_k φ_jk = 1`, so the regional integral is preserved to
  rounding, for any choice of `area1`, `area2`. A uniform field on a fully
  covered target returns `v·(1+ε_k)`, `ε_k` = local metric mismatch: that is
  the price of exact mass, and why intensive fields use `frac` or `dst` without
  areas.
- `area1` without `area2` (or the reverse) with `norm="dst"`, any
  `norm /= "frac"` with `stat /= "mean"`, and `norm /= "frac"` on a
  non-conservative map are errors (`fesm_error`).
- `frac2` on a non-conservative map (bilin, nn, distance kernels) is the
  indicator `merge(1,0,mask2)`.
- Pass-throughs: `mask_pack` restricts as today (`frac2 = 0` outside);
  `fill_method`/`filt_*` run after and do not touch `frac2` (filled values are
  not coverage, and are not conservative).

Bit-identity: `norm` absent or `"frac"` runs today's fast path unchanged.
`frac2` present with `norm="frac"` is accumulated in a separate register in the
same loop, so `var2` is unchanged bit for bit.

`map_field`: the 8 real specifics (grid/points × grid/points × dp/sp) gain
`norm`, `frac2(:,:)` or `frac2(:)` (target shape, always `real(dp)`),
`area1`, `area2` (source/target shape, `real(dp)`). The 4 `int` specifics do
not get them (categorical data has no amount).

## 2. Coupler: model areas, merge, split, integral

### 2.1 Model area per registered grid

```fortran
subroutine coupler_set_area(cpl, name, area)
    type(coupler_class), intent(inout) :: cpl
    character(len=*),    intent(in)    :: name
    real(dp),            intent(in)    :: area(:,:)   ! [m²]
```

Default: `grid%area·xy_conv²` (m²), set at registration. A host overrides it to
match what its model integrates: Yelmo `dx·dy` (plane), VILMA2 the quadrature
`R²Δλ·w_j` (automatic once Gauss bounds exist, §3). The coupler passes these as
`area1`/`area2` whenever it applies `norm="dst"`/`"sum"`, so mass accounting
and `coupler_integral` share one metric per grid.

### 2.2 `coupler_integral`

```fortran
function coupler_integral(cpl, name, var, mask, missing_value) result(total)
    type(coupler_class), intent(in) :: cpl
    character(len=*),    intent(in) :: name
    real(dp),            intent(in) :: var(:,:)
    logical,  optional,  intent(in) :: mask(:,:)
    real(dp), optional,  intent(in) :: missing_value   ! cells equal to it are skipped
    real(dp) :: total                                   ! Σ area·var, Neumaier-compensated
```

Consistent check: for a regional -> global exchange with `norm="dst"`,
`coupler_integral(cpl, glob, v2) == coupler_integral(cpl, reg, v1)` to
~1e-13 relative, because both sides use the registered areas and §1.3's
identity holds in that metric. Nothing is cached.

### 2.3 `coupler_merge`: N regional fields + background -> global

Accumulator form, so that N arrays of different shapes need no pointer
wrappers and every domain is added where it is computed:

```fortran
! coupler_merge_t: opaque; per-domain names + cov(:,:,d), amt(:,:,d) on the global grid
call coupler_merge_begin(acc, cpl, dst [, kind])     ! kind: "extensive" (default) | "intensive"
call coupler_merge_add(acc, cpl, src, var_src [, owner, missing_value])
call coupler_merge_end(acc, cpl, background, var_dst [, frac_tot, overlap_tol])
!   acc: intent(out) / (inout) / (inout); cpl: intent(inout) (builds/caches src->dst "con")
!   var_src(:,:), owner(:,:) [0,1] on the src grid (default 1), background(:,:): real(dp), in
!   var_dst: real(dp), allocatable, inout; frac_tot(:,:) = F_k unclipped, out; overlap_tol default 1e-2
```

Per domain `d` (source cells `j` with valid value, weighted by owner `o_dj`):

```
F_dk = r_k Σ_j w_jk o_dj / wfull_dst_k                                 coverage
A_dk = Σ_j φ_jk o_dj area_dj v_dj / area_k      (extensive)            amount / area
A_dk = r_k Σ_j w_jk o_dj v_dj / wfull_dst_k     (intensive)            frac·mean
F_k  = Σ_d F_dk
extensive:  var_k = (1 - min(F_k,1))·bg_k + Σ_d A_dk
intensive:  var_k = (1 - min(F_k,1))·bg_k + Σ_d A_dk / max(F_k,1)
```

- Conservation (extensive): `∫var = ∫(1-min(F,1))·bg + Σ_d ∫ o_d·v_d`, exact in
  the registered areas. Regional amounts are never rescaled; only the
  background is capped. Intensive merges preserve a constant field to rounding.
- Owner weights: `Σ_d o_d = 1` wherever domains overlap is the caller's
  (kryos') responsibility. Absent `owner` means `o = 1`. If `max(F_k) >
  1 + overlap_tol`, `merge_end` stops, naming the domains that cover the
  offending cells: overlapping domains need owner weights. Small excess within
  `overlap_tol` (sub-cell mismatch of owner boundaries) is accepted, logged as
  the double-counted amount.
- Helper for priority-ordered owners (kryos calls it in priority order, once at
  init; uses a cached regional -> regional `con` map):

  ```fortran
  subroutine coupler_owner_reduce(cpl, src_hi, owner_hi, dst, owner)
      ! owner := max(0, owner - remap(owner_hi, src_hi -> dst, "con", norm="dst"))
  ```

- Order independence: `merge_end` sums `F_dk`, `A_dk` in the order of the
  domain names (sorted), so the result is bit-identical for any order of
  `merge_add` calls. Adding the same `src` twice is an error.
- Cached: the `src -> dst` maps, and per map slot the vectors
  `area_src/wfull_src·q` and `1/area_dst` (rebuilt when `coupler_set_area` is
  called). The accumulator holds `2·N·n_global` dp values (T64: 0.3 MB per
  domain).

### 2.4 `coupler_split`: global -> one regional grid

```fortran
subroutine coupler_split(cpl, src, var_src, dst, var_dst, method, ref_src, ref_dst, &
                         missing_value, fill_method, mask2)
    type(coupler_class),   intent(inout) :: cpl
    character(len=*),      intent(in)    :: src, dst
    real(dp),              intent(in)    :: var_src(:,:)            ! global
    real(dp), allocatable, intent(inout) :: var_dst(:,:)            ! regional
    character(len=*), optional, intent(in) :: method                ! default "bilin"
    real(dp), optional,    intent(in)    :: ref_src(:,:), ref_dst(:,:)
```

Without references: `remap`. With both: `var_dst = ref_dst + remap(var_src -
ref_src)`, where a missing anomaly counts as 0 (the reference is kept). One
reference without the other is an error. The caller loops over domains. For
`geo` under decision 4, VILMA2 already exports anomalies (`dz_bed`, `dz_ss`), so
kryos calls `split` without references and adds the domain reference itself;
the anomaly form serves absolute global products (climate, ocean). Cached: the
map only.

## 3. Gaussian grids: bounds from Gauss weights

Problem: cell edges are axis midpoints, in three copies (`grid_area`,
`conservative::axis_edges`, `cell_corners_grid`). For a Gauss grid this neither
reproduces the quadrature weights nor reaches the poles (nlat=96: the polar row
stops near ±89.5°), hence VILMA2's global rescale.

Change, in this order:

1. `grid_axis_class` gains `real(dp), allocatable :: x_bnds(:), y_bnds(:)`
   (`0:nx`, `0:ny`), always filled at `grid_init`. Default fill uses today's
   midpoint/extrapolation arithmetic (`(a+b)/2` and `0.5*(a+b)` are the same
   floating-point operation), and all three consumers read the bounds instead
   of recomputing them. Pure refactor, bit-identical. Per-path clamps (±90 in
   `conservative_latlon`) and `unwrap_lon` stay where they are. A flag
   `G%explicit_bnds` records whether bounds were given. Caveat: for lon180
   grids built from 0..360 centres, `grid_area` and `cell_corners_grid` take raw
   midpoints across the seam (wrong seam cell) while `conservative_latlon`
   unwraps first; step 1 keeps both, and the seam fix is a separate,
   non-bit-identical change.
2. `grid_init_from_opts` gains `x_bnds=`, `y_bnds=` (validated: monotonic,
   each centre inside its cell).
3. New constructor:

   ```fortran
   subroutine grid_init_gaussian(grid, name, nlon, nlat, planet, lon180, x0)
       ! lat from gaussian_latitudes_calc(nlat, lat, wts); y_bnds from
       ! sin φ_{j+1/2} = sin φ_{j-1/2} + w_j, starting at -1, last edge set to +1
   ```

   Gauss-Legendre nodes interlace with the cumulative weights, so each node
   lies inside its band.
4. `grid_area` for grids with `explicit_bnds` and a lon-lat/Gauss `mtype` uses
   the latitude-band area (sphere: `R²Δλ(sinφ₂ - sinφ₁)`; ellipsoid: the
   authalic zone formula). On a sphere this equals `R²Δλ·w_j`, the SHTns
   quadrature weight. Grids without explicit bounds keep the Karney 4-corner
   area (bit-identical; §10, Q6).
5. Polar cells now have two corners at the pole: `build_sphere_poly` (now
   `axis_cell_poly` and the target strips, §11) must drop repeated vertices (degenerate edges break the great-circle clip).
6. `grid_cdo`: the writer emits cdo `ybounds` (2·ny values) only when
   `explicit_bnds`; the reader parses `ybounds` into `y_bnds`. Existing
   `grid_*.txt` files have no `ybounds` and read exactly as today. Also,
   cdo then generates the same cells when `gen="cdo"`.

Result for VILMA2: `ll -> Gauss` uses `conservative_latlon`, whose metric
`Δλ·Δsinφ` on Gauss bounds is proportional to the quadrature, so the SHTns
integral of the remapped field equals the band-metric integral of the source,
with no rescale. Use the `"Spherical Earth"` planet for the Gauss grid (SHTns is
spherical); the ratio is radius-independent.

## 4. Map cache identity

Today the key is `scrip-<method>_<src>_<dst>.nc`. Equal names with different
axes collide silently (VILMA2 works around it by tagging names with
dimensions).

- **Grid signature**: `cs` (mtype, units, planet a/f, lon180, xy_conv, all
  projection parameters) + `nx, ny` + first/last/sum of `x`, `y`, `x_bnds`,
  `y_bnds`. Compared with tolerance (`1e-9·cell size`), because a grid read
  back from `grid_<name>.txt` (`f0.12`) differs from the in-memory grid in the
  last bits, and an exact hash would raise false mismatches.
- **Short hash**: FNV-1a 64 over the signature with reals rounded to 6
  significant digits; the first 8 hex characters (`h8`) go into file names only
  (below).
- **Metadata**: new cache files carry global attributes `coords_src_sig`,
  `coords_dst_sig`, `coords_method`, `coords_gen`, `coords_weight_cutoff`,
  `coords_format`. cdo-generated files get them added after cdo returns.
- **On load** (`map_load_wm`, so plain `map_init` callers are covered too):
  1. always check `n_src`, `n_dst` against the grids (stops: a wrong-size map
     is never valid; today it indexes out of bounds);
  2. signature present and different -> policy `on_mismatch = "error"`
     (default; message names both signatures and the file) or `"rebuild"`
     (regenerate and overwrite, logged);
  3. no signature (legacy or foreign file) -> load as today, with a one-time
     "unverified" log line. Bit-identical.
- **File names** stay as today by default. Opt-in `coupler_init(...,
  hash_names=.true.)` / `map_init(..., hash_names=.true.)` names files
  `scrip-<method>_<src>-<h8>_<dst>-<h8>.nc`, so two grids with one name can
  coexist in one folder. Not the default: it would orphan every existing cache
  and force regeneration (and a cdo-made cache would be replaced by a
  coords-made one, which is not bit-identical).

Known hole: a legacy file of the same size but other axes is accepted
unverified. This is exactly the VILMA2 migration (midpoint -> Gauss bounds at
equal `nlat`), so that migration must clear its map folder or rename its grid
(§9).

## 5. Coupler housekeeping

- **Identity shortcut**: `remap` with `src == dst` copies (`mask_pack`,
  `reset`, `fill_method`, `filt_*` still apply; `mask2 = var /= missing`,
  `frac2 = merge(1,0,mask2)`); no map, no registration needed.
  `coupler_prime` with `src == dst` is a no-op. Not bit-identical for a caller
  that today remaps a grid onto itself through the coupler (`(w·v)/w` can
  differ from `v` by 1 ulp). yelmox already short-circuits; other callers must
  be checked by grep before landing.
- **Pass-through**: every `remap` specific gains `missing_value, norm, mask2,
  frac2, reset, mask_pack, fill_method, filt_method, filt_par`, forwarded
  unchanged to `map_field`. With `norm="dst"`/`"sum"` the coupler supplies the
  registered `area1`/`area2`. For 3D, `mask2`/`frac2` are 3D and the map is
  applied per layer. With `reset=.false.`, `var_dst` must already be allocated
  at the target shape (error otherwise). The `int` specifics get only
  `missing_value, mask2, reset, mask_pack`.
- **Containers without caps or deep copies**:

  ```fortran
  type grid_slot ; character(256) :: name ; type(grid_class), allocatable :: grid
                   real(dp), allocatable :: area(:,:) ; logical :: from_disk = .false. ; end type
  type map_slot  ; character(256) :: src, dst ; character(32) :: method
                   type(map_class), allocatable :: map
                   real(dp), allocatable :: c_src(:), rdst(:) ; end type      ! §2.3 cache
  ! coupler_class: grids(:), maps(:) allocatable slot arrays + ngrids, nmaps,
  !                map_fldr, gen, on_mismatch, hash_names
  ```

  Growth doubles the slot array and moves each `grid`/`map` component with
  `move_alloc`, so the large arrays are never copied, and indices stay valid.
  Slots are never assigned with `=`. Rejected: a linked list (pointer
  ownership, explicit cleanup). `GRID_MAX`/`MAP_MAX` go.
- **Deep copies**: `get_map` passes `cpl%grids(ig)%grid` straight to
  `map_init` (no `resolve_grid` copy). A grid resolved from `grid_<name>.txt`
  is read once into a `from_disk` slot and reused. `coupler_add_grid` keeps
  `intent(in)` and copies once (init-time). New `coupler_take_grid(cpl, name,
  grid)` takes `type(grid_class), allocatable, intent(inout)` and moves it,
  for grids where a second copy costs GBs (1 km Antarctica: ~1.4 GB per copy).
  Adding a name that is already a `from_disk` slot: no-op if the signatures
  agree, error otherwise.

## 6. Precision

All accumulation, weights, coverage and areas stay `real(dp)` (already true).
`remap`/`map_field` keep their concrete `dp` and `sp` specifics (not `wp`;
`coupler.f90` note and `sustainability-plan.md` Step 2). New outputs `frac2`
and inputs `area1`/`area2`/`owner` are `real(dp)` in every specific. `merge`,
`split` and `integral` are `dp`-only, consistent with kryos decision 2 (backends
convert at their own boundary).

## 7. Performance

- Regional projected -> Gauss/lon-lat goes through `conservative_spherical`:
  a k-d tree query with radius `chord(srcdeg + 1.5·tgtdeg)`, then a
  great-circle clip per candidate pair. For a fine source on a coarse target
  the candidate disk is ~7× the overlapping sources, so most clips are
  wasted. Benchmarked in `bench_conservative.f90` (ANT-16KM -> T64); since §11
  the loop is source-major (timings there).
- Improvements, in order: (i) a spherical-cap pre-reject before the clip
  (rejects only pairs that cannot overlap, so the weights do not change);
  (ii) opt-in adaptive NSUB for source cells smaller than ~0.1° (fewer
  vertices; changes the weights of newly generated maps, so not the default).
- OpenMP: source-parallel with a `static` schedule (§11; it was
  target-parallel, bounded by the ~10³ covered targets of one ice sheet on T64).
- Everything is built once and cached (disk + memory). Per coupling step the
  cost is one `O(n_links)` apply per domain plus `O(N·n_global)` in
  `merge_end`, negligible next to the models. The `q`, `r`, `c_src` vectors
  cost `O(n_links)` once per load.

## 8. Tests (test/coords/)

| File | Checks |
|---|---|
| `test_coupler_global.f90` | Polar-stereographic disk domain (load inside a disk, missing outside) onto global 1° lon-lat and onto `grid_init_gaussian` T64. (a) `norm="dst"` with registered areas: `coupler_integral` equal on both grids to 1e-12 relative. (b) `frac2 ∈ [0,1]`; `Σ wfull_dst·frac2 = Σ_valid wcov_src` to 1e-10; 1 for targets inside the disk, 0 outside. (c) Two overlapping regional grids (different projections), owners from `coupler_owner_reduce`: merged integral = background part + Σ owned regional integrals; `merge(A,B)` bit-identical to `merge(B,A)` (`==`); without owners `max(frac_tot) > 1 + tol`. (d) `kind="intensive"` merge of a constant returns the constant to 1e-14. (e) `split` with `var_src == ref_src` returns `ref_dst` exactly. (f) A pruned and an unpruned map give the same integrals to 1e-13. |
| `test_gauss_bounds.f90` | `grid_init_gaussian(nlat=96, nlon=192, "Spherical Earth")`: bounds monotonic, ±90 at the ends, each node inside its band; `area = R²Δλ·w_j` to 1e-13; `Σ area = 4πR²`; ll -> Gauss `con` of 1 is 1; `Σ_k Δλ w_j v_k` equals the band-metric source integral to 1e-13 (no rescale); `grid_cdo` write/read keeps `y_bnds`; a polar-cell spherical clip returns no NaN. |
| `test_grid_signature.f90` | Two grids with one name and different `dy` in one folder: `on_mismatch="error"` stops (run as a separate process, exit code checked), `"rebuild"` regenerates; an in-memory grid and its `grid_<name>.txt` round trip have equal signatures; a legacy file without attributes loads; size mismatch stops; `hash_names=.true.` gives two files that coexist. |
| `test_coupler_bitident.f90` | Fixed call set (con/bilin/nn, dp/sp/int, `stat` variants, `reset=.false.`, `mask_pack`, `fill_method`, 3D) written to a stream file; `bitident.sh <ref-commit>` builds the same test at the pre-change commit in a worktree and compares with `cmp`. In-binary: `frac2` present vs absent gives identical `var2`; `remap` vs direct `map_field` identical. |

Plus extend `test_coupler.f90`: more than 16 grids / 64 maps (growth), identity
shortcut, pass-through of `mask2`.

## 9. Migration order

0. `test_coupler_bitident` + `bitident.sh` against the current head, before
   any change.
1. Bounds refactor (§3.1): one source of cell edges. Bit-identical.
2. Coupler containers, deep-copy removal, pass-through, identity shortcut (§5).
   Bit-identical except self-remaps (checked by grep first).
3. Map signature + size check in `map_load_wm` (§4). Bit-identical.
4. Prune criterion fix + `wfull`/`wcov` in generator and cache + `norm`/`frac2`
   (§1). Default paths untouched; new caches carry the new variables.
5. `coupler_set_area`, `coupler_integral`, `merge`, `split`,
   `coupler_owner_reduce` (§2).
6. Gauss bounds: `grid_init_gaussian`, explicit-bounds `grid_area`, polar
   polygon fix, `grid_cdo` `ybounds` (§3.2–3.6).
7. Performance pre-reject (§7 i).

Downstream:
- **VILMA2** (`vilma_remap.f90`): build the Gauss grid with
  `grid_init_gaussian` (or `y_bnds` from the SHTns weights, `mtype="gaussian"`
  instead of today's `"latlon"`); delete `conserve_mass` and its rescale; clear
  the map folder once (§4 hole). The dimension-tagged names may stay
  (coexistence), but no longer carry correctness. If a check is kept, it must
  use band areas, not the ellipsoidal `ll%area` the current check uses.
- **yelmox** (`libs/kryos.f90:~680-751`): `remap_2D`/`remap_3D` reduce to
  `remap`; `prime_map` reduces to `coupler_prime`. `remap_method_smooth` stays
  in kryos unless Q9 decides otherwise. Re-run the phase 1 bit-identity gate.
- **climber-x**: the `reset=.false.` overlay in `ice_to_geo` can become
  `merge` when it adopts the coupler (not required).

## 10. Open questions

1. Default model area for projected grids: `grid%area` (true ellipsoidal) or
   the plane `dx·dy` that Yelmo integrates? It decides which volume a merge
   conserves; the API supports both via `coupler_set_area`.
2. `norm="sum"` is the source-partition sum (metric-free), not SCRIP's raw
   `Σ w v` (units depend on the generator). Is the raw sum needed anywhere?
3. Prune criterion: `cutoff·max(w)` changes newly generated maps everywhere
   slivers are dropped (existing caches, hence the bit-identity gate, are
   unaffected). Acceptable, or restrict the change to targets with
   `> 1/cutoff` links?
4. Legacy cache without a signature: accept unverified (proposed default) or
   rebuild?
5. Support `norm /= "frac"` on cdo SCRIP `con` maps (would need
   `dst_grid_frac` and the file's normalization option), or coords-only?
6. Use the band area for every lon-lat grid (correct cell definition, but
   changes `grid%area` for existing callers), or only for explicit bounds?
7. Merge API: accumulator (proposed) vs one-shot `coupler_merge(cpl, dst,
   srcs(:), fields(:))` with a pointer wrapper type.
8. `overlap_tol` default 1e-2, and stop vs. warn above it.
9. Absorb `remap_method_smooth` as `method="auto"` (bilin when refining, con
   when coarsening, by `sqrt` of median model area), or keep it as kryos policy?
10. `from_disk` auto-registration: OK that a later `coupler_add_grid` of the
    same name is an error when the signatures differ?
11. `frac2` always `real(dp)` (proposed), or the kind of `var2`?

## 11. Spherical clip geometry (follow-up to step 7)

Two problems of `conservative_spherical` found in step 7, both fixed. They
change newly generated projected -> lon-lat/Gauss `con` maps only (in the
bit-identity harness: `bi-ps` onto `bi-ll2`, `bi-ls`, `bi-gauss`, 49 of 409
streams in `gen` mode; `cached` 409 of 409 identical).

1. **Source polygons.** The 4 inverse-projected corners were interpolated
   linearly in lon/lat, so a cell across the 0/360 meridian or around the pole
   ran the long way round (one garbage meridian strip for every pole-covering
   ANT/GRL map). Now each edge is sampled NSUB = 8 times in the grid's own
   plane and every sample is inverse-projected (`axis_cell_poly`); the pole
   cell is a polygon around the pole. NSUB stays 8.
2. **Concave target cells.** The poleward parallel of a lon-lat cell, traced
   by great-circle chords, is concave, and Sutherland-Hodgman clipped away the
   slivers behind it. Now the source polygon is clipped against the cell's
   NSUB strips in longitude (two meridians, one chord of each parallel; convex;
   a triangle at a pole) and the piece areas are summed
   (`target_strips_overlap`). Strips of adjacent cells share edges, so the
   targets tile the sphere. Target longitudes come from `lon_edges` (the
   unwrapped axis, which also fixes the wrong default seam cell of a lon180
   target built from 0..360 centres) and latitude bounds are clamped to ±90,
   as in `conservative_latlon`.

Measured residuals, ANT-16KM polar stereographic (WGS84, phi = -71,
pole-centred) onto global targets that cover it; "before" is 6271c75, "fix 1"
the source polygons alone (concave clip only):

| Residual | Target | before | fix 1 | after |
|---|---|---|---|---|
| target row holding 70°S, `1 - wcov_dst/wfull_dst` (lon 10–350) | 1° lon-lat [-71,-70] | 1.3e-3 | 1.3e-3 | 8.6e-10 |
| same | T64 [-70.89,-69.02] | 2.4e-3 | 2.4e-3 | 4.6e-9 |
| source cells at 69.5–70.5°S, max `|1 - Σ_k φ_jk|` | 1° | 0.98 (seam) | 1.7e-2 | 1.7e-9 |
| all source cells, max `|1 - Σ_k φ_jk|` | 1° | 0.98 | 3.0e-2 | 2.4e-9 |
| same | T64 | 0.97 | 0.11 | 1.9e-9 |
| source cells with `|1 - Σ_k φ_jk|` > 1e-2 | 1° / T64 | 7748 / 11238 | 7592 / 11076 | 0 / 0 |

On a spherical planet `Σ_k φ_jk = 1` to rounding (ANT-100KM -> spherical T64:
7.7e-13; 1.8e-10 for the cell holding the pole, which collects the round-off of
~1500 geodesic piece areas). On WGS84 it closes to ~1e-8 (100 km cells: 3.4e-8):
the pieces have great-circle edges, and their geodesic areas on the ellipsoid
add up to the cell's only to O(f·L²/R). A spherical-excess area metric would
close to rounding and be cheaper than the Karney area, but changes the map
metric; not done. `bench_conservative`'s `max|co-cdo|` for ANT-16KM -> T64
(1.7, in the polar row) comes from cdo: on the target cells inside the domain
(86–50°S) coords is within 4.4e-3 of the cell-centre value of the linear test
field, cdo off by up to 0.8.

Loop and cost. Each source polygon is built once (source-major loop; candidate
targets from a k-d tree over the target centres; the cap pre-reject of §7 i
kept). A strip with every source vertex outside one of its edges is skipped,
and a source polygon inside one strip takes its own area (computed once) as the
overlap. The static schedule gives every thread one contiguous source range,
so each target's links are in ascending source order for any thread count
(weights identical serial / 1 / 4 / 10 threads). ANT-16KM -> T64 (WGS84,
`bench_conservative`, mean of 3): before 3.77 s serial, 2.69 s on 4 threads,
1.60 s on 10; after 4.32 s, 1.31 s, 0.81 s. The serial +15 % is the NSUB
inverse projections per edge (~150 ns each).
