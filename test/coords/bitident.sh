#!/usr/bin/env bash
# Bit-identity gate for the coupler / map_field apply surface.
#
# Usage (from the fesm-utils root, after `configme config ...` generated the
# root Makefile):
#     test/coords/bitident.sh <ref-commit> [openmp=0]
#
# Builds test/coords/test_coupler_bitident.f90 (the copy in the working tree)
# against libfesmutils at <ref-commit> (in a temporary git worktree) and at the
# working tree, runs each binary, and compares every stream file in
# bitident_out/ with cmp. Two comparisons:
#   gen    - both runs start from an empty map folder (maps generated);
#   cached - both runs load the map cache that the reference run generated
#            (the apply path on identical map files).
# Exit status 0 when every stream is identical in the compared modes.
# BITIDENT_MODES selects them (default "gen cached"). Across a change that
# alters newly generated maps by design (the prune criterion fix, see
# docs/dev/coupler-global-design.md 1.2), only "cached" is expected identical
# against a reference before it: BITIDENT_MODES=cached.
# Set BITIDENT_KEEP=1 to keep the scratch directory.

set -euo pipefail

ref=${1:?usage: test/coords/bitident.sh <ref-commit> [openmp=0]}
omp=${2:-openmp=0}

root=$(git rev-parse --show-toplevel)
src=test/coords/test_coupler_bitident.f90
[ -f "$root/Makefile" ] || { echo "bitident: no root Makefile (run configme config first)"; exit 2; }

scratch=$(mktemp -d "${TMPDIR:-/tmp}/bitident.XXXXXX")
refwt="$scratch/ref"
cleanup() {
    git -C "$root" worktree remove --force "$refwt" >/dev/null 2>&1 || true
    if [ "${BITIDENT_KEEP:-0}" = "1" ]; then echo "bitident: kept $scratch"; else rm -rf "$scratch"; fi
}
trap cleanup EXIT

# Root Makefile for the reference tree: its own config/Makefile template with
# the compiler block of the current root Makefile (everything between the
# template's <COMPILER_CONFIGURATION> placeholder and the verbatim tail).
make_ref_makefile() {
    local tmpl="$root/config/Makefile" gen="$root/Makefile" out=$1 reftmpl=$2
    local p T G
    p=$(grep -n '^<COMPILER_CONFIGURATION>$' "$tmpl" | cut -d: -f1)
    T=$(wc -l < "$tmpl"); G=$(wc -l < "$gen")
    sed -n "${p},$((G - T + p))p" "$gen" > "$scratch/compiler_block.mk"
    sed -e "/^<COMPILER_CONFIGURATION>$/{r $scratch/compiler_block.mk" -e 'd;}' "$reftmpl" > "$out"
}

# Extra make target, read after the root Makefile, that links one program.
cat > "$scratch/bitident.mk" <<'EOF'
bitident-exe : fesmutils-static
	$(FC) $(DFLAGS) $(FFLAGS) -o $(BITIDENT_EXE) $(BITIDENT_SRC) -L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
EOF

build() {   # build <tree> <exe>
    local tree=$1 exe=$2
    make -C "$tree" -f Makefile -f "$scratch/bitident.mk" bitident-exe "$omp" \
         BITIDENT_SRC="$root/$src" BITIDENT_EXE="$exe" > "$exe.build.log" 2>&1 \
        || { echo "bitident: build failed in $tree (log: $exe.build.log)"; tail -n 30 "$exe.build.log"; exit 2; }
}

run() {     # run <exe> <dir> [map folder to seed]
    local exe=$1 dir=$2 seed=${3:-}
    mkdir -p "$dir"
    [ -n "$seed" ] && cp -R "$seed" "$dir/bitident_maps"
    (cd "$dir" && "$exe" > run.log 2>&1) || { echo "bitident: $exe failed in $dir"; tail -n 20 "$dir/run.log"; exit 2; }
}

compare() { # compare <label> <dirA> <dirB>
    local label=$1 a=$2/bitident_out b=$3/bitident_out n=0 nd=0 f
    for f in "$a"/*.bin; do
        n=$((n+1))
        if [ ! -f "$b/$(basename "$f")" ]; then
            echo "  MISSING [$label] $(basename "$f" .bin)"; nd=$((nd+1))
        elif ! cmp -s "$f" "$b/$(basename "$f")"; then
            echo "  DIFF    [$label] $(basename "$f" .bin)"; nd=$((nd+1))
        fi
    done
    echo "bitident [$label]: $((n-nd)) of $n streams identical"
    [ "$nd" -eq 0 ]
}

echo "bitident: reference $ref ($(git -C "$root" rev-parse --short "$ref")) vs working tree"
git -C "$root" worktree add --detach "$refwt" "$ref" > /dev/null 2>&1
make_ref_makefile "$refwt/Makefile" "$refwt/config/Makefile"
mkdir -p "$refwt/bin"

build "$refwt" "$scratch/ref.x"
build "$root"  "$scratch/head.x"

run "$scratch/ref.x"  "$scratch/ref_gen"
run "$scratch/ref.x"  "$scratch/ref_cached"  "$scratch/ref_gen/bitident_maps"
run "$scratch/head.x" "$scratch/head_gen"
run "$scratch/head.x" "$scratch/head_cached" "$scratch/ref_gen/bitident_maps"

rc=0
for mode in ${BITIDENT_MODES:-gen cached}; do
    compare "$mode" "$scratch/ref_$mode" "$scratch/head_$mode" || rc=1
done
[ $rc -eq 0 ] && echo "bitident: PASS" || echo "bitident: FAIL"
exit $rc
