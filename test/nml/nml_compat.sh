#!/usr/bin/env bash
# Parse compatibility of nml against a reference commit.
#
# Usage (from the fesm-utils root, after `configme config ...` generated the
# root Makefile):
#     test/nml/nml_compat.sh <ref-commit> <file.nml> [<file.nml> ...]
#
# Lists every key of the files with test/nml/nml_keys.f90 built at the working
# tree (nml_group_names, nml_group_keys), then reads each key as a string
# vector with test/nml/test_nml_compat.f90 (the working-tree copy) built
# against libfesmutils at <ref-commit> (in a temporary git worktree) and at the
# working tree, and diffs the two outputs (one line per key:
# `file group key : n : v1 | v2 | ...`, or `: absent`). Exit status 0 when
# every key parses the same. Set NML_COMPAT_KEEP=1 to keep the scratch
# directory (key list and both outputs).

set -euo pipefail

ref=${1:?usage: test/nml/nml_compat.sh <ref-commit> <file.nml> [...]}
shift
[ $# -gt 0 ] || { echo "usage: test/nml/nml_compat.sh <ref-commit> <file.nml> [...]"; exit 2; }

root=$(git rev-parse --show-toplevel)
[ -f "$root/Makefile" ] || { echo "nml_compat: no root Makefile (run configme config first)"; exit 2; }

scratch=$(mktemp -d "${TMPDIR:-/tmp}/nml_compat.XXXXXX")
refwt="$scratch/ref"
cleanup() {
    git -C "$root" worktree remove --force "$refwt" >/dev/null 2>&1 || true
    if [ "${NML_COMPAT_KEEP:-0}" = "1" ]; then echo "nml_compat: kept $scratch"; else rm -rf "$scratch"; fi
}
trap cleanup EXIT

# Root Makefile for the reference tree: its own config/Makefile template with
# the compiler block of the current root Makefile (as test/coords/bitident.sh).
make_ref_makefile() {
    local tmpl="$root/config/Makefile" gen="$root/Makefile" out=$1 reftmpl=$2
    local p T G
    p=$(grep -n '^<COMPILER_CONFIGURATION>$' "$tmpl" | cut -d: -f1)
    T=$(wc -l < "$tmpl"); G=$(wc -l < "$gen")
    sed -n "${p},$((G - T + p))p" "$gen" > "$scratch/compiler_block.mk"
    sed -e "/^<COMPILER_CONFIGURATION>$/{r $scratch/compiler_block.mk" -e 'd;}' "$reftmpl" > "$out"
}

# Extra make target, read after the root Makefile, that links one program.
cat > "$scratch/nml_compat.mk" <<'MK'
nml-compat-exe : fesmutils-static
	$(FC) $(DFLAGS) $(FFLAGS) -o $(NML_COMPAT_EXE) $(NML_COMPAT_SRC) -L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
MK

build() {   # build <tree> <source> <exe>
    local tree=$1 src=$2 exe=$3
    make -C "$tree" -f Makefile -f "$scratch/nml_compat.mk" nml-compat-exe openmp=0 \
         NML_COMPAT_SRC="$root/$src" NML_COMPAT_EXE="$exe" > "$exe.build.log" 2>&1 \
        || { echo "nml_compat: build failed in $tree (log: $exe.build.log)"; tail -n 30 "$exe.build.log"; exit 2; }
}

files=()
for f in "$@"; do
    [ -f "$f" ] || { echo "nml_compat: no file $f"; exit 2; }
    files+=("$(cd "$(dirname "$f")" && pwd)/$(basename "$f")")
done

echo "nml_compat: reference $ref ($(git -C "$root" rev-parse --short "$ref")) vs working tree, ${#files[@]} files"
git -C "$root" worktree add --detach "$refwt" "$ref" > /dev/null 2>&1
make_ref_makefile "$refwt/Makefile" "$refwt/config/Makefile"

build "$root"  test/nml/nml_keys.f90        "$scratch/keys.x"
build "$refwt" test/nml/test_nml_compat.f90 "$scratch/ref.x"
build "$root"  test/nml/test_nml_compat.f90 "$scratch/head.x"

"$scratch/keys.x" "${files[@]}" > "$scratch/keys.txt"
"$scratch/ref.x"  "$scratch/keys.txt" > "$scratch/ref.out"
"$scratch/head.x" "$scratch/keys.txt" > "$scratch/head.out"

nkeys=$(grep -c . "$scratch/keys.txt" || true)
if diff "$scratch/ref.out" "$scratch/head.out" > "$scratch/diff.txt"; then
    echo "nml_compat: $nkeys keys, all parse the same"
    echo "nml_compat: PASS"
else
    cat "$scratch/diff.txt"
    echo "nml_compat: $nkeys keys, $(grep -c '^<' "$scratch/diff.txt" || true) parse differently (< $ref, > working tree)"
    echo "nml_compat: FAIL"
    exit 1
fi
