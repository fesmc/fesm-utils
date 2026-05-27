#!/usr/bin/env python3
"""
Build driver for the fesm-utils autotools dependencies (fftw, lis).

Builds either component, for any machine defined under machines/, with or
without OpenMP, from one command.

The `utils` library is NOT built here: it is configured and built by `configme`
(which owns the central machine/compiler fragments + netCDF auto-detection) as
the `fesm-utils/utils` subpackage. Use:
    configme install fesm-utils         # clone + configure + build (fftw/lis + utils)
    configme config  fesm-utils/utils   # (re)generate the utils Makefile only

Examples
--------
    ./build.py -m dkrz_levante -c ifx                       # fftw + lis, both variants
    ./build.py -m macbook -c gfortran --component lis --variant serial
    ./build.py -m pik_hpc2024 -c ifx --component fftw --variant omp
    ./build.py --list-machines

Configuration for each machine lives in machines/<name>.toml. This driver does
not replace the existing install_*.sh scripts; it wraps the fftw/lis autotools
builds behind one consistent interface.
"""

import argparse
import glob
import os
import re
import shlex
import subprocess
import sys
from pathlib import Path

try:
    import tomllib  # Python >= 3.11
except ModuleNotFoundError:
    try:
        import tomli as tomllib  # pip install tomli (older Python)
    except ModuleNotFoundError:
        sys.exit(
            "error: need Python >= 3.11 (for tomllib) or `pip install tomli`."
        )

ROOT = Path(__file__).resolve().parent
MACHINE_DIR = ROOT / "machines"
FFTW_SRC = "fftw-3.3.10"
LIS_SRC = "lis-2.1.6"
COMPONENTS = ["fftw", "lis"]


def die(msg):
    sys.exit(f"error: {msg}")


def list_machines():
    # TEMPLATE.toml is a documented skeleton to copy, not a real machine.
    return sorted(
        p.stem for p in MACHINE_DIR.glob("*.toml") if p.stem != "TEMPLATE"
    )


def load_machine(name):
    path = MACHINE_DIR / f"{name}.toml"
    if not path.exists():
        die(f"unknown machine '{name}'. Available: {', '.join(list_machines())}")
    with open(path, "rb") as f:
        return tomllib.load(f)


def configure_vars(machine, component, compiler):
    """Resolve the autotools configure variables for (component, compiler).

    A per-component compiler table fully replaces the machine-level one.
    """
    comp = machine.get("components", {}).get(component, {})
    table = comp.get("compilers", {}).get(compiler)
    if table is None:
        table = machine.get("compilers", {}).get(compiler)
    if table is None:
        avail = ", ".join(sorted(machine.get("compilers", {}))) or "(none)"
        die(
            f"compiler '{compiler}' not defined for machine "
            f"'{machine['machine']['name']}'. Available: {avail}"
        )
    return table


def homebrew_gcc():
    """Newest Homebrew `gcc-N` on macOS (basename), or None if none installed.

    Apple Clang -- the default `cc` -- has no OpenMP support, so FFTW/LIS
    configure fail their `--enable-openmp` check. Homebrew never installs an
    unversioned `gcc` (that would shadow Apple's), only versioned binaries like
    `gcc-15`, so we glob the standard prefixes and pick the highest version.
    """
    best = None
    for prefix in ("/opt/homebrew/bin", "/usr/local/bin"):
        for path in glob.glob(f"{prefix}/gcc-[0-9]*"):
            m = re.fullmatch(r"gcc-(\d+)", os.path.basename(path))
            if m:
                ver = int(m.group(1))
                if best is None or ver > best[0]:
                    best = (ver, os.path.basename(path))
    return best[1] if best else None


def apply_default_cc(component, table, omp):
    """On macOS, default CC to Homebrew gcc when the config leaves it unset.

    Mutates `table` in place. Warns (only for OpenMP builds, where it actually
    breaks) if no Homebrew gcc is found and the build will fall back to Apple
    Clang.
    """
    if sys.platform != "darwin" or "CC" in table:
        return
    cc = homebrew_gcc()
    if cc:
        table["CC"] = cc
        print(f"  {component}: macOS -- using Homebrew {cc} as C compiler (CC={cc})")
    elif omp:
        print(
            "  \033[33mwarning:\033[0m no Homebrew gcc found; the C compiler will be "
            "Apple Clang, which lacks OpenMP support, so this --enable-openmp build "
            "will likely fail.\n"
            "           Install one with: brew install gcc"
        )


def modules_for(machine, component):
    mods = list(machine.get("modules", {}).get("load", []))
    comp_mods = (
        machine.get("components", {}).get(component, {}).get("modules", {}).get("load")
    )
    if comp_mods is not None:
        mods += comp_mods
    return mods


def module_prefix(mods):
    return "".join(f"module load {m} && " for m in mods)


def var_args(table):
    """Format configure variables as KEY='value' shell words."""
    return " ".join(f"{k}={shlex.quote(str(v))}" for k, v in table.items())


def intel_runtime_preamble(table):
    """Handle the `intel_runtime_libs` flag for Intel/ifx autotools builds.

    Older configs hard-coded `ac_cv_c_libs` as a long list of machine-specific
    spack paths, which only existed on the one cluster they were copied from.
    Instead, ask the compiler where its runtime libs live -- robust and
    machine-independent.

    If the flag is set, mutate `table` in place (drop the flag and any static
    `ac_cv_c_libs`) and return (preamble, extra_var):
      - preamble: shell that derives the lib dir from `$CC -print-file-name`.
      - extra_var: an `ac_cv_c_libs=...` configure word referencing it.
    Otherwise return ("", "").
    """
    if not table.pop("intel_runtime_libs", False):
        return "", ""
    cc = str(table.get("CC", "icx"))
    # A static ac_cv_c_libs would override the derived one; drop it.
    table.pop("ac_cv_c_libs", None)
    preamble = (
        f'_oneapi_libdir="$(dirname "$({shlex.quote(cc)} '
        '-print-file-name=libimf.so)")" && '
    )
    extra_var = (
        'ac_cv_c_libs="-L$_oneapi_libdir -lsvml -lirng -limf -lirc -ldl -lm"'
    )
    return preamble, extra_var


DRY_RUN = False


def run(cmd, cwd, login=False):
    print(f"\n\033[1m>>> [{cwd.relative_to(ROOT)}] {cmd}\033[0m\n", flush=True)
    if DRY_RUN:
        return
    # A login shell (`bash -lc`) is used ONLY when the machine config loads
    # modules, because `module` is a shell function that needs the lmod init a
    # login shell sources. When no modules are loaded (e.g. the generic `linux`
    # machine), a login shell would re-source the user's rc files and reload a
    # default module set, silently overriding the environment they already set
    # up (a sourced module script, a hand-loaded gcc, etc.). So with no modules
    # we use a plain non-login `bash -c` that inherits that environment verbatim.
    # On failure raise a CalledProcessError carrying the human-readable command
    # (`cmd`) rather than the ['bash', '-lc', 'set -e; ...'] wrapper list, so the
    # summary shows a one-line command you can copy-paste and re-run.
    shell_flag = "-lc" if login else "-c"
    proc = subprocess.run(["bash", shell_flag, f"set -e; {cmd}"], cwd=cwd)
    if proc.returncode != 0:
        raise subprocess.CalledProcessError(proc.returncode, cmd)


# --- Per-component builders -------------------------------------------------
#
# Each builder takes the same signature (machine, compiler, omp) and is
# registered in BUILDERS below. To add a new package, write one build_<name>
# function and add it to BUILDERS (and to COMPONENTS) -- nothing else changes.
# `build_autotools` covers the common configure/make/install case shared by any
# autotools package.


def build_autotools(machine, component, compiler, omp, src, base_opts, omp_opt):
    """Configure + make + install an autotools package into <name>-{omp,serial}."""
    suffix = "omp" if omp else "serial"
    prefix = ROOT / f"{component}-{suffix}"
    opts = [f"--prefix={shlex.quote(str(prefix))}"] + list(base_opts)
    if omp:
        opts.append(omp_opt)
    table = dict(configure_vars(machine, component, compiler))
    apply_default_cc(component, table, omp)
    preamble, extra_var = intel_runtime_preamble(table)
    vargs = " ".join(filter(None, [var_args(table), extra_var]))
    mods = modules_for(machine, component)
    cmd = (
        module_prefix(mods)
        + preamble
        + f"./configure {' '.join(opts)} {vargs} && make clean && make && make install"
    )
    run(cmd, ROOT / src, login=bool(mods))


def build_fftw(machine, compiler, omp):
    build_autotools(
        machine, "fftw", compiler, omp, FFTW_SRC,
        base_opts=["--disable-doc"], omp_opt="--enable-openmp",
    )


def build_lis(machine, compiler, omp):
    build_autotools(
        machine, "lis", compiler, omp, LIS_SRC,
        base_opts=["--enable-f90", "--enable-saamg"], omp_opt="--enable-omp",
    )


BUILDERS = {
    "fftw": build_fftw,
    "lis": build_lis,
}


def main():
    ap = argparse.ArgumentParser(
        description="Build driver for fesm-utils autotools deps (fftw, lis).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument("-m", "--machine", help="machine name (see --list-machines)")
    ap.add_argument("-c", "--compiler", help="compiler key, e.g. ifx, ifort, gfortran")
    ap.add_argument(
        "--component",
        default="all",
        choices=COMPONENTS + ["all"],
        help="which component to build (default: all)",
    )
    ap.add_argument(
        "--variant",
        default="both",
        choices=["omp", "serial", "both"],
        help="build the OpenMP-enabled variant, the serial variant, or both (default: both)",
    )
    ap.add_argument(
        "--list-machines", action="store_true", help="list available machines and exit"
    )
    ap.add_argument(
        "--dry-run",
        action="store_true",
        help="print the build commands without executing them",
    )
    ap.add_argument(
        "--check",
        action="store_true",
        help="validate the machine config (compilers defined) and print the "
        "resolved commands, without building",
    )
    args = ap.parse_args()

    global DRY_RUN
    # --check is a no-build validation pass: it reuses the dry-run print path,
    # while the resolution code (configure_vars) raises clear errors for any
    # missing compiler table.
    DRY_RUN = args.dry_run or args.check

    if args.list_machines:
        print("Available machines:")
        for name in list_machines():
            print(f"  - {name}")
        return

    if not args.machine or not args.compiler:
        ap.error("-m/--machine and -c/--compiler are required")

    machine = load_machine(args.machine)

    components = COMPONENTS if args.component == "all" else [args.component]
    omp_variants = {"omp": [True], "serial": [False], "both": [False, True]}[
        args.variant
    ]

    # Each (component, variant) is built independently: a failure in one is
    # recorded and the rest still run, so a single broken target doesn't abort
    # the whole build. Results are summarized at the end.
    results = []  # (label, ok, error_message_or_None)
    for component in components:
        for omp in omp_variants:
            label = f"{component}-{'omp' if omp else 'serial'}"
            print(f"\n=== Building {label} "
                  f"[{args.machine}/{args.compiler}] ===")
            try:
                BUILDERS[component](machine, args.compiler, omp)
                results.append((label, True, None))
            except (subprocess.CalledProcessError, SystemExit) as e:
                print(f"\n\033[31m✗ {label} failed:\033[0m {e}")
                results.append((label, False, str(e)))

    succeeded = [label for label, ok, _ in results if ok]
    failed = [(label, err) for label, ok, err in results if not ok]

    if args.check:
        title, note = "Config check", " (nothing was built)"
    elif DRY_RUN:
        title, note = "Dry run", " (nothing was built)"
    else:
        title, note = "Build summary", ""
    print(f"\n\033[1m{title}:\033[0m {args.machine}/{args.compiler}{note}")
    if succeeded:
        print(f"  \033[32m✓\033[0m {', '.join(succeeded)}")
    if failed:
        print(f"  \033[31m✗\033[0m {', '.join(label for label, _ in failed)}")
        for label, err in failed:
            first = (err or "").strip().splitlines()
            print(f"      {label}: {first[0] if first else 'failed'}")
        sys.exit(1)


if __name__ == "__main__":
    main()
