#!/usr/bin/env python3
"""
Unified build driver for fesm-utils.

Builds any combination of the three components (fftw, lis, utils), for any
machine defined under machines/, with or without OpenMP, from one command.

Examples
--------
    ./build.py -m dkrz_levante -c ifx                       # everything, both variants
    ./build.py -m macbook -c gfortran --component utils --variant serial
    ./build.py -m pik_hpc2024 -c ifx --component lis --variant omp
    ./build.py --list-machines

Configuration for each machine lives in machines/<name>.toml. This driver does
not replace the existing install_*.sh scripts or utils/config.py + Makefile --
those still work by hand. It just wraps them behind one consistent interface.
"""

import argparse
import glob
import os
import re
import shlex
import shutil
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
COMPONENTS = ["fftw", "lis", "utils"]


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


def run(cmd, cwd):
    print(f"\n\033[1m>>> [{cwd.relative_to(ROOT)}] {cmd}\033[0m\n", flush=True)
    if DRY_RUN:
        return
    subprocess.run(["bash", "-lc", f"set -e; {cmd}"], cwd=cwd, check=True)


# --- Per-component builders -------------------------------------------------
#
# Each builder takes the same signature (machine, compiler, omp, debug) and is
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
    cmd = (
        module_prefix(modules_for(machine, component))
        + preamble
        + f"./configure {' '.join(opts)} {vargs} && make clean && make && make install"
    )
    run(cmd, ROOT / src)


def build_fftw(machine, compiler, omp, debug):
    build_autotools(
        machine, "fftw", compiler, omp, FFTW_SRC,
        base_opts=["--disable-doc"], omp_opt="--enable-openmp",
    )


def build_lis(machine, compiler, omp, debug):
    build_autotools(
        machine, "lis", compiler, omp, LIS_SRC,
        base_opts=["--enable-f90", "--enable-saamg"], omp_opt="--enable-omp",
    )


def build_utils(machine, compiler, omp, debug):
    openmp = 1 if omp else 0

    preamble = module_prefix(modules_for(machine, "utils"))

    # The utils Makefile is owned by `configme` (central machine/compiler
    # fragments + auto-detected netCDF). If it already exists (e.g. generated by
    # an earlier `configme` run), reuse it; otherwise generate it now via the
    # config-only `configme <pkg>` command -- which does NOT invoke build.py, so
    # there is no recursion.
    makefile = ROOT / "utils" / "Makefile"
    if makefile.exists():
        print("  utils: using existing Makefile (skipping configure)")
        run(preamble + f"make openmp={openmp} debug={debug} fesmutils-static",
            ROOT / "utils")
        return

    if shutil.which("configme") is None:
        die(
            "utils/Makefile not found and `configme` is not installed.\n"
            "  configme generates the utils Makefile from central machine/"
            "compiler settings.\n"
            "  Install it:\n"
            "      pip install git+https://github.com/fesmc/configme\n"
            "  (if the `configme` command is then not found, add the Python user\n"
            "   bin to PATH: export PATH=\"${PATH}:${HOME}/.local/bin\")\n"
            "  Then re-run this build. To bypass configme entirely you can still\n"
            "  configure by hand: cd utils && python3 config.py config/<machine>_<compiler>"
        )

    mname = machine["machine"]["name"]
    # Generate utils/Makefile (run from the repo root so `configme fesm-utils`
    # resolves this checkout), with modules loaded so nf-config/nc-config are
    # found, then build.
    run(preamble + f"configme fesm-utils -m {shlex.quote(mname)} "
        f"-c {shlex.quote(compiler)}", ROOT)
    run(preamble + f"make openmp={openmp} debug={debug} fesmutils-static",
        ROOT / "utils")


BUILDERS = {
    "fftw": build_fftw,
    "lis": build_lis,
    "utils": build_utils,
}


def main():
    ap = argparse.ArgumentParser(
        description="Unified build driver for fesm-utils (fftw, lis, utils).",
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
        "--debug",
        type=int,
        default=0,
        choices=[0, 1, 2],
        help="utils debug level: 0 normal, 1 debug, 2 profile (default: 0)",
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
        help="validate the machine config (compilers defined, utils config "
        "files present) and print the resolved commands, without building",
    )
    args = ap.parse_args()

    global DRY_RUN
    # --check is a no-build validation pass: it reuses the dry-run print path,
    # while the resolution code (configure_vars, build_utils) raises clear
    # errors for any missing compiler table or utils config file.
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
                BUILDERS[component](machine, args.compiler, omp, args.debug)
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
