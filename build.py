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
    vargs = var_args(configure_vars(machine, component, compiler))
    cmd = (
        module_prefix(modules_for(machine, component))
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
    cfg = machine.get("utils", {}).get("config", {}).get(compiler)
    if cfg is None:
        die(
            f"no utils config mapped for compiler '{compiler}' on machine "
            f"'{machine['machine']['name']}'. Add one under utils/config/ and "
            f"reference it in [utils.config]."
        )
    cfg_path = ROOT / "utils" / "config" / cfg
    if not cfg_path.exists():
        die(
            f"utils config file '{cfg_path.relative_to(ROOT)}' (mapped from "
            f"[utils.config].{compiler}) does not exist. Copy "
            f"utils/config/TEMPLATE to create it."
        )
    openmp = 1 if omp else 0
    cmd = (
        module_prefix(modules_for(machine, "utils"))
        + f"python3 config.py config/{shlex.quote(cfg)} && "
        + f"make openmp={openmp} debug={debug} fesmutils-static"
    )
    run(cmd, ROOT / "utils")


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

    done = []
    for component in components:
        for omp in omp_variants:
            label = f"{component}-{'omp' if omp else 'serial'}"
            print(f"\n=== Building {label} "
                  f"[{args.machine}/{args.compiler}] ===")
            BUILDERS[component](machine, args.compiler, omp, args.debug)
            done.append(label)

    if args.check:
        print(
            f"\n\033[1m✓ config valid:\033[0m "
            f"{args.machine}/{args.compiler} resolves for "
            + ", ".join(done)
            + " (nothing was built)."
        )
    elif DRY_RUN:
        print("\n\033[1mDry run:\033[0m " + ", ".join(done) + " (nothing was built).")
    else:
        print("\n\033[1mBuilt:\033[0m " + ", ".join(done))


if __name__ == "__main__":
    main()
