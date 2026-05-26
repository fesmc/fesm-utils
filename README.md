# fesm-utils

Convenience repository holding external libraries needed to run fast Earth system models, like [CLIMBER-X](https://github.com/cxesmc/climber-x).

## Directory structure

`fesm-utils` manages the build of three components:

- `fftw` — FFTW, built with autotools.
- `lis` — the LIS solver library, built with autotools.
- `utils` — a self-contained collection of Fortran modules compiled into a static library.

All three are built through a single driver, `build.py`, configured per machine via `machines/<name>.toml`.

## Building

`build.py` builds any combination of the three components, for any machine, with or without OpenMP. Each (component × variant) is installed into its own subfolder, so an OpenMP-enabled and a serial version coexist:

```bash
fftw-omp     fftw-serial
lis-omp      lis-serial
utils/include-omp     utils/include-serial
```

### Usage

```bash
./build.py --list-machines                                  # show known machines

./build.py -m dkrz_levante -c ifx                           # everything, omp + serial
./build.py -m macbook -c gfortran --component utils --variant serial
./build.py -m pik_hpc2024 -c ifx --component lis --variant omp
./build.py -m linux -c gfortran --dry-run                   # print commands, don't run
```

Options:

- `-m/--machine` — a file under `machines/` (e.g. `dkrz_levante`, `awi_albedo`, `pik_hpc2024`, `macbook`, `linux`).
- `-c/--compiler` — `ifx`, `ifort`, or `gfortran` (whichever the machine defines).
- `--component` — `fftw`, `lis`, `utils`, or `all` (default `all`).
- `--variant` — `omp`, `serial`, or `both` (default `both`).
- `--debug` — `utils` debug level: `0` normal, `1` debug, `2` profile (default `0`).
- `--dry-run` — print the commands that would run, without executing.
- `--check` — validate a machine config (compilers defined, `utils` config files present) and print the resolved commands, without building. Run this first when setting up a new machine.

### Adding a machine

Everything machine-specific lives in **config files**, never in `build.py` — you should not need to edit the driver to support a new machine. A machine is described by:

1. `machines/<name>.toml` — compilers, flags, and (on clusters) modules.
2. *(only if you build `utils`)* one `utils/config/<name>_<compiler>` file per compiler, describing the Fortran compiler and NetCDF paths for the `utils` library.

Two annotated skeletons are provided to copy: `machines/TEMPLATE.toml` and `utils/config/TEMPLATE`. Every field is documented inline in those files.

#### Step by step

Suppose your machine is `my_cluster` and you want `gfortran`.

1. **Copy the machine template** and rename it:

   ```bash
   cp machines/TEMPLATE.toml machines/my_cluster.toml
   ```

   Edit it: set `[machine].name`, keep a `[compilers.gfortran]` table (an empty table means "use the compiler's autotools defaults"), and — on an HPC cluster — list the modules to load under `[modules].load`. Find module names with `module avail` or your cluster's documentation.

2. **If you build the `utils` component**, copy the utils template and rename it to match:

   ```bash
   cp utils/config/TEMPLATE utils/config/my_cluster_gfortran
   ```

   Then map it in `machines/my_cluster.toml`:

   ```toml
   [utils.config]
   gfortran = "my_cluster_gfortran"
   ```

   In the `utils` config file, the main thing to get right is the NetCDF location. After loading any NetCDF module, discover the paths with:

   ```bash
   nf-config --prefix     # NetCDF-Fortran install root  -> NETCDFFI_ROOT
   nc-config --prefix     # NetCDF-C install root         -> NETCDFC_ROOT
   ```

   (If you don't build `utils` on this machine, skip this step and delete the `[utils.config]` section from the `.toml`.)

3. **Validate before building** — this catches typos, undefined compilers, and missing `utils` config files, and prints the exact commands that would run:

   ```bash
   ./build.py -m my_cluster -c gfortran --check
   ```

4. **Build for real:**

   ```bash
   ./build.py -m my_cluster -c gfortran
   ```

#### What the `[compilers.<name>]` fields mean

Each `[compilers.<name>]` table is a set of autotools configure variables appended verbatim to the `fftw`/`lis` `./configure` call (the same role as the `COMPILER_OPTS` string in the old `install_*.sh` scripts). The common ones:

| field          | meaning                          | example                  |
| -------------- | -------------------------------- | ------------------------ |
| `CC`           | C compiler command               | `gcc`, `icx`, `icc`      |
| `F77`          | Fortran compiler command         | `gfortran`, `ifx`, `ifort` |
| `FFLAGS`       | Fortran compile flags            | `-Ofast -march=core-avx2` |
| `CFLAGS`       | C compile flags                  | `-Ofast -march=core-avx2` |
| `ac_cv_c_libs` | explicit C runtime link paths    | *(advanced, see below)*  |

#### Per-component overrides (advanced)

If one component needs a different toolchain than the rest, override the modules and/or compiler table for that component only — the override **replaces** the machine-level table for that component. `machines/pik_hpc2024.toml` does this for `lis`, which must build with an older Intel module than `fftw`/`utils`.

#### Troubleshooting

| symptom                                                                 | fix                                                                                                                                                                                              |
| ----------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `unknown machine '...'`                                                 | the `.toml` file name must match the `-m` argument. Run `./build.py --list-machines`.                                                                                                            |
| `compiler '...' not defined for machine`                                | add (or fix the name of) a `[compilers.<name>]` table in the `.toml`.                                                                                                                            |
| `no utils config mapped` / `utils config file ... does not exist`       | add the `[utils.config]` mapping and create the `utils/config/<name>_<compiler>` file (copy `utils/config/TEMPLATE`).                                                                            |
| `utils` build can't find `netcdf.mod` or `-lnetcdff`                    | fix `NETCDFFI_ROOT` / `NETCDFC_ROOT` (and the module load) in the `utils` config file. Verify with `nf-config --prefix`.                                                                          |
| `ifx` build of `fftw`/`lis` fails at the **link** stage with missing Intel runtime symbols (`svml`, `irc`, `imf`, …) | set `ac_cv_c_libs` in the `[compilers.ifx]` table to the explicit Intel/GCC runtime library paths for your cluster. See `machines/dkrz_levante.toml` for a working example — but note those paths are Levante-specific and must be replaced with your cluster's equivalents. |

### Building `utils` by hand

The `utils` static library can also be built directly with its own Makefile, independently of `build.py`:

```bash
cd utils
python config.py config/dkrz_levante_ifx   # replace with config file for your system
make clean
make fesmutils-static
```

## Use the libraries

A symlink can be made to these libraries for use within, e.g., CLIMBER-X:

```bash
cd climber-x/utils/
ln -s /path/to/fesm-utils ./
```

## To get the original libraries

This step is mainly only relevant for the maintainers of this repository,
or in the case that a user wants to try out different versions.

```bash
### FFTW
wget https://www.fftw.org/fftw-3.3.10.tar.gz
tar -xvf fftw-3.3.10.tar.gz
rm fftw-3.3.10.tar.gz

### LIS
wget https://www.ssisc.org/lis/dl/lis-2.1.6.zip
unzip lis-2.1.6.zip
rm lis-2.1.6.zip
```

Then build with `build.py` as described above.
