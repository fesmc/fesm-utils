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
./build.py -m generic -c gfortran --dry-run                 # print commands, don't run
```

Options:

- `-m/--machine` — a file under `machines/` (e.g. `dkrz_levante`, `awi_albedo`, `pik_hpc2024`, `macbook`, `generic`).
- `-c/--compiler` — `ifx`, `ifort`, or `gfortran` (whichever the machine defines).
- `--component` — `fftw`, `lis`, `utils`, or `all` (default `all`).
- `--variant` — `omp`, `serial`, or `both` (default `both`).
- `--debug` — `utils` debug level: `0` normal, `1` debug, `2` profile (default `0`).
- `--dry-run` — print the commands that would run, without executing.

### Adding a machine

Copy an existing file in `machines/` and edit the compiler flags, modules, and (for `utils`) the config-file mapping. Each `[compilers.<name>]` table is a set of autotools configure variables passed verbatim to the `fftw`/`lis` `./configure` call. See `machines/pik_hpc2024.toml` for an example of per-component module and compiler overrides.

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
