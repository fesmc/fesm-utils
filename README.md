# fesm-utils

Convenience repository holding external libraries needed to run fast Earth system models, like [CLIMBER-X](https://github.com/cxesmc/climber-x).

## Directory structure

Currently, `fesm-utils` manages installation of `fftw`, `lis` and a collection of useful modules `utils`.

The `utils` subdirectory is self-contained and the modules can be compiled into a static library with a Makefile in the directory. `lis` and `fftw` are installed using the install scripts in the main directory (see below).

## Unified build (recommended): `build.py`

`build.py` is a single entry point that builds any combination of the three
components (`fftw`, `lis`, `utils`), for any machine, with or without OpenMP.
Per-machine configuration lives in `machines/<name>.toml`.

```bash
./build.py --list-machines                                  # show known machines

./build.py -m dkrz_levante -c ifx                           # build everything, omp + serial
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

Each component keeps its native build: `fftw` and `lis` are configured/compiled
with autotools and installed into `fftw-{omp,serial}` / `lis-{omp,serial}`;
`utils` is built via its own `config.py` + `Makefile`. To add a new machine,
copy an existing `machines/*.toml` and edit the compiler/module settings.

The component-specific scripts below still work and are not affected by `build.py`.

## Configure and compile `lis` and `fftw` (legacy scripts)

To compile, run the install script and specify your compiler (currently `ifx`, `ifort` or `gfortran`):

```bash
./install.sh ifx
```

There are additional install scripts specific to certain HPC systems. If you are using one of these systems, these scripts should work:

```bash
./install_pik.sh ifx  # PIK HPC2024 (foote) cluster
./install_dkrz.sh ifx # DKRZ levante cluster
./install_awi.sh ifx  # AWI albedo cluster
```

Running the `install.sh` script will compile and "install" the following
library versions into separate subfolders:

```bash
fftw-omp
fftw-serial
lis-omp
lis-serial
```

These can then be linked to in any external program. See the internals of `install.sh` if you would like to customize any installation options further.

## Configure and compile `utils`

To make a static library of the `utils` modules, configure the Makefile for your system and compile.

```bash
cd utils
python config.py config/dkrz_levante_ifx  # replace with config file for your system
make clean
make fesmutils-static
```

## Use the libraries

Now a symlink can be made to these libraries for use within, e.g., CLIMBER-X:

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

Then, proceed with configure/compile instructions above.
