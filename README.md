# fesm-utils

Convenience repository holding external libraries needed to run fast Earth system models, like [CLIMBER-X](https://github.com/cxesmc/climber-x).

## Directory structure

Currently, `fesm-utils` manages installation of `fftw`, `lis` and a collection of useful modules `modules`.

## Configure and compile `lis` and `fftw`

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
