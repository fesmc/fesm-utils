# climber-x-exlib

Convenience repository holding external libraries needed to run CLIMBER-X.

## Configure and compile each library

To compile with your default compiler (most likely gcc/gfortran) run:

```bash
./install.sh
```

To compile with an `ifx` or `ifort` compiler, run with argument `ifx` or `ifort`:

```bash
./install.sh ifx
./install.sh ifort
```

Specifically to deal with an issue on the PIK HPC2024 (Foote) cluster,
run the PIK specific script:

```bash
./install_pik.sh ifx
```

Running the `install.sh` script will compile and "install" the following
library versions:

```bash
exlib/fftw-omp
exlib/fftw-serial
exlib/lis-omp
exlib/lis-serial
```

See the internals of `install.sh` if you would like to customize any
installation options further.

## Use the libraries

Now a symlink can be made to these libraries for use within CLIMBER-X:

```bash
cd climber-x/utils
ln -s /path/to/climber-x-exlib/exlib
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
