#!/bin/bash

## This script will work on the DKRZ levante system, and 
## the AWI albedo system. 

# Make sure we stop on errors
set -e 

# Set compiler options
if [[ $1 = "ifx" ]]; then
    #COMPILER_OPTS="F77=ifx 'FFLAGS=-Ofast -march=core-avx2 -mtune=core-avx2 -traceback' 'CFLAGS=-Ofast -march=core-avx2 -mtune=core-avx2 -traceback'"
    COMPILER_OPTS="CC=icx F77=ifx ac_cv_c_libs='-L/sw/spack-levante/intel-oneapi-compilers-2023.2.1-kfv7xx/compiler/2023.2.1/linux/compiler/lib/intel64_lin -L/sw/spack-levante/intel-oneapi-compilers-2023.2.1-kfv7xx/compiler/2023.2.1/linux/bin-llvm/../lib -L/sw/spack-levante/gcc-11.2.0-bcn7mb/lib/gcc/x86_64-pc-linux-gnu/11.2.0 -L/sw/spack-levante/gcc-11.2.0-bcn7mb/lib/gcc/x86_64-pc-linux-gnu/11.2.0/../../../../lib64 -L/lib/../lib64 -L/usr/lib/../lib64 -L/sw/spack-levante/gcc-11.2.0-bcn7mb/lib/gcc/x86_64-pc-linux-gnu/11.2.0/../../.. -L/lib -L/usr/lib -L/sw/spack-levante/intel-oneapi-compilers-2023.2.1-kfv7xx/compiler/2023.2.1/linux/lib -lsvml -lirng -limf -lm -lgcc_s -lirc -ldl -lirc_s'"
elif [[ $1 = "ifort" ]]; then
    COMPILER_OPTS="CC=icc F77=ifort 'FFLAGS=-Ofast -march=core-avx2 -mtune=core-avx2 -traceback' 'CFLAGS=-Ofast -march=core-avx2 -mtune=core-avx2 -traceback'"
elif [[ $1 = "gfortran" ]]; then
    COMPILER_OPTS=""
else
    echo "Compiler not recognized: $1"
    exit 1
fi

echo "COMPILER_OPTS = $COMPILER_OPTS"
echo ""

### FFTW ###

SRCDIR=$PWD
FFTWSRC=fftw-3.3.10

# with omp enabled
cd $FFTWSRC
eval "./configure --disable-doc --prefix=$SRCDIR/fftw-omp --enable-openmp $COMPILER_OPTS"
make clean
make
make install
cd $SRCDIR

# serial (without omp) 
cd $FFTWSRC
eval "./configure --disable-doc --prefix=$SRCDIR/fftw-serial $COMPILER_OPTS"
make clean
make
make install
cd $SRCDIR

### LIS ###

SRCDIR=$PWD
LISSRC=lis-2.1.6

# with omp enabled
cd $LISSRC
eval "./configure --prefix=$SRCDIR/lis-omp --enable-f90 --enable-omp $COMPILER_OPTS"
make clean
make
make install
cd $SRCDIR

# serial (without omp) 
cd $LISSRC
eval "./configure --prefix=$SRCDIR/lis-serial --enable-f90 $COMPILER_OPTS" 
make clean
make
make install
cd $SRCDIR

echo ""
echo "" 
./check.sh
echo ""
