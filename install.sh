#!/bin/bash

# Make sure we stop on errors
set -e 

# Set compiler options
if [[ $1 = "ifx" ]]; then
    COMPILER_OPTS="CC=icx F77=ifx 'FFLAGS=-Ofast -march=core-avx2 -mtune=core-avx2 -traceback' 'CFLAGS=-Ofast -march=core-avx2 -mtune=core-avx2 -traceback'"
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
