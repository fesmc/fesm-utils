#!/bin/bash

SRCDIR=$PWD

NOW=$SRCDIR/fftw-omp
if test -d $NOW; then
  CHECK="Path exists"
else
  CHECK="Path does not exist"
fi
echo "$CHECK: $NOW"

NOW=$SRCDIR/fftw-serial
if test -d $NOW; then
  CHECK="Path exists"
else
  CHECK="Path does not exist"
fi
echo "$CHECK: $NOW"

NOW=$SRCDIR/lis-omp
if test -d $NOW; then
  CHECK="Path exists"
else
  CHECK="Path does not exist"
fi
echo "$CHECK: $NOW"

NOW=$SRCDIR/lis-serial
if test -d $NOW; then
  CHECK="Path exists"
else
  CHECK="Path does not exist"
fi
echo "$CHECK: $NOW"

