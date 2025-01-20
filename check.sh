#!/bin/bash

SRCDIR=$PWD

NOW=$SRCDIR/exlib/fftw-omp
if test -d $NOW; then
  CHECK="Path exists"
else
  CHECK="Path does not exist"
fi
echo "$CHECK: $NOW"

NOW=$SRCDIR/exlib/fftw-serial
if test -d $NOW; then
  CHECK="Path exists"
else
  CHECK="Path does not exist"
fi
echo "$CHECK: $NOW"

NOW=$SRCDIR/exlib/lis-omp
if test -d $NOW; then
  CHECK="Path exists"
else
  CHECK="Path does not exist"
fi
echo "$CHECK: $NOW"

NOW=$SRCDIR/exlib/lis-serial
if test -d $NOW; then
  CHECK="Path exists"
else
  CHECK="Path does not exist"
fi
echo "$CHECK: $NOW"

