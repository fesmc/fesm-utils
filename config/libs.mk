# config/libs.mk — machine/compiler specifics for building the vendored C/Fortran
# libraries FFTW, LIS and SHTns via their own autotools (configure/make/install).
#
# The root Makefile's `fftw`, `lis` and `shtns` targets consume the variables
# set here. Selection is by $(MACHINE) and $(COMPILER), which `configme` writes
# into the generated Makefile. Override on the command line for a one-off, e.g.
#     make fftw openmp=1 MACHINE=dkrz_levante COMPILER=ifx
#
# This file (plus the target recipes) replaces the old machines/*.toml + build.py:
# it is the single in-repo home for the autotools build quirks, sitting next to
# the Fortran build config (the configme-generated compiler block + common.mk).
#
# The library builds assume the compiler/netCDF modules are already loaded in
# your shell — the same assumption the utils build makes for nf-config/nc-config.
# On clusters, `module load` the toolchain before running `make fftw|lis|shtns`.

# ---- Variant: install-prefix suffix + each package's OpenMP configure flag ----
LIB_VARIANT   = serial
FFTW_OMP_OPT  =
LIS_OMP_OPT   =
SHTNS_OMP_OPT =
ifeq ($(openmp),1)
    LIB_VARIANT   = omp
    FFTW_OMP_OPT  = --enable-openmp
    LIS_OMP_OPT   = --enable-omp
    SHTNS_OMP_OPT = --enable-openmp
endif

# ---- Defaults: rely on autotools autodetection (generic linux / gfortran) -----
# Empty => the variable is omitted from ./configure, so autotools uses its own
# default (finds gcc/gfortran on PATH). Machine blocks below override as needed.
CC_LIB     ?=
F77_LIB    ?=
CFLAGS_LIB ?=
FFLAGS_LIB ?=
MARCH_LIB  ?=
SHTNS_CFLAGS ?=
INTEL_RUNTIME_LIBS ?= 0

# autoconf >= 2.70 probes for the newest C standard the compiler supports (C23,
# then C11, then C99) and appends the matching -std flag to CC. FFTW/LIS are
# pre-C23 C: they use old-style (void (*)()) casts that are then called with
# arguments, which C23 turns from a warning into a hard error. Pre-seed the
# probe's cache vars to "no" so no -std flag is added and the compiler keeps its
# own (pre-C23) default. Unknown to configure scripts that predate the probe
# (e.g. FFTW's), so this is always safe. Applied on every machine.
CONF_CSTD = ac_cv_prog_cc_c23=no ac_cv_prog_cc_c11=no ac_cv_prog_cc_c99=no

# =============================================================================
#  Per-machine / per-compiler quirks
# =============================================================================

# ---- macbook (local macOS, gfortran) ----------------------------------------
# Apple Clang (the default `cc`) has no OpenMP, so FFTW/LIS --enable-openmp fails
# under it. Pick the newest Homebrew gcc-N and pin it to -std=gnu17: GCC 15
# defaults to gnu23 (empty () means "no params"), which rejects FFTW/LIS's
# old-style casts. Embedding the flag in CC (not CFLAGS) keeps the packages' own
# CFLAGS from clobbering it. If no Homebrew gcc is found, fall back to autotools
# default (a serial build still works; an OpenMP build will fail — brew install gcc).
ifeq ($(MACHINE),macbook)
  _BREW_GCC := $(shell ls /opt/homebrew/bin/gcc-[0-9]* /usr/local/bin/gcc-[0-9]* 2>/dev/null \
                 | sed -E 's;.*/;;' | sort -t- -k2 -n | tail -1)
  ifneq ($(_BREW_GCC),)
    CC_LIB := $(_BREW_GCC) -std=gnu17
  endif
endif

# ---- dkrz_levante (AMD EPYC 7763 / Zen3: AVX2, no AVX-512) -------------------
ifeq ($(MACHINE),dkrz_levante)
  MARCH_LIB := znver3
  ifeq ($(COMPILER),ifx)
    CC_LIB := icx
    F77_LIB := ifx
    INTEL_RUNTIME_LIBS := 1
  endif
  ifeq ($(COMPILER),ifort)
    CC_LIB := icc
    F77_LIB := ifort
  endif
endif

# ---- awi_albedo (AMD EPYC 7702 / Zen2: AVX2, no AVX-512) ---------------------
ifeq ($(MACHINE),awi_albedo)
  MARCH_LIB := znver2
  ifeq ($(COMPILER),ifx)
    CC_LIB := icx
    F77_LIB := ifx
    INTEL_RUNTIME_LIBS := 1
  endif
  ifeq ($(COMPILER),ifort)
    CC_LIB := icc
    F77_LIB := ifort
  endif
endif

# ---- pik_hpc2024 -------------------------------------------------------------
# NB: lis must be built with an OLDER Intel toolchain than fftw/utils (oneAPI
# 2024.0 errors out on lis). Since fftw/lis are now separate `make` targets,
# `module swap` to intel/oneAPI/2023.2.0 before `make lis`, then swap back.
# Under 2023.2.0, the lis ifx build must drop CC (icx) — see LIS_CC_LIB below.
ifeq ($(MACHINE),pik_hpc2024)
  MARCH_LIB := core-avx2
  ifeq ($(COMPILER),ifx)
    CC_LIB := icx
    F77_LIB := ifx
  endif
  ifeq ($(COMPILER),ifort)
    CC_LIB := icc
    F77_LIB := ifort
  endif
  ifeq ($(COMPILER),gfortran)
    CC_LIB := gcc
    F77_LIB := gfortran
  endif
endif

# ---- generic linux (no HPC file applies) ------------------------------------
ifeq ($(MACHINE),linux)
  MARCH_LIB := core-avx2
  ifeq ($(COMPILER),ifx)
    CC_LIB := icx
    F77_LIB := ifx
  endif
  ifeq ($(COMPILER),ifort)
    CC_LIB := icc
    F77_LIB := ifort
  endif
endif

# =============================================================================
#  Derived flags — composed from the machine/compiler selections above
# =============================================================================

# Intel builds carry an -Ofast/-march optimization string (and -traceback, which
# SHTns keeps after opt-flag stripping). core-avx2 targets also add -mtune.
ifneq ($(filter $(COMPILER),ifx ifort),)
  ifneq ($(MARCH_LIB),)
    ifeq ($(MARCH_LIB),core-avx2)
      CFLAGS_LIB := -Ofast -march=core-avx2 -mtune=core-avx2 -traceback
      FFLAGS_LIB := -Ofast -march=core-avx2 -mtune=core-avx2 -traceback
    else
      CFLAGS_LIB := -Ofast -march=$(MARCH_LIB) -traceback
      FFLAGS_LIB := -Ofast -march=$(MARCH_LIB) -traceback
    endif
  endif
  # SHTns owns its own -O/-fp-model scheme, so its target strips the opt/arch
  # flags and passes --enable-march=$(MARCH_LIB) instead; -traceback is the only
  # non-opt residual worth forwarding.
  SHTNS_CFLAGS := -traceback
endif

# lis compiler override (default: same as the rest of the build). pik/ifx must
# drop CC (icx) for the lis build under the 2023.2.0 toolchain.
LIS_CC_LIB ?= $(CC_LIB)
ifeq ($(MACHINE)/$(COMPILER),pik_hpc2024/ifx)
  LIS_CC_LIB :=
endif

# Intel/ifx-built fftw/lis fail at the link stage with missing Intel runtime
# symbols (svml, irc, imf, ...) unless the Intel runtime lib dir is on the
# linker path. Derive it from the compiler itself (machine-independent) and pass
# a matching ac_cv_c_libs to configure.
CONF_INTEL_LIBS =
ifeq ($(INTEL_RUNTIME_LIBS),1)
  _ONEAPI_LIBDIR := $(shell dirname "$$($(firstword $(CC_LIB)) -print-file-name=libimf.so)" 2>/dev/null)
  CONF_INTEL_LIBS = ac_cv_c_libs="-L$(_ONEAPI_LIBDIR) -lsvml -lirng -limf -lirc -ldl -lm"
endif

# ---- Assembled ./configure variable strings ---------------------------------
# _confvars(cc): CC/F77/CFLAGS/FFLAGS (each omitted when empty) + C-std pins +
# the Intel runtime-lib fix. Single-quoted so multi-word values survive the shell.
_confvars = $(if $(1),CC='$(1)') $(if $(F77_LIB),F77='$(F77_LIB)') \
            $(if $(CFLAGS_LIB),CFLAGS='$(CFLAGS_LIB)') $(if $(FFLAGS_LIB),FFLAGS='$(FFLAGS_LIB)') \
            $(CONF_CSTD) $(CONF_INTEL_LIBS)

FFTW_CONF_VARS = $(call _confvars,$(CC_LIB))
LIS_CONF_VARS  = $(call _confvars,$(LIS_CC_LIB))

# SHTns uses AC_PROG_FC (FC, not F77), builds its CPU kernels with a separate
# "kernel compiler" (--enable-kernel-compiler, defaults to gcc regardless of CC),
# and pins its target arch via --enable-march rather than -march in CFLAGS.
SHTNS_CONF_VARS = $(if $(F77_LIB),FC='$(F77_LIB)') \
                  $(if $(CC_LIB),CC='$(CC_LIB)' --enable-kernel-compiler='$(CC_LIB)') \
                  $(if $(MARCH_LIB),--enable-march=$(MARCH_LIB)) \
                  $(if $(SHTNS_CFLAGS),CFLAGS='$(SHTNS_CFLAGS)') \
                  $(CONF_CSTD)
