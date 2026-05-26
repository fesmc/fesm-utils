# Shared build configuration for fesm-utils/utils (dependency wiring).
#
# Loaded after the compiler and machine fragments (configme assembles them in
# the order: compiler -> machine -> netCDF -> common). utils depends only on
# netCDF; it references FFLAGS / FFLAGS_OPENMP (compiler) and LIB_NC (machine or
# auto-detected netCDF).

# OpenMP build (make openmp=1): append the compiler's OpenMP flag. (The objdir
# serial/omp swap is handled by the template.)
ifeq ($(openmp), 1)
    FFLAGS += $(FFLAGS_OPENMP)
endif

LFLAGS_EXTRA ?=
LFLAGS = $(LIB_NC) $(LFLAGS_EXTRA)
