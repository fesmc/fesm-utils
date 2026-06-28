.SUFFIXES: .f .F .F90 .f90 .o .mod
.SHELL: /bin/sh

# Command-line options at make call
debug  ?= 0
openmp ?= 0 

# PATH options
srcdir = src
testdir = test
bindir = bin
libdir = libs
objdir_serial = include-serial
objdir_omp = include-omp

objdir = $(objdir_serial)
ifeq ($(openmp), 1)
    objdir = $(objdir_omp)
endif

## COMPILER CONFIGURATION ##
# Loaded from a config/<host> fragment by `python config.py config/<host>`,
# or assembled by `configme`. The fragment defines the compiler (FC,
# FFLAGS, FFLAGS_OPENMP, DFLAGS_*), netCDF wiring (INC_NC, LIB_NC), and
# any host-specific overrides (e.g. LFLAGS_EXTRA).

<COMPILER_CONFIGURATION>

## SHARED DEPENDENCY WIRING ##
# Pulled in after the compiler/machine fragment so it can reference FFLAGS,
# FFLAGS_OPENMP, LIB_NC, and LFLAGS_EXTRA set there. Applies the openmp=
# build switch and sets LFLAGS.

include config/common.mk

# Determine whether to use normal flags or debugging flags
DFLAGS   = $(DFLAGS_NODEBUG)
ifeq ($(debug), 1)
	DFLAGS   = $(DFLAGS_DEBUG)
endif

# Debugging flags with profiling output enabled
ifeq ($(debug), 2)
	DFLAGS   = $(DFLAGS_PROFILE)
endif

###############################################
##
## List of rules and source files
##
###############################################

$(objdir)/derivatives.o: $(srcdir)/derivatives.f90 $(objdir)/precision.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/distances.o: $(srcdir)/distances.f90 $(objdir)/precision.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/grid_to_cdo.o: $(srcdir)/coords/grid_to_cdo.f90 $(objdir)/precision.o $(objdir)/constants.o $(objdir)/coordinates.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/gaussian_filter.o: $(srcdir)/coords/gaussian_filter.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/gaussian_quadrature.o: $(srcdir)/gaussian_quadrature.f90 $(objdir)/precision.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/index.o: $(srcdir)/coords/index.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/interp1D.o: $(srcdir)/coords/interp1D.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/interp2D.o: $(srcdir)/coords/interp2D.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/mapping_scrip.o: $(srcdir)/coords/mapping_scrip.f90 $(objdir)/ncio.o $(objdir)/interp2D.o \
								$(objdir)/gaussian_filter.o $(objdir)/index.o $(objdir)/grid_to_cdo.o \
								$(objdir)/weight_map.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/geodesic.o: $(srcdir)/coords/geodesic.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/planet.o: $(srcdir)/coords/planet.f90 $(objdir)/precision.o $(objdir)/constants.o $(objdir)/geodesic.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/projection_oblimap2.o: $(srcdir)/coords/projection_oblimap2.f90 $(objdir)/precision.o $(objdir)/constants.o $(objdir)/planet.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/coordinates.o: $(srcdir)/coords/coordinates.f90 $(objdir)/precision.o $(objdir)/constants.o \
								$(objdir)/projection_oblimap2.o $(objdir)/planet.o $(objdir)/ncio.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_NC) -c -o $@ $<

$(objdir)/gaussian_latitudes.o: $(srcdir)/coords/gaussian_latitudes.f90 $(objdir)/precision.o $(objdir)/constants.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/kdtree.o: $(srcdir)/coords/kdtree.f90 $(objdir)/precision.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/polygons.o: $(srcdir)/coords/polygons.f90 $(objdir)/precision.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/weight_map.o: $(srcdir)/coords/weight_map.f90 $(objdir)/precision.o $(objdir)/constants.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/mapping.o: $(srcdir)/coords/mapping.f90 $(objdir)/precision.o $(objdir)/constants.o \
								$(objdir)/coordinates.o $(objdir)/planet.o $(objdir)/kdtree.o $(objdir)/weight_map.o \
								$(objdir)/mapping_scrip.o $(objdir)/conservative.o $(objdir)/interp2D.o $(objdir)/gaussian_filter.o \
								$(objdir)/map_io.o $(objdir)/grid_to_cdo.o $(objdir)/ncio.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_NC) -c -o $@ $<

$(objdir)/map_io.o: $(srcdir)/coords/map_io.f90 $(objdir)/precision.o $(objdir)/weight_map.o $(objdir)/ncio.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_NC) -c -o $@ $<

$(objdir)/ncio_interp.o: $(srcdir)/coords/ncio_interp.f90 $(objdir)/precision.o $(objdir)/constants.o \
								$(objdir)/ncio.o $(objdir)/mapping.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_NC) -c -o $@ $<

$(objdir)/conservative.o: $(srcdir)/coords/conservative.f90 $(objdir)/precision.o $(objdir)/constants.o \
								$(objdir)/coordinates.o $(objdir)/projection_oblimap2.o $(objdir)/planet.o \
								$(objdir)/kdtree.o $(objdir)/weight_map.o $(objdir)/polygons.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/coords.o: $(srcdir)/coords/coords.f90 $(objdir)/precision.o $(objdir)/constants.o $(objdir)/geodesic.o \
								$(objdir)/planet.o $(objdir)/projection_oblimap2.o $(objdir)/coordinates.o \
								$(objdir)/gaussian_latitudes.o $(objdir)/kdtree.o $(objdir)/weight_map.o \
								$(objdir)/mapping.o $(objdir)/map_io.o $(objdir)/ncio_interp.o $(objdir)/conservative.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/ncio.o: $(srcdir)/ncio.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_NC) -c -o $@ $<

$(objdir)/ncio_transpose.o: $(srcdir)/ncio_transpose.f90 $(objdir)/ncio.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_NC) -c -o $@ $<

$(objdir)/nml.o: $(srcdir)/nml.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/precision.o: $(srcdir)/precision.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/constants.o: $(srcdir)/constants.f90 $(objdir)/precision.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/staggering.o: $(srcdir)/staggering.f90 $(objdir)/precision.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/subgrid.o: $(srcdir)/subgrid.f90 $(objdir)/precision.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/timer.o: $(srcdir)/timer.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/timeout.o: $(srcdir)/timeout.f90 $(objdir)/nml.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/timestepping.o: $(srcdir)/timestepping.f90 $(objdir)/nml.o $(objdir)/ncio.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/variable_io.o: $(srcdir)/variable_io.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/varslice.o: $(srcdir)/varslice.f90 $(objdir)/precision.o $(objdir)/nml.o $(objdir)/ncio.o \
					$(objdir)/mapping.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

fesmutils_libs = 	$(objdir)/constants.o \
					$(objdir)/geodesic.o \
					$(objdir)/planet.o \
					$(objdir)/projection_oblimap2.o \
					$(objdir)/coordinates.o \
					$(objdir)/gaussian_latitudes.o \
					$(objdir)/kdtree.o \
					$(objdir)/polygons.o \
					$(objdir)/weight_map.o \
					$(objdir)/mapping.o \
					$(objdir)/map_io.o \
					$(objdir)/ncio_interp.o \
					$(objdir)/conservative.o \
					$(objdir)/coords.o \
					$(objdir)/derivatives.o \
					$(objdir)/distances.o \
					$(objdir)/gaussian_filter.o \
					$(objdir)/gaussian_quadrature.o \
					$(objdir)/grid_to_cdo.o \
					$(objdir)/index.o \
					$(objdir)/interp1D.o \
					$(objdir)/interp2D.o \
					$(objdir)/mapping_scrip.o \
					$(objdir)/ncio.o \
					$(objdir)/ncio_transpose.o \
					$(objdir)/nml.o \
					$(objdir)/precision.o \
					$(objdir)/subgrid.o \
					$(objdir)/staggering.o \
					$(objdir)/timer.o \
					$(objdir)/timeout.o \
					$(objdir)/timestepping.o \
					$(objdir)/variable_io.o \
					$(objdir)/varslice.o

###############################################
##
## Compilation of complete programs
##
###############################################

libs : $(fesmutils_libs)
		@echo " "
		@echo "    fesm-utils is ready."
		@echo " "

# Static library compilation
fesmutils-static: $(fesmutils_libs)
	ar rc $(objdir)/libfesmutils.a $(fesmutils_libs)
	ranlib $(objdir)/libfesmutils.a
	@echo " "
	@echo "    $(objdir)/libfesmutils.a is ready."
	@echo " "

# coords library test programs
test-coords : fesmutils-static
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_geodinverse.x $(testdir)/coords/test_geodinverse.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_proj.x $(testdir)/coords/test_proj.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_gauss.x $(testdir)/coords/test_gauss.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_rotated.x $(testdir)/coords/test_rotated.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_kdtree.x $(testdir)/coords/test_kdtree.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_weightmap.x $(testdir)/coords/test_weightmap.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_mapinit.x $(testdir)/coords/test_mapinit.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_bilinear.x $(testdir)/coords/test_bilinear.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_ccsm3_grl.x $(testdir)/coords/test_ccsm3_grl.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_mapio.x $(testdir)/coords/test_mapio.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_polygons.x $(testdir)/coords/test_polygons.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_conservative.x $(testdir)/coords/test_conservative.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_conservative_xsys.x $(testdir)/coords/test_conservative_xsys.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_conservative_sphere.x $(testdir)/coords/test_conservative_sphere.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_map_init_con.x $(testdir)/coords/test_map_init_con.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_map_cache.x $(testdir)/coords/test_map_cache.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_map_cdo.x $(testdir)/coords/test_map_cdo.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_varslice_map.x $(testdir)/coords/test_varslice_map.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_ncio_interp.x $(testdir)/coords/test_ncio_interp.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_parity.x $(testdir)/coords/test_parity.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		@echo " "
		@echo "    coords test programs ready in $(bindir)/."
		@echo " "

bench-coords : fesmutils-static
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/bench_conservative.x $(testdir)/coords/bench_conservative.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		@echo " "
		@echo "    coords benchmark ready: $(bindir)/bench_conservative.x (needs cdo on PATH)."
		@echo " "

.PHONY : usage
usage:
	@echo ""
	@echo "    * USAGE * "
	@echo ""
	@echo " make test-coords : compiles the coords test programs"
	@echo " make clean       : cleans object files"
	@echo ""

clean:
	rm -f $(bindir)/*.x
	rm -f  *.x gmon.out $(objdir_serial)/*.o $(objdir_serial)/*.mod $(objdir_serial)/*.a $(objdir_serial)/*.so
	rm -f  *.x gmon.out $(objdir_omp)/*.o $(objdir_omp)/*.mod $(objdir_omp)/*.a $(objdir_omp)/*.so
	rm -rf *.x.dSYM
