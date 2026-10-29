.SUFFIXES: .f .F .F90 .f90 .o .mod
.SHELL: /bin/sh

# Command-line options at make call
debug  ?= 0
openmp ?= 0 
pic    ?= 0

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
# FFLAGS_OPENMP, LIB_NC, and LFLAGS_EXTRA set there. Applies the openmp= /
# pic= build switches and sets LFLAGS.

include config/common.mk

## EXTERNAL LIBRARY BUILD CONFIG ##
# Machine/compiler specifics for the autotools builds of FFTW, LIS and SHTns
# (the `fftw`, `lis`, `shtns` targets below). Reads $(MACHINE)/$(COMPILER) set
# in the compiler block above. Kept separate from common.mk so the utils build
# does not carry the C-library quirks.

include config/libs.mk

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

$(objdir)/grid_cdo.o: $(srcdir)/coords/grid_cdo.f90 $(objdir)/precision.o $(objdir)/constants.o $(objdir)/coordinates.o \
								$(objdir)/fesm_log.o
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
								$(objdir)/gaussian_filter.o $(objdir)/index.o $(objdir)/grid_cdo.o \
								$(objdir)/weight_map.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/geodesic.o: $(srcdir)/coords/geodesic.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/planet.o: $(srcdir)/coords/planet.f90 $(objdir)/precision.o $(objdir)/constants.o $(objdir)/geodesic.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/projection_oblimap2.o: $(srcdir)/coords/projection_oblimap2.f90 $(objdir)/precision.o $(objdir)/constants.o $(objdir)/planet.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/coordinates.o: $(srcdir)/coords/coordinates.f90 $(objdir)/precision.o $(objdir)/constants.o \
								$(objdir)/projection_oblimap2.o $(objdir)/planet.o $(objdir)/ncio.o \
								$(objdir)/gaussian_latitudes.o $(objdir)/fesm_log.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_NC) -c -o $@ $<

$(objdir)/gaussian_latitudes.o: $(srcdir)/coords/gaussian_latitudes.f90 $(objdir)/precision.o $(objdir)/constants.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/kdtree.o: $(srcdir)/coords/kdtree.f90 $(objdir)/precision.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/polygons.o: $(srcdir)/coords/polygons.f90 $(objdir)/precision.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/weight_map.o: $(srcdir)/coords/weight_map.f90 $(objdir)/precision.o $(objdir)/constants.o \
								$(objdir)/fesm_log.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/grid_signature.o: $(srcdir)/coords/grid_signature.f90 $(objdir)/precision.o $(objdir)/coordinates.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/mapping.o: $(srcdir)/coords/mapping.f90 $(objdir)/precision.o $(objdir)/constants.o \
								$(objdir)/coordinates.o $(objdir)/planet.o $(objdir)/kdtree.o $(objdir)/weight_map.o \
								$(objdir)/mapping_scrip.o $(objdir)/conservative.o $(objdir)/interp2D.o $(objdir)/gaussian_filter.o \
								$(objdir)/map_io.o $(objdir)/grid_cdo.o $(objdir)/ncio.o $(objdir)/grid_signature.o \
								$(objdir)/fesm_log.o
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
								$(objdir)/mapping.o $(objdir)/map_io.o $(objdir)/ncio_interp.o $(objdir)/conservative.o \
								$(objdir)/grid_cdo.o $(objdir)/grid_signature.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/coupler.o: $(srcdir)/coupler.f90 $(objdir)/precision.o $(objdir)/coordinates.o \
								$(objdir)/mapping.o $(objdir)/coords.o $(objdir)/fesm_log.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/ncio.o: $(srcdir)/ncio.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_NC) -c -o $@ $<

# Regenerate src/ncio.f90 from the fypp template src/ncio.fypp. This is a
# developer step run only after editing the template; the generated ncio.f90
# is committed, so the normal build never needs fypp (or Python). fypp is
# vendored at tools/fypp.py.
.PHONY : ncio-gen
ncio-gen:
	python3 tools/fypp.py $(srcdir)/ncio.fypp $(srcdir)/ncio.f90
	@echo "Regenerated $(srcdir)/ncio.f90 from $(srcdir)/ncio.fypp"

$(objdir)/ncio_transpose.o: $(srcdir)/ncio_transpose.f90 $(objdir)/ncio.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_NC) -c -o $@ $<

$(objdir)/nml.o: $(srcdir)/nml.f90 $(objdir)/fesm_log.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/precision.o: $(srcdir)/precision.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/fesm_log.o: $(srcdir)/fesm_log.f90 $(objdir)/precision.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/fesm_fields.o: $(srcdir)/fesm_fields.f90 $(objdir)/precision.o $(objdir)/fesm_log.o \
					$(objdir)/ncio.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/fesm_units.o: $(srcdir)/fesm_units.f90 $(objdir)/precision.o $(objdir)/fesm_log.o \
					$(objdir)/fesm_fields.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/constants.o: $(srcdir)/constants.f90 $(objdir)/precision.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/root_finder.o: $(srcdir)/root_finder.f90 $(objdir)/precision.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/phys_constants.o: $(srcdir)/phys_constants.f90 $(objdir)/precision.o \
					$(objdir)/nml.o $(objdir)/ncio.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/staggering.o: $(srcdir)/staggering.f90 $(objdir)/precision.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/subgrid.o: $(srcdir)/subgrid.f90 $(objdir)/precision.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/multitopo.o: $(srcdir)/multitopo.f90 $(objdir)/precision.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/timer.o: $(srcdir)/timer.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/timeout.o: $(srcdir)/timeout.f90 $(objdir)/nml.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/timestepping.o: $(srcdir)/timestepping.f90 $(objdir)/nml.o $(objdir)/ncio.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/fesm_clock.o: $(srcdir)/fesm_clock.f90 $(objdir)/precision.o $(objdir)/fesm_log.o \
					$(objdir)/nml.o $(objdir)/ncio.o $(objdir)/phys_constants.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/variable_io.o: $(srcdir)/variable_io.f90 $(objdir)/ncio.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/varslice.o: $(srcdir)/varslice.f90 $(objdir)/precision.o $(objdir)/constants.o $(objdir)/nml.o $(objdir)/ncio.o \
					$(objdir)/mapping.o $(objdir)/interp2D.o $(objdir)/fesmdata.o $(objdir)/fesm_log.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/fesmdata.o: $(srcdir)/fesmdata.f90 $(objdir)/precision.o $(objdir)/constants.o $(objdir)/nml.o \
					$(objdir)/ncio.o $(objdir)/coordinates.o $(objdir)/grid_cdo.o $(objdir)/mapping.o \
					$(objdir)/fesm_log.o $(objdir)/fesm_units.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/regions.o: $(srcdir)/regions.f90 $(objdir)/precision.o $(objdir)/nml.o $(objdir)/ncio.o \
					$(objdir)/coordinates.o $(objdir)/mapping.o $(objdir)/fesmdata.o $(objdir)/fesm_log.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/topodata.o: $(srcdir)/topodata.f90 $(objdir)/precision.o $(objdir)/nml.o $(objdir)/ncio.o \
					$(objdir)/coordinates.o $(objdir)/mapping.o $(objdir)/fesmdata.o $(objdir)/fesm_log.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/topoprep.o: $(srcdir)/topoprep.f90 $(objdir)/precision.o $(objdir)/constants.o $(objdir)/fesm_log.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/htopo.o: $(srcdir)/htopo.f90 $(objdir)/precision.o $(objdir)/constants.o $(objdir)/nml.o \
					$(objdir)/ncio.o $(objdir)/coordinates.o $(objdir)/grid_cdo.o $(objdir)/interp2D.o \
					$(objdir)/phys_constants.o $(objdir)/topodata.o $(objdir)/regions.o $(objdir)/fesm_log.o \
					$(objdir)/topoprep.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/series.o: $(srcdir)/series.f90 $(objdir)/precision.o $(objdir)/constants.o $(objdir)/nml.o $(objdir)/ncio.o \
					$(objdir)/interp1D.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/tsgen.o: $(srcdir)/tsgen.f90 $(objdir)/precision.o $(objdir)/constants.o $(objdir)/nml.o $(objdir)/ncio.o \
					$(objdir)/series.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/zremap.o: $(srcdir)/zremap.f90 $(objdir)/precision.o $(objdir)/constants.o $(objdir)/nml.o $(objdir)/ncio.o \
					$(objdir)/fesm_log.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/fesm_datasrc.o: $(srcdir)/fesm_datasrc.f90 $(objdir)/precision.o $(objdir)/constants.o \
					$(objdir)/fesm_log.o $(objdir)/fesm_fields.o $(objdir)/fesm_units.o \
					$(objdir)/nml.o $(objdir)/ncio.o $(objdir)/coords.o $(objdir)/varslice.o $(objdir)/series.o
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
					$(objdir)/grid_signature.o \
					$(objdir)/mapping.o \
					$(objdir)/map_io.o \
					$(objdir)/ncio_interp.o \
					$(objdir)/conservative.o \
					$(objdir)/coords.o \
					$(objdir)/coupler.o \
					$(objdir)/derivatives.o \
					$(objdir)/distances.o \
					$(objdir)/gaussian_filter.o \
					$(objdir)/gaussian_quadrature.o \
					$(objdir)/grid_cdo.o \
					$(objdir)/index.o \
					$(objdir)/interp1D.o \
					$(objdir)/interp2D.o \
					$(objdir)/mapping_scrip.o \
					$(objdir)/ncio.o \
					$(objdir)/ncio_transpose.o \
					$(objdir)/nml.o \
					$(objdir)/precision.o \
					$(objdir)/fesm_log.o \
					$(objdir)/fesm_fields.o \
					$(objdir)/fesm_units.o \
					$(objdir)/phys_constants.o \
					$(objdir)/root_finder.o \
					$(objdir)/subgrid.o \
					$(objdir)/multitopo.o \
					$(objdir)/staggering.o \
					$(objdir)/timer.o \
					$(objdir)/timeout.o \
					$(objdir)/timestepping.o \
					$(objdir)/fesm_clock.o \
					$(objdir)/variable_io.o \
					$(objdir)/varslice.o \
					$(objdir)/fesmdata.o \
					$(objdir)/regions.o \
					$(objdir)/topodata.o \
					$(objdir)/topoprep.o \
					$(objdir)/htopo.o \
					$(objdir)/series.o \
					$(objdir)/tsgen.o \
					$(objdir)/zremap.o \
					$(objdir)/fesm_datasrc.o

###############################################
##
## Compilation of complete programs
##
###############################################

# Static library compilation. The archive is a file target, so its mtime only
# changes when an object does (hosts such as yelmo depend on it to rebuild
# after a type change).
fesmutils-static: $(objdir)/libfesmutils.a

$(objdir)/libfesmutils.a: $(fesmutils_libs)
	ar rc $@ $(fesmutils_libs)
	ranlib $@
	@echo " "
	@echo "    $@ is ready."
	@echo " "

# coords library test programs
test-coords : fesmutils-static
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_geodinverse.x $(testdir)/coords/test_geodinverse.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_proj.x $(testdir)/coords/test_proj.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_gauss.x $(testdir)/coords/test_gauss.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_gaussian_filter.x $(testdir)/coords/test_gaussian_filter.f90 \
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
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_conservative_desclat.x $(testdir)/coords/test_conservative_desclat.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_conservative_rotated.x $(testdir)/coords/test_conservative_rotated.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_map_init_con.x $(testdir)/coords/test_map_init_con.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_map_cache.x $(testdir)/coords/test_map_cache.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_map_cdo.x $(testdir)/coords/test_map_cdo.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_map_con_proj.x $(testdir)/coords/test_map_con_proj.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_lon180_wrap.x $(testdir)/coords/test_lon180_wrap.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_varslice_map.x $(testdir)/coords/test_varslice_map.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_ncio_interp.x $(testdir)/coords/test_ncio_interp.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_parity.x $(testdir)/coords/test_parity.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_grid_cdo_rw.x $(testdir)/coords/test_grid_cdo_rw.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_map_init_names.x $(testdir)/coords/test_map_init_names.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_pole_fill.x $(testdir)/coords/test_pole_fill.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_regional_drop.x $(testdir)/coords/test_regional_drop.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_coupler.x $(testdir)/coords/test_coupler.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_grid_cf_read.x $(testdir)/coords/test_grid_cf_read.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_coupler_bitident.x $(testdir)/coords/test_coupler_bitident.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_grid_signature.x $(testdir)/coords/test_grid_signature.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_coupler_global.x $(testdir)/coords/test_coupler_global.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_gauss_bounds.x $(testdir)/coords/test_gauss_bounds.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_latlon_area.x $(testdir)/coords/test_latlon_area.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		@echo " "
		@echo "    coords test programs ready in $(bindir)/."
		@echo " "

# root_finder test program
# phys_constants test program. Run it from the repo root so that the relative
# path to par/phys_const_earth.nml resolves:
#     ./bin/test_phys_constants.x
test-phys-constants : fesmutils-static
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_phys_constants.x \
			$(testdir)/phys_constants/test_phys_constants.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		@echo " "
		@echo "    phys_constants test program ready: $(bindir)/test_phys_constants.x"
		@echo " "

test-root : fesmutils-static
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_roots.x $(testdir)/root/test_roots.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		@echo " "
		@echo "    root_finder test program ready: $(bindir)/test_roots.x"
		@echo " "

bench-coords : fesmutils-static
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/bench_conservative.x $(testdir)/coords/bench_conservative.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		@echo " "
		@echo "    coords benchmark ready: $(bindir)/bench_conservative.x (needs cdo on PATH)."
		@echo " "

# ncio test programs: precision (bit rounding in nc_write), attributes on an open file (ncid)
test-ncio : fesmutils-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_NC) -o $(bindir)/test_ncio_precision.x $(testdir)/ncio/test_ncio_precision.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_NC) -o $(bindir)/test_ncio_attr_ncid.x $(testdir)/ncio/test_ncio_attr_ncid.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		@echo " "
		@echo "    ncio test programs ready: $(bindir)/test_ncio_precision.x, $(bindir)/test_ncio_attr_ncid.x"
		@echo " "

# variable_io (variable tables, var_io_write) test program
test-var-io : fesmutils-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_NC) -o $(bindir)/test_var_io.x $(testdir)/ncio/test_var_io.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		@echo " "
		@echo "    variable_io test program ready: $(bindir)/test_var_io.x"
		@echo " "

# series (tabulated time-series reader) test program
test-series : fesmutils-static
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_series.x $(testdir)/series/test_series.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		@echo " "
		@echo "    series test program ready: $(bindir)/test_series.x"
		@echo " "

# varslice (time-slicing of netcdf forcing) test program
test-varslice : fesmutils-static
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_varslice.x $(testdir)/series/test_varslice.f90 \
			$(INC_NC) -L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		@echo " "
		@echo "    varslice test program ready: $(bindir)/test_varslice.x"
		@echo " "

# multitopo (multi-resolution topography) test program
test-multitopo : fesmutils-static
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_multitopo.x $(testdir)/multitopo/test_multitopo.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		@echo " "
		@echo "    multitopo test program ready: $(bindir)/test_multitopo.x"
		@echo " "

# regions (FesmData regions, zones and basins) test program
test-regions : fesmutils-static
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_regions.x $(testdir)/fesmdata/test_regions.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		@echo " "
		@echo "    regions test program ready: $(bindir)/test_regions.x"
		@echo " "

# topodata (FesmData topography) test program
test-topodata : fesmutils-static
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_topodata.x $(testdir)/fesmdata/test_topodata.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		@echo " "
		@echo "    topodata test program ready: $(bindir)/test_topodata.x"
		@echo " "

# topoprep and the topography preparation of htopo (synthetic data, no files needed)
test-topoprep : fesmutils-static
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_topoprep.x $(testdir)/fesmdata/test_topoprep.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		@echo " "
		@echo "    topoprep test program ready: $(bindir)/test_topoprep.x"
		@echo " "

# htopo (hi-res geometry hub) test program
test-htopo : fesmutils-static
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_htopo.x $(testdir)/fesmdata/test_htopo.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		@echo " "
		@echo "    htopo test program ready: $(bindir)/test_htopo.x"
		@echo " "

# units test programs: fesmdata units (water-flux factor, reference
# temperature) and fesm_units (units_conv)
test-units : fesmutils-static
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_units.x $(testdir)/fesmdata/test_units.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_fesm_units.x $(testdir)/units/test_fesm_units.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		@echo " "
		@echo "    units test programs ready: $(bindir)/test_units.x, $(bindir)/test_fesm_units.x"
		@echo " "

# varslice climatology (sub-annual cycle without years) test program
test-varslice-clim : fesmutils-static
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_varslice_clim.x $(testdir)/varslice/test_varslice_clim.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		@echo " "
		@echo "    varslice climatology test program ready: $(bindir)/test_varslice_clim.x"
		@echo " "

# tsgen (time-series generator) test program
test-tsgen : fesmutils-static
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_tsgen.x $(testdir)/tsgen/test_tsgen.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		@echo " "
		@echo "    tsgen test program ready: $(bindir)/test_tsgen.x"
		@echo " "

# zremap (elevation remapping) test program
test-zremap : fesmutils-static
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_zremap.x $(testdir)/zremap/test_zremap.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		@echo " "
		@echo "    zremap test program ready: $(bindir)/test_zremap.x"
		@echo " "

# fesm_log (shared logging and error reporting) test program
test-log : fesmutils-static
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_log.x $(testdir)/log/test_log.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		@echo " "
		@echo "    fesm_log test program ready: $(bindir)/test_log.x"
		@echo " "

# fesm_clock (run clock, cadences, schedules) test programs
test-clock : fesmutils-static
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_clock_exact.x $(testdir)/clock/test_clock_exact.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_clock_chunks.x $(testdir)/clock/test_clock_chunks.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_clock_reject.x $(testdir)/clock/test_clock_reject.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_clock_restart.x $(testdir)/clock/test_clock_restart.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_clock_forcing.x $(testdir)/clock/test_clock_forcing.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_clock_negative.x $(testdir)/clock/test_clock_negative.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_clock_schedule.x $(testdir)/clock/test_clock_schedule.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_clock_calendar.x $(testdir)/clock/test_clock_calendar.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_clock_convert.x $(testdir)/clock/test_clock_convert.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		@echo " "
		@echo "    fesm_clock test programs ready: $(bindir)/test_clock_*.x"
		@echo " "

# fesm_fields (exchange field containers) test programs: $(CHECK_FIELDS), and
# the expected-abort programs $(CHECK_FIELDS_ABORT) in test/fields/abort/
test-fields : fesmutils-static
		for t in $(CHECK_FIELDS); do \
			$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/$$t.x $(testdir)/fields/$$t.f90 \
				-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS) || exit 1; \
		done
		for t in $(CHECK_FIELDS_ABORT); do \
			$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/$$t.x $(testdir)/fields/abort/$$t.f90 \
				-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS) || exit 1; \
		done
		@echo " "
		@echo "    fesm_fields test programs ready: $(bindir)/test_fields_*.x, $(bindir)/fields_abort_*.x"
		@echo " "

# nml test programs: $(CHECK_NML), and the expected-abort programs
# $(CHECK_NML_ABORT) in test/nml/abort/
test-nml : fesmutils-static
		for t in $(CHECK_NML); do \
			$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/$$t.x $(testdir)/nml/$$t.f90 \
				-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS) || exit 1; \
		done
		for t in $(CHECK_NML_ABORT); do \
			$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/$$t.x $(testdir)/nml/abort/$$t.f90 \
				-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS) || exit 1; \
		done
		@echo " "
		@echo "    nml test programs ready: $(bindir)/test_nml_*.x, $(bindir)/nml_abort_*.x"
		@echo " "

# fesm_datasrc (field data sources) and its preparatory nml/varslice steps:
# test programs $(CHECK_DATASRC), and the expected-abort programs
# $(CHECK_DATASRC_ABORT) in test/datasrc/abort/
test-datasrc : fesmutils-static
		for t in $(CHECK_DATASRC); do \
			$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/$$t.x $(testdir)/datasrc/$$t.f90 \
				-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS) || exit 1; \
		done
		for t in $(CHECK_DATASRC_ABORT); do \
			$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/$$t.x $(testdir)/datasrc/abort/$$t.f90 \
				-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS) || exit 1; \
		done
		@echo " "
		@echo "    fesm_datasrc test programs ready in $(bindir)/."
		@echo " "

# zremap documentation demo (Greenland SMB remapped to the LGM; needs ../ice_data)
demo-zremap : fesmutils-static
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/demo_zremap_grl.x $(testdir)/zremap/demo_zremap_grl.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		@echo " "
		@echo "    zremap demo program ready: $(bindir)/demo_zremap_grl.x"
		@echo " "

# map_field apply-throughput benchmarks (per-timestep remap cost, no cdo needed)
bench-apply : fesmutils-static
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/bench_apply.x $(testdir)/coords/bench_apply.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/bench_alloc.x $(testdir)/coords/bench_alloc.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		@echo " "
		@echo "    apply benchmarks ready: $(bindir)/bench_apply.x, $(bindir)/bench_alloc.x"
		@echo " "

###############################################
##
## Test suite: make check
##
## `make check` builds each self-contained test program (via the test-*
## targets above) and runs it from the repo root, failing on a non-zero exit.
## One PASS/FAIL/SKIP line per test; full output in logs/check/<test>.log.
## A group runs all its tests before failing; `make -k check` also continues
## past a failing group. The cdo-based tests (test_map_cdo, test_map_con_proj)
## report SKIP when cdo is not on PATH. Tests that need ../ice_data run under
## `make check-data` instead. Benchmarks and demos are not run.
##
###############################################

CHECK_COORDS = test_geodinverse test_proj test_gauss test_gaussian_filter test_rotated \
               test_kdtree test_weightmap test_mapinit test_bilinear test_ccsm3_grl \
               test_mapio test_polygons test_conservative test_conservative_xsys \
               test_conservative_sphere test_conservative_desclat test_conservative_rotated \
               test_map_init_con test_map_cache test_map_cdo test_map_con_proj \
               test_lon180_wrap test_varslice_map test_ncio_interp test_parity \
               test_grid_cdo_rw test_map_init_names test_pole_fill test_regional_drop \
               test_coupler test_grid_cf_read test_coupler_bitident test_grid_signature \
               test_coupler_global test_gauss_bounds test_latlon_area
CHECK_DATA   = test_regions test_topodata test_htopo
CHECK_CLOCK  = test_clock_exact test_clock_chunks test_clock_reject test_clock_restart \
               test_clock_forcing test_clock_negative test_clock_schedule test_clock_calendar \
               test_clock_convert
CHECK_FIELDS = test_fields_declare test_fields_put_get test_fields_mean test_fields_accum \
               test_fields_mask test_fields_units test_fields_link test_fields_restart \
               test_fields_exchange test_fields_nonfinite
CHECK_FIELDS_ABORT = $(sort $(basename $(notdir $(wildcard $(testdir)/fields/abort/fields_abort_*.f90))))
CHECK_DATASRC = test_find_group test_varslice_legacy test_varslice_schema test_time_axis \
                test_nml_quoted test_varslice_units_attr test_varslice_get test_varslice_interval \
                test_static test_transient test_clim test_scalar test_anomaly test_bind test_examples
CHECK_DATASRC_ABORT = $(sort $(basename $(notdir $(wildcard $(testdir)/datasrc/abort/datasrc_abort_*.f90))))
CHECK_NML    = test_nml_validate
CHECK_NML_ABORT = $(sort $(basename $(notdir $(wildcard $(testdir)/nml/abort/nml_abort_*.f90))))
ICE_DATA     = ../ice_data/v2/NorthTest

# $(call check_run,<tests>): run each $(bindir)/<test>.x, log its output and
# print one summary line; exit non-zero if any test did. A test that exits 0
# after printing "SKIP: <test> (<reason>)" (a missing optional tool) is
# reported as skipped.
check_run = mkdir -p logs/check; nfail=0; \
	for t in $(1); do \
		log=logs/check/$$t.log; ./$(bindir)/$$t.x > $$log 2>&1; rc=$$?; \
		if [ $$rc -ne 0 ]; then \
			echo "  FAIL  $$t (exit $$rc, log: $$log)"; tail -n 20 $$log; nfail=$$((nfail+1)); \
		elif grep -q 'SKIP:' $$log; then \
			echo "  SKIP  $$(grep -m1 'SKIP:' $$log | sed 's/.*SKIP: *//')"; \
		else \
			echo "  PASS  $$t"; \
		fi; \
	done; \
	[ $$nfail -eq 0 ]

# $(call check_abort,<tests>): as check_run, for programs that must stop
# through fesm_error: each passes only if it exits non-zero, the log holds
# "stopped by fesm_error" and the message line announced by its
# "EXPECT: <msg>" line, and names the offending field `probe`.
check_abort = mkdir -p logs/check; nfail=0; \
	for t in $(1); do \
		log=logs/check/$$t.log; ./$(bindir)/$$t.x > $$log 2>&1; rc=$$?; \
		msg=$$(grep -m1 'EXPECT:' $$log | sed 's/.*EXPECT: *//'); \
		if [ $$rc -eq 0 ]; then \
			echo "  FAIL  $$t (exit 0, expected an abort; log: $$log)"; tail -n 20 $$log; nfail=$$((nfail+1)); \
		elif grep -q 'stopped by fesm_error' $$log && grep -qxF "    $$msg" $$log && \
		     grep -v 'EXPECT:' $$log | grep -q 'probe'; then \
			echo "  PASS  $$t (aborted: $$msg)"; \
		else \
			echo "  FAIL  $$t (exit $$rc without the expected error '$$msg'; log: $$log)"; tail -n 20 $$log; \
			nfail=$$((nfail+1)); \
		fi; \
	done; \
	[ $$nfail -eq 0 ]

.PHONY : check check-data check-coords check-phys-constants check-root check-ncio \
         check-var-io check-series check-varslice check-multitopo check-units \
         check-varslice-clim check-tsgen check-zremap check-log check-clock check-fields \
         check-datasrc check-nml check-nml-compat check-topoprep

check : check-coords check-phys-constants check-root check-ncio check-var-io \
        check-series check-varslice check-multitopo check-units check-varslice-clim \
        check-tsgen check-zremap check-log check-clock check-fields check-datasrc \
        check-nml check-topoprep
		@echo " "
		@echo "    make check: all tests passed (data tests: make check-data)."
		@echo " "

check-coords : test-coords
		@$(call check_run,$(CHECK_COORDS))
check-phys-constants : test-phys-constants
		@$(call check_run,test_phys_constants)
check-root : test-root
		@$(call check_run,test_roots)
check-ncio : test-ncio
		@$(call check_run,test_ncio_precision test_ncio_attr_ncid)
check-var-io : test-var-io
		@$(call check_run,test_var_io)
check-series : test-series
		@$(call check_run,test_series)
check-varslice : test-varslice
		@$(call check_run,test_varslice)
check-multitopo : test-multitopo
		@$(call check_run,test_multitopo)
check-units : test-units
		@$(call check_run,test_units test_fesm_units)
check-varslice-clim : test-varslice-clim
		@$(call check_run,test_varslice_clim)
check-tsgen : test-tsgen
		@$(call check_run,test_tsgen)
check-zremap : test-zremap
		@$(call check_run,test_zremap)
check-log : test-log
		@$(call check_run,test_log)
check-topoprep : test-topoprep
		@$(call check_run,test_topoprep)
check-clock : test-clock
		@$(call check_run,$(CHECK_CLOCK))
check-fields : test-fields
		@$(call check_run,$(CHECK_FIELDS)); rc_run=$$?; \
		$(call check_abort,$(CHECK_FIELDS_ABORT)) && [ $$rc_run -eq 0 ]
check-datasrc : test-datasrc
		@$(call check_run,$(CHECK_DATASRC)); rc_run=$$?; \
		$(call check_abort,$(CHECK_DATASRC_ABORT)) && [ $$rc_run -eq 0 ]
check-nml : test-nml
		@$(call check_run,$(CHECK_NML)); rc_run=$$?; \
		$(call check_abort,$(CHECK_NML_ABORT)) && [ $$rc_run -eq 0 ]

# nml parse compatibility against a reference commit (not part of `check`):
#     make check-nml-compat ref=v1.4 [files="a.nml b.nml"]
# Every key of the files, read as a string vector, must parse the same at
# `ref` and in the working tree (test/nml/nml_compat.sh). files defaults to
# the yelmox, yelmo, chion and FastIsostasy namelists next to fesm-utils.
NML_COMPAT_FILES = $(wildcard ../yelmox/yelmox/*.nml ../yelmox/input/*.nml \
                   ../yelmo/input/yelmo_defaults.nml ../chion/input/*.nml ../FastIsostasy/par/*.nml)
check-nml-compat :
		@if [ -z "$(ref)" ]; then echo "  usage: make check-nml-compat ref=<commit> [files=...]"; exit 2; fi
		@if [ -z "$(or $(files),$(NML_COMPAT_FILES))" ]; then \
			echo "  SKIP  check-nml-compat (no namelist files: pass files=...)"; \
		else \
			test/nml/nml_compat.sh $(ref) $(or $(files),$(NML_COMPAT_FILES)); \
		fi

# FesmData tests on the NorthTest grids; need ice_data checked out next to
# fesm-utils (not available in CI).
check-data :
		@if [ -d $(ICE_DATA) ]; then \
			$(MAKE) --no-print-directory test-regions test-topodata test-htopo && \
			$(call check_run,$(CHECK_DATA)); \
		else \
			echo "  SKIP  $(CHECK_DATA) (needs $(ICE_DATA), not found)"; \
		fi

###############################################
##
## External libraries: FFTW, LIS, SHTns
##
## Each target runs the vendored package's own autotools build
## (configure/make/make install) into <lib>/<lib>-{serial,omp}. The
## machine/compiler configure flags come from config/libs.mk. Pick the variant
## with openmp=0 (serial, default) or openmp=1; pic=1 builds FFTW/LIS
## position-independent (for linking into a shared library):
##   make fftw openmp=0      make lis openmp=1      make shtns openmp=0
##   make libs openmp=1      # fftw + lis + shtns, one variant
##
## Build serially (a one-time autotools build; do not pass -j to `make libs`).
## SHTns links the FFTW built here, so fftw must be installed before shtns.
###############################################

FFTW_SRC     = fftw/fftw-3.3.10
LIS_SRC      = lis/lis-2.1.11
SHTNS_SRC    = SHTns/shtns-3.7.5
FFTW_PREFIX  = $(CURDIR)/fftw/fftw-$(LIB_VARIANT)
LIS_PREFIX   = $(CURDIR)/lis/lis-$(LIB_VARIANT)
SHTNS_PREFIX = $(CURDIR)/SHTns/shtns-$(LIB_VARIANT)

.PHONY : fftw lis shtns libs

fftw:
	@echo "=== fftw ($(LIB_VARIANT)) [machine=$(MACHINE) compiler=$(COMPILER)] ==="
	cd $(FFTW_SRC) && ./configure --prefix=$(FFTW_PREFIX) --disable-doc $(FFTW_OMP_OPT) $(CONF_PIC) \
		$(FFTW_CONF_VARS) && $(MAKE) clean && $(MAKE) && $(MAKE) install
	@echo "    installed: $(FFTW_PREFIX)"

lis:
	@echo "=== lis ($(LIB_VARIANT)) [machine=$(MACHINE) compiler=$(COMPILER)] ==="
	cd $(LIS_SRC) && ./configure --prefix=$(LIS_PREFIX) --enable-f90 --enable-saamg $(LIS_OMP_OPT) $(CONF_PIC) \
		$(LIS_CONF_VARS) && $(MAKE) clean && $(MAKE) && $(MAKE) install
	@echo "    installed: $(LIS_PREFIX)"

shtns:
	@echo "=== shtns ($(LIB_VARIANT)) [machine=$(MACHINE) compiler=$(COMPILER)] ==="
	cd $(SHTNS_SRC) && ./configure --prefix=$(SHTNS_PREFIX) $(SHTNS_OMP_OPT) \
		CPPFLAGS='-I$(FFTW_PREFIX)/include' LDFLAGS='-L$(FFTW_PREFIX)/lib' \
		$(SHTNS_CONF_VARS) && $(MAKE) clean && $(MAKE) && $(MAKE) install
	@echo "    installed: $(SHTNS_PREFIX)"

libs: fftw lis shtns
	@echo ""
	@echo "    fftw, lis, shtns ($(LIB_VARIANT)) built and installed."
	@echo ""

# Everything for one variant: the external libraries + the utils static library.
# This is what `configme install` runs (once per variant, openmp=0 then =1).
.PHONY : all
all: libs fesmutils-static

.PHONY : usage
usage:
	@echo ""
	@echo "    * USAGE * "
	@echo ""
	@echo " make fesmutils-static : builds libfesmutils.a (openmp=0|1, pic=0|1)"
	@echo " make test-phys-constants : builds the phys_constants test program"
	@echo " make fftw|lis|shtns   : builds one external library (openmp=0|1; pic=0|1 fftw/lis only)"
	@echo " make libs             : builds fftw + lis + shtns (openmp=0|1)"
	@echo " make test-coords      : compiles the coords test programs"
	@echo " make test-root        : compiles the root_finder test program"
	@echo " make test-ncio        : compiles the ncio test programs"
	@echo " make test-var-io      : compiles the variable_io test program"
	@echo " make test-regions     : compiles the regions test program"
	@echo " make test-topodata    : compiles the topodata test program"
	@echo " make test-htopo       : compiles the htopo test program"
	@echo " make test-topoprep   : compiles the topography preparation test program"
	@echo " make test-units       : compiles the units test programs (fesmdata, fesm_units)"
	@echo " make test-varslice-clim : compiles the varslice climatology test program"
	@echo " make test-datasrc     : compiles the fesm_datasrc test programs (with its nml/varslice steps)"
	@echo " make check            : builds and runs the self-contained tests (fails on any failure)"
	@echo " make check-<group>    : one group, e.g. check-coords, check-ncio, check-zremap"
	@echo " make check-data       : runs the FesmData tests that need ../ice_data"
	@echo " make check-nml-compat ref=<commit> [files=...] : nml parsing of namelist files vs a reference commit"
	@echo " make ncio-gen         : regenerates src/ncio.f90 from src/ncio.fypp (devs only)"
	@echo " make clean            : cleans object files"
	@echo ""

clean:
	rm -f $(bindir)/*.x
	rm -f  *.x gmon.out $(objdir_serial)/*.o $(objdir_serial)/*.mod $(objdir_serial)/*.a $(objdir_serial)/*.so
	rm -f  *.x gmon.out $(objdir_omp)/*.o $(objdir_omp)/*.mod $(objdir_omp)/*.a $(objdir_omp)/*.so
	rm -rf *.x.dSYM
