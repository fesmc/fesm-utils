.SUFFIXES: .f .F .F90 .f90 .o .mod
.SHELL: /bin/sh

# Command-line options at make call
debug  ?= 0
openmp ?= 0 

# PATH options
srcdir = src
bindir = bin
libdir = libs
objdir_serial = include-serial
objdir_omp = include-omp

objdir = $(objdir_serial)
ifeq ($(openmp), 1)
    objdir = $(objdir_omp)
endif

## COMPILER CONFIGURATION ##
# (should be loaded from config directory)

<COMPILER_CONFIGURATION>

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

$(objdir)/grid_to_cdo.o: $(srcdir)/coordinates-light/grid_to_cdo.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/gaussian_filter.o: $(srcdir)/coordinates-light/gaussian_filter.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/index.o: $(srcdir)/coordinates-light/index.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/interp1D.o: $(srcdir)/coordinates-light/interp1D.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/interp2D.o: $(srcdir)/coordinates-light/interp2D.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/mapping_scrip.o: $(srcdir)/coordinates-light/mapping_scrip.f90 $(objdir)/ncio.o $(objdir)/interp2D.o \
								$(objdir)/gaussian_filter.o $(objdir)/index.o $(objdir)/grid_to_cdo.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/ncio.o: $(srcdir)/ncio.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_NC) -c -o $@ $<

$(objdir)/nml.o: $(srcdir)/nml.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/timer.o: $(srcdir)/timer.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/timeout.o: $(srcdir)/timeout.f90 $(objdir)/nml.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/timestepping.o: $(srcdir)/timestepping.f90 $(objdir)/nml.o $(objdir)/ncio.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/variable_io.o: $(srcdir)/variable_io.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/varslice.o: $(srcdir)/varslice.f90 $(objdir)/nml.o $(objdir)/ncio.o $(objdir)/mapping_scrip.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

fesmutils_libs = 	$(objdir)/gaussian_filter.o \
					$(objdir)/grid_to_cdo.o \
					$(objdir)/index.o \
					$(objdir)/interp1D.o \
					$(objdir)/interp2D.o \
					$(objdir)/mapping_scrip.o \
					$(objdir)/ncio.o \
					$(objdir)/nml.o \
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

test : fesmutils-static
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test.x $(srcdir)/test.f90 \
			-L${CURDIR}/$(objdir) -lfesmutils $(LFLAGS)
		@echo " "
		@echo "    $(bindir)/test.x is ready."
		@echo " "

.PHONY : usage
usage:
	@echo ""
	@echo "    * USAGE * "
	@echo ""
	@echo " make test       : compiles test.x"
	@echo " make clean      : cleans object files"
	@echo ""

clean:
	rm -f $(bindir)/*.x
	rm -f  *.x gmon.out $(objdir_serial)/*.o $(objdir_serial)/*.mod $(objdir_serial)/*.a $(objdir_serial)/*.so
	rm -f  *.x gmon.out $(objdir_omp)/*.o $(objdir_omp)/*.mod $(objdir_omp)/*.a $(objdir_omp)/*.so
	rm -rf *.x.dSYM
