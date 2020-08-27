#
# Meant to be called as a sub-make; relies on variables set in ../../Makefile
#
COMPILE = $(FC) $(FCFLAGS) $(FCINCLUDE) -c
%.o: %.F90
	$(COMPILE) $<
%: %.o
		$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)

VPATH = ./:$(RRTMGP_ROOT)/extensions/:$(RRTMGP_ROOT)/examples:/$(TEST_ROOT)/test/util/src/

#
# Extra sources -- extensions to RRTMGP classes, shared infrastructure, local sources
#
ADDITIONS  = mo_load_coefficients.o mo_simple_netcdf.o mo_garand_atmos_io.o
ADDITIONS += mo_test_files_io.o mo_heating_rates.o mo_fluxes_byband.o mo_fluxes_byband_kernels.o mo_rrtmgp_clr_all_sky.o

#
# Targets
#
all: test_flux_compute_topography

test_flux_compute_topography:  test_flux_compute_topography.o $(ADDITIONS)

test_flux_compute_topography.o: test_flux_compute_topography.F90 $(ADDITIONS)

# This program has a different name

# test_flux_compute_topography:   $(ADDITIONS) test_flux_compute_topography.o

mo_fluxes_byband.o: mo_fluxes_byband.F90 mo_fluxes_byband_kernels.o

mo_load_coefficients.o: mo_load_coefficients.F90 mo_simple_netcdf.o
mo_test_files_io.o:     mo_test_files_io.F90     mo_simple_netcdf.o

clean:
	-rm -f test_flux_compute_topography *.o  *.mod
