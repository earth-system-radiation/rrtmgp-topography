#
# Environment variable RRTMGP_ROOT needs to be set before building anything
#
export RRTMGP_BUILD = $(RRTMGP_ROOT)/build
export TEST_ROOT = ${CURDIR}/..

-include $(RRTMGP_BUILD)/Makefile.rules
-include $(RRTMGP_BUILD)/Makefile.conf
-include Makefile.libs
export LDFLAGS += -L$(RRTMGP_BUILD) -lrte -lrrtmgp
export FCINCLUDE += -I$(RRTMGP_BUILD)

export FCFLAGS += -I${RRTMGP_BUILD} -I${NFHOME}/include -I${NCHOME}/include
export LDFLAGS += -L${RRTMGP_BUILD} -L${NFHOME}/lib     -L${NCHOME}/lib -lrte -lrrtmgp -lnetcdff -lnetcdf
# compile tests (all)
all: test_flux_compute_top


# tests flux_compute including topography
test_flux_compute_top:
	@echo --- building test flux_compute ---
	( cd ../test/flux_compute_topography/build ; make )


# clean
clean:
	( cd ../test/flux_compute_topography/build ; make clean )

