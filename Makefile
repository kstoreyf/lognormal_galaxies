# C++ compiler
CXX = g++
# Fortran90 compiler
FC = gfortran
# Location of the GSL library
GSL_HOME = /usr/local
# Location of the FFTW-3
FFTW_HOME = /usr/local
# Location of "omp.h"
OMP_INC = /usr/local/lib/gcc/5/gcc/x86_64-apple-darwin15.2.0/5.3.0/include/
# Flag to link openmp
OMP_FLAG = -fopenmp
#OMP_FLAG = 
# Additional links to library [like -static if needed]
ADD_LDFLAGS =

prog	= eisensteinhubaonu compute_xi compute_pkG generate_Poisson calculate_pk calculate_cross lognormal aux_codes

all: 
	@for p in $(prog); do \
	cd $$p; \
	$(MAKE) CXX=${CXX} FC=${FC} GSL_HOME=${GSL_HOME} FFTW_HOME=${FFTW_HOME} OMP_INC=${OMP_INC} OMP_FLAG=${OMP_FLAG} ADD_LDFLAGS=${ADD_LDFLAGS}; \
	cd ../; \
	done
clean:
	@for p in $(prog); do \
	cd $$p; \
	$(MAKE) clean; \
	cd ../; \
	done

tidy:
	@for p in $(prog); do \
	cd $$p; \
	$(MAKE) tidy; \
	cd ../; \
	done
