# C++ compiler
CXX = icpc
# Fortran90 compiler
FC = ifort
# Flag to link openmp
OMP_FLAG = -fopenmp
# Additional links to library [like -static if needed]
ADD_LDFLAGS = -static
OMP_INC = /usr/local/lib

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
