DEBUG_OR_OPTIMIZE=-O3
GSL_LIB=${GSL_HOME}/lib/
GSL_INC=${GSL_HOME}/include/
FFTW_LIB=${FFTW_HOME}/lib/
FFTW_INC=${FFTW_HOME}/include/
CXXFLAGS=$(DEBUG_OR_OPTIMIZE) -Wall -msse3 $(OMP_FLAG)
LDFLAGS=
LINK  =-L${FFTW_LIB} -L${GSL_LIB} -lgsl -lgslcblas
INC   =-I${FFTW_INC} -I${GSL_INC} -I${OMP_INC}
TARGET= calc_cpk_const_los_v2
.SUFFIXES: .cpp
.cpp.o:
	$(CXX) -c $(CXXFLAGS) $< -o $@ $(INC)
default: $(TARGET)
calc_cpk_const_los_v2: matrix_inversion.o calc_cpk_const_los_v2.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LINK) -lfftw3_threads -lpthread -lfftw3 -lm $(LDFLAGS)
clean:
	rm -f *.o *~
tidy: clean
	rm -f $(TARGET)
