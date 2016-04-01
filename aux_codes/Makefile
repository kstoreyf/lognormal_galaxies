DEBUG_OR_OPTIMIZE=-O3
GSL_LIB=${GSL_HOME}/lib/
GSL_INC=${GSL_HOME}/include/
CXXFLAGS=$(DEBUG_OR_OPTIMIZE) -Wall -msse3
LINK  =-L${GSL_LIB} -lgsl -lgslcblas
INC   =-I${GSL_INC}
TARGET= calc_xi_gm discretize_pk kaiser_pk calc_pk
.SUFFIXES: .cpp
.cpp.o:
	$(CXX) -c $(CXXFLAGS) $< -o $@ $(INC)
default: $(TARGET)

calc_xi_gm: ../compute_pkG/integration_modified.o calc_xi_gm.o
	$(CXX) -o $@ $^ $(LINK)
discretize_pk: discretize_pk.o
	$(CXX) -o $@ $^ $(LINK)
kaiser_pk: ../generate_Poisson/spline.o kaiser_pk.o 
	$(CXX) -o $@ $^ $(LINK)
calc_pk: ../compute_pkG/integration_modified.o calc_pk.o
	$(CXX) -o $@ $^ $(LINK)

clean:
	rm -f *.o *~
tidy: clean
	rm -f $(TARGET)
