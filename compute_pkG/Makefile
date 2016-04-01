DEBUG_OR_OPTIMIZE=-O3
GSL_LIB=${GSL_HOME}/lib/
GSL_INC=${GSL_HOME}/include/
CXXFLAGS=$(DEBUG_OR_OPTIMIZE) -Wall -msse3
LINK  =-L${GSL_LIB} -lgsl -lgslcblas
INC   =-I${GSL_INC}
TARGET= calc_pkG 
.SUFFIXES: .cpp
.cpp.o:
	$(CXX) -c $(CXXFLAGS) $< -o $@ $(INC)
default: $(TARGET)
calc_pkG: integration_modified.o calc_pkG.o
	$(CXX) -o $@ $^ $(LINK)
clean:
	rm -f *.o *~
tidy: clean
	rm -f $(TARGET)
