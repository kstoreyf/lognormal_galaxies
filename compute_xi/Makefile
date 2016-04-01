FFLAGS	= -O3 
OBJS	= linearpk.o spline.o hunt.o
.SUFFIXES: .f90
.f90.o: 
	$(FC) $(FFLAGS) -c $<
default:  compute_xi
compute_xi: $(OBJS) compute_xi.o
	$(FC) $(FFLAGS) -o $@ $(OBJS) $@.o $(LDFLAGS)
clean:
	rm -f *.o *.mod
tidy: clean
	rm -f compute_xi
