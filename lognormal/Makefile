FFLAGS	= -O3 
.SUFFIXES: .f90
.f90.o: 
	$(FC) $(FFLAGS) -c $<
default:  read_lognormal
read_lognormal: $(OBJS) read_lognormal.o
	$(FC) $(FFLAGS) -o $@ $(OBJS) $@.o $(LDFLAGS)
clean:
	rm -f *.o *.mod
tidy: clean
	rm -f read_lognormal
