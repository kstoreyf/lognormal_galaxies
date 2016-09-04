FFLAGS	= -O3
OBJS	= cosmo.o growth.o spline.o hunt.o qromb.o dverk.o tf_fit.o tf_fit_nu.o eisensteinhu.o
.SUFFIXES: .f90
.f90.o: 
	$(FC) $(FFLAGS) -c $<
default: compute_pk
compute_pk: $(OBJS) compute_pk.o
	$(FC) $(FFLAGS) -o $@ $(OBJS) $@.o $(LDFLAGS)
clean:
	rm -f *.o *.mod *~
tidy: clean
	rm -f compute_pk
