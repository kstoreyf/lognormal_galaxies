*** Two-point Correlation Function ***
September 18, 2008: E.Komatsu

Here we provide a program for computing the two-point correlation function,
xi(r), from the power spectrum.

The program generates a spherically-averaged two-point correlation function
as a function of R (in units of h^-1 Mpc) given by

xi(R) = int (k^2 dk)/(2 pi^2) P(k) sin(kR)/(kR)

where P(k) is a spherically-averaged power spectrum given by

P(k) = int_0^1 dmu P(k,mu)

- To compile and use the  program, edit Makefile
and simply "make"
- It will generate executables called "compute_xi".
- Running compute_xi will generate the data file called "Rh_xi.txt".
