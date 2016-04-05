*** Codes for generating log-normal realisations, and computing the monopole, quadrupole, and hexadecapole power spectra ***

<< HISTORY >>
- Originally developed by Donghui Jeong (Penn State Univ.) in 2011
- Enhanced by Chi-Ting Chiang (Stony Brook Univ.) through 2015
  - arXiv:1306.4157 is based on this code
- Packaged by Eiichiro Komatsu on December 28, 2015
- (v2) includes a fix in "calc_pk_const_los_ngp" by Issha Kayo, Jan 12, 2016
- (v3) Added cross-power spectrum code by Donghui Jeong, February 8, 2016
- (v4) Generate velocities from the matter density field, instead of the galaxy density field divided by the linear bias, by Aniket Agrawal, March 4, 2016
- (v4.1) New Makefile [by Ryu Makiya], making it easier to change compilers etc, and re-packaged with cleaned python scripts, March 30, 2016

This package consists of the following steps:

1) eisensteinhubaonu/
- Generate the input linear power spectrum using Eisenstein&Hu's fitting formula including the effect of massive neutrinos. Based on Eisenstein & Hu, ApJ, 511, 5 (1999), but extended to include the Baryon Acoustic Oscillation (BAO)

2) compute_xi/
- Fourier transform the power spectrum to generate the two-point correlation function, xi(r)

3) compute_pkG/
- Use xi(r) to compute the power spectrum of the log-normal field, G, as defined in Appendix B of arXiv:1306.4157

4) generate_Poisson/
- Generate log-normal realisations of locations and velocities of galaxies, and store them in "lognormal/"

5) calculate_pk/
- Calculate the monopole, quadrupole, and hexadecapole power spectra and store them in "pk/"

6) calculate_cross/
- Calculate the monopole, quadrupole, and hexadecapole cross-power spectra and store them in "cross_pk/"

7) lognormal/
- If necessary, you can read the simulated binary data using a sample code given in this directory

8) aux_codes/										
- This directory contains some useful codes, such as:
	- Calculate the predicted galaxy-matter cross power spectrum and the cross correlation function
	- Calculate the monopole and quadrupole power spectra using the Kaiser prediction, taking into account the galaxy-matter cross
	- Calculate the theoretical power spectrum averaged over the same Fourier grids as calc_pk

=============
User's Manual
=============

1) Edit Makefile in the top directory "codes_lognormal/". The necessary information includes:

CXX	 = [C++ compiler]
FC 	 = [fortran 90 compiler]
GSL_HOME = [Location of the GSL library]
FFTW_HOME= [Location of the FFTW-3]
OMP_INC	 = [Location of "omp.h"]
OMP_FLAG = [Option to link OpenMP, e.g., -fopenmp]
ADD_LDFLAGS=[Other links to libraries, e.g., -static]


2) In the top directory, execute "make". This will compile all the codes. Make sure to execute "make" twice to see if there was a failure in compiling the codes. It may generate many warnings such as "unused variable", but you can ignore them. Upon successful compilation, you should see:
----
make[1]: Nothing to be done for `default'.
make[1]: Nothing to be done for `default'.
make[1]: Nothing to be done for `default'.
make[1]: Nothing to be done for `default'.
make[1]: Nothing to be done for `default'.
make[1]: Nothing to be done for `default'.
make[1]: Nothing to be done for `default'.
----
or something similar.

3) To generate the input power spectrum [wavenumber_pk.txt], correlation function [Rh_xi.txt], the log-normal galaxy power spectrum [pkG_rmaxNNNN_bYYYY.dat] and the log-normal matter power spectrum [pkG_rmaxNNNN_b1.0.dat] (used to generate velocities), run the shell script "generate_inputs.sh". In this file, edit:

----
compute_pk2 1.3         # specify the redshift
compute_pkG 1.455       # specify the linear galaxy bias
----

to use the desired values of the redshift and galaxy bias. To run the shell script, type

sh generate_inputs.sh

This may take a few minutes, so be patient. The output files contain:

wavenumber_pk.txt
- 1st column: wavenumber [h/Mpc]
- 2nd column: power spectrum [Mpc^3/h^3]
[you may compare it with eisensteinhubaonu/wavenumber_pk_mnu0.2_z1.3.txt]

wavenumber_fnu.txt
- 1st column: wavenumber [h Mpc^-1]
- 2nd column: logarithmic growth rate - fnu - for massive (massless) neutrinos which varies (does not vary) as a function of 1st column

Rh_xi.txt
- 1st column: comoving separation [Mpc/h]
- 2nd column: correlation function [dimensionless]
[you may compare it with compute_xi/Rh_xi_mnu0.2_z1.3.txt]

pkG_rmaxNNNN_bYYYY.dat [NNNN refers to the maximum separation used, and YYYY refers to the value of the galaxy bias]
- 1st column: wavenumber [h/Mpc]
- 2nd column: galaxy power spectrum [Mpc^3/h^3]
[you may compare it with compute_pkG/pkG_rmax10000_b1.455_mnu0.2_z1.3.dat]

pkG_rmaxNNNN_b1.0.dat [NNNN refers to the maximum separation used]
- 1st column: wavenumber [h/Mpc]
- 2nd column: matter power spectrum [Mpc^3/h^3]

4) Open "pkG_rmaxNNNN_bYYYY.dat" and look for negative values in the 2nd column. Remove them by hand if you find any

5) To generate the log-normal realisation, use a Python wrapper as "python run_genPoissonmock_iseed.py". It will ask for the number of realisations and the initial random number seed. This code is Open-MP parallelised. How to use the parallelisation depends on platforms, but a sample shell script is given in "run_parallel_example.sh". 

The important parameters to specify in "run_genPoissonmock_iseed.py" are:
----
# cosmological parameters
Omegam=0.272            # Omega matter for computing the velocity field
Omegade=1-Omegam
zz=1.3                  # redshift for computing the velocity field
aHz=100.*pow(Omegam*pow(1.+zz,3)+Omegade,0.5)/(1.+zz)

# parameters for lognormal realization
Ngalaxies = 8345000     # number of galaxy in integer
bias = 1.455            # linear bias for computing the matter fluctuation and then velocity field (need to be consistent with input P(k))

Lx = 3666       # box size in x [h^-1 Mpc]
Ly = 1486       # box size in y [h^-1 Mpc]
Lz = 734        # box size in z [h^-1 Mpc]
Pnmax = 2500   # mesh number of the longest dimension for lognormal realization
----

Make sure to use the same cosmological parameters and the bias parameter that were used to generate the input power spectra.

The larger "Pnmax", the larger the memory that the code requires. Start with a smaller value such as 512 when testing the code.

The results and the used parameter files will be stored in "lognormal/". The positions and velocities of galaxies are stored in binary files. To read them, use a sample code "read_lognormal.f90" in that directory. Edit Makefile and compile it.

6) To compute the power spectra, use a Python wrapper as "python run_calc_pk_ngp.py". It will ask for the number of realisations and the initial random number seed used for generating log-normal realisations. Use the consistent values. This code does not generate random numbers, but these numbers are used to identify the right filenames. This code is Open-MP parallelised also.

You should choose the parameters in "run_calc_pk_ngp.py" consistently with those in "run_genPoissonmock_iseed.py". The important new parameters are:
----
losx = 0        # line-of-sight direction of x (this is a vector, and enter (0,0,0) for real-space calculation)
losy = 0        # line-of-sight direction of y
losz = 0        # line-of-sight direction of z
----

When all of them are equal to 0, the real-space power spectra are calculated. When one of them is equal to 1, the positions of galaxies are shifted by the velocity in that direction. AT MOST one of them should be equal to 1 - do not set more than one components to be 1!

The results and the used parameter files will be stored in "pk/". The format is:

pk_ngp_isN_rlzM_zXXXX_bYYYY_tophat_Rs10_los0.0.0.dat [or los1.0.0 etc depending on which component was set to 1]
1st column: wavenumber [h/Mpc]
2nd column: monopole power spectrum [Mpc^3/h^3]
3rd column: quadrupole power spectrum [Mpc^3/h^3]
4th column: hexadecapole power spectrum [Mpc^3/h^3]
5th column: the number of Fourier modes averaged

To check the results, compute some number of real-space monopole power spectra, average them, and compare the average with the input matter power spectrum "wavenumber_pk.txt" times the bias squared. You may find a deviation on small scales due to the resolution effect of Fourier meshes, but it should reproduce the input power spectrum very precisely on large scales.

The first few bins may return "nan", but do not worry about them.

7) To compute the predictions for the monopole and quadrupole of the redshift space power spectrum in Kaiser limit, the cross correlation function and the cross power spectrum between generated galaxy and matter density fields are needed. We also need to average the theoretical power spectrum over the same Fourier grids as in the measurements. To do this, there are Python wrappers such as "run_calc_pk.py", "run_calc_xi_gm.py", "run_discretize_pk.py", and "run_kaiser_pk.py". The details of how to run these codes are described in a readme file in "aux_codes".
A script to generate the Kaiser prediction using all of the codes in aux_codes is provided as "generate_kaiser.sh". The script produces a cross correlation function file(xi_gm_kmax*), a cross power spectrum file(pk_rmax*_b1.dat), a Kaiser P(k) file (pk_kaiser.dat) containing 3 columns - k, P_0(k), P_2(k) and discretized versions of the Kaiser P(k) - pkt_kaiser_2.dat(discretized monopole) and pkt_kaiser_3.dat(discretized quadrupole).
NOTE : The script assumes filenames and other parameters consistent with the other provided scripts. Please make necessary changes if changing the default scripts.  
