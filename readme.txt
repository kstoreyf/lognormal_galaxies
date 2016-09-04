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
- (v5) By Ryu Makiya, Sep 04, 2016.
	-- Add new python script, run.py, which executes the all steps of the simulation
	-- Includes a new option of estimating the power spectrum, in which the power spectrum is estimated in the cubic box
	   which is large enough to emcompass the whole survey region. (set calc_mode_pk = 1 in .ini file to use it)

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

3)  To run the code, type:

> ./run.py example.ini

where exmaple.ini is the exmaple of input file.
It contains the input parameters of simulation e.g. cosmological parameters, survey geometries, etc.

You can also run each code independently, for example:

> ./calculate_pk/calc_pk_const_los_ngp

In this case code will ask you the input parameters.

4) The resulting files are as follows:

	a) Files in data/inputs
	If your set gen_inputs = False, calculation of following files will be skipped.

	exmaple_pk.txt
	- Input power spectrum
	- 1st column: wavenumber [h/Mpc]
	- 2nd column: power spectrum [Mpc^3/h^3]
	- Note: you can also use the power spectrum computed by external code (e.g. CAMB) by specifying the file name in .ini file

	example_fnu.txt
	- logarithmic growth rate fnu for massive (massless) neutrinos which varies (does not vary) as a function of k
	- 1st column: wavenumber [h Mpc^-1]
	- 2nd column: fnu

	example_Rh_xi.txt
	- correlation function calculated from input Pk
	- 1st column: comoving separation [Mpc/h]
	- 2nd column: correlation function [dimensionless]

	example_pkG_bYYYY.dat [YYYY refers to the value of the galaxy bias]
	- power spectrum of the log-normal field, G
	- 1st column: wavenumber [h/Mpc]
	- 2nd column: galaxy power spectrum [Mpc^3/h^3]

	example_pkG_b1.0.dat
	- pkG for matter density field
	- 1st column: wavenumber [h/Mpc]
	- 2nd column: matter power spectrum [Mpc^3/h^3]

	b) Files in data/lognormal
	For each realization, the code generates the followng two binary files:
 
	example_lognormal_rlz{# of realization}.bin (for galaxy field) 
	example_density_lognormal_rlz{# of realization}.bin (for density field)

	The positions and velocities of galaxies are stored in those files. To read them, use a sample code "lognormal/read_lognormal.f90". 
	The larger "Pnmax", the larger the memory that the code requires. Start with a smaller value such as 512 when testing the code.

	c) Files in data/pk
	exmaple_pk_rlz{# of realization}.dat
	- The power spectrum calculated from each mock galaxy catalog
	- 1st column: wavenumber [h/Mpc]
	- 2nd column: monopole power spectrum [Mpc^3/h^3]
	- 3rd column: quadrupole power spectrum [Mpc^3/h^3]
	- 4th column: hexadecapole power spectrum [Mpc^3/h^3]
	- 5th column: the number of Fourier modes averaged

	When all of los parameters, losx, losy and losz, are equal to 0, the real-space power spectra are calculated. 
	When one of them is equal to 1, the positions of galaxies are shifted by the velocity in that direction. 
	AT MOST one of them should be equal to 1 - do not set more than one components to be 1!

5) To check the results, compute some number of real-space monopole power spectra, average them, and compare the average with the input matter power spectrum times the bias squared. You may find a deviation on small scales due to the resolution effect of Fourier meshes, but it should reproduce the input power spectrum very precisely on large scales.
The first few bins may return "nan", but do not worry about them.

6) To compute the predictions for the monopole and quadrupole of the redshift space power spectrum in Kaiser limit, the cross correlation function and the cross power spectrum between generated galaxy and matter density fields are needed. We also need to average the theoretical power spectrum over the same Fourier grids as in the measurements. To do this, there are Python wrappers such as "run_calc_pk.py", "run_calc_xi_gm.py", "run_discretize_pk.py", and "run_kaiser_pk.py". The details of how to run these codes are described in a readme file in "aux_codes".
A script to generate the Kaiser prediction using all of the codes in aux_codes is provided as "generate_kaiser.sh". The script produces a cross correlation function file(xi_gm_kmax*), a cross power spectrum file(pk_rmax*_b1.dat), a Kaiser P(k) file (pk_kaiser.dat) containing 3 columns - k, P_0(k), P_2(k) and discretized versions of the Kaiser P(k) - pkt_kaiser_2.dat(discretized monopole) and pkt_kaiser_3.dat(discretized quadrupole).
NOTE : The script assumes filenames and other parameters consistent with the other provided scripts. Please make necessary changes if changing the default scripts.  
