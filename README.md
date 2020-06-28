# Lognormal mock catalogs

This code generates lognormal mock galaxy catalogs. 
It is via Eiichiro Komatsu on bitbucket: https://bitbucket.org/komatsu5147/lognormal_galaxies/src/master/.
I have just added some tweaks for ease of generating many catalogs without overwriting.

This repo is also keeping track of the configs used to generate the mocks used in my https://github.com/kstoreyf/continuous-estimator project.

The original README follows below.

# Log-normal galaxies
Codes for generating log-normal realisations of galaxies in redshift-space, and computing the monopole, quadrupole, and hexadecapole power spectra.

Reference: A. Agrawal, R. Makiya, C.-T. Chiang, D. Jeong, S. Saito, and E. Komatsu, arXiv:1706.09195

## History

- Originally developed by Donghui Jeong (Penn State Univ.) in 2011
- Enhanced by Chi-Ting Chiang (Stony Brook Univ.) through 2015 (arXiv:1306.4157 is based on this code)
- Packaged by Eiichiro Komatsu on December 28, 2015
- (v2) includes a fix in "calc\_pk\_const\_los\_ngp" by Issha Kayo, Jan 12, 2016
- (v3) Added cross-power spectrum code by Donghui Jeong, February 8, 2016
- (v4) Generate velocities from the matter density field, instead of the galaxy density field divided by the linear bias, by Aniket Agrawal, March 4, 2016
- (v4.1) New Makefile [by Ryu Makiya], making it easier to change compilers etc, and re-packaged with cleaned python scripts, March 30, 2016
- (v5) By Ryu Makiya, Sep 04, 2016.
    - Added new python script, run.py, enabling to execute the all steps of the simulation all at once 
    - Includes a new option for the estimation of the power spectrum, in which the Pk is estimated in the cubic box which is large enough to encompass the whole survey region. (set calc\_mode\_pk = 1 in .ini file to use it)

## Overview
This package consists of the following steps:

1. eisensteinhubaonu/
	- Generate the input linear power spectrum using Eisenstein & Hu's fitting formula including the effect of massive neutrinos. Based on Eisenstein & Hu, ApJ, 511, 5 (1999), but extended to include the Baryon Acoustic Oscillation (BAO)

1. compute_xi/
	- Fourier transform the power spectrum to generate the two-point correlation function, $\xi$(r)

1. compute_pkG/
	- Compute the power spectrum of the log-normal field, G, as defined in Appendix B of arXiv:1306.4157, using $\xi(r)$

1. generate_Poisson/
	- Generate log-normal realisations of locations and velocities of galaxies, and store them in "lognormal/"

1. calculate_pk/
	- Calculate the monopole, quadrupole, and hexadecapole power spectra and store them in "pk/"

1. calculate_cross/
	- Calculate the monopole, quadrupole, and hexadecapole galaxy-matter cross-power spectra and store them in "cross\_pk/"

1. lognormal/
	- If necessary, you can read the simulated binary data using a sample code given in this directory

1. aux_codes/										
	- This directory contains some useful codes, such as:
		- Calculate the predicted galaxy-matter cross power spectrum and the cross correlation function
		- Calculate the monopole and quadrupole power spectra using the Kaiser prediction, taking into account the galaxy-matter cross
		- Calculate the theoretical power spectrum averaged over the same Fourier grids as calc_pk


## User's Manual
1. First you need to install [GSL - GNU Scientific Library](https://www.gnu.org/software/gsl/) and [FFTW](http://www.fftw.org/)

1. Edit Makefile in the top directory "codes\_lognormal/". The necessary information includes:

        CXX	 = [C++ compiler]  
        FC 	= [fortran 90 compiler]  
        GSL\_HOME = [Location of the GSL library]  
        FFT_HOME = [Location of the FFTW-3]  
        OMP\_INC = [Location of "omp.h"]  
        OMP\_FLAG = [Option to link OpenMP, e.g., -fopenmp]  
        ADD\_LDFLAGS = [Other links to libraries, e.g., -static]


1. In the top directory, execute "make". This will compile all the codes. Make sure to execute "make" twice to see if there was a failure in compiling the codes. It may generate many warnings such as "unused variable", but you can ignore them. Upon successful compilation, you should see:

        make[1]: Nothing to be done for `default'.  
        make[1]: Nothing to be done for `default'.  
        make[1]: Nothing to be done for `default'.  
        make[1]: Nothing to be done for `default'.  
        make[1]: Nothing to be done for `default'.  
        make[1]: Nothing to be done for `default'.  
        make[1]: Nothing to be done for `default'.  

    or something similar.

1. To run the code, type:

        ./run.py example.ini

    where exmaple.ini is the exmaple of input file.
    It contains the input parameters of simulation e.g. cosmological parameters, survey geometries, etc.  
    You can also run each code independently, for example:  

        ./calculate_pk/calc_pk_const_los_ngp   

    In this case code will ask you the input parameters.

1. The resulting files are as follows:
    + Files in data/inputs  
        _Note_: If you set gen_inputs = False, calculation of following files will be skipped.
        + exmaple_pk.txt
            + Input power spectrum
            + 1st column: wavenumber [h/Mpc]
            + 2nd column: power spectrum [Mpc^3/h^3]
            + _Note_: you can also use the power spectrum computed by external code (e.g. CAMB) by specifying the file name in .ini file
        + example_fnu.txt
            + logarithmic growth rate fnu for massive (massless) neutrinos which varies (does not vary) as a function of k
            + 1st column: wavenumber [h Mpc^-1]
            + 2nd column: fnu
        + example_Rh_xi.txt
            + correlation function calculated from input Pk
            + 1st column: comoving separation [Mpc/h]
            + 2nd column: correlation function [dimensionless]
        + example_pkG.dat
            + power spectrum of the log-normal field, G
            + 1st column: wavenumber [h/Mpc]
            + 2nd column: galaxy power spectrum [Mpc^3/h^3]
        + example_mpkG.dat
            + pkG for matter density field
            + 1st column: wavenumber [h/Mpc]
            + 2nd column: matter power spectrum [Mpc^3/h^3]
        + example_cpkG.dat
            + galaxy-matter cross pkG
            + 1st column: wavenumber [h/Mpc]
            + 2nd column: matter power spectrum [Mpc^3/h^3]

    + Files in data/lognormal  
    For each realization, the code generates the followng two binary files:
        + example\_lognormal\_rlz{# of realization}.bin
            Those files contain positions and velocities of galaxies. The larger "Pnmax", the larger the memory that the code requires. Start with a smaller value such as 512 when testing the code. To read them, use a sample code "lognormal/read_lognormal.f90". The byte-by-byte description of the file is:
			
			|Bytes|Format|Explanations                |
			|:---:|:----:|:--------------------------:|
			|1-8  |double|box length in x [Mpc/h]     |
			|9-16 |double|box length in y [Mpc/h]     |
			|17-24|double|box length in z [Mpc/h]     |
			|25-28|int   |number of galaxies          |
			|29-32|float |galaxy position in x [Mpc/h]|
			|33-36|float |galaxy position in y [Mpc/h]|
			|37-40|float |galaxy position in z [Mpc/h]|
			|41-44|float |galaxy velocity in x [km/s] |
			|45-48|float |galaxy velocity in y [km/s] |
			|49-52|float |galaxy velocity in z [km/s] |
			|53-  |float |position and velocity of $i$-th galaxy|

        + example\_density\_lognormal\_rlz{# of realization}.bin  
            The matter density in each cell is stored in those files. The byte-by-byte description of the file is:
					
			|Bytes|Format|Explanations                |
			|:---:|:----:|:--------------------------:|
			|1-8  |double|box length in x [Mpc/h]     |
			|9-16 |double|box length in y [Mpc/h]     |
			|17-24|double|box length in z [Mpc/h]     |
			|25-28|int   |mesh size in x              |
			|29-32|int   |mesh size in y              |
			|33-36|int   |mesh size in z              |
			|37-  |double|matter density in each cell |

    + Files in data/pk
        + exmaple\_pk\_rlz{# of realization}.dat  
            The power spectrum calculated from each mock galaxy catalog.
            + 1st column: wavenumber [h/Mpc]
            + 2nd column: monopole power spectrum [Mpc^3/h^3]
            + 3rd column: quadrupole power spectrum [Mpc^3/h^3]
            + 4th column: hexadecapole power spectrum [Mpc^3/h^3]
            + 5th column: the number of Fourier modes averaged

            When all of los parameters, losx, losy and losz, are equal to 0, the real-space power spectra are calculated. When one of them is equal to 1, the positions of galaxies are shifted by the velocity in that direction. __AT MOST__ one of them should be equal to 1 - do not set more than one components to be 1! 


1. To check the results, compute some number of real-space monopole power spectra, average them, and compare the average with the input matter power spectrum times the bias squared. 
You may find a deviation on small scales due to the resolution effect of Fourier meshes, but it should reproduce the input power spectrum very precisely on large scales.The first few bins may return "nan", but do not worry about them.  
_NOTE_: In the cubic mode (calc\_mode\_pk == 1), the estimated power spectrum is automatically convolved with the survey geometry function Wk, thus you should take into account the effect of Wk when comparing it with the input matter power spectrum. 


1. To compute the predictions for the monopole and quadrupole of the redshift space power spectrum in Kaiser limit, the cross correlation function and the cross power spectrum between generated galaxy and matter density fields are needed. We also need to average the theoretical power spectrum over the same Fourier grids as in the measurements. To do this, there are Python wrappers such as "run\_calc\_pk.py", "run\_calc\_xi\_gm.py", "run\_discretize\_pk.py", and "run\_kaiser\_pk.py". The details of how to run these codes are described in a readme file in "aux\_codes".  
A script to generate the Kaiser prediction using all of the codes in aux_codes is provided as "generate\_kaiser.sh". The script produces a cross correlation function file (xi\_gm\_kmax\*), a cross power spectrum file (pk\_rmax\*\_b1.dat), a Kaiser P(k) file (pk\_kaiser.dat) containing 3 columns -- k, P\_0(k), P\_2(k) and discretized versions of the Kaiser P(k) - pkt\_kaiser\_2.dat (discretized monopole) and pkt\_kaiser\_3.dat (discretized quadrupole).  
_NOTE_: The script assumes filenames and other parameters consistent with the other provided scripts. Please make necessary changes if changing the default scripts.
