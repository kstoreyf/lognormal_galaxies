*** Linear matter density power spectrum from Eisenstein&Hu's transfer function with the baryonic oscillation and massive neutrinos ***
November 23, 2015: E.Komatsu => compute_pk.f90
March 4, 2016: modified by Aniket Agrawal to output the scale-dependent growth rates => compute_pk2.f90


The linear matter power spectrum is given by [E.g., Eq.(74) of Komatsu et al., ApJS, 180, 330 (2009)]

P(k,z) = Delta_R^2*(2*k^2/5/H0^2/Omega_M)^2
        *D(k,z)^2*T(k)^2*(k/k0)^(ns-1)
        *(2*pi^2)/k^3
where

- Delta_R^2 is the amplitude of the curvature perturbation in comoving gauge (or uniform curvature gauge) on superhorizon scale (e.g., Delta_R^2=2.46d-9)
- D(k,z) is the linear growth factor normalized such that (1+z)D(z)->1 during the matter era. It depends on k because of massive neutrinos, as computed by Hu & Eisenstein, 498, 497 (1998)
- T(k) is the linear transfer function. It combines the effects of massive neutrinos and baryon acoustic oscillation as follows:

  T(k) = [trans_bao(k) - trans_nowiggle(k)] + trans_nu(k)

where "trans_bao(k)" and "trans_nowiggle(k)" are the transfer functions with MASSLESS neutrinos. The former incluldes BAO, whereas the latter does not. Both are given in Eisenstein & Hu, ApJ, 496, 605 (1998). "trans_nu(k)" is the transfer function with massive neutrinos, but does not include BAO. It is given in Eisenstein & Hu, ApJ, 511, 5 (1999). 

- k0 is the wavenumber at which the initial power spectrum is 
normalized (e.g., k0=0.002/Mpc)

*** OUTPUT Of COMPUTE_PK (wavenumber_pk.txt) ***
1st column: wavenumber [h Mpc^-1]
2nd column: power spectrum [h^-3 Mpc^3]

- To compile and use the program, edit Makefile and simply "make"
- It will generate "compute_pk"

A second code "compute_pk2.f90" is also provided, which generates 4 output files - wavenumber_pk.txt, wavenumber_gnu.txt, wavenumber_dgnudlna.txt, wavenumber_fnu.txt. The format of wavenumber_pk.txt is the same as above. 

*** OUTPUT Of COMPUTE_PK2 (wavenumber_gnu.txt) ***
1st column: wavenumber [h Mpc^-1]
2nd column: growth rate (gnu) for massless neutrinos (independent of wavenumber)
3rd column: growth rate (gnu) for massive neutrinos

*** OUTPUT Of COMPUTE_PK2 (wavenumber_dgnudlna.txt) ***
1st column: wavenumber [h Mpc^-1]
2nd column: derivative of growth rate (dgnu/dlna) for massless neutrinos (independent of wavenumber)
3rd column: growth rate (dgnu/dlna) for massive neutrinos

*** OUTPUT Of COMPUTE_PK2 (wavenumber_fnu.txt) ***
1st column: wavenumber [h Mpc^-1]
2nd column: logarithmic growth rate - fnu - for massive (massless) neutrinos which varies (does not vary) as a function of 1st column
