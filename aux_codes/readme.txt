***Some codes to calculate various quantities***
March 3, 2016: Aniket Agrawal

Some codes in this folder use spline_array.h and integration_modified.o which are present in the generate_Poisson and compute_pkG folders. The following codes are provided in this folder : 

1) calc_xi_gm : This code calculates the cross correlation function of matter and galaxy for the fields generated in the lognormal code. The cross power spectrum is useful in computing the cross-correlation coefficient and the Kaiser prediction for the power spectrum.

***INPUT***
The code takes the Gaussian power spectra of galaxy and matter (which should have been calculated at the time of generating the fields) and a value of kmax for the integration. The provided version assumes that the input power spectra are specified at the exact same k values.

***OUTPUT (xi_gm_kmax*.dat)***
1st column : separation [Mpc/h]
2nd column : cross correlation function

2) calc_pk : This code calculates the power spectrum from an input correlation function by calculating

P(k) = 4*pi*\int_0_rmax xi(r)*r^2*j0(kr) dr 

where rmax can be set by the user. This code can be used to calculate, for example, the cross power spectrum from the cross correlation function. 

3) kaiser_pk : This code calculates the Kaiser prediction for the monopole and quadrupole power spectrum. 

***INPUT***
The code takes as input the matter power spectrum, galaxy-matter cross power spectrum, linear bias, f(k) file

***OUTPUT (pk_kaiser.dat)***
1st column : wavenumber [h/Mpc]
2nd column : monopole Kaiser power spectrum [h^-3 Mpc^-3]
3rd column : quadrupole Kaiser power spectrum [h^-3 Mpc^-3]

4) discretize_pk : This code calculates the discretized version of an input power spectrum on the standard Fourier grid (set by choosing kbin=0.0 while measuring P(k)) for a given box. It can be modified to also include non-standard Fourier grids. The smoothing scale in this code refers to the scale used in the input power spectrum to generate lognormal density fields (please set to 0.0 if no smoothing used), P(k)=P_0(k)*exp(-k^2*Rs^2).

Sample codes (run_*.py) to run these codes are provided as Python scripts the top directory in codes_lognormal. To execute use
python run_*.py
