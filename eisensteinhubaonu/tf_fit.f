c   The following routines implement all of the fitting formulae in 
c   Eisenstein \& Hu (1997) 
c
c   Program TF_fit: sample driver
c   Subroutine TFset_parameters(): sets all the scalar parameters
c   Subroutine TFtransfer_function(): calculates various transfer functions
c   Functions TF_zerobaryon, TF_nowiggles, sound_horizon_fit, kpeak,
c       alpha_gamma: implement various scaling approximations of \S 4.2
c
c
c ------------------------ DRIVER ROUTINE --------------------------- 
c The following is a driver routines you might use. 
c Basically, the driver routine asks for Omega_0, the baryon fraction,
c the hubble constant, and the CMB temperature, calls TFset_parameters() to
c set all the parameters of the fit.  A loop over wavenumbers kmin to an
c inputed kmax sampled at numk per decade calls TFtransfer_function. 
c
c IMPORTANT: TFtransfer_function asks for wavenumbers in Mpc^{-1} so
c	     multiply by hubble to convert from h Mpc^{-1}
c
c The latter returns values of the various transfer functions at the given 
c wavenumber which are output to the file "trans.dat" 
c
c Also included is an example of how to call the functions
c
c	TF_nowiggles:	    shape approximation
c	TF_zerobaryon:      zero baryon TF 
c	k_peak:		    approximate first peak location
c	sound_horizon_fit:  approximate sound horizon
c	alpha_gamma:  	    small scale modification to Gamma
c
c
c
c INPUT:  omega0 -- the matter density (baryons+CDM) in units of critical 
c	  f_baryon -- the ratio of baryon density to matter density 
c	  hubble -- the Hubble constant, in units of 100 km/s/Mpc
c	  Tcmb -- the CMB temperature in Kelvin.  2.728(4) is COBE and is the default
c	          reached by setting Tcmb=0.
c	  kmax -- maximum k in  h Mpc^{-1}
c	  numk -- number of k per decade 
c
c OUTPUT: The file "trans.dat" with columns
c
c	 (1)  k -- wavenumber in h Mpc^{-1}
c	 (2)  tf_full -- The full fitting formula, eq. (16), for the matter
c			transfer function. 
c	 (3)  tf_baryon -- The baryonic piece of the full fitting formula, eq. 21.
c	 (4)  tf_cdm -- The CDM piece of the full fitting formula, eq. 17.
c	 (5)  tf_nowiggles -- An approximate form, eqs. (30)-(31), that fits
c			only the non-oscillatory part of the transfer 
c			function.  Appropriate only for low baryon fractions.
c	 (6)  tf_zerobaryon -- The transfer function of the zero-baryon case,
c			eq. (29); i.e. what would have occured were the
c			baryons CDM instead. 			
c
c	and the approximate values of k_peak,sound_horizon,alpha_gamma 
c	to stdout, a more accurate form of sound_horizon lives in
c
c	GLOBALVARIABLES: Various intermediate fit parameters are stored in 
c                        this common block for easy access.


c ----------------------------- DRIVER ------------------------------- 


c$$$       program TFfit
c$$$
c$$$
c$$$        real    omega0,f_baryon,hubble,Tcmb,kmax,kmin
c$$$	real	k,tf_full,tf_baryon,tf_cdm,tf_nowiggles,tf_zerobaryon,
c$$$     *		k_peak
c$$$	integer numk
c$$$
c$$$c  cosmological parameters
c$$$
c$$$	write(6,*) 'Omega_0,f_baryon,h,T_cmb?'
c$$$	read*,     omega0,f_baryon,hubble,Tcmb
c$$$	omhh = omega0*hubble*hubble
c$$$
c$$$c  call routine to set fitting parameters
c$$$
c$$$        call TFset_parameters(omhh, f_baryon, Tcmb)
c$$$
c$$$
c$$$c  loop over k, call subroutine and functions to calc TFs
c$$$ 
c$$$	open(10,file='trans.dat')
c$$$
c$$$	write(6,*) 'k_max (h Mpc^{-1}),#pts per decade (10,50)'
c$$$	read*,kmax,numk
c$$$
c$$$	if (kmax.le.0) kmax=10.
c$$$	if (numk.le.0) numk=50
c$$$
c$$$	kmin = 0.0001
c$$$	numk = numk*log10(kmax/kmin)
c$$$
c$$$	do i=1,numk
c$$$
c$$$	 k=10.**(i*(log10(kmax/kmin)/numk))*kmin
c$$$         call TFtransfer_function(k*hubble,omhh,f_baryon,
c$$$     &     tf_full,tf_baryon,tf_cdm)
c$$$	 write(10,50) k,tf_full,tf_baryon,tf_cdm,
c$$$     &               TF_nowiggles(k*hubble,omhh,f_baryon,Tcmb),
c$$$     &               TF_zerobaryon(k*hubble/omhh*(Tcmb/2.7)**2)
c$$$
c$$$
c$$$	end do
c$$$
c$$$c  example of how to use the scaling functions
c$$$
c$$$	write(6,*) 'Some useful approximate scalings:'
c$$$
c$$$	write(6,10) k_peak(omhh,f_baryon)/hubble
c$$$10	FORMAT(1X,' First peak location (h Mpc^{-1}):  ',E13.4)
c$$$	write(6,20) sound_horizon_fit(omhh,f_baryon)*hubble
c$$$20	FORMAT(1X,' Approx. sound horizon (h^{-1} Mpc):',E13.4)
c$$$	write(6,30) alpha_gamma(omhh,f_baryon)
c$$$30	FORMAT(1X,' alpha_gamma:                       ',E13.4)
c$$$
c$$$
c$$$50	FORMAT(1X,7E13.5)
c$$$
c$$$      end

c
c
c PART I:------------------- FITTING FORMULAE ROUTINES ----------------- 
c
c There are two routines and a set of functions.  
c   TFset_parameters() sets all the scalar parameters, while 
c   TFtransfer_function() calculates various transfer functions 
c
c Global variables -- We've left many of the intermediate results as
c global variables in case you wish to access them, e.g. by declaring
c them as a common block in your main program. 
c
c Note that all internal scales are in Mpc, without any Hubble constants! 
c

	subroutine TFset_parameters(omhh,f_baryon,Tcmb)
        implicit none
	real*8 y,omhh,obhh,Tcmb
	real*8 theta_cmb,z_equality,k_equality,z_drag,R_drag,R_equality,
     &       sound_horizon,k_silk,alpha_c,beta_c,alpha_b,beta_b,
     &	     f_baryon,beta_node
	common/GLOBALVARIABLES/theta_cmb,z_equality,k_equality,z_drag,
     &       R_drag,R_equality,sound_horizon,k_silk,alpha_c,beta_c,
     &       alpha_b,beta_b,beta_node

c Set all the scalars quantities for Eisenstein & Hu 1997 fitting formula */
c Input omhh -- The density of CDM and baryons, in units of critical dens,
c                multiplied by the square of the Hubble constant, in units
c                of 100 km/s/Mpc */
c       f_baryon -- The fraction of baryons to CDM */
c       Tcmb -- The temperature of the CMB in Kelvin, 2.728(4) is COBE and is
c		the default reached by inputing Tcmb=0 -- reset on output. */
c Output nothing, but set many global variables in common block 
c       GLOBALVARIABLES. You can access them yourself, if you want:
c
c	theta_cmb,	/* Tcmb in units of 2.7 K */ 
c	z_equality,	/* Redshift of matter-radiation equality, really 1+z */
c	k_equality,	/* Scale of equality, in Mpc^-1 */
c	z_drag,		/* Redshift of drag epoch */
c	R_drag,		/* Photon-baryon ratio at drag epoch */
c	R_equality,	/* Photon-baryon ratio at equality epoch */
c	sound_horizon,	/* Sound horizon at drag epoch, in Mpc */
c	k_silk,		/* Silk damping scale, in Mpc^-1 */
c	alpha_c,	/* CDM suppression */
c	beta_c,		/* CDM log shift */
c	alpha_b,	/* Baryon suppression */
c	beta_b,		/* Baryon envelope shift */



c Are inputs reasonable?
	if (f_baryon.le.0) f_baryon=1.e-5
	if (Tcmb.le.0) Tcmb=2.728
        if (omhh.le.0.0) then
	   write(6,*) 'TFset_parameters(): Illegal input'  
	   pause
	end if

c        if (hubble.gt.10.0) then
c	   write(6,*) 'TFset_parameters(): WARNING, Hubble constant in 
c     &	               100km/s/Mpc desired'
c	end if

c Auxiliary variables
        obhh = omhh*f_baryon
        theta_cmb = Tcmb/2.7

c Main variables
        z_equality = 2.50e4*omhh*theta_cmb**(-4.) - 1.D0
        k_equality = 0.0746*omhh*theta_cmb**(-2.) 

          z_drag = 0.313*omhh**(-0.419)*(1.+0.607*omhh**(0.674))
          z_drag = 1e0 + z_drag*obhh**(0.238*omhh**(0.223))
        z_drag = 1291e0 * omhh**(0.251)/
     &           (1e0 + 0.659*omhh**(0.828)) * z_drag
 
        R_drag = 31.5*obhh*theta_cmb**(-4.)*1000e0/(1e0 + z_drag) 
        R_equality = 31.5*obhh*theta_cmb**(-4.) 
     &    	     *1000e0/(1e0 + z_equality) 

        sound_horizon = 2./3./k_equality*sqrt(6./R_equality)*
     &	    log(( sqrt(1.+R_drag)+sqrt(R_drag+R_equality) )
     &       /(1.+sqrt(R_equality)))

        k_silk = 1.6*obhh**(0.52)*omhh**(0.73)* 
     &           (1e0 + (10.4*omhh)**(-0.95))

          alpha_c = ((46.9*omhh)**(0.670)*(1e0+(32.1*omhh)**(-0.532)))
          alpha_c = alpha_c**(-f_baryon) 
	alpha_c = alpha_c*((12.0*omhh)**(0.424)*(1e0 + 
     &             (45.0*omhh)**(-0.582)))**(-f_baryon**3.)

    
          beta_c = 0.944/(1+(458.*omhh)**(-0.708))
          beta_c = 1.+beta_c*((1.-f_baryon)**((0.395*omhh)**(-0.0266)) 
     &    	- 1e0)
	  beta_c = 1./beta_c

          y = (1e0+z_equality)/(1e0+z_drag)
          alpha_b = y*(-6.*sqrt(1.+y)+(2.+3.*y)*log((sqrt(1.+y)+1.)
     &   	    /(sqrt(1.+y)-1.)))
        alpha_b = 2.07*k_equality*sound_horizon*
     &            (1.+R_drag)**(-0.75)*alpha_b


        beta_b = 0.5+f_baryon+(3.-2.*f_baryon)*
     &           sqrt((17.2*omhh)**2.+1e0)

        beta_node = 8.41*omhh**(0.435)

        return

        end


        subroutine TFtransfer_function(k,omhh,f_baryon,tf_full,
     &                tf_baryon,tf_cdm)
        implicit none
c  Calculate transfer function from the fitting parameters stored in
c  GLOBALVARIABLES.
c
c  Input: 
c	 k -- wavenumber in Mpc^{-1}  
c        omhh -- The density of CDM and baryons, in units of critical dens,
c                multiplied by the square of the Hubble constant, in units
c                of 100 km/s/Mpc */
c        f_baryon -- The fraction of baryons to CDM */
c	
c  Output:
c	 tf_full -- The full fitting formula, eq. (16), for the matter
c	            transfer function. 
c	 tf_baryon -- The baryonic piece of the full fitting formula, eq. 21.
c	 tf_cdm -- The CDM piece of the full fitting formula, eq. 17.
c

        real*8 TF_pressureless,s_tilde
	real*8 k,tf_full,tf_baryon,tf_cdm,q,ks
	real*8 theta_cmb,z_equality,k_equality,z_drag,R_drag,R_equality,
     &       sound_horizon,k_silk,alpha_c,beta_c,alpha_b,beta_b,
     &	     f_baryon,beta_node,omhh
	common/GLOBALVARIABLES/theta_cmb,z_equality,k_equality,z_drag,
     &       R_drag,R_equality,sound_horizon,k_silk,alpha_c,beta_c,
     &       alpha_b,beta_b,beta_node


c  Reasonable k?

        if (k.le.0) then
           write(6,*) 'TFtransfer_function(): Illegal k'
           pause
        end if 


c  Auxiliary Variables

	    q = k/13.41/k_equality
	    ks = k*sound_horizon

c  Main Variables

              tf_cdm = 1./(1.+(ks/5.4)**4.)
	    tf_cdm = tf_cdm*TF_pressureless(q,1d0,beta_c) + 
     &           (1.-tf_cdm)*TF_pressureless(q,alpha_c,beta_c)


	      s_tilde = sound_horizon/(1.+(beta_node/ks)**3.)**(1./3.) 
	      tf_baryon = TF_pressureless(q,1d0,1d0)/(1d0+(ks/5.2d0)**2d0)
	      tf_baryon = tf_baryon + alpha_b/(1.+(beta_b/ks)**3)
     &                       *exp(-(k/k_silk)**(1.4))
	      tf_baryon = tf_baryon *(sin(k*s_tilde)/(k*s_tilde))
	    tf_full = f_baryon*tf_baryon + (1-f_baryon)*tf_cdm

         return

	 end

c       auxiliary function: Pressureless TF

	double precision function TF_pressureless(q,a,b)
        implicit none
	  real*8 q,a,b
	
	  TF_pressureless = Log(exp(1.)+1.8*b*q)
	  TF_pressureless = TF_pressureless/(TF_pressureless + 
     &                      (14.2/a + 386/(1.+69.9*q**1.08))*q**2)

	  return

	end	
c
c
c
c
c
c
c PART II:------------------- Scaling Functions ROUTINES ----------------- 
c
c       omhh -- The density of CDM and baryons, in units of critical dens,
c                multiplied by the square of the Hubble constant, in units
c                of 100 km/s/Mpc */
c       f_baryon -- The fraction of baryons to CDM */
c
c
c	TF_zerobaryon:     
c	  Input:  q = k/omhh * (Tcmb/2.7)**2    (k in Mpc^{-1})
c	  Output: zero baryon TF Eq(29)
c	TF_nowiggles:      
c	  Input:  k = wavenumber in Mpc^{-1}, omhh, f_baryon, Tcmb
c	  Output: shape approximation TF  Eq(30-31)
c	  Calls: TF_zerobaryon,sound_horizon_fit,alpha_gamma
c 	sound_horizon_fit: 
c         Input:  omhh,f_baryon	
c	  Output: approximate sound horizon in Mpc	
c	kpeak:		   
c	  Input:  omhh,f_baryon
c         Output: first peak location in Mpc
c	  Calls:  sound_horizon_fit
c	alpha_gamma:	   
c	  Input: omhh,f_baryon
c	  Output: effective small scale suppression


	double precision function TF_zerobaryon(q)
        implicit none
	  real*8 q
	  TF_zerobaryon = log(2.0*exp(1.)+1.8*q)
	  TF_zerobaryon = TF_zerobaryon/(TF_zerobaryon
     &                    +(14.2 + 731.0/(1+62.5*q))*q**2)

	  return

	end

	double precision function TF_nowiggles(k,omhh,f_baryon,Tcmb)
        implicit none
	  real*8 k,omhh,f_baryon,q_eff,a,Tcmb
          real*8 alpha_gamma,sound_horizon_fit,TF_zerobaryon
	
	    if (Tcmb.le.0) Tcmb=2.728
	    a = alpha_gamma(omhh,f_baryon)
	    q_eff = k/omhh*(Tcmb/2.7)**2
	    q_eff = q_eff/(a+(1.-a)/
     &              (1.+(0.43*k*sound_horizon_fit(omhh,f_baryon))**4))

	    TF_nowiggles = TF_zerobaryon(q_eff)
	  
	  return
        end


	double precision function sound_horizon_fit(omhh,f_baryon)
        implicit none
	 real*8 omhh,obhh,f_baryon
	 obhh = f_baryon*omhh
         sound_horizon_fit = 44.5*log(9.83/omhh)
     &                      /sqrt(1.+10.0*obhh**(0.75))

	return
	end


	double precision function k_peak(omhh,f_baryon)
        implicit none
	 real*8 omhh,obhh,f_baryon
         real*8 sound_horizon_fit
	 obhh = f_baryon*omhh
         k_peak = 5.*3.14159/2.*(1.+0.217*omhh)/
     &             sound_horizon_fit(omhh,f_baryon)

	 return
	end


	double precision function alpha_gamma(omhh,f_baryon)
        implicit none
	 real*8 omhh,f_baryon

         alpha_gamma = 1.-0.328*log(431.0*omhh)*f_baryon 
     &                + 0.38*log(22.3*omhh)*(f_baryon)**2
    
	 return
	end 
