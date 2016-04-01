c	Subroutines extracted from Wayne Hu's "power.f" available at
c     http://background.uchicago.edu/~whu/transfer/power.f

       subroutine TFset_parameters_nu(omhh,f_baryon,f_nu,N_nu,Tcmb)
	implicit none
        real*8 omhh,f_baryon,f_nu,N_nu,Tcmb
        real*8 y_d,alpha_nu,
     *     beta_c,sound_horizon,theta_cmb,omega,omegal,z_equality
        real*8 obhh,k_equality,z_drag,R_drag,R_equality,p_c,p_cb,f_c,
     *     f_cb,f_nub
        common/GLOBALVARIABLESNU/y_d,alpha_nu,
     *     beta_c,sound_horizon,theta_cmb,omega,omegal,z_equality

c Auxiliary variable

        obhh = omhh*f_baryon
        theta_cmb = Tcmb/2.7

c Main variables

        z_equality = 2.50e4*omhh*theta_cmb**(-4.) - 1.D0
        k_equality = 0.0746*omhh*theta_cmb**(-2.)

          z_drag = 0.313*omhh**(-0.419)*(1.+0.607*omhh**(0.674))
          z_drag = 1e0 + z_drag*obhh**(0.238*omhh**(0.223))
        z_drag = 1291e0 * omhh**(0.251)/
     &           (1e0 + 0.659*omhh**(0.828)) * z_drag

	y_d = (1.+z_equality)/(1.+z_drag)

        R_drag = 31.5*obhh*theta_cmb**(-4.)*1000e0/(1e0 + z_drag)
        R_equality = 31.5*obhh*theta_cmb**(-4.)
     &               *1000e0/(1e0 + z_equality)

        sound_horizon = 2./3./k_equality*sqrt(6./R_equality)*
     &      log(( sqrt(1.+R_drag)+sqrt(R_drag+R_equality) )
     &       /(1.+sqrt(R_equality)))

	 p_c  = -(5.-sqrt(1.+24*(1.-f_nu-f_baryon)))/4.
	 p_cb = -(5.-sqrt(1.+24*(1.-f_nu)))/4.
	 f_c  = 1.-f_nu-f_baryon
	 f_cb = 1.-f_nu
	 f_nub= f_nu+f_baryon


	 alpha_nu= (f_c/f_cb)* (2.*(p_c+p_cb)+5.)/(4.*p_cb+5)
	 alpha_nu= alpha_nu*(1.-0.553*f_nub+0.126*f_nub**3)
	 alpha_nu= alpha_nu/(1.-0.193*sqrt(f_nu)+0.169*f_nu)
	 alpha_nu= alpha_nu*(1.+y_d)**(p_c-p_cb)
	 alpha_nu= alpha_nu*(1.+ (p_cb-p_c)/2. *
     *			(1.+1./(4.*p_c+3.)/(4.*p_cb+7.))/(1.+y_d))
	 beta_c=1./(1.-0.949*f_nub)

	 return

	end 


	real*8 function TF_master(k,omhh,f_baryon,f_nu,N_nu)
        implicit none
	real*8 omhh,f_nu,f_baryon,N_nu
        real*8 y_d,alpha_nu,
     *     beta_c,sound_horizon,theta_cmb,omega,omegal,z_equality
        real*8 k,q,gamma_eff,q_eff,q_nu
        common/GLOBALVARIABLESNU/y_d,alpha_nu,
     *     beta_c,sound_horizon,theta_cmb,omega,omegal,z_equality
	 q = k*theta_cmb**2/omhh
	 gamma_eff=(sqrt(alpha_nu) + (1.-sqrt(alpha_nu))/
     *		(1.+(0.43*k*sound_horizon)**4))

	 q_eff = q/gamma_eff
	 TF_master= dlog(dexp(1.d0)+1.84*beta_c*sqrt(alpha_nu)*q_eff)
	 TF_master = TF_master/(TF_master + q_eff**2*
     *	             (14.4 + 325./(1.+60.5*q_eff**1.11)))

	 q_nu = 3.92*q*sqrt(N_nu/f_nu)
	 TF_master = TF_master*
     *     (1.+(1.2*f_nu**(0.64)*N_nu**(0.3+0.6*f_nu))/
     *	   (q_nu**(-1.6)+q_nu**(0.8)))

	 return 

	end
