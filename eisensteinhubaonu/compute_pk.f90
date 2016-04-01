PROGRAM Compute_Pk
  USE cosmo
  USE growth
  ! A sample program for computing the linear power spectrum
  ! of density fluctuations using Eisenstein&Hu's transfer function
  ! that includes the effects of massive neutrinos and
  ! the baryonic oscillation, multiplied by the growth function
  ! that includes the effect of massive neutrinos. 
  ! The modification due to massive neutrinos is given by subroutines 
  ! extracted from Wayne Hu's "power.f" (given in tf_fit_nu.f), and 
  ! the growth  function fit derived in Hu & Eisenstein (1998).
  ! The baryon acoustic oscillation is included by computing the difference
  ! with respect to the no-wiggle power spectrum for massless neutrinos.
  ! Ref: Eisenstein & Hu, ApJ, 496, 605 (1998) for transfer functions with massless nu
  !      Eisenstein & Hu, ApJ, 511, 5 (1999) for no-wiggle transfer function with massive nu
  ! - k is in units of h Mpc^-1
  ! - P(k) is in units of h^-3 Mpc^3
  ! November 23, 2015: E.Komatsu
  IMPLICIT none
  integer :: j
  double precision :: g,z,D1,D,trans,trans_nu,trans_nowiggle
  double precision :: k_ov_h,pk,dlnk,lnk,zin
  double precision :: ob0,h0,ns,run,deltaR2
  double precision :: tf_baryon,tf_cdm
  double precision :: mnu,Nnu,onu0,fnu,pcb,yfs
  double precision :: q,zeq
  double precision :: TF_master
  external g,TF_master
  character(len=128) :: filename
  integer :: n
! Specify three cosmological parameters
! The data type has been defined in MODULE cosmo.
  ode0=0.728d0
  om0=0.272d0
  w=-1d0
! Specify four more cosmological parameters
! These are not defined in MODULE cosmo.
  ob0=0.0456d0
  h0=0.704d0
  ns=0.963d0
  run=0d0
  deltaR2=2.46d-9*(0.809/0.8140756)**2d0
! neutrino parameters
  mnu=0.2d0 ! the total neutrino mass in units of eV
  Nnu=3.046d0 ! effective number of neutrino species
  onu0=(mnu/93.1d0)/h0**2d0 ! Omega_nu of massive neutrinos today
  fnu=onu0/om0 ! neutrino fraction in mass
  pcb=(5d0-dsqrt(25d0-24d0*fnu))/4d0
  zeq=2.5d4*om0*h0**2d0*(2.7d0/2.726d0)**4d0
  print*,'cosmological parameters:'
  print*,'Omega_m [incl. baryons and neutrinos]=',om0
  print*,'Omega_l=',ode0
  print*,'Omega_b=',ob0
  print*,'h=',h0
  print*,'w=',w
  print*,'ns=',ns
  print*,'dns/dlnk=',run
  print*,'deltaR2=',deltaR2
  print*,'total mass of neutrinos=',mnu,' eV'
! tabulate g(z) by calling "setup_growth"
  call setup_growth
! ask for redshift
  print*,'redshift?'
  read*,z
  D1=g(z)*(1d0+zeq)/(1d0+z) ! linear growth factor, normalized such that (1+z)D(z)=1+zeq during the matter era
! now output P(k,z) as a function of k
  open(1,file='wavenumber_pk.txt')
  k_ov_h=1d-5 ! h/Mpc
  dlnk=2d-2
  CALL TFset_parameters(om0*h0**2d0,ob0/om0,2.726d0)
  CALL TFset_parameters_nu(om0*h0**2d0,ob0/om0,fnu,Nnu,2.726d0)
  do while (k_ov_h<1d4)
     CALL eisensteinhu(k_ov_h*h0,om0*h0**2d0,ob0/om0,trans_nowiggle) ! no-wiggle P(k) with massless neutrinos
     trans_nu=TF_master(k_ov_h*h0,om0*h0**2d0,ob0/om0,fnu,Nnu) ! no-wiggle P(k) with massive neutrinos
     CALL TFtransfer_function(k_ov_h*h0,om0*h0**2d0,ob0/om0,trans,tf_baryon,tf_cdm) ! P(k) with BAO and massless neutrinos
     trans=(trans-trans_nowiggle)+trans_nu ! add BAO to no-wiggle P(k) with massive neutrinos
     ! modify growth by neutrinos
     q=k_ov_h/(om0*h0)*(2.726d0/2.7d0)**2d0
     yfs=17.2d0*fnu*(1d0+0.488d0/fnu**(7d0/6d0))*(Nnu*q/fnu)**2d0
     D=((1d0-fnu)**(0.7d0/pcb)+(D1/(1d0+yfs))**0.7d0)**(pcb/0.7d0)*D1**(1d0-pcb)
     D=D/(1d0+zeq)
     ! Eq.(74) of Komatsu et al., ApJS, 180, 330 (2009) with kWMAP=0.002/Mpc
     ! Remember that k_ov_h is in units of h/Mpc whereas "k" in Eq.(74) is in units 
     ! of 1/Mpc.
     pk=deltaR2*(2d0*k_ov_h**2d0*2998d0**2d0/5d0/om0)**2d0*D**2d0 &
       *trans**2d0*(k_ov_h*h0/0.002d0)**(ns-1d0+0.5d0*run*dlog(k_ov_h*h0/0.002d0)) &
       *2d0*3.14159d0**2d0/k_ov_h**3d0
     write(1,'(2E18.8)')k_ov_h,pk
     lnk=dlog(k_ov_h)+dlnk
     k_ov_h=dexp(lnk)
  enddo
  close(1)
END PROGRAM Compute_Pk
