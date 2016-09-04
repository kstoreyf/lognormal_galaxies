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
  ! March 4, 2016: modified by Aniket Agrawal to output the scale-dependent growth factors
  IMPLICIT none
  integer :: j
  double precision :: g,dgdlna,dlngdlna,growth_nu,z,D1,D,trans,trans_nu,trans_nowiggle   !modified by Aniket
  double precision :: gnu,dgnudlna   !modified by Aniket
  double precision :: k_ov_h,pk,dlnk,lnk,zin
  double precision :: ob0,h0,ns,run,deltaR2
  double precision :: tf_baryon,tf_cdm
  double precision :: mnu,Nnu,onu0,fnu,pcb,yfs
  double precision :: q,zeq
  double precision :: TF_master
  external g,dgdlna,TF_master   !modified by Aniket
  character(len=128) :: ofile_prefix, arg
  integer :: n

  print *, "Calculate the linear power spectrum using Eisenstein & Hu's transfer function"

  call ReadParams

  ! set parameters
  run = 0.d0
  Nnu=3.046d0 ! effective number of neutrino species
  onu0=(mnu/93.1d0)/h0**2d0 ! Omega_nu of massive neutrinos today
  fnu=onu0/om0 ! neutrino fraction in mass
  pcb=(5d0-dsqrt(25d0-24d0*fnu))/4d0
  zeq=2.5d4*om0*h0**2d0*(2.7d0/2.726d0)**4d0

  ! tabulate g(z) by calling "setup_growth"
  call setup_growth
  ! linear growth factor, normalized such that (1+z)D(z)=1+zeq during the matter era
  D1=g(z)*(1d0+zeq)/(1d0+z) ! now output P(k,z) as a function of k

  ! open output files
  open(1,file='./data/inputs/'//trim(ofile_prefix)//'_pk.txt')
  open(2,file='./data/inputs/'//trim(ofile_prefix)//'_gnu.txt') 
  open(3,file='./data/inputs/'//trim(ofile_prefix)//'_dgnudlna.txt')
  open(4,file='./data/inputs/'//trim(ofile_prefix)//'_fnu.txt')

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
     if(fnu/=0d0)then   !modified by Aniket
        yfs=17.2d0*fnu*(1d0+0.488d0/fnu**(7d0/6d0))*(Nnu*q/fnu)**2d0   !modified by Aniket
     else          !modified by Aniket
        yfs=1d10   !modified by Aniket
     endif         !modified by Aniket
     D=((1d0-fnu)**(0.7d0/pcb)+(D1/(1d0+yfs))**0.7d0)**(pcb/0.7d0)*D1**(1d0-pcb)
     D=D/(1d0+zeq)
     ! Eq.(74) of Komatsu et al., ApJS, 180, 330 (2009) with kWMAP=0.002/Mpc
     ! Remember that k_ov_h is in units of h/Mpc whereas "k" in Eq.(74) is in units 
     ! of 1/Mpc.
     pk=deltaR2*(2d0*k_ov_h**2d0*2998d0**2d0/5d0/om0)**2d0*D**2d0 &
       *trans**2d0*(k_ov_h*h0/0.002d0)**(ns-1d0+0.5d0*run*dlog(k_ov_h*h0/0.002d0)) &
       *2d0*3.14159d0**2d0/k_ov_h**3d0
     dlngdlna=dgdlna(z)/g(z)  !modified by Aniket
     gnu=(1d0+z)*D            !modified by Aniket
     dgnudlna=-pcb+(1d0-pcb)*dlngdlna+pcb*(D1/(1d0+yfs))**0.7d0 &
          *(1d0+dlngdlna)/((1d0-fnu)**(0.7d0/pcb)+(D1/(1d0+yfs))**0.7d0)   !modified by Aniket
     dgnudlna=dgnudlna*gnu        !modified by Aniket
     growth_nu=1d0+dgnudlna/gnu   !modified by Aniket                                   
     write(1,'(2E18.8)')k_ov_h,pk
     write(2,'(3F18.8)')k_ov_h,g(z),gnu             !modified by Aniket
     write(3,'(3F18.8)')k_ov_h,dgdlna(z),dgnudlna   !modified by Aniket
     write(4,'(2E18.8)')k_ov_h,growth_nu            !modified by Aniket
     lnk=dlog(k_ov_h)+dlnk
     k_ov_h=dexp(lnk)
  enddo
  close(1)
  close(2)
  close(3)
  close(4)
  stop
!!$======================================================================
CONTAINS
    SUBROUTINE ReadParams
        IF(iargc() == 0) THEN
            print '(A,$)', '> enter the prefix of output file name: '
            read *, ofile_prefix
            print '(A,$)', '> Omega_m [incl. baryons and neutrinos] = '
            read *, om0
            print '(A,$)','> Omega_l = '
            read *, ode0
            print '(A,$)','> Omega_b = '
            read *, ob0
            print '(A,$)','> h = '
            read *, h0
            print '(A,$)','> w = '
            read *, w
            print '(A,$)','> ns = '
            read *, ns
            print '(A,$)','> dns/dlnk = '
            read *, run
            print '(A,$)','> deltaR2 = '
            read *, deltaR2
            print '(A,$)','> total mass of neutrinos in [eV] = '
            read *, mnu
            print '(A)', '> redshift = '
            read *, z
        ELSEIF(iargc() == 11) THEN
            call getarg(1, ofile_prefix)
            call getarg(2, arg); read(arg,*) om0
            call getarg(3, arg); read(arg,*) ode0
            call getarg(4, arg); read(arg,*) ob0
            call getarg(5, arg); read(arg,*) h0
            call getarg(6, arg); read(arg,*) w
            call getarg(7, arg); read(arg,*) ns
            call getarg(8, arg); read(arg,*) run
            call getarg(9, arg); read(arg,*) deltaR2
            call getarg(10, arg); read(arg,*) mnu
            call getarg(11, arg); read(arg,*) z
        ELSE
            print *, 'number of input parameters should be 11!!'
            print *, 'quit'
            stop
        ENDIF
    END SUBROUTINE ReadParams
!!$======================================================================
END PROGRAM Compute_Pk
