!-----------------------------------------------------------------------------
!
! Linear Growth Factor, g(z)=D(z)(1+z), which is constant and equal to 
! unity during the matter dominated era. Also computes its derivative 
! with respect to lna, dg/dlna.
! 
! *** Radiation density is ignored in the calculation ***
!
! << sample code >>
!
! USE cosmo
! USE growth
! double precision :: g,dgdlna,redshift
! external g,dgdlna
! ! Specify three cosmological parameters in double precision.
! ! The data type has been defined in MODULE cosmo.
! ode0=0.723d0
! om0=0.277d0
! w   = -1d0
! CALL setup_growth
! redshift=10d0
! print*,'omega matter=',om0
! print*,'omega de=',ode0
! print*,'w=',w
! print*,'g(z) at z=',redshift,'is',g(redshift)
! print*,'dg/dlna(z) at z=',redshift,'is',dgdlna(redshift)
! end
!
! August 23, 2008
! E. Komatsu
!
!-----------------------------------------------------------------------------
MODULE Growth
  INTEGER, parameter :: nw=501
  DOUBLE PRECISION, dimension(nw) :: xa,ga,g2a,dga,dg2a
contains
  SUBROUTINE Setup_Growth
    ! tabulate g(x), where x=ln(a)
    IMPLICIT none
    INTEGER :: i,ind
    DOUBLE PRECISION :: x,xini,y(2),c(24),work(2,9),tol=1.d-8
    EXTERNAL fderiv
    ind=1 
    do i=1,nw
       x = (i-nw)*0.01 ! x=ln(a)=-5.00,-4.99,...,0
       xa(i) = x
       if(i.eq.1)then
          xini = xa(i)
          y(1)=1d0 ! initial condition, g=1, at z=-5.00 (z=147.413)
          y(2)=0d0 ! initial condition, g'=1, at z=-5.00 (z=147.413)
       else
          call dverk(2,fderiv,xini,y,x,tol,ind,c,2,work)
          xini = x
       endif
       ga(i) = y(1)
       dga(i)= y(2)
    enddo
    CALL spline(xa,ga,nw,1.d30,1.d30,g2a)  ! evaluate g2a(i)=d^2[g(i)]/dy^2
    CALL spline(xa,dga,nw,1.d30,1.d30,dg2a)! evaluate g2a(i)=d^2[dg/dlna(i)]/dy^2
    return
  END SUBROUTINE Setup_Growth
END MODULE Growth
!-----------------------------------------------------------------------------
SUBROUTINE fderiv(n,x,y,yprime)
  USE cosmo
! x = ln(a) = -ln(1+z)
! y(1) = g
! y(2) = g'
! See Eq.(76) of Komatsu et al. (2008) [WMAP 5-year]
  IMPLICIT none 
  INTEGER, intent(IN) :: n
  DOUBLE PRECISION :: x,y(n),yprime(n)
  DOUBLE PRECISION :: ok,ode,ok0,h2,z,a
  DOUBLE PRECISION, parameter :: or0 = 0d0 ! ignore radiation density
  z=dexp(-x)-1d0
  a=dexp(x)
  ok0=1d0-om0-ode0
  h2=om0/a**3d0+or0/a**4d0+ok0/a**2d0+ode0/a**(3d0*(1.+w))
  ok=ok0/h2*dexp(-2d0*x)
  ode=ode0/h2*dexp(-3d0*x*(1d0+w)) 
  yprime(1)= y(2) 
  yprime(2)= -(2.5d0+0.5d0*(ok-3d0*w*ode))*y(2) &
             -(2d0*ok+1.5d0*(1d0-w)*ode)*y(1)
  return
END SUBROUTINE fderiv
!-----------------------------------------------------------------------------
DOUBLE PRECISION FUNCTION g(z)
  USE growth
  ! Returns g(z) = D(z)(1+z) by interpolating the table generated by
  ! Setup_Growth. 
  IMPLICIT none
  DOUBLE PRECISION, intent(IN) :: z
  DOUBLE PRECISION :: a,b,hh,x
  INTEGER :: jlo
  x=-dlog(1d0+z) ! x=ln(a)=-ln(1+z)
  if(x==0d0)then
     g=ga(nw)
  else
     CALL hunt(xa,nw,x,jlo)
     hh=xa(jlo+1)-xa(jlo)
     a=(xa(jlo+1)-x)/hh
     b=(x-xa(jlo))/hh
     g=a*ga(jlo)+b*ga(jlo+1)+((a**3-a)*g2a(jlo)+(b**3-b)*g2a(jlo+1))*(hh**2)/6.
  endif
  return
END FUNCTION g
!-----------------------------------------------------------------------------
DOUBLE PRECISION FUNCTION dgdlna(z)
  USE growth
  ! Returns dg/dlna(z) by interpolating the table generated by
  ! Setup_Growth.
  IMPLICIT none
  DOUBLE PRECISION, intent(IN) :: z
  DOUBLE PRECISION :: a,b,hh,x
  INTEGER :: jlo
  x=-dlog(1d0+z) ! x=ln(a)=-ln(1+z)
  if(x==0d0)then
     dgdlna=dga(nw)
  else
     CALL hunt(xa,nw,x,jlo)
     hh=xa(jlo+1)-xa(jlo)
     a=(xa(jlo+1)-x)/hh
     b=(x-xa(jlo))/hh
     dgdlna=a*dga(jlo)+b*dga(jlo+1)+((a**3-a)*dg2a(jlo)+(b**3-b)*dg2a(jlo+1))*(hh**2)/6.
  endif
  return
END FUNCTION dgdlna