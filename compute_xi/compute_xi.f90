PROGRAM Compute_Xi
  USE linearpk
  ! A sample program for computing the two-point correlation function, xi(r),
  ! from the linear power spectrum in real space.
  ! The program generates spherically-averaged correlation functions.
  ! NOTE: the power spectrum is smoothed by a Gaussian with the width of 
  ! sigma_smooth = 1h^-1 Mpc, i.e., P(k) -> P(k)*exp(-k^2*sigma_smooth^2), 
  ! in order to suppress the aliasing (oscillations).   
  ! September 18, 2008: E.Komatsu
  IMPLICIT none
  double precision :: sigma_smooth=0.1d0 ! h^-1 Mpc
  double precision, allocatable, dimension(:) :: xi,Rh ! h^-1 Mpc
  double precision :: k_ov_h,k0_ov_h,kmax_ov_h ! h Mpc^-1
  double precision :: linear_pk
  double precision :: dk,factor,dummy
  double precision :: pk
  character(len=128) :: filename
  integer :: n,i,j,nrad=750
  external linear_pk
  ! maximum wavenumber to which P(k) is integrated is set to be 3 times
  ! 1/(the smoothing length). Use a larger value for more accurate results.
  kmax_ov_h=3d0/sigma_smooth
  ! read in linear P(k)
  filename='wavenumber_pk.txt'
  n=1037    ! # of lines in the file
  CALL open_linearpk(filename,n)
  open(11,file=filename,status='old')
  read(11,*)k0_ov_h,dummy
  close(11)
  k_ov_h=k0_ov_h
  ALLOCATE(xi(nrad),Rh(nrad))
  xi=0d0
  dk=1d-4
  Rh=0d0
  do while (k_ov_h<=kmax_ov_h)
     pk=linear_pk(k_ov_h)*exp(-(k_ov_h*sigma_smooth)**2d0)
     factor=(k_ov_h)**2d0*dk/(2d0*3.1415926535d0**2d0) ! h^3 Mpc^-3
     do j=1,nrad
        if(Rh(j)<=10d0)then
           Rh(j)=0.d0+0.1d0*dble(j-1)
        elseif (Rh(j)<=500d0)then
           Rh(j)=Rh(j-1)+1d0
        elseif (Rh(j)<=1d3)then
           Rh(j)=Rh(j-1)+5d0
        else
           Rh(j)=dexp(dlog(Rh(j-1))+0.046d0)
        endif
        if(Rh(j)==0d0)then
           xi(j)=xi(j)+pk*factor
        else
           xi(j)=xi(j)+pk*factor*sin(k_ov_h*Rh(j))/(k_ov_h*Rh(j))
        endif
     enddo
     k_ov_h=k_ov_h+dk
  enddo
  CALL close_linearpk
  close(11)
  ! write out the result
  open(2,file='Rh_xi.txt',status='unknown')
  do i=1,nrad
     write(2,'(2E16.5)')Rh(i),xi(i)
  enddo
  close(2)
  DEALLOCATE(xi,Rh)
END PROGRAM Compute_Xi
