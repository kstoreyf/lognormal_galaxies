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
  double precision, allocatable, dimension(:) :: xi,Rh ! h^-1 Mpc
  double precision :: k_ov_h,k0_ov_h,kmax_ov_h ! h Mpc^-1
  double precision :: linear_pk
  double precision :: dk,factor,dummy
  double precision :: pk
  character(len=128) :: input_pk_fname, ofile_prefix, arg
  integer :: n,i,j,nrad=1000
  integer :: lenaw=10000, jj
  double precision, allocatable :: aw(:)
  double precision :: answer = 0.d0, err, r
  external linear_pk

  call ReadParams

  ! maximum wavenumber to which P(k) is integrated.
  ! Use a larger value for more accurate results.
  kmax_ov_h=30.d0

  ! read in linear P(k)
  CALL open_linearpk(input_pk_fname,n)
  open(11,file=input_pk_fname,status='old')
  read(11,*)k0_ov_h,dummy
  close(11)

  k_ov_h=k0_ov_h
  ALLOCATE(xi(nrad),Rh(nrad))
  xi=0d0
  dk=1d-4
  Rh=0d0

  allocate(aw(lenaw))
  CALL intdeoini(lenaw, 1d-307, 1d-8, aw)

  do j=1,nrad
    Rh(j) = 10.d0**(-3.d0+8.d0*(dble(j)/nrad))
  end do

  do j=1,nrad
    r = Rh(j)
    if (r > 0.d0) then
      CALL intdeo(int_pk_func, k0_ov_h, r, aw, xi(j), err)
    endif
  enddo

  CALL close_linearpk

  close(11)

  ! write out the result
  open(2,file='./data/inputs/'//trim(ofile_prefix)//'_Rh_xi.txt',status='unknown')
  do i=1,nrad
     write(2,'(2E16.5)')Rh(i),xi(i)
  enddo
  close(2)

  DEALLOCATE(xi,Rh)
!!$======================================================================
CONTAINS
!!$======================================================================
    DOUBLE PRECISION FUNCTION int_pk_func(k) result(ans)
        DOUBLE PRECISION, INTENT(IN) :: k
        DOUBLE PRECISION :: pk

        if (k < kmax_ov_h) then
           pk=linear_pk(k)
        else
           pk = 0.0
        endif

        if(r == 0.d0)then
            ans = pk*k**2d0/(2d0*3.1415926535d0**2d0) ! h^3 Mpc^-3
        else
            ans = pk*sin(k*r)/(k*r)*k**2d0/(2d0*3.1415926535d0**2d0) ! h^3 Mpc^-3
        endif
    END FUNCTION
!!$======================================================================
    SUBROUTINE ReadParams
        ! read input and output file names
        IF (iargc() == 0) THEN
            print '(A,$)', '> enter the prefix of output file name: '
            read *, ofile_prefix 
            print '(A,$)', '> enter the input Pk file name: '
            read *, input_pk_fname
            print '(A,$)', '> enter the number of lines in the input Pk file: '
            read *, n
        ELSEIF (iargc() == 3) THEN
            call getarg(1, ofile_prefix)
            call getarg(2, input_pk_fname)
            call getarg(3, arg); read(arg,*) n
        ELSE
            print *, 'number of input parameters should be 3!!'
            print *, 'quit'
            stop
        ENDIF
    END SUBROUTINE 
!!$======================================================================
END PROGRAM Compute_Xi
