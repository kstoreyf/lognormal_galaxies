PROGRAM Read_Lognormal
  implicit none
  integer::i,j
  integer::ngal,nshow
  real(kind=8)::Lx,Ly,Lz
  real(kind=4),allocatable::pos(:,:),vel(:,:)
  character(len=128)::filename

  nshow=30 ! # of lines to display
  ! filename='lognormal_is1_rlz0_Pnmax128_b1.455_tophat_Rs10.bin'
  filename = 'data/catalog_example/lognormal/example_lognormal_rlz0.bin'

  open(unit=1,file=filename,form='unformatted',action='read',access='stream')
  read(1) Lx,Ly,Lz
  read(1) ngal
  write(*,*) "(Lx,Ly,Lz) = ",Lx,Ly,Lz
  write(*,*) "number of galaxies = ",ngal
  allocate(pos(3,ngal))
  allocate(vel(3,ngal))
  print*,'Showing the first',nshow,' lines (x[Mpc/h], y[Mpc/h], z[Mpc/h], vx[km/s], vy[km/s], vz[km/s]):'
  do i=1,ngal
     do j=1,3
        read(1) pos(j,i)
     enddo
     do j=1,3
        read(1) vel(j,i)
     enddo
     if (i<=nshow) then
        write(*,'(6E13.5)') pos(1,i),pos(2,i),pos(3,i),vel(1,i),vel(2,i),vel(3,i)
     endif
  enddo
  close(1)
end program Read_Lognormal
