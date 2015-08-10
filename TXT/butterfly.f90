subroutine butterfly(newsig,sigma)
  ! Written by Jack Blake, Bath University, 2011
  ! 
  ! This is a fortran code for sharing data between all
  ! processors by using the butterfly topology structure.
  ! There is a BLAS routine for this but I enjoyed writing
  ! it nonetheless, and it's pretty much the same efficiency.

  include "mpif.h"

  real (kind=8), intent(in)   :: newsig
  real (kind=8), intent(out)  :: sigma
  integer                     :: ierr, myid, numprocs, dest

  call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, numprocs, ierr)
  
  sigma = newsig
  
  ! Butterfly Topology 

  ! loop over the number of 'rounds' of exchange needed
  do i = 1,int(log(real(numprocs,8))/log(2.0_8))
     
     ! First find out which section of the array the
     ! processor is in. eg. in the first round, the
     ! sections are (0,1) (2,3) (4,5) etc, but in the
     ! second round they are: (0,1,2,3) (4,5,6,7) etc.
     k = 1
     dest = myid
     do 
        if (mod(dest,2**i) == dest) exit
        dest = dest - 2**i
        k = k+1
     enddo

     ! destination of this processors 'send's
     dest = mod(myid + 2**(i-1),2**i) + (k-1)*2**i

     ! send first or receive first?
     if (myid < dest) then
        ! Send
        call MPI_Send(sigma,1,MPI_DOUBLE_PRECISION,dest,myid, &
             MPI_COMM_WORLD,ierr)
        ! Receive
        call MPI_RECV(newsig,1,MPI_DOUBLE_PRECISION,dest, &
             MPI_ANY_TAG,MPI_COMM_WORLD, stat,ierr)
     else
        ! Receive
        call MPI_RECV(newsig,1,MPI_DOUBLE_PRECISION,dest, &
             MPI_ANY_TAG,MPI_COMM_WORLD, stat,ierr)
        ! Send
        call MPI_Send(sigma,1,MPI_DOUBLE_PRECISION,dest,myid, &
             MPI_COMM_WORLD,ierr)
     endif

     ! update sigma
     sigma = sigma + newsig
     
  enddo

end subroutine butterfly
