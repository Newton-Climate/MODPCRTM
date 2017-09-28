    program multicall
use mpi
implicit none

logical PASS1
!common /GTBRDFC/PASS1

integer ierr,num_procs,id
call mpi_init(ierr)
write(*,*), 'Got past MPI_INIT'
      call MPI_Comm_rank ( MPI_COMM_WORLD, id, ierr )
      call MPI_Comm_size ( MPI_COMM_WORLD, num_procs, ierr )
pass1 = .true.
write(*,*) 'entering driver'
call driver(id)
write(*,*), 'got past driver'
call mpi_finalize(ierr)


end program multicall
