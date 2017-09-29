      subroutine run_modpcrtm

      use mpi

      implicit none
      logical PASS1
      common /GTBRDFC/PASS1

      integer ierr,num_procs,id
      call mpi_init(ierr)
         call MPI_Comm_rank ( MPI_COMM_WORLD, id, ierr )
         call MPI_Comm_size ( MPI_COMM_WORLD, num_procs, ierr )
        pass1 = .true.
        call driver(id)
      call mpi_finalize(ierr)

      return
      end
