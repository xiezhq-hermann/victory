program main
  
  use mpi_mod
  use global_parameter
  use math_mod
  use parquet_util
  use parquet_kernel
  use parquet_equation
  use parquet_selfenergy
  use parquet_phys
  
  implicit none

  ! ... local vars ...
  integer       :: ite
  real(dp)      :: t1, t2, t3, t4, t5, t2_1
  logical           :: Converged

  complex(dp), allocatable :: Grt(:, :, :)

  ! ------------------- initialize the mpi enviroment ----------------------
  call parallel_start

  ! ------------------- print the license message --------------------------
  master = 0
  if (id == master) call license_message(ntasks)

  call readin

  allocate(Grt(Nx, Ny, Nf))
  ! ------------------- start the self-consisent loop ----------------------
  Converged = .False.
  ite = 1
  do while (.NOT. Converged .and. ite < 50)
     if (id == master) call loop_message(ite)
     t1 = MPI_WTIME()
 
     ! calculate single-particle green's function
     call pa_Gkw_Chi0(ite, Grt)
     t2 = MPI_WTIME()
     if (id == master) write(*, "(a, f12.6)") '  time spent on pa_Gkw_Chi0 is:', t2-t1 
     
     ! determine reducible vertex function and its kernel approximation
        call reducible_vertex(ite)
        t2_1 = MPI_WTIME()
        if (id == master) write(*, "(a, f12.6)") '  time spent on reducible_vertex is:', t2_1-t2
        call get_kernel_function(ite)
     t3 = MPI_WTIME()
     if (id == master) write(*, "(a, f12.6)") '  time spent on kernel calculation is:', t3-t2_1

     ! solve the parquet equation
     call solve_parquet_equation(ite)
     t4 = MPI_WTIME()
     if (id == master) write(*, "(a, f12.6)") '  time spent on solve_parquet_equation is:', t4-t3

     ! calculate the self-energy 
     call self_energy(ite, Grt, converged)
     t5 = MPI_WTIME()
     if (id == master) write(*, "(a, f12.6)") '  time spent on self_energy is:', t5-t4

     ! --- update irreducible vertex in each channel ---
     G_d = F_d - G_d   
     G_m = F_m - G_m
     G_s = F_s - G_s
     G_t = F_t - G_t
        
     if (id == master) write(*, "(a, f12.6)") 'total time cost:', t5-t1
     ite = ite + 1
  end do

  ! ------------------- calculate eigen values in each channel -------------------------------
  call solve_eigen_equation    

  ! ------------------- clean the memory -----------------------------------
  call MPI_barrier(MPI_COMM_WORLD, rc)
  call Memory_Release

  ! ------------------- finalize the mpi enviroment ------------------------
  call parallel_end

  deallocate(Grt)
end program main
