### Preparation 

- Unified memory, successful compiled by PGI compiler and run correctly (2*2 cluster size)

### Computation Kernel

No MPI, Small GPU Kernels, No log writing, implement device function first and then design computation kernels.

#### solve_parquet_equation

- index_operation, list_index, kernel

#### reducible_vertex

- index_operation, ZGEMM (cuBLAS), list_index

#### self_energy

- fftd2b (cuFFT), FDfit (math module), nfourier, index_operation, kernel

We are GeekPie_HPC team and we are working on an application wrote in Fortran named victory, it's about calculating the energy and some other states of eletrons.
In this afternoon, we compiled our program by PGI compiler and fixed some incompatible issues. There are plenty global parameters are shared and used by many functions and subroutines in our program, so we're still working on the data managment and unified memory. And now we planed to reform our data manage structure to enable some non-trivial directives.