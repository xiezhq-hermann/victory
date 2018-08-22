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

Today we have a clear understanding to our data, and successfully manage the data transfer both in unified memory and manually data region. Therefore, we make some directives workable but stuck in terrible busy memory copy, which leads to a 30 times faster, no, slower. Then we deploy the cuBlas library and more directives to help us accelerating the linear operation, which also reduce the memory transfer. It's somehow troublesome to deploy the cuBlas library in our project, we tried it out for much time and finally successful with the help of Matt and . Thank you for your kind help. Now the directives achieve a reasonable performance, 30% faster in one subroutine, we'll keep working on optimizing the directives. Besides we're writing some more CUDA kernels to be auxiliary part in further code. That's all, we think it could be better tomorrow.