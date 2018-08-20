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

