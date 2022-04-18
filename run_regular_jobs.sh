#!/bin/bash

module load mpi/OpenMPI/4.0.5-GCC-10.2.0 &&
make &&
sbatch mpi_regular_2.job &&
sbatch mpi_regular_4.job &&
sbatch mpi_regular_8.job &&
sbatch mpi_regular_16.job