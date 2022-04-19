#!/bin/bash

module load compiler/GCC/11.2.0 mpi/OpenMPI/4.1.1-GCC-11.2.0 && make &&
sbatch mpi_regular_2.job &&
sbatch mpi_regular_4.job &&
sbatch mpi_regular_8.job &&
sbatch mpi_regular_16.job