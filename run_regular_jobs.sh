#!/bin/bash
sbatch mpi_regular_2.job &&
sbatch mpi_regular_4.job &&
sbatch mpi_regular_8.job &&
sbatch mpi_regular_16.job