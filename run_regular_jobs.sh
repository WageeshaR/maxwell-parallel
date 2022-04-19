#!/bin/bash

module load compiler/GCC/11.2.0 && make &&
sbatch omp_regular_2.job &&
sbatch omp_regular_4.job &&
sbatch omp_regular_8.job &&
sbatch omp_regular_16.job