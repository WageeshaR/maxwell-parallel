#!/bin/bash

module load compiler/GCC/11.2.0 system/CUDA/11.0.2-GCC-9.3.0 && make &&
sbatch cuda_regular_2_5.job &&
sbatch cuda_regular_2_10.job &&
sbatch cuda_regular_4_5.job &&
sbatch cuda_regular_5_10.job &&
sbatch cuda_regular_8_25.job &&
sbatch cuda_regular_10_20.job &&
sbatch cuda_regular_10_25.job &&
sbatch cuda_regular_10_50.job