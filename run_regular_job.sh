#!/bin/bash

module load compiler/GCC/11.2.0 && make &&
sbatch original_regular.job