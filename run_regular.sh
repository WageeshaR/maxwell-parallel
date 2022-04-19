#!/bin/bash

make;

./maxwell -x 100 -y 100 -f 1 -o omp/my_sim -n 1000 -l 2 > omp_regular_2.log && printf "\n\n" >> omp_regular_2.log &&
./maxwell -x 400 -y 400 -f 1 -o omp/my_sim -n 1000 -l 2 >> omp_regular_2.log && printf "\n\n" >> omp_regular_2.log &&
./maxwell -x 1000 -y 1000 -f 1 -o omp/my_sim -n 1000 -l 2 >> omp_regular_2.log;

./maxwell -x 100 -y 100 -f 1 -o omp/my_sim -n 1000 -l 4 > omp_regular_4.log && printf "\n\n" >> omp_regular_4.log &&
./maxwell -x 400 -y 400 -f 1 -o omp/my_sim -n 1000 -l 4 >> omp_regular_4.log && printf "\n\n" >> omp_regular_4.log &&
./maxwell -x 1000 -y 1000 -f 1 -o omp/my_sim -n 1000 -l 4 >> omp_regular_4.log;

./maxwell -x 100 -y 100 -f 1 -o omp/my_sim -n 1000 -l 8 > omp_regular_8.log && printf "\n\n" >> omp_regular_8.log &&
./maxwell -x 400 -y 400 -f 1 -o omp/my_sim -n 1000 -l 8 >> omp_regular_8.log && printf "\n\n" >> omp_regular_8.log &&
./maxwell -x 1000 -y 1000 -f 1 -o omp/my_sim -n 1000 -l 8 >> omp_regular_8.log;

./maxwell -x 100 -y 100 -f 1 -o omp/my_sim -n 1000 -l 16 > omp_regular_16.log && printf "\n\n" >> omp_regular_16.log &&
./maxwell -x 400 -y 400 -f 1 -o omp/my_sim -n 1000 -l 16 >> omp_regular_16.log && printf "\n\n" >> omp_regular_16.log &&
./maxwell -x 1000 -y 1000 -f 1 -o omp/my_sim -n 1000 -l 16 >> omp_regular_16.log;