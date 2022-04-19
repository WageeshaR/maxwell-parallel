#!/bin/bash

make;

./maxwell -x 100 -y 100 -f 1 -o cuda/my_sim -n 1000 -p 2 -q 2 -r 5 -s 5 > cuda_regular_2_5.log && printf "\n\n" >> cuda_regular_2_5.log &&
./maxwell -x 400 -y 400 -f 1 -o cuda/my_sim -n 1000 -p 2 -q 2 -r 5 -s 5 >> cuda_regular_2_5.log;

./maxwell -x 100 -y 100 -f 1 -o cuda/my_sim -n 1000 -p 2 -q 2 -r 10 -s 10 > cuda_regular_2_10.log && printf "\n\n" >> cuda_regular_2_10.log &&
./maxwell -x 400 -y 400 -f 1 -o cuda/my_sim -n 1000 -p 2 -q 2 -r 10 -s 10 >> cuda_regular_2_10.log;

./maxwell -x 100 -y 100 -f 1 -o cuda/my_sim -n 1000 -p 4 -q 4 -r 5 -s 5 > cuda_regular_4_5.log && printf "\n\n" >> cuda_regular_4_5.log &&
./maxwell -x 400 -y 400 -f 1 -o cuda/my_sim -n 1000 -p 4 -q 4 -r 5 -s 5 >> cuda_regular_4_5.log && printf "\n\n" >> cuda_regular_4_5.log &&
./maxwell -x 1000 -y 1000 -f 1 -o cuda/my_sim -n 1000 -p 4 -q 4 -r 5 -s 5 >> cuda_regular_4_5.log;

./maxwell -x 100 -y 100 -f 1 -o cuda/my_sim -n 1000 -p 5 -q 5 -r 10 -s 10 > cuda_regular_5_10.log && printf "\n\n" >> cuda_regular_5_10.log &&
./maxwell -x 400 -y 400 -f 1 -o cuda/my_sim -n 1000 -p 5 -q 5 -r 10 -s 10 >> cuda_regular_5_10.log && printf "\n\n" >> cuda_regular_5_10.log &&
./maxwell -x 1000 -y 1000 -f 1 -o cuda/my_sim -n 1000 -p 5 -q 5 -r 10 -s 10 >> cuda_regular_5_10.log;

./maxwell -x 400 -y 400 -f 1 -o cuda/my_sim -n 1000 -p 8 -q 8 -r 25 -s 25 > cuda_regular_8_25.log && printf "\n\n" >> cuda_regular_8_25.log &&
./maxwell -x 1000 -y 1000 -f 1 -o cuda/my_sim -n 1000 -p 8 -q 8 -r 25 -s 25 >> cuda_regular_8_25.log;

./maxwell -x 400 -y 400 -f 1 -o cuda/my_sim -n 1000 -p 10 -q 10 -r 20 -s 20 > cuda_regular_10_20.log && printf "\n\n" >> cuda_regular_10_20.log &&
./maxwell -x 1000 -y 1000 -f 1 -o cuda/my_sim -n 1000 -p 10 -q 10 -r 20 -s 20 >> cuda_regular_10_20.log;

./maxwell -x 1000 -y 1000 -f 1 -o cuda/my_sim -n 1000 -p 10 -q 10 -r 25 -s 25 > cuda_regular_10_25.log;