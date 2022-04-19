#!/bin/bash

make;

./maxwell -x 100 -y 100 -f 1 -o out/my_sim -n 1000 > original_regular.log && printf "\n\n" >> original_regular.log &&
./maxwell -x 400 -y 400 -f 1 -o out/my_sim -n 1000 >> original_regular.log && printf "\n\n" >> original_regular.log &&
./maxwell -x 1000 -y 1000 -f 1 -o out/my_sim -n 1000 >> original_regular.log;