#!/bin/bash

make;

mpiexec -n 2 ./maxwell -x 100 -y 100 -f 1 -o mpi/my_sim -n 1000 > mpi_regular_2.log && printf "\n\n" >> mpi_regular_2.log &&
mpiexec -n 2 ./maxwell -x 400 -y 400 -f 1 -o mpi/my_sim -n 1000 >> mpi_regular_2.log && printf "\n\n" >> mpi_regular_2.log &&
mpiexec -n 2 ./maxwell -x 1000 -y 1000 -f 1 -o mpi/my_sim -n 1000 >> mpi_regular_2.log;

mpiexec -n 4 ./maxwell -x 100 -y 100 -f 1 -o mpi/my_sim -n 1000 > mpi_regular_4.log && printf "\n\n" >> mpi_regular_4.log &&
mpiexec -n 4 ./maxwell -x 400 -y 400 -f 1 -o mpi/my_sim -n 1000 >> mpi_regular_4.log && printf "\n\n" >> mpi_regular_4.log &&
mpiexec -n 4 ./maxwell -x 1000 -y 1000 -f 1 -o mpi/my_sim -n 1000 >> mpi_regular_4.log;

mpiexec -n 8 ./maxwell -x 100 -y 100 -f 1 -o mpi/my_sim -n 1000 > mpi_regular_8.log && printf "\n\n" >> mpi_regular_8.log &&
mpiexec -n 8 ./maxwell -x 400 -y 400 -f 1 -o mpi/my_sim -n 1000 >> mpi_regular_8.log && printf "\n\n" >> mpi_regular_8.log &&
mpiexec -n 8 ./maxwell -x 1000 -y 1000 -f 1 -o mpi/my_sim -n 1000 >> mpi_regular_8.log;