# A very simple 2D Maxwell Solver using Yee FDTD

This is a simple solver for Maxwell's equations in two dimensions. The application uses the method introduced by Kane Yee in 1966 [1].

The problem is set up with an electric current rotating around the centre of a 1.0m x 1.0m box. By default, the problem is solved on a 4000 x 4000 grid.

At the end of execution, a VTK file is produced for visualisation purposes. This can be loaded into VisIt for analysis.

## Building

The application can be built with the provided Makefile. e.g.

```
$ make
```

This will build a maxwell binary.

## Running

The application can be ran in its default configuration with:

```
$ ./maxwell
```

This will run with 4000 x 4000 grid and will output status every 100 iterations. At the end of execution a VTK file will be produced. 

There are numourous other options available. These can be queried with:

```
$ ./maxwell --help
```

**Beware:** In the default configuration, the VTK file will be in the region of 2 GB for a 4000 x 4000 problem. For visualisation purposes you should decrease the grid resolution. If you would like to visualise the entire execution, you should enable checkpointing and increase the output frequency. You can do this like so:

```
$ mkdir mpi
$ mpiexec -n 8 ./maxwell -x 100 -y 100 -c -f 1 -o mpi/my_sim -n 1000
```

This will run a 100 x 100 problem (perfect for manageable visualisation), enable checkpointing every timestep, and will output the VTK files to the out directory. It will run for 1000 steps. *-n* flag indicates the number of parallel nodes to run the application.

## Validation of output

There are two comparison modes built into the application to perform validation in contrast to the original application. Comparision mode (comp_mode) can be passed as an argument using *-e* flag followed by values 1 or 2

* comp_mode == 1: Perform comparison only with total maginitudes of E and B fields at the end of each iteration
* comp_mode == 2: Perform comparison with total output files at the end of each iteration

At the end of either of comparision modes, the application with print to the sysout *total_error* in scientific format.

**Important:** To run in comparison mode, the original application has to be run in prior. Follow the procedure in the ***main*** branch to create required output files.

```
$ mpiexec -n 8 ./maxwell -x 100 -y 100 -f 1 -o -n 1000 -e 1
```

Executing above line will run the application with comp_mode == 1. Notice the absence of *-c* flag. It is not required as it only compares total magnitudes.

```
$ mpiexec -n 8 ./maxwell -x 100 -y 100 -c -f 1 -o -n 1000 -e 2
```

Executing above line will run the application with comp_mode == 2. *-c* flag is required to run in this mode as the comparison performs at file writing stage.
