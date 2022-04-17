# CUDA parallelisation of 2D Maxwell Solver

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

The application cannot be run in its default configuration as there are additional flags required to run: Available options can be queried with:

```
$ ./maxwell --help
```

```
$ mkdir out
$ ./maxwell -x 100 -y 100 -c -f 1 -o cuda/my_sim -n 1000 -p 2 -q 2 -r 5 -s 5
```

This will run a 100 x 100 problem in 2x2 CUDA grid resolution and 5x5 CUDA block resolution, enable checkpointing every timestep, and will output the VTK files to the out directory. It will run for 1000 steps.

## Validation of output

There are two comparison modes built into the application to perform validation in contrast to the original application. Comparision mode (comp_mode) can be passed as an argument using *-e* flag followed by values 1 or 2

* comp_mode == 1: Perform comparison only with total maginitudes of E and B fields at the end of each iteration
* comp_mode == 2: Perform comparison with total output files at the end of each iteration

At the end of either of comparision modes, the application with print to the sysout *total_error* in scientific format.

**Important:** To run in comparison mode, the original application has to be run in prior. Follow the procedure in the ***main*** branch to create required output files.

```
$ ./maxwell -x 100 -y 100 -f 1 -n 1000 -p 2 -q 2 -r 5 -s 5 -e 1
```

Executing above line will run the application with comp_mode == 1. Notice the absence of *-c* flag. It is not required as it only compares total magnitudes.

```
$ ./maxwell -x 100 -y 100 -c -f 1 -n 1000 -p 2 -q 2 -r 5 -s 5 -e 2
```

Executing above line will run the application with comp_mode == 2. *-c* flag is required to run in this mode as the comparison performs at file writing stage.


