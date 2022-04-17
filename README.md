# OpenMP parallelisation of 2D Maxwell Solver

Below commands can be used to run the application in different configuration modes

## Building

The application can be built with the provided Makefile. e.g.

```
$ make
```

This will build a maxwell binary using -fopenmp extension.

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
$ mkdir out
$ ./maxwell -x 100 -y 100 -c -f 1 -o out/my_sim -n 1000 -l 8
```

This will run a 100 x 100 problem (perfect for manageable visualisation), enable checkpointing every timestep, and will output the VTK files to the out directory. It will run for 1000 steps. The
argument *l* specifies the number of OpenMP parallel threads to run in each parallel region.

## Validation of output

There are two comparison modes built into the application to perform validation in contrast to the original application. Comparision mode (comp_mode) can be passed as an argument using *-e* flag followed by values 1 or 2

* comp_mode == 1: Perform comparison only with total maginitudes of E and B fields at the end of each iteration
* comp_mode == 2: Perform comparison with total output files at the end of each iteration

At the end of either of comparision modes, the application with print to the sysout *total_error* in scientific format.

