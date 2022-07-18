# Critical path tool to collect performance model factors on the fly

## How to build 
### For MPICH-based MPI

```lang=BASH
$ mkdir BUILD && cd BUILD
$ export MPICH_CC=clang
$ export MPICH_CXX=clang++
$ cmake .. -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ \
-DMPI_C_COMPILER=mpicc -DMPI_CXX_COMPILER=mpicxx -DMPI_LINK_FLAGS=-flto
$ make
```
For IntelMPI change `mpicc` to `mpigcc` and `mpicxx` to `mpigxx` 

### For OpenMPI

```lang=BASH
$ mkdir BUILD && cd BUILD
$ export OMPI_CC=clang
$ export OMPI_CXX=clang++
$ cmake .. -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ \
-DMPI_C_COMPILER=mpicc -DMPI_CXX_COMPILER=mpicxx -DMPI_LINK_FLAGS=-flto
$ make
```

## How to use (basic mode)

Build your code for production. Then execute with the critical path tool loaded.

### OpenMP-only code

```lang=BASH
$ env OMP_TOOL_LIBRARIES=<cpt-path>/BUILD/libompt_criticalpath.so <your-exe>
```

### MPI-only code

```lang=BASH
$ mpirun -np 4 env LD_PRELOAD=<cpt-path>/BUILD/libompt_mpicriticalpath.so \
 <your-exe>
```

### hybrid MPI+OpenMP code

```lang=BASH
$ mpirun -np 4 env LD_PRELOAD=<cpt-path>/BUILD/libompt_mpicriticalpath.so \
 OMP_TOOL_LIBRARIES=<cpt-path>/BUILD/libompt_mpicriticalpath.so <your-exe>
```

## How to use (targeted mode)

In current state, the critical path tool allows to focus on a single region of code by marking start and end of the region. 
**Note**: Make sure to only track a single iteration in case start and end markers get executed multiple times!

For begin and end markers, OpenMP or MPI-specific functions can be used:

OpenMP:
```lang=C
omp_control_tool(omp_control_tool_start, 0, NULL);  // Mark start of the region
omp_control_tool(omp_control_tool_end, 0, NULL);    // Mark end of the region
```

MPI:
```lang=C
MPI_Pcontrol(1);                                    // Mark start of the region
MPI_Pcontrol(0);                                    // Mark end of the region
```

To delay the start of the tool to the start marker, add the following `ANALYSIS_OPTIONS` environmental variable option for the execution:

```lang=BASH
$ env ANALYSIS_OPTIONS=start_stopped=1 \
 OMP_TOOL_LIBRARIES=<cpt-path>/BUILD/libompt_criticalpath.so <your-exe>
```