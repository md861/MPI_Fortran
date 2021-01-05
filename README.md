# MPI_Fortran
Example send and receive Master-Slave algorithm for multiplication of (random) matrices and vectors over distributed nodes. The algorithm is a basic FIFO (First In First Out) distribution of jobs to slaves. 

## Usage:
Open a terminal in the directory of the codes. Type in following for:
* Compilation: `mpifort [-o <outputfile>] filename`. Here `filename` is either *fMatMatMul.f90* or *fMatVecMul.f90*. The `-o` is an optional flag, that could be used to create an output file of desired name.
* Run: `mpirun -np nprocs [--use-hwthread-cpus, --oversubscribe] <outputfile>`. Here `nprocs` is the number of parallel processes to be used. By default, the processes are distributed over available CPUs. If `nprocs` exceeds the maximum number of CPUs available, an error is issued. The `--use-hwthread-cpus` option allows for forcing the number of threads when `nprocs` is within the maximum number of threads available on the hardware (otherwise an error message is issued). The `--oversubscribe` option allows for ignoring the number of CPUs/threads available.

## Example 
An example compilation and run for *fMatVecMul.f90* with just a single slave is given below:

`$ mpifort -o test fMatVecMul.f90`  
`$ mpirun -np 2 --use-hwthread-cpus test`  
![](/images/SingleSlaveExample.png)

## Profiling tools
[Extrae](https://tools.bsc.es/tools_hands-on) is an open-source profiling tool that provides various event logs for MPI, OpenMP (amongst other) libraries. The log files generated during run time can be viewed with [Paraver](https://tools.bsc.es/paraver) - a graphical tool for viewing program performance and profiling. Some examples are shown below for the *fMatVecMul.f90* file when run with 4 threads.
![](/images/MPIstats.png)![](/images/MPI_calls.png)![](/images/com_size.png)
