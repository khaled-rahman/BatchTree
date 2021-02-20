# BatchPrEL

## System Requirements

Users need to have following softwares/tools installed in their PC. The source code was compiled and run successfully in both linux and macOS.
```
GCC version >= 4.9
OpenMP version >= 4.5
Python3 packages: matplotlib, scipy, networkx, numpy
```
Some helpful links for installation can be found at [GCC](https://gcc.gnu.org/install/), [OpenMP](https://clang-omp.github.io) and [Environment Setup](http://heather.cs.ucdavis.edu/~matloff/158/ToolsInstructions.html#compile_openmp).

## Compile BatchPrEL
To compile BatchPrEL type the following command on terminal:
```
$ make clean
$ make
```
This will generate an executible file in bin folder.

## Run BatchPrEL from command line

Input file must be in matrix market format ([check here for details about .mtx file](https://math.nist.gov/MatrixMarket/formats.html)). A lot of datasets can be found at [suitesparse website](https://sparse.tamu.edu). We provide few example input files in datasets/input directory. To run BatchPrEL, use the following command:
```
$ ./bin/BatchPrEL -input ./datasets/input/sf_ba6000.mtx -output ./datasets/output/ -iter 10 -batch 256 -threads 32 -algo 2
```
Here, `-input` is the full path of input file, `-output` is the directory where output file will be saved, `-iter` is the number of iterations, `-batch` is the size of minibatch which is 256 here, `-threads` is the maximum number of threads which is 32 and `-algo` is the choice of algorithm to run which is 2 represending cache blocking stochastic minibatch update algorithm. All options are described below:
```
-input <string>, full path of input file (required).
-output <string>, directory where output file will be stored.
-batch <int>, size of minibatch.
-init <int>, any of 0 or 1, 1 - random initialization, 0 - greedy initialization.
-initf <string>, this will overwrite -init and initialize coordinates from a given file.
-iter <int>, number of iteration.
-threads <int>, number of threads, default value is maximum available threads in the machine.
-algo <int>, an integer among 0, 1, 2, 3, 4, 5 and 6.
        0 - for sequential layout generation.
        1 - for naive parallel approach for layout generation.
        2 - for parallel layout generation using cache blocking stochastic minibatch update.
        3 - for parallel layout generation using linlog mode, (0,-1)-energy model.
        4 - for parallel layout generation using barnes-hut repulsive force approximation.
        5 - for parallel layout generation using greedy force approximation approach.
        6 - for parallel layout generation, 80% time in -algo 4 and last 20% time in -algo 2.
        7 - for parallel layout generation, this will calculate forces using linear algebric format.
        8 - for parallel layout generation, (1,-1)-energy model.
-bht <float>, barnes-hut threshold parameter.
-engt <float>, convergenece criteria. It should be a floating point number between 0 to 1 exclusive.
It indicates if energy value is improved less than this percenrages then optimization will stop and return layout. If this value is set, then number of iteration will be overwritten.
-weight <float>, a real number to include edge weight.
-h, show help message.

default: -weight 1.0 -bht 1.2 -engt 0.01 -batch 256 -iter 600 -threads MAX -algo 2 -init 0
```

### Contact 
If you have questions, please don't hesitate to ask me (Md. Khaledur Rahman) by sending email to `morahma@iu.edu`.
