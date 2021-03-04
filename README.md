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
$ ./bin/BatchPrEL -input datasets/raw/Graph_8.txt.mtx -output dataset/output/ -iter 200 -lr 0.8 -batch 128 -algo 2 -nsamples 20 -label datasets/raw/Graph_8.txt.labels -lrforlo 0.5 -iter2nd 200
```
Here, `-input` is the full path of input file, `-output` is the directory where output file will be saved, `-iter` is the number of iterations, `-batch` is the size of minibatch which is 128 here, and `-algo` is the choice of algorithm to run which is 2 represending cache blocking stochastic minibatch update algorithm. All options are described below:
```
-input <string>, full path of input file (required).
-output <string>, directory where output file will be stored.
-batch <int>, size of minibatch.
-init <int>, any of 0 or 1, 1 - random initialization, 0 - greedy initialization.
-iter <int>, number of iteration.
-threads <int>, number of threads, default value is maximum available threads in the machine.
-algo <int>, an integer among 2 and 3.
        2 - for parallel layout generation using cache blocking minibatch update.
        3 - for parallel layout generation using linlog mode, (0,-1)-energy model.
-h, show help message.

default: -batch 128 -iter 600 -threads MAX -algo 2 -init 0
```

### Contact 
If you have questions, please don't hesitate to ask me (Md. Khaledur Rahman) by sending email to `morahma@iu.edu`.
