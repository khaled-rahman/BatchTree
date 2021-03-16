# BatchPrEL
BatchPrEL generates edge-crossing-free and label-overlapping-free layout of a tree.

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

Input file must be in matrix market format ([check here for details about .mtx file](https://math.nist.gov/MatrixMarket/formats.html)). A lot of datasets can be found at [suitesparse website](https://sparse.tamu.edu). We provide few example input files in datasets/input directory. 

#### Layout generation ####
To generate crossing-free layout using BatchPrEL, use the following command:
```
$ ./bin/BatchPrEL -input datasets/raw/Graph_8.txt.weighted.mtx -output datasets/output/ -iter 120 -lr 0.8 -batch 256 -nsamples 20 -label datasets/raw/Graph_8.txt.labels -lrforlo 0.5 -iter2nd 150
```
Here, `-input` is the full path of input file, `-output` is the directory where output file will be saved, `-iter` is the number of iterations, `-batch` is the size of minibatch which is 128 here, and `-algo` is the choice of algorithm to run which is 2 represending cache blocking stochastic minibatch update algorithm. 

#### Post-processing step ####
After running the first step, to remove existing label overlaps and increase area utilization, we can use this post-processing step. In this step, we will initialize the layout generated in previous step as follows:
```
$ ./bin/BatchPrEL -input datasets/raw/Graph_8.txt.weighted.mtx -output datasets/output/ -label datasets/raw/Graph_8.txt.labels -scalingbox 400 -box 100 -psamples 600 -initf datasets/output/Graph_8.txt.weighted.mtxBatchPrEL128PARAOUT200.txt
```
We can run this post-processing step multiple times to get the desired output.

All options are described below:
```
-input <string>, full path of input file (required).
-output <string>, directory where output file will be stored.
-batch <int>, size of minibatch.
-initf <string>, a layout file to initializee the coordinates of all vertices.
-iter <int>, number of iteration.
-lr <float>, learning rate for edgecrossing free drawing.
-threads <int>, number of threads, default value is maximum available threads in the machine.
-lrforlo <float>, learning rate for overlap removal step.
-iter2nd <int>, iterations for overlap removal step.
-scalingbox <float>, scaling box sizee for post-processing.
-box <float>, box size for post-processing.
-psamples <int>, number of samples for post-processing.
-h, show help message.

default: -batch 128 -iter 600 -threads MAX -algo 2 -init 0
```

### Contact 
If you have questions, please don't hesitate to ask me (Md. Khaledur Rahman) by sending email to `morahma@iu.edu`.
