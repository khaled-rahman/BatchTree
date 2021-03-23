#include <omp.h>
#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <functional>
#include <fstream>
#include <iterator>
#include <ctime>
#include <cstring>
#include <cmath>
#include <string>
#include <sstream>
#include <random>

#include "../sample/algorithms.hpp"

#include "../sample/Coordinate.hpp"

using namespace std;

#define VALUETYPE double
#define INDEXTYPE int


void helpmessage(){
	printf("\n");
	printf("Usage of BatchLayout tool:\n");
	printf("-input <string>, full path of input file (required).\n");
	printf("-output <string>, directory where output file will be stored.\n");
	printf("\nFor detailed instructions, visit https://github.com/khaled-rahman/BatchTree\n\n");
}

void myTest(){
	Coordinate<VALUETYPE> p1(3,4), p2(2,2), p3(6,3);
	auto pt = p1.getProjection(p2, p3);
	printf("x = %f, y = %f\n", pt.x, pt.y);
	

}

void TestAlgorithms(int argc, char *argv[]){
	VALUETYPE lr = 1.0, bhThreshold = 1.2, lrforlo = 0.5;
	INDEXTYPE psamples = 1000, init = 0, batchsize = 128, iterations = 200, numberOfThreads = omp_get_max_threads(), algoOption = 3, nsamples=10, iter=200;
	string inputfile = "", initfile = "", outputfile = "", labelfile = "", algoname = "CACHE", initname = "GREEDY";
	VALUETYPE scalingbox = 10000, pbox = 50, expc = 0.01;
	for(int p = 0; p < argc; p++){
		if(strcmp(argv[p], "-input") == 0){
			inputfile = argv[p+1];
		}
		if(strcmp(argv[p], "-output") == 0){
			outputfile = argv[p+1];
		}
		if(strcmp(argv[p], "-batch") == 0){
			batchsize = atoi(argv[p+1]);
		}
		if(strcmp(argv[p], "-nsamples") == 0){
                        nsamples = atoi(argv[p+1]);
                }
		if(strcmp(argv[p], "-init") == 0 && init != 2){
			init = atoi(argv[p+1]);
			if(init < 0 || init > 1) init = 1;
		}
		if(strcmp(argv[p], "-iter") == 0){
			iterations = atoi(argv[p+1]);
		}
		if(strcmp(argv[p], "-iter2nd") == 0){
			iter = atoi(argv[p+1]);
		}
		if(strcmp(argv[p], "-threads") == 0){
			numberOfThreads = atoi(argv[p+1]);
		}
		if(strcmp(argv[p], "-algo") == 0){
			algoOption = atoi(argv[p+1]);
		}
		if(strcmp(argv[p], "-psamples") == 0){
                        psamples = atoi(argv[p+1]);
                }
		if(strcmp(argv[p], "-box") == 0){
                        pbox = atof(argv[p+1]);
                }
		if(strcmp(argv[p], "-expc") == 0){
                        expc = atof(argv[p+1]);
                }
		if(strcmp(argv[p], "-scalingbox") == 0){
                        scalingbox = atof(argv[p+1]);
                }
		if(strcmp(argv[p], "-bht") == 0){
			bhThreshold = atof(argv[p+1]);
		}
		if(strcmp(argv[p], "-initf") == 0){
			initfile = string(argv[p+1]);
			init = 1;
			iter = 0;
			iterations = 0;
		}
		if(strcmp(argv[p], "-label") == 0){
                        labelfile = string(argv[p+1]);
                }
		if(strcmp(argv[p], "-lr") == 0){
			lr = atof(argv[p+1]);
			//if(lr < 0.0 || lr > 0.05) lr = 0.005;
		}
		if(strcmp(argv[p], "-lrforlo") == 0){
			lrforlo = atof(argv[p+1]);
		}
		if(strcmp(argv[p], "-h") == 0){
			helpmessage();
			exit(1);
		}
	}
	if(inputfile.size() < 1 || argc < 2){
		helpmessage();
		exit(1);
	}
	if(init == 1){
		initname = "RANDOM";
	}
	if(init == 2){
		initname = "FILE";
	}
	printf("\n");
        printf("@@@@@                   @      @@@@@@        @@@@@@ @     \n");
        printf("@    @                  @      @    @        @      @     \n");
        printf("@  @   @@@@@@ @  @@@@@@ @      @    @ @ @@@  @      @     \n");
        printf("@@@         @ @@ @    @ @@@@@@ @@@@@@ @@   @ @@@@@@ @     \n");
        printf("@  @   @@@@@@ @  @      @    @ @      @      @      @     \n");
        printf("@    @ @    @ @  @      @    @ @      @      @      @     \n");
        printf("@@@@@  @@@@@@ @@ @@@@@@ @    @ @      @      @@@@@@ @@@@@@\n\n");
	CSR<INDEXTYPE, VALUETYPE> A_csr;
        SetInputMatricesAsCSR(A_csr, inputfile);
        A_csr.Sorted();
	if(algoOption < 0 || algoOption > 8) algoOption = 2;
	if(omp_get_max_threads() < numberOfThreads){
		numberOfThreads = omp_get_max_threads();
		if(A_csr.rows < 1000000 && numberOfThreads > 32){
			numberOfThreads = 32;
		}
	}
	vector<VALUETYPE> outputvec;
	algorithms algo = algorithms(A_csr, inputfile, outputfile, labelfile, init, 1.0, lr, initfile, algoOption);
	algoname = "BATCHPREDNS";
        outputvec = algo.cacheBlockingminiBatchForceDirectedAlgorithm(iterations, numberOfThreads, batchsize, nsamples, lr, lrforlo, iter, scalingbox, psamples, pbox, expc);
	string avgfile = "Results.txt";
        ofstream output;
       	output.open(avgfile, ofstream::app);
	output << algoname << "\t" << initname << "\t";
       	output << iterations << "\t" << numberOfThreads << "\t" << batchsize << "\t";
        output << outputvec[0] << "\t" << outputvec[1];
	output << endl;
	output.close();
}
int main(int argc, char* argv[]){
	
	TestAlgorithms(argc, argv);
 	//myTest();       
	return 0;
}
