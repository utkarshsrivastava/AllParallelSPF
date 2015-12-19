/*
 * mainDriver.cpp
 *
 *  Created on: Oct 8, 2015
 *      Author: afarasat, utkarsh
 */
#include <iostream>
#include <string.h>
#include <math.h>
#include <vector>
#include <cstdlib>
#include "solution.h"
#include "dataPreparation.h"
#include "DataRetriever.h"
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <fstream>

using namespace std;

// std::vector<int> observation(2600000*20000);


int main(int argc, char* args[]) {
	// *****************************
	// Parallel Setup
//	const long  MAX_ELEM_ALLOWED = 992147483640;
	int NUM_PROCS, NUM_THREADS;
	//NUM_PROCS = omp_get_num_procs(); // number of processes requested
	NUM_PROCS = atoi(args[2]);
	int P = NUM_PROCS;
	//NUM_PROCS = (int)(args[0]);
	long Q ;
	long N =17770;
	int K ;
	int runMax = atoi(args[3]);
	 K = atoi(args[4]==NULL?"20":args[4]);
	bool netFlix=atoi(args[1]);
	NUM_THREADS = NUM_PROCS;
	omp_set_num_threads(NUM_THREADS); 
	// 	USAGE : ./output netflix NUM_PROCS runMax K Q N
	// Parameters Setting
	float mu = 1.0; // This represent parameter of l2 regularization
	float init_mu = mu;
	float lambda = 1.0; // This represent parameter of l1 regularization
	float gamma = 1.0; // This represent parameter of Frobenius norm
	float init_gamma = gamma;
	static int numOfIteration = 2000; // Number of Iteration of main optimization loop
	double sparcityParameter =0.05 ;
	float stepSize; // The stepsize in the optimization process 
	float init_stepSize = 1;
	// number of processes requested

	
	// Generating Random Indexes for whole otimization
	
	//*****************************
	std::vector<int> observation;
	/*
	int *p;
	p = new (nothrow) int[Q*N];
	long memtotal=Q*N;
	while((p == nullptr) && memtotal < MAX_ELEM_ALLOWED){
		memtotal += N;
		p = new (nothrow) int[memtotal];
	}
	if(p == nullptr) {cout << "MEMORY COULD NOT BE ALLOCATED!"; return 0;}
	delete[] p;
	observation.resize(memtotal);
	*/
	
	//std::vector<int> data(Q*N);
	//int ** data; // full data
	//int ** observation; // observation = data-missing value
	std::vector<double> W;
	std::vector<double> C;
	int bestItr;
	long currObj, bestObj, initObj;
	double bestTime;
	int bestITR [runMax];
	long bestOBJ[runMax];
	double bestTIME[runMax];
	double totTIME[runMax];
	long improvObj[runMax];
	std::cout <<" USAGE : ./output 0 NUM_PROCS runMax K Q N or  USAGE : ./output 1 NUM_PROCS runMax K <dirpath> " 
			<<std::endl<<"Data Genration section started ..."<< std::endl;
			//	dataPreparation dataObj(Q,N,K);//if dataPreparation
 DataRetriever *drptr=nullptr;	
	if (netFlix)
	{
	
/*//		K = atoi(args[4]);
		char *nfpath=new char[256];
		*nfpath =*"/gpfs/courses/cse603/students/ITS/ParallelSpraseMatrixFactorization/AlisOpenMP/NewVersion1/NetFlix/";
        	 if(argc>=5 && args[5]!=NULL)
		{
		strncpy(nfpath,args[5],sizeof(nfpath));
		}
		cout<<"dir path : "<<*nfpath<<endl */;
		DataRetriever dR("/gpfs/courses/cse603/students/ITS/ParallelSpraseMatrixFactorization/AlisOpenMP/NewVersion1/NetFlix/");	
		drptr=&dR;
		Q=100480507/17770 +1;
		cout<< " dataSize :"<<Q*N<<endl;
	}
	else
	{
          Q = atol(args[5]==NULL?"3500":args[5]);
          N =atol(args[6]==NULL?"17771":args[6]);
cout<<"K Q N"<<K<<" "<<Q<<N<<std::endl;
		DataRetriever dR(Q,N,K);
		dR.initialization();
		drptr=&dR;
	}
	int index;
    	double timeStart,timeEnd;
	int num_rand = numOfIteration*NUM_PROCS;
	 // Generating Random Indexes for whole otimization
	
	std::cout << "data is ready now!!!"<<num_rand<<std::endl;
	std::vector<double> I_index(num_rand);
	std::vector<double> K_index(num_rand);
	std::vector<double> J_index(num_rand);
	double timeS_Obj, timeE_Obj;
	int numRun;
	std::ofstream myFile;
	try{
		myFile.open("Results//Q"+std::to_string(Q)+"N"+std::to_string(N)+"K"+std::to_string(K)+"P"+std::to_string(P)+".txt");
	}
	catch (exception& e){
		std::cout << "c++11 is required to use std::to_string \n Usage : icc -o3 -openmp -std=C++11 <files>" << std::endl;
	}
	int currentIteration;
	int counter;	
	cout<<"****************"<<" OBSERVATION "<<endl;
	for (numRun = 0; numRun < runMax; numRun++){
		
		mu = init_mu ;
        	gamma = init_gamma;
        	stepSize = init_stepSize; // The stepsize in the optimization process
		timeStart = omp_get_wtime();
		#pragma omp parallel shared(I_index,K_index,J_index,numOfIteration,NUM_PROCS)
		{
		srand(time(NULL)+omp_get_thread_num()); 
		#pragma omp for private(index) 	
			for (index = 0; index < num_rand; index++){
					I_index[index] =(rand() % (Q));
					K_index[index] = (rand() % (K));
					J_index[index] = (rand() % (N));
			}
		}
			if (drptr==NULL) return -1;///SANITY CHECK
		solution solutionObj(*drptr, NUM_PROCS,sparcityParameter,Q,N,K);// Q*N=711084000 records , so greater than 65000
		
	std::cout << "data is ready now _ !!"<<drptr->dataSize()<<std::endl;
//		W = solutionObj.getW();
//		C = solutionObj.getC();	
		timeS_Obj = omp_get_wtime();
		currObj =  solutionObj.objectiveFunction(lambda, mu,gamma);
		initObj = currObj;
		//myFile << "Run: " << numRun << std::endl;
		//myFile << "Iteration" << "\t" << "Objective Function" << std::endl;
		//myFile << "Initial Value" << "\t" << currObj << std::endl;
		//cout << "Objective function: " << currObj << endl;
		bestObj = currObj;
		bestItr = 0;
		bestTime = timeStart;
		timeE_Obj = omp_get_wtime();
		std::cout << "Time Objective Function: " << timeE_Obj - timeS_Obj <<  std::endl;
		std::cout<<"****************"<<" Optimization Started..."<<endl;
		currentIteration = 0;
		stepSize = 2/(currentIteration+2);
		counter = 1;
		while (currentIteration < numOfIteration){
			//NUM_PROCS = 1;
			omp_set_num_threads(NUM_PROCS);
			int tid;
			#pragma omp parallel shared(W,C)
			{ 
			int iW, kW, kC, jC;
			//#pragma omp for schedule(dynamic) private(iW, kW, kC, jC,procIterator)
			//for (procIterator = 0; procIterator < NUM_PROCS; procIterator++){
				tid = omp_get_thread_num();
				//std::cout<< "ThreadID: " << tid << " rand: " << I_index[currentIteration*NUM_PROCS+tid]<< std::endl;
				iW = I_index[currentIteration*NUM_PROCS+tid];
				kW = K_index[currentIteration*NUM_PROCS+tid];
				jC = J_index[currentIteration*NUM_PROCS+tid];
				//iW =(rand() % (Q));
				//kW = (rand() % (K));
				//printf("process id: %d Iterator: %d iW: %d kW: %d\n",tid,procIterator,iW,kW); 
			
				W[iW*K+kW] = solutionObj.updateWij(iW,kW,mu,stepSize);
				//W[iW]  = solutionObj.updateRowWij(iW,mu,stepSize);
				kC = kW; //jC = (rand() % (N));
				C[kC*N+jC] = solutionObj.updateCij(kC,jC,gamma,stepSize);
				//C[jC] = solutionObj.updateColumnCij(jC,gamma,stepSize);
				//cout << "C[i][j]: " << C[kC][jC] << endl;
			//}
			}
			//#pragma omp barrier	
			currObj = solutionObj.objectiveFunction(lambda, mu,gamma);
		//	myFile << currentIteration  << "\t" << currObj << std::endl;
			//cout << "Iteration: " << currentIteration << ", Objective function: " << currObj << endl;
			if (currObj < bestObj){
				bestObj = currObj;
				bestItr = currentIteration;
				bestTime = omp_get_wtime();
				bestITR[numRun] = bestItr;
				bestOBJ[numRun] = bestObj;
				bestTIME[numRun] = bestTime - timeStart;
				improvObj[numRun] = initObj - bestObj;
			}
			stepSize = 2/(currentIteration+2);	
		//	stepSize /= 1.1;
		//	mu = mu/1.05;
		//	gamma = gamma/1.05;
			/*
			if(currentIteration % 200 == 0 && currentIteration != 0 ){
				counter++;
				stepSize = ((float)1/(float)counter) * init_stepSize;
				mu = ((float)1/(float)counter) * init_mu;
				gamma = ((float)1/(float)counter) * init_gamma;
			} 
			*/
			currentIteration++;
			
		}
		//std::cout << "mu:" << mu << " step:" << stepSize << "gamma:" << gamma << std::endl;
		timeEnd =omp_get_wtime();
		//myFile << "Best Objective Function: " << bestObj << "\t" << "Best Iteration: " << bestItr << "\t" << "Best Time" << bestTime-timeStart << std::endl;
		//myFile << "Num of Procs: " <<  NUM_PROCS  << "\t" << " Total time: " << timeEnd-timeStart << std::endl;
		//std::cout << "Num of Procs: " <<  NUM_PROCS  <<" Total time: " << timeEnd-timeStart << std::endl;
		totTIME[numRun] = timeEnd-timeStart;
		solutionObj.~solution();
		vector<double>().swap(W);
		vector<double>().swap(C);
		//W.erase(W.begin(),W.end());
		//C.erase(C.begin(),C.end());
		
	}
	double sumBestOBJ = 0.0;
	double sumbestTime = 0.0;
	double sumTotTIME = 0.0;
	int sumBestITR = 0;
	double sumImprov = 0.0;
	myFile << "Best Objective Function" << "\t" << "Best Iteration" << "\t" << "Best Time" << "\t" << "Total Time"<< "\t" << "Improvement" << std::endl;
	for (int i = 0; i < runMax; i++){
		sumBestOBJ +=  bestOBJ[i];
		sumBestITR +=  bestITR[i];
		sumbestTime += bestTIME[i];
		sumTotTIME += totTIME[i];
		sumImprov += improvObj[i]; 
		myFile << bestOBJ[i]  << "\t" << bestITR[i] << "\t" << bestTIME[i] << "\t" << totTIME[i] << "\t" << improvObj[i] << std::endl;
	}
	myFile << "Average Objective Function" << "\t" << "Average Iteration" << "\t" << "Average Best Time" << "\t" << "Average Total Time"<< "\t"<<"Average Improvement"<< std::endl;
	myFile <<(float) sumBestOBJ/(float)runMax  << "\t" <<(float)sumBestITR/(float)runMax << "\t" << (float)sumbestTime/(float)runMax << "\t" << (float)sumTotTIME/(float)runMax << "\t"<<(float)sumImprov/(float)runMax <<std::endl;
	myFile.close();

return -1;
}
