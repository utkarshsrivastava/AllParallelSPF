/*
 * mainDriver.cpp
 *
 *  Created on: Oct 8, 2015
 *      Author: afarasat
 */
#include <iostream>
#include <string.h>
#include <math.h>
#include <vector>
#include <cstdlib>
#include "solution.h"
#include "dataPreparation.h"
#include "DataRetriever.h"

using namespace std;




int main(int arg, char* args[]) {
	// *****************************
	// Parameters Setting
	double mu (20.0); // This represent parameter of l2 regularization
	double lambda (20.0); // This represent parameter of l1 regularization
	double gamma (20.0); // This represent parameter of Frobenius norm
	static int numOfIteration (10); // Number of Iteration of main optimization loop
	double sparcityParameter = 0.1;
	double stepSize (20); // The stepsize in the optimization process 
	int numOfNodes;     // number of processes requested
	std::vector<double> I_index;
	std::vector<double> K_index;
 	std::vector<double> J_index;
	
	// Generating Random Indexes for whole otimization
	
	//*****************************
	int ** data; // full data
	int ** observation1; // observation = data-missing values
	int Q(9429),N(17770), K(10);
	// Generating Random Indexes for whole otimization
	for (int index = 0; index < numOfIteration*numOfNodes; index++){
                I_index[index] =(rand() % (Q-1));
		K_index[index] = (rand() % (K-1));
		J_index[index] = (rand() % (N-1));
	}
	std::vector<double> W;
	std::vector<double> C;
	char workinDirectory[1024]="/gpfs/user/afarasat/CSE603Project/NetFlixData/Files";
	DataRetriever d(workinDirectory);
	//dataPreparation dataObj(Q,N,K);
	//dataObj.initialization();
	//data = dataObj.getMainData();
	//observation1 = dataObj.getObservation();
	cout<<"****************"<<" OBSERVATION "<<endl;
	for (int i=0; i < Q; i++){
		for (int j =0; j < N; j++){
			//cout << observation[i][j] << " ";
		}
		//cout << '\n';
	}

	solution solutionObj(&d, sparcityParameter, Q, N, K);

	W = solutionObj.getW();
	C = solutionObj.getC();

	cout << "Objective function: " << solutionObj.objectiveFunction(lambda, mu,gamma) << endl;

	cout<<"****************"<<" Optimization "<<endl;

	for (int currentIteration = 0;currentIteration < numOfIteration; currentIteration++){
		int iW =(rand() % (Q-1));
		int kW = (rand() % (K-1));
		W[iW*K+kW] = solutionObj.updateWij(iW,kW,mu,stepSize);
		//W[iW]  = solutionObj.updateRowWij(iW,mu,stepSize);
		cout << "W[i][k]: " << W[iW*K+kW] << endl;
		int kC =(rand() % (K-1)); int jC = (rand() % (N-1));
		C[kC*N+jC] = solutionObj.updateCij(kC,jC,gamma,stepSize);
		//C[jC] = solutionObj.updateColumnCij(jC,gamma,stepSize);
		cout << "C[i][j]: " << C[kC*N+jC] << std::endl;
		//cout << "Objective function: " << solutionObj.objectiveFunction(lambda, mu,gamma) << endl;
		stepSize = stepSize/1.2;
		mu = mu/1.2;
		gamma = gamma/1.2;
	}
	for (int i=0; i < Q; i++){
		for (int j =0; j < K; j++){
		//	cout << W[i*K+j] << " ";
		}
		//cout << '\n';
	}

return -1;
}


