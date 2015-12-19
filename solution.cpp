/*
 * solution.cpp
 *
 *  Created on: Oct 10, 2015
 *      Author: afarasat, utkarsh
 */

#include "solution.h"
# include <iostream>
#include <cstdlib>
#include <math.h>
#include <omp.h>
#include <time.h>
#include "DataRetriever.h"

//#include <mkl.h>

using namespace std;
std::vector<double> _W;
std::vector<double> _C;
long _Q=1, _N=1;
 int _K=1;
//int ** _observation;
int num_Procs;  // Number of processor;
//std::vector<int> observ;
DataRetriever *dataRet= new DataRetriever(_Q,_N,_K);
solution::solution() {

}

solution::solution( DataRetriever& dR,int & num_proc, double & sparsityThreshold, long & Q, long  & N, int & K) {
	num_Procs = num_proc;
	//_tbservation = observation;
	delete dataRet;
	dataRet=&dR;
 cout<<"netflix flag : "<<dataRet->ntflx;
////	dataRet->observation()	;
	cout<<"Initializing  W and C ... " << '\n';
	_Q = Q; _N = N; _K = K;
	double seed, randt;
        long i, j, m, n;
	_W.resize(_Q*_K);
	_C.resize(_K*_N);
	omp_set_num_threads(num_Procs);
	double timeS, timeE;
	timeS = omp_get_wtime();
	#pragma omp sections //shared(_W,_C,_Q,_N,_K,sparsityThreshold) //private(i,j,m,n,randt)
	{
		#pragma omp section 
		{
		srand(time(NULL)+omp_get_thread_num());
		#pragma omp parallel for shared(_W,_Q,_K,sparsityThreshold) private(i,randt)
		for (i = 0; i < _Q*_K; i++) {
		//	for (j = 0; j < _K; j++) {
				randt = ((double) rand() / (RAND_MAX));
				if (randt < 5*sparsityThreshold) {
					_W[i]= 10 * ((double) rand() / (RAND_MAX));
				}else{
					_W[i]=0;
				}
		//	}
		}
		}
		#pragma omp section
		{
	        srand(time(NULL)+omp_get_thread_num());
		cout<<"_C started"<<endl;
	        #pragma omp parallel for shared(_C,_N,_K) private(n,m)
		for (n = 0; n < _K*_N; n++) {
//			for (m = 0; m < _N; m++) {
				_C[n]=10 * ((double) rand() / (RAND_MAX)) - 1;
//			}
	
		}
		}
	}
	timeE = omp_get_wtime();
	std::cout<<"W and C initialized in "<<timeE-timeS<<" seconds"<<std::endl;
_W.shrink_to_fit();
_C.shrink_to_fit();
}
solution::~solution() {
	vector<double>().swap(_W);
	vector<double>().swap(_C);
	delete dataRet;	//.shrink_to_fit
//	_W.erase(_W.begin(),_W.end());
//	_C.erase(_C.begin(),_C.end());
}
std::vector<double> solution::getW(){
	return _W;
}
std::vector<double> solution::getC(){
	return _C;
}

double solution::logestfunct(float WiCj)
{
if (WiCj>=14)return 1;
else if (WiCj<=-14)return 0;
return (1/(1+exp(-WiCj)));
}

long solution::objectiveFunction (/*const DataRetriever& dataRet,*/ const double & lambda, const double & mu, const double & gamma){
	long objectiveFunction(0.0);
	long sumLog (0.0); long sumWl1(0.0); long sumWl2(0.0); long sumC(0.0);
	float WiCj;
	 bool obs_ij;             //observation data
        double logesticFunc;    // = 1/1+exp(...)
	long iter_1,iter_2,iter_3; 
	long i,j,k;
	omp_set_num_threads(num_Procs);
	#pragma omp sections
	{
	#pragma omp section
	{
	#pragma omp parallel for shared(_Q,_N,_K) private(i,j,obs_ij,k,WiCj,logesticFunc) reduction(+:sumLog)
	for (i = 0; i < _Q; i++){
		iter_1 = i * _K;
		for (j = 0; j < _N; j++){
			obs_ij = dataRet->observation(i,j);
			if (obs_ij >1){
				WiCj = 0.0;
				for (k=0; k < _K; k++){
					WiCj += _W[iter_1+k]*_C[j+k*_N];
				}
				logesticFunc=logestfunct(WiCj);
				sumLog += ((obs_ij) ? log(0.01+logesticFunc) : log(1.01-logesticFunc));
			}
		}
	  }
cout<<"-"<<endl;
	}
	#pragma omp section
	{
	#pragma omp parallel for shared(_Q,_K) private(i) reduction(+:sumWl1,sumWl2)
	for (i = 0; i < _Q*_K; i++){
		
//		for (k = 0; k < _K; k++){
                        sumWl1 += (_W[i] < 0) ? -_W[i] : _W[i];
                        sumWl2 += _W[i]*_W[i];
  //              }       
	}
cout<<"--"<<endl;
	}
	#pragma omp section
	{
	#pragma omp parallel for shared(_N,_K) private(i) reduction(+:sumC)
	for (i = 0; i < _K*_N; i++){
//		for (j = 0; j < _N; j++){
			sumC += _C[i]*_C[i];
//		}
	}
	}	
	}
	/*
	//Print elements of the objective function
        std::cout << "SumLog: " << sumLog<<std::endl;
	std::cout << "SumWl1: " << sumWl1<<std::endl;
        std::cout << "SumWl2: " << sumWl2<<std::endl;
        std::cout << "SumC: " << sqrt(sumC)<<std::endl;
	*/

	objectiveFunction = -sumLog + lambda*sumWl1 + mu*sumWl2 + gamma*sqrt(sumC);
	return objectiveFunction;
}
double solution::objectiveFunctionW (double  & lambda, double & mu){
	double objectiveFunction(0.0);
	int ** _observation;
	double sumLog (0.0); double sumWl1(0.0);  double sumWl2(0.0);
	for (int i = 0; i < _Q; i++){
		for (int j = 0; j < _N; j++){
			if (dataRet->observation(i,j) != -1){
				double WiCj(0.0);
				for (int k=0; k < _K; k++){
					WiCj += _W[i*_K+k]*_C[j+k*_N];
				}
				//sumLog+= -log(pow(1/(1+exp(-WiCj)),_observation[i][j])*pow(1-(1/(1+exp(-WiCj))),(1-_observation[i][j])));
				sumLog += _observation[i][j]* log(1/(1+exp(-WiCj))) + (1-_observation[i][j])*log(1-1/(1+exp(-WiCj)));
			}
		}
		double sumTemp(0.0);
		for (int k = 0; k < _K; k++){

			sumWl1 += (_W[i*_K+k] < 0) ? -_W[i*_K+k] : _W[i*_K+k];
			sumTemp += _W[i*_K+k]*_W[i*_K+k];
		}
		sumWl2 += pow(sumTemp,0.5);
	}

	objectiveFunction = -sumLog + lambda*sumWl1 + mu*sumWl2;
	return objectiveFunction;
}
double solution::objectiveFunctionC (double & gamma){
	double objectiveFunction(0.0);
	double sumLog (0.0); double sumC(0.0);
int ** _observation;
	for (int i = 0; i < _Q; i++){
		for (int j = 0; j < _N; j++){
			if (dataRet->observation(i,j) != -1){
				double WiCj(0.0);
				for (int k=0; k < _K; k++){
					WiCj += _W[i*_K+k]*_C[j+k*_N];
				}
				sumLog+= (_observation[i][j]* log(1/(1+exp(-WiCj)))+(1-_observation[i][j])*log(1-(1/(1+exp(-WiCj)))));
				//sumLog+= -log(pow(1/(1+exp(-WiCj)),observation[i][j])   *    pow(1-(1/(1+exp(-WiCj))),(1-observation[i][j])));
			}
		}
	}
	for (int i = 0; i < _K; i++){
		for (int j = 0; j < _N; j++){
			sumC += _C[i*_N+j]*_C[i*_N+j];
		}
	}
	objectiveFunction = -sumLog + gamma*sqrt(sumC);
	return objectiveFunction;
}
double solution::updateWij(/*const DataRetriever& dataRet,*/ const long & index_i, const long  & index_k, const double & mu, const double & stepSize){
	double delF(0.0);
	double WikCkj (0.0);
	double yipi[_N];
	double sum (0.0);
	long j, k;
	 bool obs_ij;
	omp_set_num_threads(num_Procs);
	#pragma omp parallel for private(j,k,WikCkj,yipi) reduction (+:sum)	
	for (j = 0; j < _N; j++){
		WikCkj = 0.0;
		for (k = 0; k < _K; k++){
				WikCkj += _W[index_i*_K+k]*_C[j+k*_N];
			}
		obs_ij=dataRet->observation(index_i,j);
		if (obs_ij >1){
			yipi[j] =obs_ij - logestfunct(WikCkj);//logest
		}else{
			yipi[j] =(rand()%2) - logestfunct(WikCkj);//logest
		}

		sum += _C[index_k*_N+j]*yipi[j];
	}
	//cout << "Sum->W: " << sum << endl;

	delF = -sum + mu*_W[index_i*_K+index_k];
	//cout << "DeltaF_w: " << delF << endl;
	//#pragma omp critical
	 _W[index_i*_K+index_k] = (double) rand() / (RAND_MAX) < 0.1 ? 0 :  _W[index_i*_K+index_k];
	_W[index_i*_K+index_k] =(_W[index_i*_K+index_k]-stepSize*delF) <= 0 ? 0 : (_W[index_i*_K+index_k]-stepSize*delF);
	_W[index_i*_K+index_k] = _W[index_i*_K+index_k] >= 20 ? 0 : _W[index_i*_K+index_k];
	//std::cout <<"W[i,j]: " << _W[index_i*_K+index_k] << std::endl;
	return _W[index_i*_K+index_k];
}
double solution::updateCij(/*const DataRetriever& dataRet,*/ const long & index_k, const long & index_j, const double & gamma,const double & stepSize){
	double delF(0.0);
	double WikCkj (0.0);
	double yipi[_Q];
	double sum (0.0);
	long j , k;
	bool obs_ij;
	omp_set_num_threads(num_Procs);
        #pragma omp parallel for private(j,k,WikCkj,yipi) reduction (+:sum)
	for (j = 0; j < _Q; j++){
		WikCkj = 0.0;
		for (k = 0; k < _K; k++){
				WikCkj += _W[j*_K+k]*_C[index_j+k*_N];
			}
		obs_ij=dataRet->observation(j,index_j);
		if (obs_ij >1){
			yipi[j] =obs_ij-logestfunct(WikCkj);//logest
		}else{
			yipi[j] = (rand()%2)-logestfunct(WikCkj);//logest
		}
		sum += _W[j*_K+index_k]*yipi[j];
	}
	//cout << "Sum->C: " << sum << endl;
	delF = -sum + gamma*_C[index_k*_N+index_j];
	//cout << "DeltaF_w: " << delF << endl;
	
	_C[index_k*_N+index_j] = (1/(1+gamma*stepSize))*(_C[index_k*_N+index_j]-stepSize*delF);
	_C[index_k*_N+index_j] =_C[index_k*_N+index_j] > 20 ? 2 * ((double) rand() / (RAND_MAX))-1 : _C[index_k*_N+index_j];
	 _C[index_k*_N+index_j] =_C[index_k*_N+index_j] < -20 ? 2 * ((double) rand() / (RAND_MAX))-1 : _C[index_k*_N+index_j];

	//_C[index_k][index_j] = _C[index_k][index_j]-stepSize*delF;
	// std::cout <<"C[i,j]: " << _C[index_k*_N+index_j] << std::endl;

	return _C[index_k*_N+index_j];
}

/*
 *

double * solution::updateRowWij(int & index_i, double & mu,double & stepSize){
	double delF [_K];
	double WikCkj (0.0);
	double yipi[_N];
	for (int j = 0; j < _N; j++){
		WikCkj = 0.0;
		for (int k = 0; k < _K; k++){
				WikCkj += _W[index_i][k]*_C[k][j];
			}
		if (_observation[index_i][j] != -1){
			yipi[j] =_observation[index_i][j] - (1/(1+exp(-WikCkj)));
		}else{
			yipi[j] =(rand()%2) - (1/(1+exp(-WikCkj)));
		}
	}
	double sum [_K];
	for (int k = 0; k < _K; k++){
		for (int n = 0; n < _N; n++){
			sum[k] += _C[k][n]*yipi[n];
		}
		delF[k] = -sum[k] + mu*_W[index_i][k];
		_W[index_i][k] =(_W[index_i][k]-stepSize*delF[k]) < 0 ? 0 : (_W[index_i][k]-stepSize*delF[k]);
	}

	return _W[index_i];
}

double * solution::updateColumnCij(int & index_j, double & gamma,double & stepSize){
	double delF[_K];
	double WikCkj (0.0);
	double yipi[_Q];
	for (int j = 0; j < _Q; j++){
		WikCkj = 0.0;
		for (int k = 0; k < _K; k++){
				WikCkj += _W[j][k]*_C[k][index_j];
			}
		if (_observation[j][index_j] != -1){
			yipi[j] = _observation[j][index_j]-(1/(1+exp(-WikCkj)));
		}else{
			yipi[j] = (rand()%2)-(1/(1+exp(-WikCkj)));
		}
	}
	double sum [_K];
	for (int k = 0; k < _K; k++){
		for (int n = 0; n < _Q; n++){
			sum[k] += _W[n][k]*yipi[n];
		}
		delF[k] = -sum[k] + gamma*_C[k][index_j];
		_C[k][index_j] = (1/(1+gamma*stepSize))*(_C[k][index_j]-stepSize*delF[k]);
	}

	//_C[index_k][index_j] = _C[index_k][index_j]-stepSize*delF;

	return _C[index_j];
}
 */
