/*
 * dataPreparation.cpp
 *
 *  Created on: Oct 8, 2015
 *      Author: afarasat
 */
#include <iostream>
#include <math.h>
#include <cstdlib>
#include "dataPreparation.h"
#include <vector>
#include <stdlib.h>     /* srand, rand */
#include <ctime>

using namespace std;
std::vector<int> Y;
std::vector<int> observation;





dataPreparation::dataPreparation(int a, int b, int c) {
	_Q = a, _N = b, _K =c;
	_dataThresold = 0.05;
	_obsThresold = 0.8;
	long MAX_ELEM_ALLOWED = _Q*_N+1;
	int *p;
        p = new (nothrow) int[_Q*_N];
        long memtotal = _Q*_N;
        while((p == nullptr) && memtotal < MAX_ELEM_ALLOWED){
                memtotal += _N;
                p = new (nothrow) int[memtotal];
        }
        //if(p == nullptr) {cout << "MEMORY COULD NOT BE ALLOCATED!"; return 0;}
        delete[] p;
        observation.resize(memtotal);
	
	for (long k = 0; k < _Q*_N; k++){
		observation.push_back(0);
	}

}

dataPreparation::~dataPreparation() {
	// TODO Auto-generated destructor stub
}
 void dataPreparation::initialization() {

    std::vector<double> W(_Q*_K);
    std::vector<double> C(_K*_N);

	/*vector<vector<double> >  C1( _K, std::vector<double> ( _N, 0 ) );
	vector<vector<double> > W1( _Q, std::vector<double> ( _K, 0 ) );
	vector<vector<int> > Y ( _Q, std::vector<int> ( _N, 0 ) );*/


	long iterW, iterC;
	srand(time(NULL));
	//cout<< "W: "<< '\n';
	for (long i = 0; i < _Q; i++) {
		iterW = i*_K;
		for (long j = 0; j < _K; j++) {
			W[iterW+j]= 0.0;
			double randt = ((double) rand() / (RAND_MAX));
			if (randt < _dataThresold) {
				W[iterW+j] = 5 * ((double) rand() / (RAND_MAX));
			}
			//cout << W[i*_K+j] << " ";
		}
		//cout << '\n';
	}
	//cout<<"C: "<<endl;
	for (int i = 0; i < _K; i++) {
		iterC = i*_N;
		for (long j = 0; j < _N; j++) {
			C[iterC+j] = 2 * ((double) rand() / (RAND_MAX)) - 1;
			//cout << C[i*_N+j]<< " ";
		}
		//cout << '\n';
	}
	//cout<< "Main Data: "<< '\n';
	std::cout << "Obs Size: " << observation.size()<<"\n";
	//std::cout<<"Obs Max Siz: "<<observation.max_size()<<"\n";
	double sum;
	long iter_1, iter_2,iter_3;
	  for (long i = 0; i < _Q; i++){
                iter_1 = i*_K;
                iter_3 = i*_N;
                for (long j = 0; j < _N; j++){
                         sum = 0;
                         for (int k = 0; k < _K; k++) {
                                sum += W[iter_1+k]*C[j+k*_N];
                        }
                        if (sum > 0){
                                //Y[iter_3+j] = 1;
                                observation[iter_3+j] = 1; //CRS_store(i,j);
                                if ((double) rand() / (RAND_MAX)<_obsThresold/(float)3 ){
                                         observation[iter_3+j] = 3;
                                }
                        }else{
                                //Y[iter_3+j] = 0;
                                //observation[iter_3+j] = 0;
                                if ((double) rand() / (RAND_MAX)<_obsThresold ){
                                        observation[iter_3+j] = 3;
                                }
                        }
//void CRS_store(i,j);
     }
    }
	








	/*
	for (int i = 0; i < _Q; i++) {
		for (int j = 0; j < _N; j++) {
			double sum = 0;
			for (int k = 0; k < _K; k++) {
				sum += W[i*_K+k]*C[j+k*+_N];
			}
			if (sum > 0){
				_Y[i][j] = 1;
				_observation[i][j] = 1;
				if ((double) rand() / (RAND_MAX)<_obsThresold/5 ){
					_observation[i][j] = -1;
				}
			}else{
				_Y[i][j] = 0;
				_observation[i][j] = 0;
				if ((double) rand() / (RAND_MAX)<_obsThresold ){
					_observation[i][j] = -1;
				}

			}
			//cout<<_Y[i][j]<<" ";
		}
		//cout << endl;

	}

	//_Y =  & Y;
	*/
}


const std::vector<int>& dataPreparation::getMainData(){
	return  Y;
}
const std::vector<int>& dataPreparation::getObservation(){
	return  observation;
}

