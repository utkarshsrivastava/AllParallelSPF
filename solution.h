/*
 * solution.h
 *
 *  Created on: Oct 10, 2015
 *      Author: afarasat
 */

#ifndef SOLUTION_H_
#define SOLUTION_H_
#include <math.h>
#include<iostream>
#include<vector>
#include "DataRetriever.h"
class solution {
public:
	solution();
	solution( DataRetriever &dataRet,int &,double &,long &,long &,int &);
	std::vector<double> getW();
	std::vector<double> getC();
	double logestfunct(float WiCj);
	long objectiveFunction(/*const DataRetriever& dataRet,*/ const double &, const double &,const double &);
	double objectiveFunctionW(double &, double &);
	double objectiveFunctionC(double &);
	double updateWij(/*const DataRetriever& dataRet,*/ const long &,const long &,const double &,const double &);
	double updateCij(/*const DataRetriever& dataRet,*/ const long &,const long &, const double &,const double &);
	double * updateRowWij(int &,double &,double &);
	double * updateColumnCij(int &,double &,double &);

	virtual ~solution();
private:
	//int _Q, _N, _K;
	//std::vector<double> _W;
	//std::vector<double> _C;
	//int ** _observation;
};

#endif /* SOLUTION_H_ */
