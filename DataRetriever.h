#ifndef DATARETRIEVER_H
#define DATARETRIEVER_H



#include <iostream>
#include <cstring>
#include <sstream>
#include <dirent.h>
#include <fstream>
#include <stdlib.h>
#include <vector>


class DataRetriever{
public:
bool observation(long row_num,long col_num);
DataRetriever(char* dir_name);
DataRetriever(long &_Q, long &_N, int &_K);
~DataRetriever();
void CRS_store (long i,long j,bool ratn);
void initialization();
const std::vector<bool>& getObservation();
bool ntflx=1;
long dataSize();

private:
 int dataLoad(char *wd);
 void CRW_colm_rate(long row_num, long col,long col_num,bool rating);
        double _dataThresold;
        double _obsThresold;
        int ** _Y;
        int ** _observation;


};
#endif
