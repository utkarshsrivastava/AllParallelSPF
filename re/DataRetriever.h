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
bool observation(int row_num,long col_num);
DataRetriever(char* dir_name);

private:
 int dataLoad(char *wd);
 void CRW_colm_rate(int row_num, long col,long col_num,int rating);

};
#endif
