
//============================================================================
// Name        : DataRetriever.cpp
// Author      : Utkarsh Srivastava
// Version     :
// Copyright   :
// Description : CRS implementation in C++, Ansi-style
//============================================================================

#include <iostream>
#include <cstring>
//#include <sys/types>
#include <sstream>
#include <dirent.h>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <regex.h>
#include "DataRetriever.h"
using namespace std;

#ifdef WINDOWS
	#include<direct.h>
	#define getCWD _getcwd
#else
	#include<unistd.h>
	#define getCWD getcwd
#endif
std::vector<bool> observ; 
std::vector<long> rows;
std::vector<long> cols;
std::vector<bool> ratings;
long N,Q,K;
const long MAX_ENTITIES=100481000;
const long MAX_MOVIES=20000;
//bool ntflx=1;
DataRetriever::~DataRetriever(){
std::vector<bool>().swap(observ);
std::vector<long>().swap( rows);
std::vector<long>().swap( cols);
std::vector<bool>().swap( ratings);

}


DataRetriever::DataRetriever(char* dir_name)
{
ntflx=1;
	 if (dir_name!=NULL && (1==dataLoad(dir_name)))
		{
	cout<<"Data from the folder was read successfully !"<<endl;	
		}

}

DataRetriever::DataRetriever(long &_Q, long &_N, int &_K)
{
ntflx=0;
       N=_N;Q=_Q;K=_K;
	 _dataThresold = 0.05;
        _obsThresold = 0.8;
        long MAX_ELEM_ALLOWED = _Q*_N+1;
        bool *p;
        p = new (nothrow) bool[_Q*_N];
        long memtotal = _Q*_N;
	short int memchekr=0;
        while((p == nullptr)&& memchekr<50 ){
                memtotal += _N;
		memchekr++;
                p = new (nothrow) bool[memtotal];
        }
delete[] p;
        observ.resize(memtotal);
	//ratings.resize(memtotal)'
	observ.shrink_to_fit();
        for (long itrtr = 0; itrtr < _Q*_N; itrtr++){
	    observ.push_back(0);
		//ratings.push_back(0);
        }
        
//if(p == nullptr) {cout << "MEMORY COULD NOT BE ALLOCATED!"; return 0;}
        //        delete[] p;
        //                observation.resize(memtotal);
        //
        //                        for (long k = 0; k < _Q*_N; k++){
        //                                        observation.push_back(0);
        //                                                }
        //
}
int DataRetriever::dataLoad(char *wd)
{
	long itr_vector;
	for (itr_vector=0;itr_vector<MAX_ENTITIES;itr_vector++)
	{
		cols.push_back(0);
		ratings.push_back(0);
	}

	for (itr_vector=0;itr_vector<MAX_MOVIES;itr_vector++)
		{
		rows.push_back(0);
		}
	DIR *dirp;
	struct dirent *dp;
std::string line="";
ifstream tmp_file;	
			
cout<<wd<<": working dir"<<endl;
if (NULL==wd) {cout<<"dir aint need No NULLS ! !";return -1;}

dirp=opendir(wd);
if (NULL==dirp)
{
	cout<<"Error opening path"<<wd;return -1;
}


int files_lredy_parsed=0;
if (dirp<=0){cout<<"Error opening current path"<<wd;return -1;}
//update row_ptr////////////
bool rowptr_done;
//if (this is UNIX)
char cmnd[1124]="rm ";//' will come here 
strcat(cmnd,"./records.count");
cout<<"if This is the first run, ignore the below error:No file found .."<<endl;system(cmnd);
memset(&cmnd[0],0,sizeof(cmnd));
strcat(cmnd," wc -l ");
strcat(cmnd,"\"");
strcat(cmnd,wd);
strcat(cmnd,"\"");
strcat( cmnd,"/*.txt | ");
strcat( cmnd,"grep -o '[0-9]\\+ \\|mv_[0-9]*' | ");
strcat( cmnd,"sed ':a;N;$!ba;s/\\nmv_/ /g' | ");
strcat( cmnd,"head -n -1 > ");// '
strcat( cmnd,"./records.count");
cout<<cmnd<<endl;
rowptr_done=system(cmnd);
if (!dirp) {return -1;cout<<"dirp failed \n";}


tmp_file.open("./records.count");


int total_lines,rownum=0;
if (tmp_file== NULL || !tmp_file){return -1;}
 std::stringstream iss(line);
while (std::getline(tmp_file,line))
{
//int k=0;
	iss.clear();
	iss<<line;
	if (!(iss>>total_lines>>rownum)){cout<<line<<total_lines<<endl;}
	rows.at(rownum)=total_lines-1;
	iss<<"";

}

//tmp_file.close();
cout<<"rows Set \n";
long tmp=0;
unsigned long rowitr =1;
long offset=0;
while(rowitr<rows.size())
{
tmp=rows[rowitr+1];
rows[rowitr+1]=offset;
offset+=tmp;
rowitr++;
//if (tmp!=0)cout<<tmp<<"ffset : "<<offset<<"arr val :"<<rowitr<<"-"<<rows[rowitr]<<endl;
}


tmp_file.close();

//column update now///////////////////////////
//if (1) return -1; //tobremoved

long col_num;
bool rating;
char mv_day[2],mv_month[2],mv_year[4];
bool isRow;
long col;
long row_num;
char tmp_ful_filename[1024];
bool tmpbool=1;
size_t colon_pos=0;
std::string rowunfiltered;
std::string file_col_values;
//std::istringstream iss;
	while (dirp!=NULL && (dp = readdir(dirp)) != NULL)
	{
		*tmp_ful_filename=*"";
		strcpy(tmp_ful_filename,wd);

		if ( ((dp->d_name)[0]== '.') && ((dp->d_name)[1]=='\0'))
		            continue;
		if ( ((dp->d_name)[0]== '.') && ((dp->d_name)[1]=='.'))
				    continue;
		++files_lredy_parsed;
	//	cout<<dp->d_name<<++files_lredy_parsed<<endl;
				if(wd[sizeof(wd)-1]=='/')
				{strcat(tmp_ful_filename,dp->d_name);}
				else
				{
					strcat(tmp_ful_filename,"/");
					strcat(tmp_ful_filename,dp->d_name);
				}


	tmp_file.open(tmp_ful_filename);
	if (tmp_file== NULL || !tmp_file){cout<<"not able to open tmpfile : "<<tmp_ful_filename<<endl<<wd<<endl;continue;}

	    isRow=1;
	    col_num=0;
	    row_num=0;
	    //check
	    		if (!std::getline(tmp_file, line)){cout<<"Couldnot read file :"<<tmp_file.uppercase<<endl;return -1;}
	    		  iss<<line;
	    		 iss >> rowunfiltered;
	    		 	colon_pos=rowunfiltered.find(':',0);
	    		 rowunfiltered=rowunfiltered.substr(0,colon_pos);//,':','_');
	    		 sscanf(rowunfiltered.c_str(),"%ld",&row_num); //convert char to int
			iss<<"";	iss.clear();
		while (std::getline(tmp_file, line) && !tmp_file.eof())
		{
			//cout<<"reading "<<line<<endl;
			file_col_values.clear();
			//line.substr(0,line.find(',',0)
		     iss<<line;
		    std::getline(iss,file_col_values,',');
		    col=atol(file_col_values.c_str());
		    std::getline(iss,file_col_values,',');
		    rating=atol(file_col_values.c_str());
		    std::getline(iss,file_col_values,',');
//if (files_lredy_parsed<=2){
//cout<<"TESTING : "<<col<<" "<<rating<<" " <<tmp_file<<" "<<line << endl;
//}
		    //date=file_col_values.c_str());
		    //
		    if (rating>5) {break;}//error
		    //if (!(iss >> col >> rating >>mv_day>>mv_month>>mv_year)) { break; } // error
		    CRW_colm_rate(row_num,col,col_num,rating);
		    col_num++;
			iss<<"";iss.clear();
		}
		tmp_file.close();
	}
closedir(dirp);

//testing :
/*
 * long j=rows[16154];
while(j<rows[16155]){
	cout<<cols[j]<<"-"<<ratings[j]<<"-"<<observation(16154,cols[j])<<endl;
j++;
}

cout<<rows[16154]<<"<-255: 256->"<<rows[16155]<<endl<<"done !";
*/
	return 1;
}



void DataRetriever::CRW_colm_rate(long row_num, long col,long  col_num,bool rating)
{
cols[rows[row_num]+col_num]=col;
ratings[rows[row_num]+col_num]=rating/4;
//if (row_num==258)cout<<"rows: "<<rows[row_num]<<": "<<row_num<<" col "<<col<<" col_num "<<col_num<<" rating "<<rating<<" date : "<<file_col_values<<endl;
//if (row_num==258)cout<<"rows: "<<rows[row_num]<<": "<<row_num<<" col "<<col<<" col_num "<<col_num<<" rating "<<rows[row_num]+col_num<<"-"<<ratings[rows[row_num]+col_num]<<endl;

}
void DataRetriever::initialization()
{
    std::vector<double> W(Q*K);
    std::vector<double> C(K*N);
 double randt;
        
        long iterW, iterC;
        srand(time(NULL));
        for (long i = 0; i < Q*K; i++) {
//                iterW = i*K;
  //              for (long j = 0; j < K; j++) {
                        W[i]= 0.0;
                         randt = ((double) rand() / (RAND_MAX));
                        if (randt < _dataThresold) {
                                W[i] = 5 * ((double) rand() / (RAND_MAX));
                        }
    //            }
        }
        for (int i = 0; i < K*N; i++) {
        //        iterC = i*N;
          //      for (long j = 0; j < N; j++) {
                        C[i] = 2 * ((double) rand() / (RAND_MAX)) - 1;
                
        }
        std::cout << "Obs Size: " << observ.size()<<"\n";
        double sum;
        long iter_1, iter_2,iter_3;
          for (long i = 0; i < Q; i++){
                iter_1 = i*K;
                iter_3 = i*N;
                for (long j = 0; j < N; j++){
                         sum = 0;
                         for (int k = 0; k < K; k++) {
                                sum += W[iter_1+k]*C[j+k*N];
                        }
                        if (sum > 0){
                                observ[iter_3+j] = 1; 
                                if ((double) rand() / (RAND_MAX)<_obsThresold/(float)3 ){
                                         observ[iter_3+j] = 3;
                                }
                        }else{
                                if ((double) rand() / (RAND_MAX)<_obsThresold ){
                                        observ[iter_3+j] = 3;
					}
				}
//CRS_store(i,j,observ[iter_3+j]);
		}
	}
}
void DataRetriever::CRS_store (long i,long j, bool ratn){
long itrtr=0;
long col_index=0;

if (rows[i]==0)
{
        while(itrtr<N || cols[rows[i-1]+itrtr]!=0)
        itrtr++;
}
rows[i]=rows[i-1]+itrtr;

itrtr=0;
//col_index=0;

while(col_index<N || cols[rows[i]+col_index]!=0)
col_index++;

cols[rows[i]+col_index]=j;
ratings[rows[i]+col_index]=ratn;
}

const std::vector<bool>&DataRetriever::getObservation(){
        return  observ;
}

long  DataRetriever::dataSize(){
return ratings.size();
}

bool DataRetriever::observation(long row_num,long col_num) 
{
//bool returnedrating=3;
if(ntflx){
long col=3;
long index=rows[row_num];
while (cols[index]!=col_num && index<rows[row_num+1])
{index++;}

if (index>=rows[row_num+1]){return 3;}
return  ratings[index];
}
else{
return observ[row_num*N+col_num];
}
}

