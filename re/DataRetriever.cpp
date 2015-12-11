
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

//class DataRetriever{
std::vector<long> rows;
std::vector<long> cols;
std::vector<bool> ratings;
/*
public :
bool observation(int row_num,long col_num);*/
DataRetriever::DataRetriever(char* dir_name)
 {
	 if (dir_name!=NULL && dir_name)
 	dataLoad(dir_name);

}
/*

private:
 int dataLoad(char *wd);
 void CRW_colm_rate(int row_num, long col,long col_num,int rating);

};*/

/*
int main(){
DataRetriever d("/gpfs/user/afarasat/CSE603Project/NetFlixData/Files");//	/gpfs/user/afarasat/CSE603Project/NetFlixData/download/mok_train");
return -1;
}*/

int DataRetriever::dataLoad(char *wd)
{
	long itr_vector;
	for (itr_vector=0;itr_vector<100000000;itr_vector++)
	{
		cols.push_back(0);
		ratings.push_back(0);
	}

	for (itr_vector=0;itr_vector<18000;itr_vector++)
		{
		rows.push_back(0);
		}
	DIR *dirp;
	struct dirent *dp;
std::string line;
ifstream tmp_file;								//cout << "!!!Hello World!!!" << endl; // prints !!!Hello World!!!
cout<<wd<<endl;
if (NULL==wd) {cout<<"dir  aint No NULLS ! !";return -1;}

dirp=opendir(wd);
if (NULL==dirp)
{
	cout<<"Error opening path"<<wd;return -1;
}


	//cout<<" Enter data source path:";
	//cin.getline(wd,sizeof(wd));

int files_lredy_parsed=0;
//dirp=opendir(wd);
if (dirp<=0){cout<<"Error opening current path"<<wd;return -1;}
//update row_ptr////////////
bool rowptr_done;
//if (this is UNIX)
char cmnd[1124]="rm ";//' will come here 
//strcat(cmnd,wd);
strcat(cmnd,"./abcd.txt");
system(cmnd);
memset(&cmnd[0],0,sizeof(cmnd));
strcat(cmnd," wc -l ");
strcat(cmnd,"\"");
strcat(cmnd,wd);
strcat(cmnd,"\"");
strcat( cmnd,"/*.txt | ");
strcat( cmnd,"grep -o '[0-9]\\+ \\|mv_[0-9]*' | ");
strcat( cmnd,"sed ':a;N;$!ba;s/\\nmv_/ /g' | ");
strcat( cmnd,"head -n -1 > ");// '
//strcat( cmnd,wd);
strcat( cmnd,"./abcd.txt");
cout<<cmnd<<endl;
rowptr_done=system(cmnd);
//if (1) return -1;

if (!dirp) return 0;


tmp_file.open("./abcd.txt");

//if (!rowptr_done){cout<"file->total records : not populated ";return -1;}

int total_lines,rownum=0;
if (tmp_file== NULL || !tmp_file){return -1;}

while (std::getline(tmp_file,line))
{

	std::istringstream iss(line);
	if (!(iss>>total_lines>>rownum)){cout<<total_lines<<endl;}
	rows.at(rownum)=total_lines-1;
}


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
int rating;
char mv_day[2],mv_month[2],mv_year[4];
bool isRow;
long col;
long row_num;
char tmp_ful_filename[1024];
bool tmpbool=1;
size_t colon_pos=0;
std::string rowunfiltered;
std::string file_col_values;
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
	    		 std::istringstream iss(line);
	    		 iss >> rowunfiltered;
	    		 	colon_pos=rowunfiltered.find(':',0);
	    		 rowunfiltered=rowunfiltered.substr(0,colon_pos);//,':','_');
	    		 sscanf(rowunfiltered.c_str(),"%ld",&row_num); //convert char to int

		while (std::getline(tmp_file, line) && !tmp_file.eof())
		{
			//cout<<"reading "<<line<<endl;
			file_col_values.clear();
			//line.substr(0,line.find(',',0)
		    std::istringstream iss(line);
		    std::getline(iss,file_col_values,',');
		    col=atol(file_col_values.c_str());
		    std::getline(iss,file_col_values,',');
		    rating=atol(file_col_values.c_str());
		    std::getline(iss,file_col_values,',');
		    //date=file_col_values.c_str());
		    if (rating>5) {break;}//error
		    //if (!(iss >> col >> rating >>mv_day>>mv_month>>mv_year)) { break; } // error
		    CRW_colm_rate(row_num,col,col_num,rating);
		    col_num++;
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



void DataRetriever::CRW_colm_rate(int row_num, long col,long  col_num,int rating)
{
cols[rows[row_num]+col_num]=col;
ratings[rows[row_num]+col_num]=rating/4;
//if (row_num==258)cout<<"rows: "<<rows[row_num]<<": "<<row_num<<" col "<<col<<" col_num "<<col_num<<" rating "<<rating<<" date : "<<file_col_values<<endl;
//if (row_num==258)cout<<"rows: "<<rows[row_num]<<": "<<row_num<<" col "<<col<<" col_num "<<col_num<<" rating "<<rows[row_num]+col_num<<"-"<<ratings[rows[row_num]+col_num]<<endl;

}
bool DataRetriever::observation(int row_num,long col_num)
{
bool returnedrating=-1;
long col=-1;
long index=rows[row_num];
while (cols[index]!=col_num && index<rows[row_num+1])
{index++;}

if (index>=rows[row_num+1]){return -1;}
return  ratings[index];
}

