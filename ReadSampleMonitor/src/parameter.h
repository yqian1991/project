/*
 * parameter.h
 *
 *  Created on: 2013-4-22
 *      Author: qianyu
 */

#ifndef PARAMETER_H_
#define PARAMETER_H_

#include "genePreDefinition.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <string.h>
#include <sys/time.h> 	//sys/time.h是Linux系统的日期时间头文件，sys/time.h通常会包含include time.h。通常如果代码可以是平台相关的，则只需要include sys/time.h。
#include <time.h>		//编写的代码如果是平台无关的，则需要在代码里include time.h，
#include <fcntl.h>
#include <unistd.h>

#define SNP_POSSIBILITY 1000 //1/1000
#define PARA_NUM_COLS 3  //means that each sample must have 3 parameters
#define FILE_NUM_COLS 6  //means that all the cols in a sampleParameters file, including two TAB and \0
#define CYCLIZLIZE_GENOME_LENGTH 0x10000000
//#define SAMPLE_DIRECTORY    "../sample/"      // store the sample results in this directory
//#define SAMPLE_PARAMETER_PATH "parameters/sampleParameters.txt"  //this file stores the parameters set by users
//#define GLOBAL_PARAMETER_PATH "parameters/globalParameters.txt"  //this file stores the parameters set by users
#define SAMPLE_PARAMETER_PATH "/home/xuzhen/qy/ReadSampleMonitor/parameters/sampleParameters.txt"  //this file stores the parameters set by users
#define SAMPLE_DIRECTORY    "/home/xuzhen/qy/sample/"
#define GLOBAL_PARAMETER_PATH "/home/xuzhen/qy/ReadSampleMonitor/parameters/globalParameters.txt"
/*
 * global parameters used in this project
 */
typedef struct range
{
	unsigned int momLen;    //300；生成母体read数据时候的基本窗口值
	unsigned int fetal;     //比母体窗口值扩大1/(3%)倍；胎儿DNA少，所以窗口值大
	unsigned int femaleY;   //1/(3% *5%)read出现在Y染色体上相对于出现在其他染色体上的概率的比例。
	unsigned int lenFlex;   //窗口值的松弛因子。不能将每个窗口设成等值固定大小
	unsigned int a_mom;     //(momLen*0.5%);母体窗口内read数的基本值
	unsigned int a_child;   //(momLen*0.01%);母体窗口内read数的基本值
	unsigned int aFlex;     //窗口内read数的松弛因子
	unsigned int snp;       //2；出现snp的几率的基本值（约为2/1000）
	unsigned int snpFlex;   //snp常集中出现在某些区域，所以对不同的区域根据松弛因子来生成一个 possibilitySNP，这样，在不同区域得到snp的几率就会出现差异
	unsigned int readLen;   //read length

    char *read;             //a tempt array stores the reads
    char *buff;
    int *pos;

    unsigned int genomeLength;
    unsigned int sample_num; //sample nums which we must calculate from sampleParameters files
    //char fileName[4][CHROMOSOME_NUM+1][100];
    char fileName1[4][100];
}range;

/*
 * parameters that mark different types of samples
 */
typedef struct sampleParameter
{
	 char parentFlag;//parentFlag='M', means it is the sample of maternal
	 char childFlag; //childFlag='G' represent girl, ='B' represent boy
	 char modeFlag;  //modeFlag='N' represent normal , 'T' represent T21
}sParameter;
//sParameter samplePara[SAMPLE_NUM];
/*
 * jump process related parameters
 */
short jump[CYCLIZLIZE_GENOME_LENGTH];

typedef struct jumpPosType
{
	unsigned int cyclizPos;
	unsigned int chrID;
	unsigned int chrInPos;
	unsigned int jumpNum;
}jumpPos;

/*----------------prepare function definition-------------------------
--------------------------------------------------------------------*/
int setParameters(range *prange);//set the global parameters
int setSample(sParameter *sparameter,range *prange);//read the parameter files to struct list
int getGenomeLength(FILE *fp);
//char *strProcess(char *buff);//omit or the blank in a char array
//char *strProcess1(char *buff, int length);//omit the \n in string
int getMomWindowLength(sParameter *parameter, int chrNum,range *prange);//get the window length of parent
int getChildWindowLength(sParameter *sParameter, int chrNum,range *prange);//get the window length of children
int getAforMom(sParameter *sParameter, int chrNum,range *prange);//get a(read number) in a window for mom
int getAforChild(sParameter *parameter, int chrNum,range *prange);//get a(read number) in a window for child
int randValue( int x,  int y);//adjust x according to y, return a value between [x+y,x-y]


/*-----------------monitor sample needed---------------------*/
int createDirectory(range *prange);//create directory to store all the result sample files
int monitorSample(sParameter *parameter, range *prange);
char* getRead(int pos,int length, char *buff,char *read);//get read from buff, which start from pos

int* snpPosition(range *prange);//get the positions where need to process snp
int* transitionCount(int genomeLength);//get the positions where need to process snp
char transversion(char ch);//base transversion in snp process
char transition(char ch);//base transition in snp process

int* generatePos(int a, int w, int readLen);//generate the position of each pos

char* errorProcess(char* read);//generate error bases randomly
char change(char ch);//called in errorProcess()
char getBaseRandom(int row);//called in change()

unsigned int new_rand ();//use linux file to generate random files

/*------------------------------------------------------------*/
void freeMemory(range*prange);
void outUsedTime(int flag);

#endif /* PARAMETER_H_ */
