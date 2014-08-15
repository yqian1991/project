/*
 * genePreDefinition.h
 *
 *  Created on: 2013-5-22
 *      Author: qianyu
 */

#ifndef GENEPREDEFINITION_H_
#define GENEPREDEFINITION_H_

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>	//是unix/linux系统定义文件状态所在的伪标准头文件。
#include <string.h>

/*--variable definition--*/
#define CYCLIZE_GENOMES_LEN 0x10000000    //268435456(*16>32亿) enough to hold pCyclizGenome,
/*#define SOURCE_GENEFILE_PATH "../../baihy/workspace/genomes/"  //reference genome stored here, 24 chrs
#define CYCLIZE_FILE_NAME "../genefile/chr.fa"//store the base file after cyclization
#define READBASE_NUM_FILE "../genefile/readOfBase.bin"
#define GENOMES_FILE_INFO "filenames11.txt"
#define GENE_DIR "../genefile/"*/
#define GENE_DIR "/home/xuzhen/qy/genefile/"
#define SOURCE_GENEFILE_PATH "/home/xuzhen/baihy/workspace/genomes/"  //reference genome stored here, 24 chrs
#define CYCLIZE_FILE_NAME "/home/xuzhen/qy/genefile/chr.fa"//store the base file after cyclization
#define READBASE_NUM_FILE "/home/xuzhen/qy/genefile/readOfBase.bin"
#define GENOMES_FILE_INFO "/home/xuzhen/qy/ReadSampleMonitor/filenames.txt"

typedef struct chrPosFile
{
	int genomesFilesNum;
	int cyclizGenomesLen;//length of genomes after cyclization
    int readInGcStratum[301];//store the reads number in each GC stratum
	char *strFileName[32];  //数组，存所有基因组文件的文件名。
	unsigned int pos[25];//mark the position of each chromosome after cyclize
	unsigned int *pCyclizGenomes;      //去除N之后，每个碱基两位，保存成一系列的unsigned int 数。
	//unsigned int *readsOfBase;//store the reads number of each base
	int *fileNameLen;    //基因组每个文件的文件名字符串长度

}chrPos;

/*--function definition--*/

int dataPreInit(chrPos* chr);
void memFree(chrPos* chr);

int setGenomeFileName(chrPos *chr);//set and check the genomes files
int cyclizGenomes(chrPos *chr);//connect the 24 chrs to a cyclizgenomes
char* getCharValSize(int val, char* seq, int size);//return ATCG from 0123
int pCyclizGenomesFile(chrPos*chr);//generate ATCG files from pCyclizGenomes which hold bases by 0123

int getReadsOfBase(chrPos* chr);//get the reads number of each base
int setReadInGcStratum(int a);//set the read numbers in each GC stratum,called in getReadsOfBase()


#endif /* GENEPREDEFINITION_H_ */
