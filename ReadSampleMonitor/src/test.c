/*
 ============================================================================
 Name        : test.c
 Author      : qianyu
 Version     : 1.1
 Copyright   : qianyu
 Description :
 ============================================================================
 */

#include "genePreDefinition.h"
#include "parameter.h"

int processGenome()
{
	chrPos *chr = (chrPos*) malloc(sizeof(chrPos));


	if( dataPreInit(chr)!=1 ){
		printf("initialization failed!!!");
		return 0;
	}

	if( setGenomeFileName(chr)!=1 ){
		printf("function setGenomeFileName() failed!!! go on\n");
		return 0;
	}

	if( cyclizGenomes(chr)!=1 ){
		printf("cyclize genomes files failed!!! go on ");
		return 0;
	}

	if( pCyclizGenomesFile(chr)!=1 ){
		printf("get total one genomes files failed!!! go on ");
		return 0;
	}

	//-----------------------------------
	printf("begin getReadsOfBase, please wait...\n");
	if( getReadsOfBase(chr)!=1){
		printf("get reads of base error");
		return 0;
	}
	memFree(chr);
	free(chr);
	printf("readOfBase generated successfully!!!\n");
	return 1;
}
int generateMonitorRead()
{
		range *prange = (range*) malloc(sizeof(range));
		sParameter *parameter=(sParameter*) malloc(sizeof(sParameter));
		if ((prange == NULL))
		{
			printf("In main(), cannot malloc pParameter. Error!\n");
			return 0;
		} else{

			setParameters(prange);

			if(setSample(parameter,prange)!=1)
			{
				printf("set sampleParameter error!!!\n");
				return 0;
			}

			if(createDirectory(prange)==0)//this function must rely on setSample()
			{
				printf("create directory error!!!\n");
				return 0;
			}

			if(monitorSample(parameter,prange)!=1 )
			{
				printf("monitor sample error!!!\n");
				return 0;
			}
		}

		freeMemory(prange);
		free(prange);
		return 1;
}
int main(void)
{
	outUsedTime(0);

	processGenome();

	generateMonitorRead();

	outUsedTime(1);

	return 0;
}
