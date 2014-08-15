/*
 * dataPreProcess.c
 * process data first, include initialize the parameters, files, directories.
 *  Created on: 2013-4-22
 *      Author: qianyu
 */


#include "parameter.h"

/*
 * read parameters from gloableParameters.txt
 * such as:momLen, readLen...
 */
int setParameters(range *prange)
{
	//initialize parameters here
	FILE *fp_p;
	int x;
	char str[100];
	char filename[100];
	strcpy(filename,GLOBAL_PARAMETER_PATH);
	fp_p=fopen(filename,"r");
	if (fp_p == NULL) {
		printf("In setSample(), cannot open file %s, error!\n",filename);
		return 0;
	}

	fscanf(fp_p,"%s",str);fscanf(fp_p,"%d",&x);
	prange->momLen=x;

	fscanf(fp_p,"%s",str);fscanf(fp_p,"%d",&x);
	prange->lenFlex=x;

	fscanf(fp_p,"%s",str);fscanf(fp_p,"%d",&x);
	prange->a_mom=x;

	fscanf(fp_p,"%s",str);fscanf(fp_p,"%d",&x);
	prange->a_child=x;

	fscanf(fp_p,"%s",str);fscanf(fp_p,"%d",&x);
	prange->aFlex=x;

	fscanf(fp_p,"%s",str);fscanf(fp_p,"%d",&x);
	prange->snp=x;

	fscanf(fp_p,"%s",str);fscanf(fp_p,"%d",&x);
	prange->snpFlex=x;

	fscanf(fp_p,"%s",str);fscanf(fp_p,"%d",&x);
	prange->readLen=x;

	prange->fetal=(prange->momLen)*(1+100/3);
	prange->femaleY=(prange->momLen)*(1+10000/15);

    fclose(fp_p);

	strcpy(filename, CYCLIZE_FILE_NAME);
	FILE *fp= fopen(filename,"r");//get the length of the genome
	if(fp==NULL){
		printf("open ../genefile/chr.fa error.\n");
		return 0;
	}
	prange->genomeLength=getGenomeLength(fp);
	fclose(fp);

	//malloc memory
	prange->read = (char*) malloc( (prange->readLen +1) * sizeof(char));
	prange->buff = (char*) malloc( ((prange->momLen)<<2) * sizeof(char));
	int posnum=(prange->a_mom)+(prange->a_child)+(prange->aFlex) ;
	prange->pos = (int*) malloc( posnum*sizeof(int));
	if( (prange->read==NULL)&&(prange->buff==NULL)&&(prange->pos==NULL) )
	{
		printf("malloc failed");
		return 0;
	}
	return 1;
}
/*
 * read sample parameters from files, this sParameter must be global to be used all
 * it means that we must allocate memory for it first
 */
int setSample(sParameter *sparameter, range *prange)
{
	FILE *fp;
	char filename[100];
	strcpy(filename,SAMPLE_PARAMETER_PATH);
	fp=fopen(filename,"r");
	if (fp == NULL) {
			printf("In setSample(), cannot open file [ %s ], error!\n",SAMPLE_PARAMETER_PATH);
			return 0;
	}
	char ch;
	while ((ch = fgetc(fp)) != '\n')
				; // read first row, namely row names not parameters

	int paraNum=0,sampleNum=0;
	char p[PARA_NUM_COLS];
	while ((ch = fgetc(fp)) != EOF){
			if ( (ch == ' ')||(ch=='\t')){
				continue;
			} else if (ch != '\n'){
				p[paraNum] = ch;
				paraNum++;
			} else if ( (ch == '\n') && (paraNum>0) ){
				(sparameter+sampleNum)->parentFlag=p[0];
				(sparameter+sampleNum)->childFlag=p[1];
				(sparameter+sampleNum)->modeFlag=p[2];
				printf("%c%c%c  ",(sparameter+sampleNum)->parentFlag,(sparameter+sampleNum)->childFlag,(sparameter+sampleNum)->modeFlag);
				paraNum = 0;
				sampleNum++;
			}
	}
	printf("\n");
	prange->sample_num=sampleNum;
	printf("sample num=%d\n",prange->sample_num);
	fclose(fp);
	return 1;
}

/*
 *get the basic  windows length of mother,mother is always normal
 *get just need to consider different chromosomes
 */
int getMomWindowLength(sParameter *parameter, int chrNum,range *prange)
{
	int length=0;
	//如果为母体样本
	//printf("parentflag=%c\n",parameter->parentFlag);
	if(parameter->parentFlag == 'M')
	{
		    //printf("M:get mother window length here\n");
	        if(chrNum <23)
	        	length = prange->momLen;
	        if(chrNum == 23 || chrNum == 24)
	        	length = 2*prange->momLen;
	}
	//printf("mother window length=%d\n",length);
	return length;
}
/*
 *get the basic windows length of child sample,sample can be girl or boy, normal or T21
 */
int getChildWindowLength(sParameter *sParameter, int chrNum,range *prange)
{
	int length=0;
	//如果为男胎儿样本
	if(sParameter->childFlag == 'B'){
	      if(sParameter->modeFlag == 'N'){
	             if(chrNum <23)
	            	 length = prange->momLen*prange->fetal;
	             if(chrNum == 23 || chrNum == 24)
	                 length =  2*prange->momLen*prange->fetal;
	      }
	      if(sParameter->modeFlag == 'T'){
	            if(chrNum<23 && chrNum!=21)
	                 length =  prange->momLen*prange->fetal;
	            if(chrNum == 21)
	            	 length = (prange->momLen*prange->fetal*2/3);
	            if(chrNum == 23 || chrNum == 24)
	                 length = 2*prange->momLen*prange->fetal;
	      }
	}
	//如果为女胎儿样本
	if(sParameter->childFlag == 'G'){
	     if(sParameter->modeFlag == 'N'){
	         if(chrNum<24)
	            length = prange->momLen*prange->fetal;
	         if(chrNum==24)
	            length = prange->momLen*prange->fetal*prange->femaleY;
	     }
	     if(sParameter->modeFlag == 'T'){
	         if(chrNum<24 && chrNum!=21)
	            length = prange->momLen*prange->fetal;
	         if(chrNum == 21)
	            length = (prange->momLen*prange->fetal*2/3);
	         if(chrNum == 24)
	            length =  prange->momLen*prange->fetal*prange->femaleY;
	     }
	}
	return length;
}
/*
 * read number for children, according to different chr
 *@return:a (read number)
 */
int getAforChild(sParameter *sParameter, int chrNum,range *prange)
{
		int a=0;
		//如果为男胎儿样本
		if(sParameter->childFlag == 'B'){
		      if(sParameter->modeFlag == 'N'){
		             if(chrNum <23)
		            	 a = 2*prange->a_child;
		             if(chrNum == 23 || chrNum == 24)
		                 a =  prange->a_child;
		      }
		      if(sParameter->modeFlag == 'T'){
		            if(chrNum<23 && chrNum!=21)
		                 a =  2*prange->a_child;
		            if(chrNum == 21)
		            	 a =  3*prange->a_child;
		            if(chrNum == 23 || chrNum == 24)
		                 a = prange->a_child;
		      }
		}
		//如果为女胎儿样本
		if(sParameter->childFlag == 'G'){
		     if(sParameter->modeFlag == 'N'){
		         if(chrNum<24)
		            a = 2*prange->a_child;
		         if(chrNum==24)
		            a = 0;
		     }
		     if(sParameter->modeFlag == 'T'){
		         if(chrNum<24 && chrNum!=21)
		            a = 2*prange->a_child;
		         if(chrNum == 21)
		            a = 3*prange->a_child;
		         if(chrNum == 24)
		            a = 0;
		     }
		}
		return a;
}
/*
 * read number for mother, according to different chr
 *@return:a read number
 */
int getAforMom(sParameter *parameter, int chrNum,range *prange)
{
	int a;
	//如果为母体样本
	//printf("parentflag=%c\n",parameter->parentFlag);
	if(parameter->parentFlag == 'M'){
		//printf("M:get mother window length here\n");
	    if(chrNum <24)
		      a = prange->a_mom;
		if(chrNum == 24)
		      a = new_rand()%3;
	}
	//printf("mother a =%d\n",a);
	return a;
}
/*
 *return a number between x-y and x+y
 */
int randValue(int x, int y)
{
	if(y==0){
		return x;
	}else{
		int flex=0;
		x-=y;
		y+=y;
		flex = new_rand()%(y);
		//printf("flex=%d \n",flex);
		x+=flex;
		//printf("in randValue():%d\n",x);
		return x;
	}
}
/*
 *get the length of the original reference genome file
 */
int getGenomeLength(FILE *fp)
{
	unsigned int length=0;//mark the genome length
	char ch;
	while ((ch = fgetc(fp)) != EOF){
		if( (ch != '\n')&&(ch!='\t')&&(ch!=' ')&&(ch!='>') )
		{//skip all the invalid chars except A T C G
			length++;
		}
	}
	return length;
}
/*
 * free the memory allocated in setParameters()
 */
void freeMemory(range*prange)
{
	free(prange->read);
	free(prange->buff);
	free(prange->pos);
}
/*
 * calculate the time consuming of this program
 */
void outUsedTime(int flag)
{
	static struct timeval tp_start, tp_end;//struct defined in time.h, used to save the accurate system time
	double time_used;
	static unsigned long start_cpu, end_cpu; // CPU time

	if (flag == 0) {
		gettimeofday(&tp_start, NULL);
		//设置开始时间: CPU time
		start_cpu = clock();
	} else {
		//设置结束时间: CPU time
		end_cpu = clock();
		//print the system time
		gettimeofday(&tp_end, NULL);
		//print the CPU time
		printf("Total CPU time used: %.2f seconds.\n", (double) (end_cpu
				- start_cpu) / CLOCKS_PER_SEC);
		time_used = (1000000 * (tp_end.tv_sec - tp_start.tv_sec)
				+ (double) (tp_end.tv_usec - tp_start.tv_usec)) / 1000000;
		printf("Total Used Time By gettimeofday(): %.2f Seconds.\n", time_used);
		gettimeofday(&tp_start, NULL);
		//设置开始时间: CPU time
		start_cpu = clock();
	}
	return;
}
/*
 * process a char array, omit the TAB

char *strProcess(char *buff)
{
	char *p=buff;
	while(*p!=0)    //没有结束，循环
	{
		if((*p=='\t')||(*p=='\n') )  //遇到TAB处理
		{
			char *q=p;  //从TAB处开始
			while(*q!=0) //直到末尾的所有字符
			{
				*q=*(q+1); //逐次前移
				q++;      //每移一个字符，指针加1，准备移下一个字符
			}
		}
		else        //当前字符不是TAB
		{
			p++;  //指针后移，指向待检查的新字符
		}
	}
	return buff;
}*/
/*
 *

char *strProcess1(char *buff,int length)
{
	int i=0;
	for(i=0;i<length;i++)
	{
		printf("strprocess1 %c\n",buff[i]);
		if(buff[i]=='\n')
		{
			buff[i]='X';
		}
	}
	printf("%s\n",buff);
	return buff;
}*/
/*
 *read readCoverage of each bases from readOfBase.bin
 */
/*int* readLines(char*fname, int start, int finish)
{
	static int lines[1000];
	int temp;

    FILE *fp=fopen(fname,"r");
    if(fp==NULL){
    	printf("open read coverage file error");
    	return 0;
    }
    //printf("ok\n");
    int i=0;
    for(i=0;i<start-1;i++){//skip the lines before line start
    	fscanf(fp,"%d",&temp);
    	//printf("%d\n",temp);
    }
    for(i=start;i<finish+1;i++){//read the needed lines from file
    	fscanf(fp,"%d",&lines[i-start]);
    	//printf("lines[%d]=%d\n",i-start,lines[i-start]);
    }
    fclose(fp);
	return lines;
}*/
