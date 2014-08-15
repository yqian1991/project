#include "genePreDefinition.h"

#define FILENAME_LEN			256
#define FILES_NUM				100
#define CHR_GENOME_LENGTH 0x10000000
//process the variables
int dataPreInit(chrPos* chr)
{
	int i=0;

	for (i = 0; i < 32; i++) {
		chr->strFileName[i] = (char*) malloc(256 * sizeof(char));
	}
	chr->fileNameLen = (int*) malloc(FILES_NUM * sizeof(int));//fileNameLen;
	chr->pCyclizGenomes=(unsigned int*) malloc(sizeof(unsigned int) * CYCLIZE_GENOMES_LEN);
	if ( (chr->pCyclizGenomes == NULL) ) {
		printf("malloc failed!\n");
		return 0;
	}
	mkdir(GENE_DIR,0777);
	return 1;
}
/*
 *read all the genomes files names to strFileName[][];
 *@return:1(success) 0(failed)
 */
int setGenomeFileName(chrPos*chr)
{
	FILE *fp;
	char ch;
	char fileName[FILENAME_LEN];
	int row=0,i=0;

	strcpy(fileName,GENOMES_FILE_INFO);
	fp = fopen(fileName, "r");
	if (fp == NULL) {
		printf("In checkDataFile(), cannot open file [ %s ], error!\n",GENOMES_FILE_INFO);
		return 0;
	}

	int genomesFilesNum = 0;
	int curFileNameLen = 0;
	//read genomesfileinfor.txt file and store the file names in pParameter->strFileName[row][curFileNameLen]
	while ((ch = fgetc(fp)) != EOF) {
		if (ch == ' ') {
			continue;
		} else if (ch != '\n') {
			fileName[curFileNameLen] = ch;
			curFileNameLen++;
		} else if (ch == '\n' && curFileNameLen > 0) {
			row = genomesFilesNum;
			for (i = 0; i < curFileNameLen; i++) {
				chr->strFileName[row][i] = fileName[i];
			}
			chr->strFileName[row][curFileNameLen] = '\0';
			chr->fileNameLen[genomesFilesNum] = curFileNameLen;
			curFileNameLen = 0;
			genomesFilesNum++;
		}
	}
	chr->genomesFilesNum=genomesFilesNum;
	fclose(fp);

	struct stat st; //获取文件信息
	int fileID=0;
	strcpy(fileName, SOURCE_GENEFILE_PATH);//SOURCE_GENEFILE_PATH="genefile/"
	int path_len = strlen(fileName);
	for (fileID = 0; fileID < genomesFilesNum; fileID++) {
		row=fileID;
		curFileNameLen=chr->fileNameLen[fileID];
		for (i = 0; i < curFileNameLen; i++) {
			fileName[path_len + i] = chr->strFileName[row][i];
		}
		fileName[path_len + curFileNameLen] = '\0';

		if (stat(fileName, &st) == 0) {
			printf("the file [ %s ] has already exist.\n", fileName);
		} else {
			printf("the file [ %s ] has not exist.\n", fileName);
			return 0;
		}
	}
	return 1;
}
/*
 *
 */
int cyclizGenomes(chrPos *chr)
{
	FILE *fp;
	int fileID;

	char curFileName[FILENAME_LEN];
	char fileName[FILENAME_LEN];
	unsigned int *pCyclizGenomes = chr->pCyclizGenomes;
	int genomesFilesNum = chr->genomesFilesNum;
	int pathLen=strlen(SOURCE_GENEFILE_PATH);

	unsigned int cyclizPos = 0;
	unsigned int curChFlag; /*  0: is n||N, 1: is a||c||g||t||A||C||G||T  */

	unsigned int curChrPos = 0;
	//unsigned int jumpChrPos;
	unsigned int baseValue;
	unsigned int baseTmp = 0;
	unsigned int posRemainder16 = 0;
	unsigned int posQuotient16 = 0;
	//unsigned int hashTMP = 0;
	//unsigned int hashMask = 0xFFFFFF; //主键是12个碱基
	char ch;

	strcpy(curFileName, SOURCE_GENEFILE_PATH);//SOURCE_GENEFILE_PATH="genefile/"
	for (fileID = 0; fileID < genomesFilesNum; fileID++) {
		strcpy(fileName, chr->strFileName[fileID]);//read one file at one time in each for recycle.
		printf("Begin to process file [ %s ], please wait ...\n", fileName);

		strcpy(curFileName + pathLen, fileName);
		fp = fopen64(curFileName, "r");
		if (fp == NULL) {
			printf("In cyclizGenomes(), cannot open file [%s], error!\n",curFileName);
			return 0;
		}

		while ((ch = fgetc(fp)) != '\n')
					; // read first row, namely annotation

		curChrPos = 0;
		chr->pos[fileID]=posQuotient16;
		//printf("compare1:pos[%d]=%d\n",fileID,posQuotient16);
		while ((ch = fgetc(fp)) != EOF) {
			switch(ch){
			     case 'A':
						baseValue = 0;
						curChFlag = 1;
						break;
				 case 'a':
						baseValue = 0;
						curChFlag = 1;
						break;
				case 'C':
						baseValue = 1;
						curChFlag = 1;
						break;
				case 'c':
						baseValue = 1;
						curChFlag = 1;
						break;
				case 'G':
						baseValue = 2;
						curChFlag = 1;
						break;
				case 'g':
						baseValue = 2;
						curChFlag = 1;
						break;
				case 'T':
						baseValue = 3;
						curChFlag = 1;
						break;
				case 't':
						baseValue = 3;
						curChFlag = 1;
						break;
				case 'U':
				case 'u':
						baseValue = 3;
						curChFlag = 1;
						break;
				case 'N':
				case 'n':
				case '.':
				case 'B':
				case 'b':
				case 'D':
				case 'd':
				case 'H':
				case 'h':
				case 'R':
				case 'r':
				case 'V':
				case 'v':
				case 'Y':
				case 'y':
						curChFlag = 0;
						break;
				default:
						//					curChFlag = 0;
						//					break;
						continue;
			}//end switch
			if(curChFlag==1){
				//chr->pos[fileID]=curChrPos;//the start position of each chromosme
				//unsigned integer is 32byte, each base use 2 bytes, so an integer can only hold 16bases
				posRemainder16 = cyclizPos & 15;//cyclizPos increase one at one time posRemainder16 ranges in [0,...,15]
				posQuotient16 = cyclizPos >> 4;//every 16 positions share the same posQuotient16
				//printf("cyclizPos=%d.posRemainder16=%d. posQuotient16=%d\n",cyclizPos,posRemainder16,posQuotient16);
				//printf("baseValue1=%d  ",baseValue);
				baseValue <<= 30 - (posRemainder16 << 1); //左移位。(remainder0,remander1,...,remainder15)在baseTMP中的顺序。
				//printf("baseValue2=%d ",baseValue);
				baseTmp |= baseValue;
				//printf("baseTmp=%d\n",baseTmp);//baseTmp ?

				if (15 == posRemainder16) {
					pCyclizGenomes[posQuotient16] = baseTmp;//
					//printf("compare2:%d\n",posQuotient16);
					baseTmp = 0;
				}
				cyclizPos++;
			}// end if(curChFlag==1)
			curChrPos++;

		}//end while ((ch = fgetc(fp)) != EOF)
		fclose(fp);
	}/* End  for(fileID=0; fileID<genomesFilesNum; fileID++) */

	//printf("cyclizPos=%d\n",cyclizPos);
	chr->pos[0]=0;
	chr->pos[genomesFilesNum]=posQuotient16+1;//bases after pos[24] doesn't really exist.
	/* The cyclization of the genome length divisible by 256 */
	//so if the bases number exceeded 512, we will also give him 512n memory
	//so totalCyclizPosLen = (cyclizPos / 512 + 1) * 512
	//and it also means that many positions will have no bases appear, because the memory is large enough
	unsigned int totalCyclizPosLen;
	if (cyclizPos % 512 > 0) {//in the above cycle, cyclizPos was increased by bases number,at last, cyclizPos equals to the total number of bases in txts
		totalCyclizPosLen = (cyclizPos / 512 + 1) * 512;//
		printf("totalCyclizPosLen=%u \n",totalCyclizPosLen);
		while (cyclizPos < totalCyclizPosLen) {
			baseValue = 3; //最后一块尾部空缺填充3--T。
			//unsigned integer is 32byte, each base use 2 bytes, so an integer can only hold 16bases
			posRemainder16 = cyclizPos & 15;
			posQuotient16 = cyclizPos >> 4;
			baseValue <<= 30 - (posRemainder16 << 1);
			baseTmp |= baseValue;

			if (15 == posRemainder16) {
				pCyclizGenomes[posQuotient16] = baseTmp;
				baseTmp = 0;
			}//end if
			cyclizPos++;
		}//end while(cyclizPos < totalCyclizPosLen)
	}//end if (cyclizPos % 512 > 0)
	else{
		totalCyclizPosLen = cyclizPos;
	}
	chr->cyclizGenomesLen=totalCyclizPosLen;
	printf("ok\n");
	return 1;
}
/*
 * return ATCG from char 0123
 */
char* getCharValSize(int val, char* seq, int size)
{
	int i = 0;
	unsigned int tmp = 0;
	unsigned int numericVal = 3;

	while (i < size) {
		tmp = val & numericVal;
		tmp = tmp >> (2 * i);

		switch (tmp) {
		case 0:
			seq[size - i - 1] = 'A';
			break;
		case 1:
			seq[size - i - 1] = 'C';
			break;
		case 2:
			seq[size - i - 1] = 'G';
			break;
		case 3:
			seq[size - i - 1] = 'T';
			break;
		default:
			//return -1;
			break;
		}

		i++;
		numericVal = numericVal << 2;
	}
	return seq;
}
/*
 *write bases in pCyclizGenomes to a file
 */
int pCyclizGenomesFile(chrPos*chr)
{
	char fileName[FILENAME_LEN];
	unsigned int genomeLen = chr->cyclizGenomesLen;
	int genomeFileNum=chr->genomesFilesNum;
	unsigned char *seq=(unsigned char*) malloc(sizeof(unsigned char)*16);

	//strcpy(fileName, SOURCE_GENEFILE_PATH);//SOURCE_GENEFILE_PATH="genefile/"
	strcpy(fileName,CYCLIZE_FILE_NAME);

	FILE *fp=fopen64(fileName,"w");
	if (fp == NULL) {
			printf("In pCyclizGenomesFile(), cannot open file [%s], error!\n",fileName);
			return 0;
	}

	int i=0,T=0;

	for(i=0;i<genomeLen/16;i++){
		int j=0;
		for(j=0;j<genomeFileNum+1;j++){//judge whether a chr file over, if over , write > to mark
			if( (chr->pos[j]==i)&&(i!=0) ){
				fwrite( "\n", 1, 1, fp );
				fwrite( ">", 1, 1, fp );
				fwrite( "\n", 1, 1, fp );
				T=0;
			}
		}
		//printf("i=%d %d\n",i,chr->pCyclizGenomes[i]);
		seq=getCharValSize(chr->pCyclizGenomes[i],seq,16);
		fwrite(seq, sizeof(char), 16, fp);//each element in pCyclizGenomes[i] has 16 bases
		T++;
		if(T==3){//make sure a line store 48 bases
			fwrite( "\n", 1, 1, fp );
			T=0;
		}
	}//end for(i=0;i<genomeLen/16;i++)
	fwrite( ">", 1, 1, fp );
	free(seq);
	fclose(fp);
	return 1;
}
/*
 * according to the relationship between GC stratum and read coverage
 * @return:readInGcStratum[]
 */
int setReadInGcStratum(int a)
{
	int result=0;
	result = 400*a/3 - 4*a*a/9 - 8000;
	if(a<90||a==90){
		result=(900*a-5*a*a)/81;
	}
	if((a>90)&&((a<150)||(a==150))){
		result=125*a-2*a*a/5-7375;
	}
	if((a>150)&&((a<210)||(a==210))){
		result=83*a-5*a*a/18-4250;
	}
	if((a>210)&&((a<300)||(a==300))){
		result=23*a-a*a/18-1450;
	}
	return result;
}
/*
 *compute the GC content in w[0,l] of each base
 *get the read number of each base
 *@return readsOfBase[]
 */
int getReadsOfBase(chrPos* chr)
{
	char ch;
	unsigned int baseCount=0;//counter of bases
	unsigned int chrBaseCount=0;//counter of bases in each chromosome
	int gcnum=0;
	for(gcnum=0;gcnum<301;gcnum++){//calculate read coverage according to GC
		chr->readInGcStratum[gcnum]=setReadInGcStratum(gcnum);
	}

	char fileName[FILENAME_LEN];
	char filename1[100];
	strcpy(fileName,CYCLIZE_FILE_NAME);
	//printf("genomeLen=%d\n",chr->cyclizGenomesLen);
	char *baseBuff=(char*) malloc( CHR_GENOME_LENGTH * sizeof(char));
	if(baseBuff==NULL){
		printf("malloc memory for baseBuf error\n!!!");
		return 0;
	}
	strcpy(filename1,READBASE_NUM_FILE);
	FILE *fp_1=fopen64(filename1,"w");
	FILE *fp=fopen64(fileName,"r");
	if( (fp==NULL)||(fp_1==NULL) ){
		printf("read file error.\n");
		return 0;
	}

	while( (ch=fgetc(fp))!=EOF ){

		if( (ch==' ')||(ch=='\t')||(ch=='\n') ){
			continue;
		}
		baseCount++;
		//printf("ch=%c,basecount=%d\n",ch,baseCount);
		baseBuff[chrBaseCount++]=ch;

		if(ch=='>'){//implies that a chromosome has stored into baseBuff, then need to process it first
			//process baseBuff[chrCount], get readsOfBase[] here
			printf("process a chr chrCount=%d\n",chrBaseCount);

			if(chrBaseCount>300){
				unsigned int index=0;
				for(index=0;index<chrBaseCount-1-300;index++){
					int gc=0;
					int i=0;
					for(i=1;i<301;i++){
						//printf("baseBuff1=%c\n",baseBuff[index+i]);
						if( (baseBuff[index+i] =='G')||(baseBuff[index+i]=='C') ){
							gc++;
						}
					}
					//printf("gc=%d ",gc);
					//printf("read=%d\n",readsOfBase[baseCount-chrCount+index]);
					fprintf(fp_1,"%d\n",(chr->readInGcStratum[gc]));
				}//end for(index=0;index<chrCount-300;index++)
				int remain[300];//process the last 300 bases of each chromosome
				int r=0;
				for(r=0;r<300;r++){//deal with the last 300 bases
					if( (baseBuff[chrBaseCount-1-300+r]=='G')||(baseBuff[chrBaseCount-1-300+r]=='C') ){
							remain[r]=1;
					}else{
							remain[r]=0;
					}
				}
				//printf("remain=%d\n",remain[0]);
				for(r=0;r<300;r++){//deal with the last 300 bases
					int gc_last=0;
					int remaingc=0;
					for(remaingc=r+1;remaingc<300;remaingc++){
						gc_last+=remain[remaingc];
					}
					fprintf(fp_1,"%d\n",(chr->readInGcStratum[gc_last]));
				}
				//process over, then do as below
				chrBaseCount=0;
				//printf("ok\n");
				baseCount--;
			}else{
				int remain[chrBaseCount];//process the last 300 bases of each chromosome
				int r=0;
				for(r=0;r<chrBaseCount;r++){//deal with the last 300 bases
					if( (baseBuff[r]=='G')||(baseBuff[r]=='C') ){
						remain[r]=1;
				    }else{
						remain[r]=0;
					}
				}
				//printf("remain=%d\n",remain[0]);
				for(r=0;r<chrBaseCount;r++){//deal with the last 300 bases
					int gc_last=0;
					int remaingc=0;
					for(remaingc=r+1;remaingc<chrBaseCount;remaingc++){
						gc_last+=remain[remaingc];
					}
					fprintf(fp_1,"%d\n",(chr->readInGcStratum[gc_last]));
				}
			}//end if(chrBaseCount<300) else
		}//end if(ch=='>')

	}//end while( (ch=fgetc(fp))!=EOF )
	printf("baseCount=%u\n",baseCount);
	free(baseBuff);
	fclose(fp);
	fclose(fp_1);
	return 1;
}
/*
 *free the memory allocated
 */
void memFree(chrPos* chr)
{
	free(chr->pCyclizGenomes);
	//free(chr->readsOfBase);
	free(chr->fileNameLen);
	int i;
	for (i = 0; i < 32; i++) {
		free(chr->strFileName[i]);
	}
}
