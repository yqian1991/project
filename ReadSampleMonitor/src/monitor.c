/*
 * monitor.c
 *
 *  Created on: 2013-4-23
 *      Author: qianyu
 */

#include "parameter.h"
/*
 * generate reads and write to files
 */
int monitorSample(sParameter *parameter, range *prange)
{
	FILE *fp;//the file store our reads
	FILE *fp1;//the file store genomes
	FILE *fp_rb;//thje file store read coverage information of each base
	int i=0;
	int lineCount=0;//linecount is the line number in readOfBase.bin,
	char file[256];
	int lines[(prange->momLen+prange->lenFlex)*100];
	unsigned int genomeLen=prange->genomeLength;
	printf("genomelength=%d\n",genomeLen);

	//start to generate data for each sample
	for(i=0;i<prange->sample_num;i++){

		//for(j=1;j<=CHROMOSOME_NUM;j++){
			//prepare files needed
			fp1=fopen64(CYCLIZE_FILE_NAME,"r");//open the cycliezd genomes file
			if(fp1 == NULL){
				printf("cannot open file %s!!!\n",CYCLIZE_FILE_NAME);
				return 0;
			}
			fp_rb=fopen64(READBASE_NUM_FILE,"r");
			if(fp_rb == NULL){
				printf("cannot open file %s!!!\n",READBASE_NUM_FILE);
				return 0;
			}
			strcpy(file,prange->fileName1[i]);
			strcpy(file+strlen(file),"momchild.fastq");
			printf("now creating file:%s\n",file);
			fp= fopen64(file,"w");
			if(fp==NULL){
				printf("open file %s error\n",file);
			}

			printf("now generating monitor reads(mom and child)...\n");
			//------generate read data and write to file---------------
			/*
			 * input:genome file
			 * output:
			 */
			int chr_num=1;
			int windowsLength=getMomWindowLength( (parameter+i),chr_num,prange);//printf("windowslength 1=%d\n",windowsLength);
			int a_mom=getAforMom(parameter+i, chr_num,prange);
			int a_child=getAforChild(parameter+i, chr_num,prange);

			unsigned int index=0;//index used to track the number of bases in genome files
			char *buff = prange->buff;//buff store a window

			int *snpPos;
			snpPos = snpPosition(prange);
			printf("snpPos finished\n");
			int *count=transitionCount(genomeLen);
			int readCount=0;//counter of reads

			int jg_flag=0;
			char ch;//read a ch (a base) from file
			//printf("index=%d\n",index);
			while(index<genomeLen){//read all the bases in a genome

				int windowsLength_adjust=randValue(windowsLength,prange->lenFlex);//adjust length
				int a_m=randValue(a_mom,prange->aFlex);//a_mom is the number of reads in single window
				int a_c=randValue(a_child,prange->aFlex);//a_child is the number of reads in single window
				//printf("adjust windowslength=%d\n",windowsLength_adjust);

				int w=0,rw=0;//w count the number of bases in a window, when it >windowslength, we must stop.
				jg_flag=0;
				while( (rw<windowsLength_adjust)&&((ch = fgetc(fp1)) != EOF)&&(jg_flag==0) ){

					if((ch==' ')||(ch=='\t')){//skip
						continue;
					}
					if(ch=='>'){//the end markd of a chr
						printf(" chr%d finished, now index=%u\n",chr_num,index);
						fwrite( "\n>\n", 1, 3, fp );
						jg_flag=1;//must stop the current process. and start to process next chr
						chr_num++;//means that a chr read over
						continue;
					}
					if( (ch!='\n')){
						index++;;//index record the bases number
						prange->buff[w]=ch;//buff[w]=ch;
						//snp test and process
						//check the index first
						int pos=(index/SNP_POSSIBILITY);
						if(index==snpPos[pos]){//need to process snp
						//n=index/1000, in n, n*2/3 times should transversion,
							int i=0, flag=0;
							for(i=0;i<sizeof(count)/sizeof(int);i++){
								//printf("snp %d %d\n",index/SNP_POSSIBILITY,count[i]);
								if((index/SNP_POSSIBILITY)==count[i]){
									flag=1;
									break;
								}
							}
							if(flag == 1){
								prange->buff[w] = transition(prange->buff[w]);

							}else{
								prange->buff[w] = transversion(prange->buff[w]);
							}
						}
						w++;
						rw=w;
					}//end if
				}//end while( read a window length)
				//printf("process window in buff, w=%d\n",w);
				jg_flag=0;

				if(w>1024){
					int a=a_m+a_c;//a is the last read numbers in a sample genome
					//printf("amom=%d,achild=%d\n",a_mom,a_child);

					lineCount=0;//count the lines read in readOfBase.bin
					while(lineCount<w){//w is the real windows length
						fscanf(fp_rb,"%d",&lines[lineCount++]);
					}
					//printf("linecount=%d\n",lineCount);

					int *pos=prange->pos;//pos[i] save the start position of read[i]
					//printf("read number a=%d \n",a);
					int num=0;//the index of pos[]
					pos=generatePos(a,w,prange->readLen);
					//printf("here\n");
					//filter the pos
					for(num=0;num<a;num++){
						//int readCoverage = lines[pos[num]];
						if( (new_rand()%2375)<(lines[pos[num]]) ){
							//get read and process to write each read to file
							char *read=prange->read;
							//printf("pos[%d]=%d\n",num,pos[num]);
							read=getRead(pos[num],prange->readLen+1,buff,prange->read);//get the read according to each pos
							//printf("read[%d]=%s ",num,read);
							//conduct error process
							if(new_rand()%10 > 1){
								read=errorProcess(read);//error bases generated
							}
							//printf("read[%d]=%s\n",num,read);
							//write read information to file
							fprintf(fp,"@%d\n",readCount++);
							fwrite(read, sizeof(char), prange->readLen, fp);
							fwrite( "\n+\n", 1, 3, fp );
							fwrite(">>>ATCGquality\n\n",1,16,fp);
						}//end if{rand()%max<readCoverage}
					}//end for(num=0;num<a;num++)
				}//end if(w>prange->readLen)
				//else{//w<readLen, don't process, write to file directly
					//fwrite(buff, sizeof(char), w, fp);
				//}
				//printf("\n---------------a window over---------------\n");
				int b=(prange->momLen)<<1;
				memset(prange->buff,0,b*sizeof(char));//printf("after memeset:%s\n",buff);
			}
			fclose(fp);
			fclose(fp1);
			fclose(fp_rb);
		//}//end for(CHROMOSOME_NUM)
	}//end for(SAMPLE PARAMETER)
	return 1;
}
/*
 *get each read from a specific window(buff), the read start from pos
 *@return:read
 */
char* getRead(int pos, int length, char *buff,char *read)
{
	int i=0,j=0;
	for(i=pos;i<pos+length;i++)
	{
		read[j]=buff[i];
		j++;
	}
	read[j]='\0';
	return read;
}
/*
 *create all the directory needed according to sample numbers and chromosome numbers
 */
int createDirectory(range *prange)
{
	int i=0;
	char file[256],str[3];
	mkdir(SAMPLE_DIRECTORY,0777);//create the folder of sample
	for(i=0;i<prange->sample_num;i++){
		sprintf(str,"%d",i);//str='i'
		//printf("in create directory:%s\n",str);
		//create a directory for sample i
		strcpy(file,SAMPLE_DIRECTORY);//fileName=sample/
		strcpy(file+strlen(file),str);//fileName=sample/0
		strcpy(file+strlen(file),"/");//fileName=sample/0/
		mkdir(file,0777);//create the folder of sample
		//printf("create directory:%s\n",file);

		strcpy(prange->fileName1[i],file);
	}
	return 1;
}
/*
 * random algorithms
 */
unsigned int new_rand ()
{
  int fd;
  unsigned int n = 0;
  fd = open ("/dev/urandom", O_RDONLY);
  if (fd > 0)
  {
    read (fd, &n, sizeof (n));
  }
  close (fd);
  return n;
}
/*
 *generate the positions where we should process snp, store the positions into an array.
 */
int* snpPosition(range *prange)
{
	int n = (prange->genomeLength/SNP_POSSIBILITY);//we need to change n positions
	static int snpPos[CYCLIZLIZE_GENOME_LENGTH/SNP_POSSIBILITY];
	int i=0;
	snpPos[0]=new_rand()%10;
	for(i=1;i<n;i++){
		snpPos[i]=randValue(snpPos[i-1]+1000,prange->snpFlex);
		//printf("index = %d\n",snpPos[i]);
	}
	return snpPos;
}
/*
 *select positions to generate transition change in random
 */
int* transitionCount(int genomeLength)
{
	static int count[CYCLIZLIZE_GENOME_LENGTH/1000];
	int n = genomeLength/SNP_POSSIBILITY;
	int number = n/3;//the possibilty of transition in snp is 1/3
	int i=0;
	for(i=0;i<number;i++){
		int index = new_rand()%n;
		count[i]=index;
	}
	return count;
}
/*
 *generate the positions which are possible to generate reads
 *sort all the positions
 */
int* generatePos(int a, int w, int readLen)
{
	static int pos[100];
	int num=0;//the index of pos[]
	pos[0]=new_rand()%5;
	for(num=1;num<a;num++){
		pos[num]=randValue(pos[num-1]+190,20);
		if( ((pos[num])+readLen)>w ){
			pos[num] = pos[num]+readLen - 800;//make sure will not exceed window length
		}
	}
	return pos;
}
/*
 *碱基颠换(transversion): 嘌呤与嘧啶之间的替代.
 *A<->T ; C<->G
 */
char transversion(char ch)
{
	//printf("transversion\n");
	switch(ch)
	{
		case  'A':ch='T';break;
		case  'T':ch='A';break;
		case  'C':ch='G';break;
		case  'G':ch='C';break;
	}
	return ch;
}
/*
 * 碱基转换(transition):嘌呤被嘌呤,或者嘧啶被嘧啶替代。
 * A<->G ; C<->T
 */
char transition(char ch)
{
	//printf("transition\n");
	switch(ch)
	{
			case  'A':ch='G';break;
			case  'T':ch='C';break;
			case  'C':ch='T';break;
			case  'G':ch='A';break;
	}
	return ch;
}
/*
 *对每个位置上的碱基随机产生错误碱基，每个位置都有可能错误成错读成任意其他碱基。
 */
char* errorProcess(char* read)
{
	int length = strlen(read);
	int pos=new_rand()%length;
	read[pos]=change(read[pos]);
	return read;
}
/*
 * used in function: errorProcess()
 */
char change(char ch)
{
	switch(ch)
	{
		case 'A':return getBaseRandom( 3 );
		case 'T':return getBaseRandom( 2 );
		case 'C':return getBaseRandom( 1 );
		case 'G':return getBaseRandom( 0 );
	}
	return ch;
}
/*
 * used in function:change()
 */
char getBaseRandom(int row)
{
	char base[4][3]={ {'A','T','C'},
	                  {'A','T','G'},
	                  {'A','C','G'},
	                  {'T','C','G'}
	                };

    int col= new_rand()%3;
    char ch=base[row][col];
    return ch;
}

/*void read(sParameter *parameter, range*prange)
{
	printf("read sampleParameter :\n");
	int i=0,j=0;
	for(i=0;i<prange->sample_num;i++)
	{
		for(j=1;j<=CHROMOSOME_NUM;j++)
		{
			printf("%c%c%c ",(parameter+i)->parentFlag,(parameter+i)->childFlag,(parameter+i)->modeFlag);
		}

	}
	printf("\n");
}*/
/*
void init_random ()
{
  unsigned int ticks;
  struct timeval tv;
  int fd;
  gettimeofday (&tv, NULL);
  ticks = tv.tv_sec + tv.tv_usec;
  fd = open ("/dev/urandom", O_RDONLY);
  if (fd > 0)
   {
      unsigned int r;
      int i;
      for (i = 0; i < 512; i++)
       {
          read (fd, &r, sizeof (r));
          ticks += r;
        }
      close (fd);
    }

  srand (ticks);
  printf("init finished ");
} */
/*
 *not all the positions in pos[n] is valid
 *use this function to filter some positions
 */
/*int* filterPos(int *pos,int *lines, int a, int* j)
{
	int i=0;
	int n=0;
	for(i=0;i<a;i++){
		int readCoverage=lines[pos[i]];
		int rand=new_rand()%30;//30 is the maximum in readCoverage, it may changes
		//printf("rand=%d,readC= %d\n",rand,readCoverage);
		if ( rand < readCoverage){
			pos[n]=pos[i];
			n++;
		}
	}
	*j=n;
	//printf("j=%d\n",*j);
	return pos;
}*/
/*char* snpProcess(char *read, int possibilitySNP,range*prange)
{
	int length=strlen(read);
	possibilitySNP = randValue(2, 0)%1000;

	if((rand()%1000)<possibilitySNP){//A<->T ; C<->G//select two poses to substitite；
		read=transversion(read,length);
		//printf("1 pos changed:%s\n",read);
		read=transversion(read,length);
		printf("transversion:%s\n",read);
		        //substitute规则：  选第几个碱基？
	} else{ //select one pos to substitite；
		read=transition(read,length);
		printf("transition:%s\n",read);
	}
	return read;
}*/
/*int createDirectory(range *prange)
{
	int flag=0;
	int i=0,j=0;
	char file[256],file1[256],str[3],str1[3];
	for(i=0;i<prange->sample_num;i++)
	{
		sprintf(str,"%d",i);//str='i'
		//printf("in create directory:%s\n",str);
		//create a directory for sample i
		strcpy(file,SAMPLE_DIRECTORY);//fileName=sample/
		strcpy(file+strlen(file),str);//fileName=sample/0
		strcpy(file+strlen(file),"/");//fileName=sample/0/
		mkdir(file,0777);//create the folder of sample
		//printf("create directory:%s\n",file);
		for(j=1;j<=CHROMOSOME_NUM;j++)
		{
			strcpy(file1,file);
			sprintf(str1,"%d",j);
			strcpy(file1+strlen(file1),str1);//fileName=sample/0/1
			strcpy(file1+strlen(file1),"/");//fileName=sample/0/1/
			mkdir(file1,0777);//create the folder of each chromosome in each sample

			//prange->fileName[i][j]=fileName;
			strcpy(prange->fileName[i][j],file1);
			//printf("in createDirectory() filename[%d][%d]=%s\n",i,j,prange->fileName[i][j]);
			flag=1;
		}
	}
	return flag;
}*/
/*int* generatePos(int a, int w, int readLen)
{
	static int pos[100];
	int num=0;//the index of pos[]
	//pos=generatePos(a,w);
	for(num=0;num<a;num++){
		pos[num]=randPos(w);//w is the real windows length, mostly ,it equals windowsLength_adjust
		while( ((pos[num])+readLen)>w ){
		 pos[num]=randPos(w);//需要改进,不随机,重复度太高
	    }
		//verify whether the pos[num] is valid here use jumpfile to process
	}
	bubbleSortA(pos,a);
	return pos;
}*/
/*
 * generate the start position of each read on a window
 */
/*int randPos(unsigned int length)
{
	int num=0;
	srand( (unsigned int)time(0) );
	num = new_rand()%(length);

	return num;
}*/
/*
 *bubble sort algorithms
 */
/*void bubbleSortA(int array[], int n)
{
    int i=0,j=0;
    for(i=0;i<n;i++){
        for(j=0;j<n-i-1;j++){
            if(array[j]>array[j+1]) {
                swap(array,j,j+1);
            }
        }
    }

}*/
/*
 * swap two integers in array
 */
/*void swap(int array[],int i, int j)
{
    int temp=0;
    temp = array[i];
    array[i] = array[j];
    array[j] = temp;
}*/
