#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "gopt.h"

/* Since the max size of human chromosome is 300,000,000, this program cannot allocate that much memory and dies,
 * possible bug to consider when running human samples. Also, because of this memory restrictions, the size of peakArray
 * cannot be more than 200000.
 */

void printWiggle(float *, long, char [],FILE *);
float* computeVector(float *,float *,int);
float* populate(float *, long);
float* memoryAllocate(float *,long);
float* computeSmoothing(float *,int);
void findPeaks(float *);

/* Constant declarations */
long VSIZE = 15000000; /* size of the vector to initialize to */
long SSIZE = 1000;     /* size of smoothing vector to initialize to*/
long DSIZE = 15000000; /* size of the duplicate vector to initialize to */
long PSIZE = 100000;    /* size of the structure to initialize to*/
int const N = 4;        /* the parameter to contorl the spread of sigma, tells the program to go +/- 4*sigma */
const char *outfilename, *infilename,*sig,*ex;
int SIGMA;
int EXCLUSION;

/* Declaration of other global variables */
struct peaks {
				float height;
				long value;
				int flag;
		}*pA;

int revcmp(struct peaks *,struct peaks *); /*Function to be used in qsort */
void printGff(struct peaks *, char [],FILE *);
void callPeaks(struct peaks *);


/* FUNCTION DEFINITIONS START HERE */

void callPeaks(struct peaks *pM){

  long k,z;
  for(k=0;k <PSIZE-1;k++){

	  if(pM[k].flag == 0)
	  	   continue;
	  if(pM[k].height == 0)
		  break;

	 for(z=k+1;z< PSIZE-2 ;z++){
		 if(pM[z].height == 0)
			 break;
		 if(abs(pM[k].value - pM[z].value) <= EXCLUSION)
     	 		  {
		 			  pM[z].flag = 0;
     	 		  }

	 }

  }


}

void findPeaks(float *dup){
    long m,k=0;
    pA = calloc(PSIZE,sizeof(struct peaks));
   		 if (pA == 0)
  		 		{
  		 			printf("ERROR: Out of memory\n");

  		 		}

	for(m=1;m<VSIZE-2;m++){
		if(dup[m] == 0)
			continue;
		if(dup[m] > dup[m-1] && dup[m] > dup[m+1]){   /* check for equality too? */
			pA[k].height = dup[m];
			pA[k].value = m;
			pA[k].flag = 1;
			k++;
		}

	}

}

float* computeVector(float *vector1,float *kernel1,int sigma){

	long j;
	float *smoothedVec;
	smoothedVec = memoryAllocate(smoothedVec,VSIZE);

	/* looping throught the whole original vector of size VSIZE */

	for(j=N*sigma;j<VSIZE-N*sigma;j++){
		//if((j-sigma) < 0 || (j+sigma) > VSIZE)  /* checking for the margins */
		//	continue;
		if(vector1[j] >= 1)      /* start from the value, where the value in the original vector is non-zero*/
		{
			int start = j - N*sigma; /* start from the current position - sigma positions */
			int end   = j+ N*sigma;  /* end to the current position + sigma positions */
			int l = 0,k;
			float factor = vector1[j];
			for(k=start;k<=end;k++){

				smoothedVec[k] = smoothedVec[k] +  factor*kernel1[l++];

			}

		}

	}

return smoothedVec;

}

float* populate(float *vec, long pos){
     vec[pos]++;
     return vec;
}

float* memoryAllocate(float *v,long vsize){

	v = calloc(vsize,sizeof(float));          /*initializing the vector*/
	if(v == NULL){
		printf("ERROR: Out of memory\n");
			}
	return v;

}

float* computeSmoothing(float *sm, int sigma){

    int i,j=0;
    for(i= -(N*sigma);i <=N*sigma; i++){
    	sm[j++] += exp(-pow(i,2)/(2*(pow(sigma,2))));

    }
	return sm;

}

int revcmp(struct peaks *v1, struct peaks *v2){
    if(v1->height > v2->height)
    	return -1;
   	else if(v1->height < v2->height)
    	return 1;
   	else
   		return 0;


}

void printWiggle(float *arr, long size, char chr[10],FILE *op){
	long i;char temp[10];
	strcpy(temp,chr);
	fprintf(op,"variableStep chrom=%s\n",temp);
	//printf("variableStep chrom=%s\n",temp);
	for(i=0;i<size;i++)
	{
		if(arr[i] <= 0)
			continue;
		fprintf(op,"%ld %f\n",i,arr[i]);
		//printf("%ld %f\n",i,arr[i]);
	}
}

void printGff(struct peaks *pB, char chr[10],FILE *op){
	long m,start,end;
	char temp[10];
	strcpy(temp,chr);
	for(m=0;m<PSIZE;m++)
	{

		start = pB[m].value - EXCLUSION;
		end   = pB[m].value + EXCLUSION;
		if((pB[m].height <=0) || (pB[m].flag == 0))    /*  change here, if you want to print by peak height */
			continue;

		fprintf(op,"%s\t%s\t%s\t%ld\t%ld\t%f\t%s\t%s\t%s\n",temp,"genetrack",".",start,end,pA[m].height,".",".",".",".");
		//printf("%s\t%s\t%s\t%ld\t%ld\t%f\t%s\t%s\t%s\n",temp,"genetrack",".",start,end,pB[m].height,".",".",".",".");


	}
        fflush(stdout);
        fflush(op);
}


/* Start of main function */

int main (int argc, const char **argv){
                                
      void *options= gopt_sort( & argc, argv, gopt_start(
      gopt_option( 'h', 0, gopt_shorts( 'h', '?' ), gopt_longs( "help", "HELP" )),
      gopt_option( 's',GOPT_ARG, gopt_shorts('s'), gopt_longs( "sigma" )),
      gopt_option( 'e',GOPT_ARG, gopt_shorts('e'), gopt_longs( "exclusion" )),
      gopt_option( 'i', GOPT_ARG, gopt_shorts( 'i' ), gopt_longs( "input" )),
      gopt_option( 'o', GOPT_ARG, gopt_shorts( 'o' ), gopt_longs( "output" ))));

      FILE *fp, *op; /* pointer to the input and output file */
      
      if( gopt( options, 'h' ) ){
          fprintf( stdout, "help text\n" );
          exit( EXIT_SUCCESS );
      }
      
      if( gopt_arg(options, 's', &sig) && strcmp(sig, "-" )){
         SIGMA = atoi(sig);
      }else{SIGMA = 5;}
      
      if( gopt_arg(options, 'e', &ex) && strcmp(ex, "-" )){
         SIGMA = atoi(ex);
      }else{SIGMA = 20;}
      
      if( gopt_arg(options, 'i', & infilename) && strcmp(infilename, "-" )){
    
          fp = fopen(infilename,"r");
          if(fp == NULL){
          fprintf(stderr, "%s: %s: could not open file for input\n", argv[0], infilename);
          exit(EXIT_FAILURE);
          } 
    
      }
      
     if( gopt_arg(options, 'o', & outfilename ) && strcmp( outfilename, "-" )){
    
        op = fopen(outfilename,"w");
        if(op == NULL){
        fprintf(stderr, "%s: %s: could not open file for output\n", argv[0], outfilename );
        exit(EXIT_FAILURE);
        }   
    }
    else{op = stdout;}
    
        fprintf(op,"%s\n","##gff-version 3");                   /* line to print the header for gff.Where to put? */
	char line[500];   /* char array to store each line*/
	char *toks[10];   /* toks is an array of 10 elements, each of which points to a char */
	char CHR[10],STRAND[5];
	long START;
	char PREVIOUS_CHR[10] = "NULL";
	float *vector, *dupvector;   /* vector to be populated */
	float *kernel;       /* smoothing vector to be populated */

	vector = memoryAllocate(vector,VSIZE);         /* call to calloc to allocate and initialize vector */
	kernel = memoryAllocate(kernel,SSIZE);   /* ONE TIME call to calloc to allocate and initialize smoothing vector */
	kernel = computeSmoothing(kernel,SIGMA);  /* Call to populate the smoothing vector */



	while(fgets(line,500,fp) != NULL) /* fgets */
	{

		if(line[0] == '#') /* ignore all the comment lines that start with # */
			continue;

		int cols = 0;
		toks[cols] = strtok(line,"\t");
		while(cols < 7)
					{
						toks[++cols] = strtok(NULL,"\t");                    /* storing all token in toks for one line */

					}
		strcpy(CHR,toks[0]);
		strcpy(STRAND,toks[6]);
		START = atoi(toks[3]);
		//printf("%s,%ld\n",CHR,START);

		if((strcmp(CHR,PREVIOUS_CHR)) && (strcmp("NULL",PREVIOUS_CHR))){

				dupvector = computeVector(vector,kernel,SIGMA);
				//printWiggle(dupvector,VSIZE,PREVIOUS_CHR,op);
				findPeaks(dupvector);
				qsort(pA,PSIZE,sizeof(struct peaks),revcmp);
				callPeaks(pA);
				printGff(pA,PREVIOUS_CHR,op);
				/* making memory free here */
				free(vector);
				vector = memoryAllocate(vector,VSIZE);
				free(dupvector);
				free(pA);

		}

		 vector = populate(vector,START);
		 strcpy(PREVIOUS_CHR,CHR);

	 }

	/* Calling on the compute function here again to make sure the last line is executed too */
	dupvector = computeVector(vector,kernel,SIGMA);
	//printWiggle(dupvector,VSIZE,CHR,op);
	findPeaks(dupvector);
	qsort(pA,PSIZE,sizeof(struct peaks),revcmp);
	callPeaks(pA);
	printGff(pA,CHR,op);
	fclose(fp);
	fclose(op);

}
