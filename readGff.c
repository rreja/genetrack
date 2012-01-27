
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


void printArray(float *, long);
float* computeVector(float *,float *,int);
//float* computeVector(float *,float *,float *,int);
float* populate(float *, long);
float* memoryAllocate(float *,long);
float* computeSmoothing(float *,int);

long VSIZE = 50000000; /* size of the vector to initialize to */
long SSIZE = 1000;     /* size of smoothing vector to initialize to*/
long DSIZE = 50000000; /* size of the duplicate vector to initialize to */



void main ()
{
	FILE *fp; /* pointer to the file */
	fp = fopen("/Users/rohitreja/Desktop/test2.txt","r");   /* reading the file and assigning it to the pointer*/

	char line[500];   /* char array to store each line*/
	char *toks[10];   /* toks is an array of 10 elements, each of which points to a char */
	char CHR[10],STRAND[5];
	long START;
	char PREVIOUS_CHR[10] = "NULL";
	static float *vector, *dupvector;   /* vector to be populated */
	static float *kernel;       /* smoothing vector to be populated */
	int SIGMA = 5;         /* SIGMA value to be taken by user */



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


		if(strcmp(CHR,PREVIOUS_CHR)){
            //dupvector = memoryAllocate(dupvector,VSIZE);
			 dupvector = computeVector(vector,kernel,SIGMA);
            //dupvector = computeVector(vector,kernel,dupvector,SIGMA);
			printArray(dupvector,VSIZE);


			/*
			 * Call the function to compute on vector and in the compute function
			 *
			 * Compute on the vector
			 */
			free(vector);                                      /* making memory free here */
			free(dupvector);
			vector = memoryAllocate(vector,VSIZE);


		}

		vector = populate(vector,START);
		//printf("%ld,%s\n",START,CHR);
		strcpy(PREVIOUS_CHR,CHR);

	 }

	/* Calling on the compute function here again to make sure the last line is executed too */
	dupvector = computeVector(vector,kernel,SIGMA);
	//dupvector = memoryAllocate(dupvector,VSIZE);
	//dupvector = computeVector(vector,kernel,dupvector,SIGMA);
	//free(dupvector);
	printArray(dupvector,VSIZE);




}


float* computeVector(float *vector1,float *kernel1,int sigma){

	long j;
	float *smoothedVec;
	smoothedVec = memoryAllocate(smoothedVec,VSIZE);

	/* looping throught the whole original vector of size VSIZE */

	for(j=sigma;j<VSIZE-sigma;j++){
		//if((j-sigma) < 0 || (j+sigma) > VSIZE)  /* checking for the margins */
		//	continue;
		if(vector1[j] >= 1)      /* start from the value, where the value in the original vector is non-zero*/
		{
			int start = j - sigma; /* start from the current position - sigma positions */
			int end   = j+ sigma;  /* end to the current position + sigma positions */
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
		printf("Insufficient memory\n");
			}
	return v;

}


float* computeSmoothing(float *sm, int sigma){

    int i,j=0;
    for(i= -(sigma);i <=sigma; i++){
    	sm[j] += exp(-pow(i,2)/(2*(pow(sigma,2))));
    	j++;
    }
	return sm;

}

void printArray(float *arr, long size)
{
	long i;
	for(i=0;i<size;i++)
	{
		if(arr[i] <= 0)
			continue;
		printf("%f,%ld\n",arr[i],i);
	}
}

