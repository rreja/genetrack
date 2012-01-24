
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


void printArray(float *, int);
void computeVector();
float* populate(float *, long);
float* memoryAllocate(float *,long);
float* computeSmoothing(float *,int);





void main ()
{
	FILE *fp; /* pointer to the file */
	fp = fopen("/Users/rohitreja/Desktop/test.txt","r");   /* reading the file and assigning it to the pointer*/

	char line[500];   /* char array to store each line*/
	char *toks[7];   /* toks is an array of 7 elements, each of which points to a char */
	char CHR[10],STRAND[5];
	long START;
	char PREVIOUS_CHR[10] = "NULL";
	static float *vector;   /* vector to be populated */
	float *smoothing;       /* smoothing vector to be populated */
	long VSIZE = 10000000; /* size of the vector to initialize to */
	long SSIZE = 1000;     /* size of smoothing vector to initialize to*/
    int SIGMA = 5;         /* SIGMA value to be take by user */


	vector = memoryAllocate(vector,VSIZE);         /* call to calloc to allocate and initialize vector */
	smoothing = memoryAllocate(smoothing,SSIZE);   /* ONE TIME call to calloc to allocate and initialize smoothing vector */
    smoothing = computeSmoothing(smoothing,SIGMA);  /* Call to populate the smoothing vector */

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

		if(strcmp(CHR,PREVIOUS_CHR)){

			computeVector();
			/*
			 * Call the function to compute on vector and in the compute function
			 * check in the vector in empty.
			 * Compute on the vector
			 */
			free(vector);                                      /* making memory free here */
			vector = memoryAllocate(vector,VSIZE);

		}

		populate(vector,START);
		printf("%f,%s\n",vector[START],CHR);
		strcpy(PREVIOUS_CHR,CHR);


	 }
	computeVector();

	/* Call on the compute function here again to make sure the last line is executed too */

}


void computeVector(){
printf("this is function to compute\n");

}

float* populate(float *vec, long pos){
     vec[pos]++;
     return vec;
}

float* memoryAllocate(float *v,long vsize){

	v = calloc(vsize,sizeof(float));          /*initializing the vector*/
	if(v == NULL){
		printf("Sufficient memory not available\n");
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

void printArray(float *arr, int size)
{
	int i;
	for(i=0;i<size;i++)
	{
		printf("%f\n",*(arr+i));
	}
}

