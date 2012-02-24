#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "gopt.h"
#include "khash.h"

/* Since the max size of human chromosome is 300,000,000, this program cannot allocate that much memory and dies,
 * possible bug to consider when running human samples. Also, because of this memory restrictions, the size of peakArray
 * cannot be more than 200000.
 */

KHASH_MAP_INIT_STR(str, int)

void printWiggle(float *, long, char [], char *);
float* computeVector(float *,float *,int);
float* populate(float *, long);
float* memoryAllocate(float *,long);
float* computeSmoothing(float *,int);
void findPeaks(float *);
void computeAll(float *,float *, char [],char *);

/* Constant declarations */
long VSIZE, PSIZE, DEBUG=1; /* vsize is the sizeof the vector to initialize to  and psize is the size of the structure to initialize to*/
long SSIZE = 1000;     /* size of smoothing vector to initialize to*/
int const N = 4;        /* the parameter to contorl the spread of sigma, tells the program to go +/- 4*sigma */
const char *outfilename, *infilename,*wiggleout, *bedout, *sig,*ex, *maxsize,*maxpeak,*tg;
int SIGMA, EXCLUSION,BED=0,BEDOUT = 0,printwig = 0,GFFOUT=0,CHR_IDX = 0,FWD_START_IDX=3,REV_START_IDX=4,STRAND_IDX=6;
float TAGCOUNT;
FILE *fp, *op, *fwg,*rwg,*ob; /* pointer to the input and output file */
char msg[1000],fwdwg[100] = "fwd_",revwg[100] = "rev_";


/* Declaration of other global variables */
struct peaks {
				float height;
				long value;
				int flag;
		}*pA;

int revcmp(struct peaks *,struct peaks *); /*Function to be used in qsort */
void printGff(struct peaks *, char [],FILE *, float *, char *);
void printBed(struct peaks *, char [],FILE *, float *, char *);
void callPeaks(struct peaks *);

void logging(char *msg){
    if (DEBUG){
        printf("*** %s\n", msg);
    }
}

void error(char *msg){
    fprintf(stderr,"(!) error %s\n", msg);
    exit(-1);
}

/* FUNCTION DEFINITIONS START HERE */

void computeAll(float *vec,float *kernel, char chr[10],char *strand){
     float *dupvector;
     
     dupvector = computeVector(vec,kernel,SIGMA);
     if(printwig)
        printWiggle(dupvector,VSIZE,chr,strand);
     findPeaks(dupvector);
     callPeaks(pA);
     if(BEDOUT)
         printBed(pA,chr,ob,vec,strand);
     if(GFFOUT)   
         printGff(pA,chr,op,vec,strand);
   
     /* making memory free here */
     free(vec);
     free(dupvector);
     free(pA);
                                
                                
}

void callPeaks(struct peaks *pM){

  long k,z;
  qsort(pM,PSIZE,sizeof(struct peaks),(void *)revcmp);
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
    float reqmem;
    reqmem = (2*VSIZE+PSIZE)/(pow(1024,3));
    pA = calloc(PSIZE,sizeof(struct peaks));
   		 if (pA == NULL)
                     {
                          fprintf(stderr,"ERROR: Out of memory, you need a minimum of %f\n",reqmem);
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
     if(pos >= VSIZE)
     {
          fprintf(stderr,"ERROR: You need to increase the size of the parameter 'N', current size=%ld\n",VSIZE);
          exit(0);
     }
     else
     {
         vec[pos]++;
         return vec;                       
     }
     
}

float* memoryAllocate(float *v,long vsize){
    long i;
    float r;
    
    v = calloc(vsize, sizeof(float));          /*initializing the vector*/
	r = vsize/pow(1024,3) * sizeof(float);
    
    if(v == NULL){
		fprintf(stderr,"ERROR: Out of memory, while allocating %fGb of RAM\n",r );
    }
    
    sprintf(msg, "allocating %.4f Gb of RAM", r); 
    logging(msg);
    
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

void printWiggle(float *arr, long size, char chr[10],char *strand){
	long i;char temp[10];
	strcpy(temp,chr);
        
        if(!(strcmp(strand,"-"))){
           fprintf(rwg,"variableStep chrom=%s\n",temp);
           for(i=0;i<size;i++)
	{
		if(arr[i] <= 1)
			continue;
		fprintf(rwg,"%ld %f\n",i,arr[i]);
	}
        fflush(rwg);
           
        }
        else{
            fprintf(fwg,"variableStep chrom=%s\n",temp);
	    for(i=0;i<size;i++){
		if(arr[i] <= 1)
			continue;
		fprintf(fwg,"%ld %f\n",i,arr[i]);
	}
        fflush(fwg);                               
        }
	
}

void printGff(struct peaks *pB, char chr[10],FILE *op, float *vec,char *strand){
	long m,start,end,z,sum;
        float var, sd;
	char temp[10];
	strcpy(temp,chr);
	for(m=0;m<PSIZE;m++)
	{
                sum = 0;
                var = 0;
		if((pB[m].height <=TAGCOUNT) || (pB[m].flag == 0))    /*  change here, if you want to print by peak height */
			continue;
                start = pB[m].value - EXCLUSION/2;
		end   = pB[m].value + EXCLUSION/2;
                for(z=start;z<=end;z++){
                      sum += vec[z];
                      var+= pow((pB[m].value - z),2)*vec[z];
                               
                }
                sd = sqrt(var/sum);
		fprintf(op,"%s\t%s\t%s\t%ld\t%ld\t%f\t%s\t%ld\treadsum=%ld;fuzziness=%f;\n",temp,"genetrack",".",start,end,pA[m].height,strand,sum,sum,sd);
		
	}
        fflush(stdout);
        fflush(op);
}

void printBed(struct peaks *pB, char chr[10],FILE *ob, float *vec,char *strand){
            
        long m,start,end;
	char temp[10];
	strcpy(temp,chr);
	for(m=0;m<PSIZE;m++)
	{
                
		if((pB[m].height <=TAGCOUNT) || (pB[m].flag == 0))    /*  change here, if you want to print by peak height */
			continue;
                start = (pB[m].value - EXCLUSION/2) - 1;
		end   = pB[m].value + EXCLUSION/2;
		fprintf(ob,"%s\t%ld\t%ld\t%s\t%s\t%s\n",temp,start,end,".",".",strand);
		
	}
        fflush(ob);                       
}

/* Start of main function */

int main (int argc, const char **argv){
    
    int ret, is_missing;
    char *dup;
    
    khash_t(str) *h;
    h = kh_init(str);
    khiter_t k;
    
               
      void *options= gopt_sort( & argc, argv, gopt_start(
      gopt_option( 'h', 0, gopt_shorts( 'h', '?' ), gopt_longs( "help", "HELP" )),
      gopt_option( 's',GOPT_ARG, gopt_shorts('s'), gopt_longs( "sigma" )),
      gopt_option( 'e',GOPT_ARG, gopt_shorts('e'), gopt_longs( "exclusion" )),
      gopt_option( 'c',GOPT_ARG, gopt_shorts('c'), gopt_longs( "tagcount" )),
      gopt_option( 'i', GOPT_ARG, gopt_shorts( 'i' ), gopt_longs( "input" )),
      gopt_option( 'b', 0, gopt_shorts( 'b' ), gopt_longs( "bed" )),
      gopt_option( 'w', GOPT_ARG, gopt_shorts( 'w' ), gopt_longs( "wiggleout" )),
      gopt_option( 't', GOPT_ARG, gopt_shorts( 't' ), gopt_longs( "bedout" )),
      gopt_option( 'N', GOPT_ARG, gopt_shorts( 'N' ), gopt_longs( "maxsize" )),
      gopt_option( 'P', GOPT_ARG, gopt_shorts( 'P' ), gopt_longs( "maxpeak" )),
      gopt_option( 'o', GOPT_ARG, gopt_shorts( 'o' ), gopt_longs( "output" ))));

      //FILE *fp, *op, *wg; /* pointer to the input and output file */
      
      if( gopt( options, 'h' ) ){
          fprintf( stdout, "\nUsage: ./genetrack -i <inputfile.gff> -s <sigma> -e <exclusion> -o <outputfile>\n");
          fprintf(stdout,"Arguments:\n-i <inputfile.gff>, in gff3 format\n");
          fprintf(stdout,"-b: if the input file is in BED format\n");
          fprintf(stdout,"-o <outputfile>,   output in GFF format\n");
          fprintf(stdout,"-t <outputfile>,  output in BED format\n");
          fprintf(stdout,"-s <smoothing>, Smoothin parameter, default = 5\n");
          fprintf(stdout,"-e: <exclusion>, Exclusion zone, default = 20\n");
          fprintf(stdout,"-c: <tag_count_threshold>, Peaks containing more than 'c' tags, default = 3\n");
          fprintf(stdout,"-w: <wigglefilename>, output in wiggle format\n");
          fprintf(stdout,"-N: The size of the largest chromosome, in millions, ex: 300 for 300,000,000 bp, default = 15\n");
          fprintf(stdout,"-P: Max Expected number of peaks in the largest chromosome, in millions,ex: 1 for 1,000,000, default = 0.1\n");
          exit( EXIT_SUCCESS );
      }
      
      if( gopt_arg(options, 's', &sig) && strcmp(sig, "-" )){
         SIGMA = atoi(sig);
      }else{SIGMA = 5;}
      
      if( gopt_arg(options, 'e', &ex) && strcmp(ex, "-" )){
         EXCLUSION = atoi(ex);
      }else{EXCLUSION = 20;}
      
      if( gopt_arg(options, 'c', &tg) && strcmp(tg, "-" )){
         TAGCOUNT = atoi(tg);
      }else{TAGCOUNT = 2.8;}
      
      if( gopt_arg(options, 'N', &maxsize) && strcmp(maxsize, "-" )){
         VSIZE = atoi(maxsize) * pow(10,6);
      } else { VSIZE = 15000000;}
      
      if( gopt_arg(options, 'P', &maxpeak) && strcmp(maxpeak, "-" )){
         PSIZE = atoi(maxpeak) * pow(10,6);
      }else{PSIZE = 100000;}
      
      if( gopt( options, 'b' ) ){
            CHR_IDX = 0;
            FWD_START_IDX=1;
            REV_START_IDX=2;
            STRAND_IDX=5;
            BED = 1;
            
      }
      
      if( gopt_arg(options, 'i', & infilename) && strcmp(infilename, "-" )){
    
          fp = fopen(infilename,"r");
          if(fp == NULL){
          fprintf(stderr, "%s: %s: could not open file for input\n", argv[0], infilename);
          exit(EXIT_FAILURE);
          } 
      }
      else{fp = stdin;}
      
      if( gopt_arg(options, 't', &bedout) && strcmp(bedout, "-" )){
    
          ob = fopen(bedout,"w");
          BEDOUT = 1;
          if(ob == NULL){
          fprintf(stderr, "%s: %s: could not open file for input\n", argv[0], bedout);
          exit(EXIT_FAILURE);
          } 
      }
      
      
      if( gopt_arg(options, 'w', &wiggleout ) && strcmp(wiggleout, "-" )){
                                
          strcat(fwdwg,wiggleout);
          strcat(revwg,wiggleout);
          fwg = fopen(fwdwg,"w");
          rwg = fopen(revwg,"w");
          printwig = 1;
          if(fwg == NULL || rwg == NULL){
          fprintf(stderr, "%s: %s: could not open file for wiggle output\n", argv[0], wiggleout);
          exit(EXIT_FAILURE);
          } 
    
      }
        
     if( gopt_arg(options, 'o', & outfilename ) && strcmp( outfilename, "-" )){
    
        op = fopen(outfilename,"w");
        fprintf(op,"%s\n","##gff-version 3");                   /* line to print the header for gff*/
        GFFOUT = 1;
        if(op == NULL){
        fprintf(stderr, "%s: %s: could not open file for output\n", argv[0], outfilename );
        exit(EXIT_FAILURE);
        }   
    }
    else{op = stdout;}          /*  Uncomment if you want to print to screen*/
    
	char line[500];   /* char array to store each line*/
	char *toks[10];   /* toks is an array of 10 elements, each of which points to a char */
	char CHR[10],STRAND[5], *ptr;;
	long START;
	char PREVIOUS_CHR[10] = "NULL";
	float *fwdvec,*revvec;       /* vectors to be populated */
	float *kernel;       /* smoothing vector to be populated */

	fwdvec = memoryAllocate(fwdvec,VSIZE);         /* call to calloc to allocate and initialize vector */
      revvec = memoryAllocate(revvec,VSIZE);         /* call to calloc to allocate and initialize vector */
	kernel = memoryAllocate(kernel,SSIZE);        /* ONE TIME call to calloc to allocate and initialize smoothing vector */
	kernel = computeSmoothing(kernel,SIGMA);      /* Call to populate the smoothing vector */
        
	while(fgets(line,500,fp) != NULL)             /* fgets */
	{

		if(line[0] == '#')                    /* ignore all the comment lines that start with # */
			continue;

		int cols = 0;
		toks[cols] = strtok(line,"\t");
		while(cols < 7)
				{
				    toks[++cols] = strtok(NULL,"\t");                    /* storing all token in toks for one line */

				}
		strcpy(CHR,toks[CHR_IDX]);
		if((strcmp(CHR,PREVIOUS_CHR)) && (strcmp("NULL",PREVIOUS_CHR))){
            
                      sprintf(msg, "processsing chromosome %s", PREVIOUS_CHR);
                      logging(msg);
                      dup = strdup(PREVIOUS_CHR);
                       k = kh_put(str, h, strdup(PREVIOUS_CHR), &ret);
                
                       if (!ret) {
                           sprintf(msg, "File not sorted,please sort your file by chromosome");
                             error(msg);
                          }
               
                                computeAll(fwdvec,kernel,PREVIOUS_CHR,"+");
                                computeAll(revvec,kernel,PREVIOUS_CHR,"-");
                                fwdvec = memoryAllocate(fwdvec,VSIZE);
                                revvec = memoryAllocate(revvec,VSIZE);

		}
                
                if( (ptr = strchr(toks[STRAND_IDX], '\n')) != NULL)   /* remove the newline characer from the strand for bed files */
                    *ptr = '\0';

                if(!(strcmp(toks[STRAND_IDX],"+"))){
                        START = atoi(toks[FWD_START_IDX]);
                        fwdvec = populate(fwdvec,START+BED);  /* when -b flag is on, add 1 to the bed file start co-ordinate */
                }else{
                        START = atoi(toks[REV_START_IDX]);    
                        revvec = populate(revvec,START);     /* we do not add +1 here as this is the reverse strand and we are considerinf the 5' end */
                } 

		 strcpy(PREVIOUS_CHR,CHR);

	 }

	/* Calling on the compute function here again to make sure the last line is executed too */
        sprintf(msg, "processsing chromosome %s", PREVIOUS_CHR);
        logging(msg);
        computeAll(fwdvec,kernel,CHR,"+");
        computeAll(revvec,kernel,CHR,"-");
	fclose(fp);
	fclose(op);

}
