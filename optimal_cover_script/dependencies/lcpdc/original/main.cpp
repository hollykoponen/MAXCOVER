
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"
#include "dc_lcp_oracle.h"

int main(int argc, char **argv){
   if(argc != 3){
      fprintf(stderr,"Usage: %s {inputfile} {logv}\n",argv[0]);
      exit(1);
   }	
   double time, totalTime;
   char filename[100],alg;
   FILE *fp;
   int n,m;
   unsigned char *x; //the text
   int *SA; //the suffix array
   int logv; //the logarithm of the difference cover to use

   //read in the text and SA
   fp = fopen(argv[1],"r");
    	
   if(!fp){
      fprintf(stderr,"Could not open file %s.\n",argv[1]);
   }
   fseek(fp,0L,SEEK_END);
   n = ftell(fp);
   //printf("%d \n", n);
   x = (unsigned char *)malloc(sizeof(unsigned char)*(n+1));
   rewind(fp);
   fread(x,sizeof(unsigned char),n,fp);
   x[n] = 0;
   //printf("%s", x);
   fclose(fp);

  SA = (int *)malloc(sizeof(int)*n);
   strcpy(filename,argv[1]);
   strcat(filename,".sa");
   fp = fopen(filename,"r");
   if(!fp){
      fprintf(stderr,"Could not open file %s.\n",filename);
      free(x);
      free(SA);
      exit(1);
   }
   fread(SA,sizeof(int),n,fp);
   fclose(fp);
   /*for(int j=0; j < n; j++){
	printf("%d=%d\n", j, SA[j]);
   }*/
   logv = atoi(argv[2]);

   time = getTime();
   dc_lcp_initialize(x,SA,n,logv);
   totalTime = getTime() - time;
   fprintf(stderr,"Total setup time: %.2f secs\n",totalTime);
   //exit(1);
   time = getTime();
   dc_lcp_construct();
   totalTime += (getTime() - time);
   fprintf(stderr,"Total lcp construction time: %.2f secs\n",totalTime);

   //int max = -1;
   //for(int i=0;i<n;i++){
   //   if(SA[i] > max){
   //      max = SA[i];
   //   }
   //}
   //fprintf(stderr,"Max LCP = %d\n",max);

   //#ifdef WRITE_LCP_TO_DISK
   strcpy(filename,argv[1]);
   strcat(filename,".lcp");
   fp = fopen(filename,"w");
   fwrite(SA,4,n-1,fp);
   //#endif
	
   dc_lcp_free();
   free(x);
   free(SA);
   return 0;
}
