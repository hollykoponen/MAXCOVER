
#include <stdio.h>
#include <stdlib.h>

#include "RMQ_succinct.hpp"
#include "secondary_dc_lcp.h"

#define MAX(i,j) (((i) > (j)) ? (i) : (j))

const int max_precomputed_cover = 8;
const int coversizes[max_precomputed_cover+1] = {1,2,3,4,5,7,9,13,20};
const unsigned cover0[] = {0};
const unsigned cover1[] = {0,1};
const unsigned cover2[] = {0,1,2};
const unsigned cover3[] = {0,1,2,4};
const unsigned cover4[] = {0,1,2,5,8};
const unsigned cover5[] = {0,1,2,3,7,11,19};   //{0,7,8,10,14,19,23};
const unsigned cover6[] = {0,1,2,5,14,16,34,42,59};
const unsigned cover7[] = {0,1,3,7,17,40,55,64,75,85,104,109,117};
const unsigned cover8[] = {0,1,3,7,12,20,30,44,65,80,89,96,114,122,
                           128,150,196,197,201,219};
const unsigned * _covers[] = { cover0, cover1, cover2, cover3, cover4,
                              cover5, cover6, cover7, cover8 };

const int _cover_sizes[] = {1,2,3,4,5,7,9,13,20}; 

unsigned *_cover;
int *_rev_cover;
int *_delta;
int _cover_size;
int _logv;
int _v;
int _mask;
int _n;
FILE *_sa_fp;
int *_sa_buffer;
int _sa_buf_size;
unsigned char *_x;
int _m;
int *_ess;
int *_rev_ess;
int *_ell;
RMQ_succinct *_rmq;

int delta_sec(unsigned int i, unsigned int j){
   return ((_delta[(j-i)%_v]-i)%_v);
}

#define BUF_CAP 1048576

void dc_lcp_initialize(unsigned char *a_x, FILE *a_sa_fp, int a_n, int a_logv){
   if(a_logv > max_precomputed_cover){
      fprintf(stderr,"Specified DC (%d) greater than max (%d)\n",a_logv,max_precomputed_cover);
      exit(1);
   }
   _logv = a_logv;
   _v = (1 << _logv);
   _mask = _v - 1;
   _n = a_n;
   _x = a_x;
   _sa_fp = a_sa_fp;
   _sa_buffer = (int *)malloc(sizeof(int)*BUF_CAP);
   _sa_buf_size = 0;
   _cover = (unsigned *)(_covers[_logv]);
   _cover_size = _cover_sizes[_logv];
   _m = _n/_v;
   _m = _m * _cover_size;
  
   //fprintf(stderr,"_m = %d\n",_m);
   //fprintf(stderr,"_v = %d\n",_v);
 
   //compute _rev_cover
   _rev_cover = (int *)malloc(sizeof(int) * _v);
   int j = 0;
   for(int i = 0; i < _v; i++){
      _rev_cover[i] = -1;
      if(_cover[j] == i){
         _rev_cover[i] = j;
         if(j < (_n & _mask)){
            _m++;
         }
         j++;
      }
   }

   //compute _delta
   _delta = (int *)malloc(sizeof(int) * _v);
   for (int i = _cover_size-1; i >= 0; i--) {
      for (j = 0; j < _cover_size; j++) {
         _delta[(_cover[j]-_cover[i])%_v] = _cover[i];
      }
   }

   //compute arrays _ess, _rev_ess
   _ess = (int *)malloc(sizeof(int) * _m);
   _rev_ess = (int *)malloc(sizeof(int) * _m);
   j = 0;
   for(int i = 0; i < _n; ){
      _sa_buf_size = fread(_sa_buffer,sizeof(int),BUF_CAP,_sa_fp);
      for(int k = 0; k < _sa_buf_size; k++){
         int si = _sa_buffer[k];
         if(_rev_cover[si&_mask] != -1){
            //this is a sample suffix
            _ess[j] = si;
            _rev_ess[(_cover_size*(si>>_logv)) + _rev_cover[si&_mask]] = j;
            j++;
         }
      }
      i += _sa_buf_size;
   }
   //fprintf(stderr,"j = %d\n",j);

   //compute _ell using _ess, _rev_ess, _rev_cover
   _ell = (int *)calloc(_m, sizeof(int));
   _ell[0] = 0;
   int len = 0;
   int computed = 0;
   int *lengths = (int *)malloc(sizeof(int) * _cover_size);
   for(int i = 0; i < _cover_size; i++){
      lengths[i] = 0;
   }
   int compares_saved = 0;
   for(int i = 0; i < _m; i++){
      //if(i % 1000000 == 0){
      //   fprintf(stderr,"Processed %d entires\n",i);
      //}
      int ihat = _rev_ess[i];
      len = lengths[i%_cover_size];
      if(len < 0) len = 0;
      compares_saved += len;
      if(ihat > 0){
         int j = _ess[ihat-1];
         while(_ess[ihat]+len < _n && j+len < _n){
            if(_x[_ess[ihat]+len] != _x[j+len]){
               break;
            }
            len++;
         }
      }
      _ell[ihat] = len;
      int a = _rev_cover[(_ess[ihat])&_mask];
      lengths[a] = len - _v;
      //#define CHECK_SAMPLE_LCP
      #ifdef CHECK_SAMPLE_LCP
      int l = 0;
      while(_x[_ess[ihat]+l] == _x[_ess[ihat-1]+l]){
         l++;
      }
      if(l != len){
         fprintf(stderr,"mismatch at i=%d, (%d,%d)\n",i,l,len);
         exit(1);
      }else{
      }
      #endif
      computed++;
   }
   //fprintf(stderr,"computed %d lcp values\n",computed);
   //fprintf(stderr,"compares_saved = %d\n",compares_saved);

   //we no longer need _ess, free it
   free(_ess);

   //preprocess _ell for RMQs
   _rmq = new RMQ_succinct(_ell,_m);
}

void dc_lcp_construct(){
   rewind(_sa_fp);
   int longcount = 0;
   int shortcount = 0;
   int s0;
   int start,mark;
   for(int i = 0; i < _n; ){
      //if(i % 10000000 == 0){
      //   fprintf(stderr,"Processed %d entires\n",i);
      //}
      mark = ftell(_sa_fp);
      _sa_buf_size = fread(_sa_buffer,sizeof(int),BUF_CAP,_sa_fp);
      if(i != 0){
         start = 0;
         //s0 is already set
      }else{
         start = 1;
         s0 = _sa_buffer[0];
      }
      for(int k = start; k < _sa_buf_size; k++){
         int s1 = _sa_buffer[k];
         //check if lcp(SA[i-1],SA[i]) < v
         int j = 0;
         while(_x[s0+j] == _x[s1+j] && j < _v){
            j++;
         }
         if(j < _v){
            shortcount++;
            _sa_buffer[k-start] = j;
         }else{
            longcount++;
            int ds0s1 = delta_sec((unsigned int)s0,(unsigned int)s1);
            int a0 = s0 + ds0s1;
            int a1 = s1 + ds0s1;
            //fprintf(stderr,"s0 = %d; s1 = %d; ds0s1 = %d; a0 = %d; a1 = %d\n",s0,s1,ds0s1,a0,a1);
            int r0 = _rev_ess[(_cover_size*(a0>>_logv)) + _rev_cover[a0&_mask]];
            int r1 = _rev_ess[(_cover_size*(a1>>_logv)) + _rev_cover[a1&_mask]];
            int rmin_auto = _ell[(_rmq->query(r0+1,r1))];
            _sa_buffer[k-start] = ds0s1 + rmin_auto;
            //#define CHECK_LCP
            #ifdef CHECK_LCP
            int j = 0;
            while(_x[s0+j] == _x[s1+j]){
               j++;
            }
            if(j != _sa_buffer[k-start]){
               fprintf(stderr,"mismatch at position %d (%d %d); r0 = %d, r1 = %d\n",i,j,_sa_buffer[k-start],r0,r1);
               int j = 0;
               while(_x[a0+j] == _x[a1+j]){
                  j++;
               }
               fprintf(stderr,"ds0s1 = %d; a0 = %d; a1 = %d; j = %d\n",ds0s1,a0,a1,j);
            }
            #endif
         }
         s0 = s1;
      }
      fseek(_sa_fp,mark,SEEK_SET);
      fwrite(_sa_buffer,sizeof(int),_sa_buf_size,_sa_fp);
      i += _sa_buf_size;
   }
   //fprintf(stderr,"longcount = %d; shortcount = %d\n",longcount,shortcount);
   fclose(_sa_fp);
}

void dc_lcp_free(){
   free(_rev_cover);
   free(_ell);
   free(_rev_ess);
   //free(_lcp);
   delete _rmq;
}


