
#include <stdio.h>
#include <stdlib.h>

#include "RMQ/RMQ_succinct.hpp"
#include "dc_lcp_oracle.h"

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
int *_sa;
int *_lcp;
unsigned char *_x;
int _m;
int *_ess;
int *_rev_ess;
int *_ell;
RMQ_succinct *_rmq;

int delta(unsigned int i, unsigned int j){
   return ((_delta[(j-i)%_v]-i)%_v);
}

void dc_lcp_initialize(unsigned char *a_x, int *a_sa, int a_n, int a_logv){
   if(a_logv > max_precomputed_cover){
      fprintf(stderr,"Specified DC (%d) greater than max (%d)\n",a_logv,max_precomputed_cover);
      exit(1);
   }
   _logv = a_logv;
   _v = (1 << _logv);
   _mask = _v - 1;
   _n = a_n;
   _x = a_x;
   _sa = a_sa;
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
   for(int i = 0; i < _n; i++){
      int si = _sa[i];
      if(_rev_cover[si&_mask] != -1){
         //this is a sample suffix
         _ess[j] = si;
         _rev_ess[(_cover_size*(si>>_logv)) + _rev_cover[si&_mask]] = j;
         j++;
      }
   }
   //fprintf(stderr,"j = %d\n",j);

   //compute _ell using _ess, _rev_ess, _rev_cover
   _ell = (int *)malloc(_m * sizeof(int));
    for (int i = 0; i < _m; ++i) {
        _ell[i] = 0;
    }
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
    // ADDED CODE
    int backup_m = _m;
    if (_m < 113) {
        backup_m = (_m + (113 - _m));
        _ell = (int*)realloc(_ell, sizeof(int) * backup_m);
        for (int i = _m; i < backup_m; ++i) {
            _ell[i] = 0;
        }
    }
    // ADDED CODE
    _rmq = new RMQ_succinct(_ell, backup_m);
//    fprintf(stderr,"Peak memory = %.2fn\n",(double)(4*3*_m)/(double)_n);
}

void dc_lcp_construct(){
   _lcp = _sa;
   int longcount = 0;
   int shortcount = 0;
   int maxlcp = 0;
   for(int i = 1; i < _n; i++){
      //if(i % 10000000 == 0){
      //   fprintf(stderr,"Processed %d entires\n",i);
      //}
      int s0 = _sa[i-1];
      int s1 = _sa[i];
      //check if lcp(SA[i-1],SA[i]) < v
      int j = 0;
      while(_x[s0+j] == _x[s1+j] && j < _v){
         j++;
      }
      if(j < _v){
         shortcount++;
         _lcp[i-1] = j;
      }else{
         longcount++;
         int ds0s1 = delta((unsigned int)s0,(unsigned int)s1);
         int a0 = s0 + ds0s1;
         int a1 = s1 + ds0s1;
         //fprintf(stderr,"s0 = %d; s1 = %d; ds0s1 = %d; a0 = %d; a1 = %d\n",s0,s1,ds0s1,a0,a1);
         int r0 = _rev_ess[(_cover_size*(a0>>_logv)) + _rev_cover[a0&_mask]];
         int r1 = _rev_ess[(_cover_size*(a1>>_logv)) + _rev_cover[a1&_mask]];
         int rmin_auto = _ell[(_rmq->query(r0+1,r1))];
         _lcp[i-1] = ds0s1 + rmin_auto;
         if(_lcp[i-1] > maxlcp){
            maxlcp = _lcp[i-1];
         }
         //#define CHECK_LCP
         #ifdef CHECK_LCP
         int j = 0;
         while(_x[s0+j] == _x[s1+j]){
            j++;
         }
         if(j != _lcp[i-1]){
            fprintf(stderr,"mismatch at position %d (%d %d); r0 = %d, r1 = %d\n",i,j,_lcp[i-1],r0,r1);
            int j = 0;
            while(_x[a0+j] == _x[a1+j]){
               j++;
            }
            fprintf(stderr,"ds0s1 = %d; a0 = %d; a1 = %d; j = %d\n",ds0s1,a0,a1,j);
         }
         #endif
      }
   }
//   fprintf(stderr,"maxlcp = %d\n",maxlcp);
   //fprintf(stderr,"longcount = %d; shortcount = %d\n",longcount,shortcount); 	
}

void dc_lcp_free(){
   free(_rev_cover);
   free(_ell);
   free(_rev_ess);
   //free(_lcp);
   delete _rmq;
}


