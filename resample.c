/**********************************************************************

  resample function ammended from testresample.c

  Real-time library interface by Dominic Mazzoni

  Based on resample-1.7:
    http://www-ccrma.stanford.edu/~jos/resample/

  License: LGPL - see the file LICENSE.txt for more information

**********************************************************************/
//#include <stdio.h>
//#include <stdlib.h>

#include "libresample.h"
#include <malloc.h>
//extern char* hoc_Emalloc(unsigned long);

int resample(double *source, double *dest, int srclen, double factor)
{
   void *handle;
   float *dst,*src;
   int i, out, o, srcused, errcount, rangecount;
   int statlen, srcpos, lendiff;
   int fwidth,srcBlock,lastFlag;
   int srcblocksize = srclen;
   int dstblocksize = (int)(srclen * factor + 10);
   int expectedlen = (int)(srclen * factor);
   int dstlen = expectedlen + 1000;
#ifdef DEBUG
   printf("-- srclen: %d  factor: %g srcblk: %d dstblk: %d expected %d\n",srclen,  factor, srcblocksize, dstblocksize,expectedlen);
#end
   src = (float*)malloc((unsigned) (srclen*sizeof(float))); // makevector(srclen);//
   for (i=0;i<srclen;i++) src[i] =  source[i];
   dst = (float*)malloc((unsigned) ((dstlen+100)*sizeof(float)));//makevector(dstlen+100);

   //   printf(" source  %x\t src %x\t dst %x\n",&source[0],&src[0],&dst[0]);
   handle = resample_open(1, factor, factor);
   fwidth = resample_get_filter_width(handle);
#ifdef DEBUG
   printf("lresample:  starting loop\n");
#end 
   out = 0;
   srcpos = 0;
   for(;;) {
      srcBlock = __min(srclen-srcpos, srcblocksize);
      lastFlag = (srcBlock == srclen-srcpos);

      o = resample_process(handle, factor,
                           &src[srcpos], srcBlock,
                           lastFlag, &srcused,
        &dst[out], __min(dstlen-out, dstblocksize));
      srcpos += srcused;
      if (o >= 0)
         out += o;
      if (o < 0 || (o == 0 && srcpos == srclen))
         break;
   }
   resample_close(handle);

   if (o < 0) {
      hoc_execerror("lresample: resample_process returned an error ", 0);
      return 0;
   }else {
#ifdef DEBUG
printf("lresample: resample done o %d out %d\n",o);
#end 
   }	   

   if (out <= 0) {
     hoc_execerror(" resample_process returned less than zero samples",0);
     free(src);
     free(dst);
     return 0;
   }

   lendiff = abs(out - expectedlen);

   if (lendiff > (int)(2*factor + 1.0)) {
      printf("lresample:   Expected ~%d, got %d samples out\n", expectedlen, out);
   }
#ifdef DEBUG   
   printf("lresample  1: dst[0] %g\t dst[end] %g\t dstlen %d\n",dst[0],dst[dstlen+99], dstlen);
#end
   int len =(int)(srclen * factor);
 for (i=0;i<len;i++) {
   if ( isnan(dst[i]) ) {len = i-1; break;}
   dest[i] = (double) dst[i];
   // printf("dst[%d]\t %g\n",i,dst[i]);
 }

#ifdef DEBUG
 printf("lresample  3:  len %d \tdest[0] %g\t dest[len] %g src %x dst %x\n",len,dest[0],dest[len-1], &src,&dst);
#end
   if (src!=NULL) free(src); //(char*)
   if (dst!=NULL) free(dst);
   
 
return len;
}


/*  resample apapted from src/ivoc/ivocvect.cpp */
 /*

int resample(double* v1, double* ans, int capacity, double f)
{
    //  double * temp;
    int i;

    printf("resample: &v1 %x\t capacity %d  factor %g destlen %d\n", &v1[0], capacity, f, (int)(capacity * f));

    if (f < 0 || f > capacity / 2) {
        printf("ivoc_resample: nyquist error with factor\n");
        //return 0;
    }
    int n = (int)(capacity * f);

    //    ans = (double*)makevector(n); zero_vector(ans, n);
    for (i = 0; i < n; i++) {
        if (isnan(v1[(int)(i/f)])) {
            printf("nan in v1: ans address %x", &ans[0]);
            return n-1;
        }
        ans[i] = v1[(int)round(i/f)];
    }

    return n;
}

 */
