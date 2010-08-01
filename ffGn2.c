 /* FFGN Fast (exact) fractional Gaussian noise and Brownian motion
  * generator.
  *	
  *	Y = FFGN(N, Hinput, MU, SIGMA) returns a vector containing a sequence of fractional Gaussian 
  *	noise or fractional Brownian motion.  The generation process uses an FFT which makes it 
  *	very fast.  The input arguments are:
  *
  *		N			is the length of the output sequence.
  *		H			is the "Hurst" index of the resultant noise (0 < H <= 2).  For 0 < H <= 1, 
  *					  the output will be fractional Gaussian noise with Hurst index H.  For 
  *					  1 < H <= 2, the output will be fractional Brownian motion with Hurst
  *					  index H-1.  Either way, the power spectral density of the output will
  *					  be nominally proportional to 1/f^(2H-1).
  *		mu			is the mean of the noise. [default = 0]
  *		sigma		is the standard deviation of the noise. [default = 1]
  *
  *	FFGN(N, H) returns a sequence of fractional Gaussian noise with a mean of zero
  *	and a standard deviation of one or fractional Brownian motion derived from such
  *	fractional Gaussian noise.
  *
  * 	References: Davies & Harte (1987); Beran (1994); Bardet et al., 2002
  *	This method is based on an embedding of the covariance matrix in a circulant matrix.	
  *
  *   Copyright © 2003-2005 by B. Scott Jackson
  *   Revision: 1.3    Date: Aug 28, 2008 by M. S. A. Zilany
  *                    Sigma is deifined for diff. sponts (mu) and Resampling has been introduced to be compatible with the AN model 
  *   Revision: 1.2    Date: March 14, 2005
  *   History:
  *       Rev. 1.2 - 3/14/05 - Added some additional documentation and input argument checking.
  *       Rev. 1.1 - 9/15/04 - Added the persistent variables and associated "if"-statement.
  *       Rev. 1.0 - 2/11/03 - Original version.
  */

/* Michael Eager: 28 June 2010 
   Converted to C code for NMODL conversion 
   Scoplib math functions replace math.h and mex library functions
*/



/* fft from Meschach Library */
/* get/resize vector to given dimension */
#include "matrix.h"
#include "matrix2.h"


#ifndef __min
#define __min(a,b) (((a) < (b))? (a): (b))
#endif


/* scoplib.h random functions allows control of default random functions in NEURON*/
double randn(){ return normrand(0.0,1.0);}


/* 
 static double *Zmag=NULL;
 static int Nlast=0;
 static int Hlast=0;
 static int Nfft=0;
*/


double ffGn(double *yffGn,int N, double tdres, double Hinput, double mu, double sigma)
{
/* persistent matlab variables converted to static C variables */
//static double * Zmag;
//Nfft=0;
//Nlast=0;
//Hlast=0;

 int Nfft = N;
 double Nlast=0;
 double Hlast=0;


 int i,nsize,nop,resamp,NfftHalf;	
 double k,H,fBn;
 double *y,*ytmp;
 VEC *Z_real=VNULL,*Z_im=VNULL,*Zmag=VNULL,*iZmag=VNULL;
  /*---- Check input arguments ---------- */
  
 if (N <= 0){
    hoc_execerror("Length of the return vector must be positive.",0);
    return 0.;
 }
  if (tdres > 1){
    hoc_execerror("Original sampling rate should be checked.",0);
    return 0.;
  }
 
  if ((Hinput < 0) || (Hinput >= 1)){
    hoc_execerror("The Hurst parameter must be in the interval (0,1].",0);
    return 0.;
  }
	
  /* See last statement regarding default sigma value
     if (sigma <= 0)
       hoc_error("Standard deviation must be greater than zero.",0);
  */
#ifdef DEBUG
  printf("ffGn 1: N  %d\ttdres  %f\tHinput %f\tmu %f\t sigma %f\n", N,tdres,Hinput,mu,sigma);
#endif


  //  Downsampling No. of points to match with those of Scott jackson (tau 1e-1)
  resamp = (int) ceil(1e-1 / (tdres));
  nop = N; 
  N = (int) ceil(N/resamp)+1; 
  if (N<10)
    N = 10;
#ifdef DEBUG
  printf("ffGn 2: N  %d\t tdres  %f\t Hinput %f\t mu %f\t sigma %f\t Nfft %d  resamp %d\n", N,tdres,Hinput,mu,sigma,Nfft,resamp);
#endif


  //Zmag = v_get(N); //Meschach library


  //  Determine whether fGn or fBn should be produced.
  if ( Hinput <= 1 ){
    H = Hinput;
    fBn = 0;
  }
  else{
    H = Hinput - 1;
    fBn = 1;
  }

  //  Calculate the fGn.
  if (H == 0.5){
    y = makevector(N);
    for (i=0;i<N;i++) y[i] = randn();  //  If H=0.5, then fGn is equivalent to white Gaussian noise.
  } else {
    //  If this function was already in memory before being called this time,
    //  AND the values for N and H are the same as the last time it was
    //  called, then the following (persistent) variables do not need to be
    //  recalculated.  This was done to improve the speed of this function,
    //  especially when many samples of a single fGn (or fBn) process are
    //  needed by the calling function.

    //    if (Nfft==0 || Nlast ==0 || Hlast==0 || N != Nlast || H != Hlast){
      //  The persistent variables must be (re-)calculated.
    Nfft = (int) pow(2,ceil(log10(2*(N-1))/log10(2)));
    NfftHalf = (int) round(Nfft/2);
#ifdef DEBUG
    printf("ffGn 2: N %d\t Nfft %d\t NfftHalf %d\n",N, Nfft, NfftHalf);
#endif
    if ( N==0 || Nfft == 0) return 0; 
    //k=[0:NfftHalf, (NfftHalf-1):-1:1];
    if(Zmag == VNULL) V_FREE(Zmag);
    Zmag = v_get(Nfft);// v_zero(Zmag); 
    if ( N==0 || Nfft == 0) return 0;      
    for (i =0;i<NfftHalf*2;i++){
      k= (i<NfftHalf) ? i : 2*NfftHalf-1-i;
      Zmag->ve[i]= 0.5 * (  pow(k+1.0,2.0*H) -  2.0*pow(k,2.0*H) + pow(abs(k-1.0),2.0*H) );
    }
#ifdef DEBUG
    printf("ffGn 3: N %d\tNfft %d\t NfftHalf %d\n", N,Nfft, NfftHalf);
#endif
    if ( N==0 || Nfft == 0) return 0; 
    iZmag = v_get(Nfft); //v_zero(iZmag);
    fft(Zmag,iZmag);
    for(i=0;i<Nfft-1;i++){
      	printf("Zmag %g\t",Zmag->ve[i]);
      if ( Zmag->ve[i] < 0 ) {
	hoc_execerror("The fast Fourier transform of the circulant covariance had negative values.",0); 
	return 0;
      }
      Zmag->ve[i] = sqrt(Zmag->ve[i]);
    }
      //  Store N and H values in persistent variables for use during subsequent calls to this function.
    Nlast = N;
    Hlast = H;
      //}
#ifdef DEBUG
    printf("ffGn 4a: N %d\tNfft %d, NfftHalf %d\n",N, Nfft, NfftHalf);
#endif
    if ( (N == 0) || (Nfft == 0)) return 0; 
    Z_real = v_get(Nfft); 
    v_zero(Z_real);
    Z_im = v_get(Nfft); 
    v_zero(Z_im);
    
    for(i=0;i<Nfft;i++) {
      /* Z_real->ve[i] = Zmag->ve[i]*randn();  */
      /* Z_im->ve[i] = Zmag->ve[i]*randn();	 */
       Z_real->ve[i] = Zmag->ve[i];  
       Z_im->ve[i] = Zmag->ve[i];	 
    }

#ifdef DEBUG
    for(i=0;i<Nfft;i++) {
      printf("Z[%d]   %g\ti%g\t Zmag %g\n",i,Z_real->ve[i],Z_im->ve[i],Zmag->ve[i]);
    }
#endif

    ifft(Z_real,Z_im);

#ifdef DEBUG
    printf("After ifft, %g",sqrt(Nfft));
    for(i=0;i<Nfft;i++) {
      printf("Z[%d]   %g\ti%g\t Zmag %g\n",i,Z_real->ve[i],Z_im->ve[i],Zmag->ve[i]);
    }

#endif


    y = (double*)makevector(N);//fft);

    for(i=0;i<N;i++){//fft-1;i++) {
    y[i] =  (Z_real->ve[i]) * sqrt(Nfft);
    }

#ifdef DEBUG
    for(i=0;i<Nfft;i++) {
      printf("Z[%d] %g\t y %g\n",i,Z_real->ve[i],y[i]);
    }
    printf("ffGn 6:  N %d\tNfft %d, NfftHalf %d\n", N,Nfft, NfftHalf);
#endif

    V_FREE(Z_im);
    V_FREE(Z_real);
    V_FREE(Zmag);
    V_FREE(iZmag);	

#ifdef DEBUG
    printf("ffGn 7:  N %d\tNfft %d, NfftHalf %d\n",N, Nfft, NfftHalf);
#endif

   //	y((N+1):end) = [];
  }


  //  Convert the fGn to fBn, if necessary.
  if (fBn){ for (i=1;i<N;i++) y[i] = y[i]+y[i-1]; }
  
#ifdef DEBUG
  printf("ffGn 8: Resampling back to original (1/tdres): match with the AN model  N %d\t resamp %d\t nop %d\t Nfft %d\n",N,resamp,nop,Nfft);
#endif
  if ( (N == 0) || (Nfft == 0)) return 0; 

  if(resamp == 1){
    for (i=0;i<__min(nop,Nfft);i++)
      yffGn[i] = y[i];
    freevector(y);
    return __min(nop,Nfft);
  }else {
    ytmp = makevector((int)(N*resamp));
    printf("ffGn 9: calling resample(%x,%x,%d,%d)\n",&y,&ytmp,N,resamp);

    resample(y,ytmp,N,resamp);  //  Resampling to match with the AN model, 
  //N= size of y
  }


//  define standard deviation

if (sigma <= 0){
  if (mu<0.5){
    sigma = 5;  
  } else {
    if (mu<18){
      sigma = 50;   //  7 when added after powerlaw
    } else {
      sigma = 200;  //  40 when added after powerlaw        
    }
  }
 }

    int ystd=0;

    /*    //correction to small values, assume mean y = 0
    for (i=0;i<nop;i++) ystd = pow(yffGn[i],2.0);
    ystd = sqrt(ystd/nop);
    */
    printf("ffGn 9: ytmp size %d  nop %d\n",(int)(N*resamp),nop);
    for (i=0;i<nop;i++){ 
      if (ytmp[i] == 0) ystd++;
      yffGn[i] = ytmp[i]*sigma;
    }
    freevector(y);
    freevector(ytmp);


#ifdef DEBUG
    printf("ffGn: done! ystd %d\t %x\n",ystd, &yffGn);
#endif

    return nop;
}

