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



/* FOURIER functions from src/ivoc/fourier.cpp */
/* 
  four1()  -- complex discrete FFT and inverse FFT
  N.R.C  p. 411


#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void four1(double data[], unsigned long nn, int isign)
{
	unsigned long n,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta;
	double tempr,tempi;

	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) {
		if (j > i) {
			SWAP(data[j],data[i]);
			SWAP(data[j+1],data[i+1]);
		}
		m=n >> 1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax=2;
	while (n > mmax) {
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) {
			for (i=m;i<=n;i+=istep) {
				j=i+mmax;
				tempr=wr*data[j]-wi*data[j+1];
				tempi=wr*data[j+1]+wi*data[j];
				data[j]=data[i]-tempr;
				data[j+1]=data[i+1]-tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}
#undef SWAP
*/
/* 
  twofft()  -- discrete FFT of two real functions simultaneously
  N.R.C  p. 414


void twofft(double data1[], double data2[], double fft1[], double fft2[],
	unsigned long n)
{
	void four1(double data[], unsigned long nn, int isign);
	unsigned long nn3,nn2,jj,j;
	double rep,rem,aip,aim;

	nn3=1+(nn2=2+n+n);
	for (j=1,jj=2;j<=n;j++,jj+=2) {
		fft1[jj-1]=data1[j];
		fft1[jj]=data2[j];
	}
	four1(fft1,n,1);
	fft2[1]=fft1[2];
	fft1[2]=fft2[2]=0.0;
	for (j=3;j<=n+1;j+=2) {
		rep=0.5*(fft1[j]+fft1[nn2-j]);
		rem=0.5*(fft1[j]-fft1[nn2-j]);
		aip=0.5*(fft1[j+1]+fft1[nn3-j]);
		aim=0.5*(fft1[j+1]-fft1[nn3-j]);
		fft1[j]=rep;
		fft1[j+1]=aim;
		fft1[nn2-j]=rep;
		fft1[nn3-j] = -aim;
		fft2[j]=aip;
		fft2[j+1] = -rem;
		fft2[nn2-j]=aip;
		fft2[nn3-j]=rem;
	}
}
*/

/* 
  realft()  -- discrete FFT of a real function with 2n data pts
  N.R.C  p. 417


void realft(double data[], unsigned long n, int isign)
{
	void four1(double data[], unsigned long nn, int isign);
	unsigned long i,i1,i2,i3,i4,np3;
	double c1=0.5,c2,h1r,h1i,h2r,h2i;
	double wr,wi,wpr,wpi,wtemp,theta;

	theta=3.141592653589793/(double) (n>>1);
	if (isign == 1) {
		c2 = -0.5;
		four1(data,n>>1,1);
	} else {
		c2=0.5;
		theta = -theta;
	}
	wtemp=sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi=sin(theta);
	wr=1.0+wpr;
	wi=wpi;
	np3=n+3;
	for (i=2;i<=(n>>2);i++) {
		i4=1+(i3=np3-(i2=1+(i1=i+i-1)));
		h1r=c1*(data[i1]+data[i3]);
		h1i=c1*(data[i2]-data[i4]);
		h2r = -c2*(data[i2]+data[i4]);
		h2i=c2*(data[i1]-data[i3]);
		data[i1]=h1r+wr*h2r-wi*h2i;
		data[i2]=h1i+wr*h2i+wi*h2r;
		data[i3]=h1r-wr*h2r+wi*h2i;
		data[i4] = -h1i+wr*h2i+wi*h2r;
		wr=(wtemp=wr)*wpr-wi*wpi+wr;
		wi=wi*wpr+wtemp*wpi+wi;
	}
	if (isign == 1) {
		data[1] = (h1r=data[1])+data[2];
		data[2] = h1r-data[2];
	} else {
		data[1]=c1*((h1r=data[1])+data[2]);
		data[2]=c1*(h1r-data[2]);
		four1(data,n>>1,-1);
	}
}

*/


/*double * yffGn=0;
 double * Zmag=0;
 double Nfft=0.;
 double Nlast=0.;
 double Hlast=0.;
*/


/* fft from Meschach Library */
/* get/resize vector to given dimension */
#include "matrix.h"
#include "matrix2.h"


 /*static double * yffGn;
static double Nfft=0;
static double Nlast=0;
static double Hlast=0;
 */

/* scoplib.h random functions allows control of default random functions in NEURON*/
double randn(){ return normrand(0,1);}


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
 double *y;
 VEC *Z_real,*Z_im,*Zmag,*iZmag;
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
  resamp = (int) ceil(1e-1 / (tdres * 1000));
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
      if(Zmag) V_FREE(Zmag);
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
      // fft(Zmag,iZmag);
      for(i=0;i<Nfft-1;i++){
	//	printf("Zmag[%d]   %g\n",i,Zmag->ve[i]);
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
#ifdef DEBUG
      printf("ffGn 4b: N %d\tNfft %d, NfftHalf %d\n",N, Nfft, NfftHalf);
#endif
      if ( (N == 0) || (Nfft == 0)) return 0; 
      //      v_zero(Z_real);
#ifdef DEBUG
      printf("ffGn 4c: N %d\tNfft %d, NfftHalf %d\n",N, Nfft, NfftHalf);
#endif
      if ( (N == 0) || (Nfft == 0)) return 0; 
      Z_im = v_get(Nfft); 
 // v_zero(Z_im);
    
      for(i=0;i<Nfft;i++) {
	Z_real->ve[i] = Zmag->ve[i]*randn(); 
	Z_im->ve[i] = Zmag->ve[i]*randn();	
  //  printf("Z[%d]   %g\ti%g\n",i,Z_real->ve[i],Z_im->ve[i]);
      }

#ifdef DEBUG
      printf("ffGn 5a: N %d\tNfft %d, NfftHalf %d\n",N, Nfft, NfftHalf);
#endif
 if ( (N == 0) || (Nfft == 0)) return 0; 
      //    ifft(Z_real,Z_im);
#ifdef DEBUG
      printf("ffGn 5b:  N %d\tNfft %d, NfftHalf %d\n",N, Nfft, NfftHalf);
#endif
        if ( (N == 0) || (Nfft == 0)) return 0; 
	y = makevector(Nfft);
#ifdef DEBUG
	printf("ffGn 5c: N  %d Nfft %d, NfftHalf %d\n", N,Nfft, NfftHalf);
#endif
    if ( (N == 0) || (Nfft == 0)) return 0; 
	for(i=0;i<Nfft-1;i++) {
	  y[i] =  Z_real->ve[i] * sqrt(Nfft);
      //      printf("y[%d]\t %g\n",i,y[i]);
    }
#ifdef DEBUG
    printf("ffGn 6:  N %d\tNfft %d, NfftHalf %d\n", N,Nfft, NfftHalf);
#endif
        if ( (N == 0) || (Nfft == 0)) return 0; 
	V_FREE(Z_im);
	V_FREE(Z_real);
	V_FREE(Zmag);
	V_FREE(iZmag);	

#ifdef DEBUG
	printf("ffGn 7:  N %d\tNfft %d, NfftHalf %d\n",N, Nfft, NfftHalf);
#endif
	if ( (N == 0) || (Nfft == 0)) return 0; 
   //	y((N+1):end) = [];
  }


  //  Convert the fGn to fBn, if necessary.
  // if (fBn){ for (i=1;i<N;i++) y[i] = y[i]+y[i-1]; }
  
#ifdef DEBUG
  printf("ffGn 8: Resampling back to original (1/tdres): match with the AN model  N %d\t resamp %d\t nop %d\t Nfft %d\n",N,resamp,nop,Nfft);
#endif
        if ( (N == 0) || (Nfft == 0)) return 0; 

	if( resamp == 1){
	  for (i=0;i<__min(nop,Nfft);i++)
	    yffGn[i] = y[i];
	}else {
	  printf("ffGn 9: resample N %d  resamp %g", N, resamp);
	  resample(y,yffGn,N,resamp);  //  Resampling to match with the AN model, 
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

    double ystd=0;
    //correction to small values, assume mean y = 0
    for (i=0;i<nop;i++) ystd = pow(yffGn[i],2.0);
    ystd = sqrt(ystd/nop);
    for (i=0;i<nop;i++)
      yffGn[i] = yffGn[i]*sigma/ystd;

    freevector(y);
    return nop;
}

