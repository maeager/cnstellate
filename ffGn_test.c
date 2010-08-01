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
*/

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
*/

void realft(double data[], unsigned long n, int isign)
{
	unsigned long i,i1,i2,i3,i4,np3;
	double c1=0.5,c2,h1r,h1i,h2r,h2i;
	double wr,wi,wpr,wpi,wtemp,theta;
	//	printf("realft: (%x,%x,%d)\n",&data[0],&n,n);
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


int ivoc_fft(double *v1, double *v2,int v1size,int v2size,int inv) 
{
  int n = 1;
  while(n < v2size) n*=2;
  printf("\t\tivoc_fft: recv (%x,%x,%d,%d)\n",&v1[0],&v2[0],v1size,v2size);
  double *data;
  data = makevector(n); zero_vector(data,n);
  int i;
  for (i=0;i<v2size;++i) data[i] = v2[i];
  printf("\t\tivoc_fft: data %x\n", &data[0]);
  realft(&data[0]-1,n,inv);

  if (v1size != n) {printf("\t\tivoc_fft: v1 must be %d",n);}
  for (i=0;i<n;++i) v1[i]=data[i];
  printf("\t\tivoc_fft: done n %d\n", n);
  freevector(data);
 return n;
}



int FFT(int inv,double *v1,double *v2,double *v3,int Nsize)
{
  int n, x,i,j;
  if (inv == 1) { // forward  
    printf("\tFFT: forward (%d,%x,%x,%x,%d)\n",inv,&v1[0],&v2[0],&v3[0],Nsize);
    n =ivoc_fft(v2,v1,Nsize,Nsize,1);
    for(x=0;x<n;x++) v2[x] = v2[x]/((double)n/2.0);
    v2[0] = v2[0] / 2.0;	// makes the spectrum appear discontinuous
    v2[1] = v2[1]/2.0;	// but the amplitudes are intuitive

    for (i=1, j=0; i<=Nsize-1; i+=2, j+=1) v3[j] = v2[i];    //v3.copy(v2, 0, 1, -1, 1, 2)   // odd elements
    for (i=0, j=0; i<=Nsize-1; i+=2, j+=1) v2[j] = v2[i];    //v2.copy(v2, 0, 0, -1, 1, 2)   // even elements
    //v2.resize(n/2+1);
    //v3.resize(n/2+1);
    v2[(int)(n/2)] = v3[0];   //highest cos started in o3.x[1
    v3[0] = v3[n/2] = 0;       // weights for sin(0*i)and sin(PI*i)
    printf("\tFFT: done (%d,%x,%x,%x,%d)\n",inv,&v1[0],&v2[0],&v3[0],Nsize);
    return (int)(n/2)+1;

  }else{ // inverse
    printf("\tFFT: inverse (%d,%x,%x,%x,%d)\n",inv,&v1[0],&v2[0],&v3[0],Nsize);
    // shuffle o3 and o4 into o2
    //n = v2.size()
    n = Nsize/2 + 1 ;
    for (i=0; i<Nsize; i+=1) v1[i] = v3[i]; //v1=v3
    for (i=0, j=0; i<=n-2; i+=1, j+=2) v3[j] = v2[i]; //v1.copy(v2, 0, 0, n-2, 2, 1)
    v3[1] = v2[n-1];
    for (i=1, j=3; i<=n-2; i+=1, j+=2) v3[j] = v1[i]; //v1.copy(v3, 3, 1, n-2, 2, 1)
    v3[0] *= 2;
    v3[1] *= 2;     
    for (i=0; i<Nsize; i+=1) printf("%g\t",v3[i]); //v1=v3
    n = ivoc_fft(v1,v3,Nsize,n,-1);
    printf("\n\tFFT: done (%d,%x,%x,%x,%d)\n",inv,&v1[0],&v2[0],&v3[0],Nsize);
    return n;
  }
}

int rFFT(int inv,double *v1,double *v2,double *v3,int Nsize)
{
  int n, x,i,j;
  if (inv == 1) { // forward  
    printf("\tFFT: forward (%d,%x,%x,%x,%d)\n",inv,&v1[0],&v2[0],&v3[0],Nsize);
    n =ivoc_fft(v2,v1,Nsize,Nsize,1);
    //for(x=0;x<n;x++) v2[x] = v2[x]/((double)n/2.0);
    //v2[0] = v2[0] / 2.0;	// makes the spectrum appear discontinuous
    //v2[1] = v2[1]/2.0;	// but the amplitudes are intuitive

    for (i=1, j=0; i<=Nsize-1; i+=2, j+=1) v3[j] = v2[i];    //v3.copy(v2, 0, 1, -1, 1, 2)   // odd elements
    for (i=0, j=0; i<=Nsize-1; i+=2, j+=1) v2[j] = v2[i];    //v2.copy(v2, 0, 0, -1, 1, 2)   // even elements
    //v2.resize(n/2+1);
    //v3.resize(n/2+1);

    v2[(int)(n/2)] = v3[0];   //highest cos started in o3.x[1

    for (i=n/2 -1,j=n/2 + 1; j<=Nsize-1; i+=-1, j+=1) v2[j] = v2[i];

    v3[0] = v3[n/2] = 0;       // weights for sin(0*i)and sin(PI*i)
    printf("\tFFT: done (%d,%x,%x,%x,%d)\n",inv,&v1[0],&v2[0],&v3[0],Nsize);
    return (int)(n/2)+1;

  }else{ // inverse
    printf("\tFFT: inverse (%d,%x,%x,%x,%d)\n",inv,&v1[0],&v2[0],&v3[0],Nsize);
    // shuffle o3 and o4 into o2
    //n = v2.size()
    n = Nsize/2 + 1 ;
    for (i=0; i<=n-1; i+=1) v1[i] = v3[i]; //v1=v3
    for (i=0, j=0; i<=n-2; i+=1, j+=2) v3[j] = v2[i]; //v1.copy(v2, 0, 0, n-2, 2, 1)
    v3[1] = v2[n-1];
    for (i=1, j=3; i<=n-1; i+=1, j+=2) v3[j] = v1[i]; //v1.copy(v3, 3, 1, n-2, 2, 1)
    v3[0] *= 2;
    v3[1] *= 2; 

    n = ivoc_fft(v1,v3,Nsize,n, -1);
    printf("\tFFT: done (%d,%x,%x,%x,%d)\n",inv,&v1[0],&v2[0],&v3[0],Nsize);
    return n;
  }
}



/* 
 static double *Zmag=NULL;
 static int Nlast=0;
 static int Hlast=0;
 static int Nfft=0;
*/

/* scoplib.h random functions allows control of default random functions in NEURON*/
double randn(){ 
  //  double nr=scop_random();//normrand(0,1);
  //printf("%g\n",nr);
  //	 return nr;
	 return normrand(0.,1.);
}


double ffGn(double *yffGn,int N, double tdres, double Hinput, double mu, double sigma)
{
  double *Zmag=NULL;
  int Nlast=0;
  int Hlast=0;
  int Nfft=0;
 int i,n,nsize,nop,resamp,k;	
 double H,fBn,NfftHalf;
 double *y,*ytmp;
 double *Ztemp=NULL,*Z_real=NULL,*Z_im=NULL;
  /*---- Check input arguments ---------- */
 Nfft=N;
 
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



  //  Downsampling No. of points to match with those of Scott jackson (tau 1e-1)
  resamp = ceil(1e-1/tdres);
  nop = N; 
  N = ceil(N/resamp)+1; 
  if (N<10)
    N = 10;

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

    //    if (Nfft==0 || Nlast ==0 || Hlast==0 || N != Nlast || H != Hlast || Zmag==NULL){
      //  The persistent variables must be (re-)calculated.
    Nfft = pow(2,ceil(log10(2*(N-1))/log10(2)));
    NfftHalf = round(Nfft/2);
    printf("Nfft  %d  NfftHalf %d\n", Nfft,NfftHalf);	

      //k=[0:NfftHalf, (NfftHalf-1):-1:1];

      Zmag = makevector(Nfft);
      Z_im = makevector(Nfft);
      Z_real = makevector(Nfft);
      
      for (i =0;i<Nfft;i++){
	if(i<=NfftHalf) { k=i;}else{ k=2*NfftHalf-i;}
	Zmag[i]=0.5 * ( pow(k+1.0,2.0*H) - 2.0*pow(k,2.0*H) + pow(fabs(k-1),2.0*H) );
#ifdef DEBUG
	printf("\tk %d Z %g\n",k, Zmag[i]);
#endif      
      }

      n=rFFT(1,Zmag,Z_real,Z_im,Nfft);
      

#ifdef DEBUG
      for (i =0;i<n;i++)
	printf("\tZmag %g Zreal %g Zim %g\n",Zmag[i],Z_real[i],Z_im[i]);
#endif   
      zero_vector(Z_im,Nfft); 
      zero_vector(Zmag,Nfft); 


      /*      for(i=0;i<n;i++){ */
      for(i=0;i<Nfft;i++){
	if ( Z_real[i] < 0 || isnan(Z_real[i]) ) {
	  printf("error %d\n",i);
	  printf("The fast Fourier transform of the circulant covariance had negative values.");
	  freevector(Zmag);
	  freevector(Z_real);
	  freevector(Z_im);
	  return 0;
	}
	Zmag[i] = sqrt(Z_real[i]);

      }

#ifdef DEBUG
      for (i =0;i<n;i++)
	printf("\tZmag %g Zreal %g Zim %g\n",Zmag[i],Z_real[i],Z_im[i]);
#endif   

      //  Store N and H values in persistent variables for use during subsequent calls to this function.
      Nlast = N;
      Hlast = H;
      //}
      zero_vector(Z_real,Nfft);
      zero_vector(Z_im,Nfft);
/* #ifdef DEBUG */
/*       for (i =0;i<n;i++) */
/* 	printf("\tZmag %g Zreal %g Zim %g\n",Zmag[i],Z_real[i],Z_im[i]); */
/* #endif    */

     for(i=0;i<Nfft;i++) { 
/*       //for(i=0;i<n;i++){ */
       /* Z_real[i] = randn() * Zmag[i];   */
       /* Z_im[i] =  randn() * Zmag[i];	  */
#ifdef DEBUG 
       Z_real[i] = Zmag[i];
       Z_im[i] = -Zmag[i];
       printf("Z[%d]  %g + i %g\t Zmag %g\n",i,Z_real[i],Z_im[i],Zmag[i]); 
#endif 
     }

     Ztemp = makevector(Nfft*2); 
     zero_vector(Ztemp,Nfft*2);
    // Inverse FFT
//    rFFT(-1,Ztemp,Z_real,Z_im, Nfft);      
     /*TEST*/    FFT(1,Z_real,Ztemp,Z_im,Nfft);      
    y = makevector(NfftHalf); 
    for(i=0;i<NfftHalf;i++)
      y[i] =  (Ztemp[i]) * sqrt(Nfft)  / 2.0;
    y[0]*=2.0;

#ifdef DEBUG
    for(i=0;i<NfftHalf;i++) printf("Ztemp %g\n",Ztemp[i]);
    for(i=0;i<NfftHalf;i++) printf("y %.15g\n",y[i]);
#endif    

     zero_vector(Ztemp,Nfft*2);
     zero_vector(Z_im,Nfft);for (i = 0; i < Nfft; ++i)  Z_im[i] = -Zmag[i]/(i+1);
     FFT(1,Z_real,Ztemp,Z_im,Nfft);      
    for(i=0;i<(Nfft*2);i++) printf("FFT(1) %g %g %g\n",Ztemp[i],Z_real[i],Z_im[i]);

     zero_vector(Ztemp,Nfft*2);
     zero_vector(Z_im,Nfft);for (i = 0; i < Nfft; ++i){Z_real[i] = Zmag[i];  Z_im[i] = -Zmag[i]/(i+1);}
     rFFT(-1,Ztemp,Z_real,Z_im,Nfft);      
    for(i=0;i<(Nfft*2);i++) printf("rFFT(-1) %g %g %g\n",Ztemp[i],Z_real[i],Z_im[i]);

    y = makevector(Nfft); 
    for(i=0;i<Nfft;i++)
      y[i] =  (Ztemp[i]) * sqrt(Nfft)  / 2.0;
    y[0]*=2.0;
    for(i=0;i<Nfft;i++) printf("y %.15g\n",y[i]);


   //	y((N+1):end) = [];
  }


  //  Convert the fGn to fBn, if necessary.
  // if (fBn){ for (i=1;i<N;i++) y[i] = y[i]+y[i-1]; }
  ytmp = makevector((int)(NfftHalf*resamp));  
  printf("ffGn: Resampling %d\t%d\t%d\n",resamp,N,nop);
  resample(y,ytmp,NfftHalf,resamp);  //  Resampling to match with the AN model, 
  //N= size of y

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

  
    //  for (i=0;i<nop;i++)
      for (i=0;i<NfftHalf*resamp;i++)
	yffGn[i] = ytmp[i]*sigma;
#ifdef DEBUG
  printf("\tffGn: Freeing vectors");
#endif    



  if(y) freevector(y); 
  if(ytmp) freevector(ytmp);

    if(Ztemp) freevector(Ztemp);
    if(Z_im) freevector(Z_im);
    if(Z_real) freevector(Z_real);
    if(Zmag) freevector(Zmag);
    
  return nop;
}
