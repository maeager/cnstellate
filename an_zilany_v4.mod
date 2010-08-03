NEURON {
   SUFFIX nothing
}
VERBATIM
extern int vector_instance_px();
extern void vector_resize();
extern double* vector_vec();
extern void* vector_arg();
ENDVERBATIM

:* anmodel Vector Method - sout.an_zilany_v4( stim , tdres,  cf, spontrate, model, species, ifspike)
VERBATIM
//#include "complex.c"  /* uncomment if an_zbcatmodel is not present in folder */

#define DEBUG

#include "zilanycarneyv4.c"
#include "resample.c"
#include "ffGn.c"


//sout.an_zilay_v4(stim, tdres,cf,spont,model, species,ifspike,wavefile)

static double an_zilany_v4(void *vv)
{

   double *stim;      //Input stimulus vector in pascals
   double *ihcout,*sout;   //Output vector containing inst. rate for channel
   double tdres,cf,spont,out;
   double cohc,cihc;
   int species,fibertype,implnt;
   int ifspike;
   int nstim, nsout,nihcout;
   int nrep;
   cohc = 1.0;
   cihc = 1.0;

//Get Instance Sout vector for synapse data
   nstim = vector_instance_px(vv, &sout);
   nsout = vector_arg_px(1, &stim);

//Get Input arguments
/*   if(ifarg(9)!=1){  //Must be changed if more input arguments added
      hoc_execerror("ERROR: input syntax must be sout.an_zilany_v4( stim, tdres,  cf, fibertype,implnt,cihc,cohc, species,nrep,ihcout[optional])", 0);
      return 0;
   }
*/
 //TDRES  resolution of stim vector
   //Bruce model uses seconds rather than msec
   tdres = (double)(*getarg(2));
   if (tdres > 0.01e-3 || tdres < 0.002e-3){
      //printf("Note: ZilanyBruceV4 resolution should be between 0.01ms (Fs = 100kHz) and 0.002ms (500kHz) for normal usage.\n");
      //tdres = 0.002e-3;
   }
   //CF of fiber
   cf = (double)(*getarg(3));

   if ((cf<80)||(cf>70e3))
   {
      hoc_execerror("ERROR: cf must be between 80 Hz and 70 kHz\n",0);
      return 0;
   }
   //Spontaneous rate of ANF
   fibertype = (long)(*getarg(4));
   if ((fibertype<=0)||(fibertype>3))
   {
      printf("an_zilany_v4: fibertype  must be 1, 2, or 3 \n Default fibertype = 3 (HSR)");
      fibertype=3;
   }
   //New variable in version 4
   implnt = (long)(*getarg(5));
   if ((implnt<0)||(implnt>1))
   {
      printf("an_zbcatmodel: implnt  must be 0 (actual) or 1 (approx)\n");
      implnt=1;
   }
   
   //Model
   cihc = (double)(*getarg(6));
   if ((cihc<0)||(cihc>1))
   {
      printf("an_zilany_v4: cihc  must be between 0 and 1\n");
      cihc=1;
   }
   cohc = (double)(*getarg(7));
   if ((cohc<0)||(cohc>1))
   {
      printf("an_zilany_v4: cohc  must be between 0 and 1\n");
      cohc=1;
   }

   //Species
   species = (int)(*getarg(8));
   if (species != 1){
     hoc_execerror("an_zilany_v4: species other than cat (1) are not implemented",0);
     return 0;
   }
   //Reps
   nrep = (int)(*getarg(9));
   
if(ifarg(10)) {
  nihcout = vector_arg_px(10, &ihcout);
  if(nihcout != nstim){
    printf("ihcout must be the same size as stim, sout\n"); return 0;}
 } else {
  ihcout = makevector(nstim);
 }

   printf("AN model: Zilany, Carney, Bruce, Nelson and others  (version 4 c2010)\n");
   printf("IHCAN(stim,%.0f,%d,%g,%d,%g,%g,ihcout,%d)\n", cf, nrep, tdres, nstim, cohc, cihc, species);
   IHCAN(stim, cf, nrep, tdres, nstim, cohc, cihc, ihcout,species);
   printf("SingleAN(ihcout,%.0f,%d,%g,%d,%d,%d,sout,%d)\n",cf,nrep,tdres,nstim,fibertype,implnt,species);
   out= SingleAN_v4(ihcout,cf,nrep,tdres,nstim,fibertype,implnt,sout,species);

   if( !ifarg(10) ) freevector(ihcout);
   return out; 
}

static double an_zilany_v4(void *vv)
{

   double *stim;      //Input stimulus vector in pascals
   double *ihcout,*sout;   //Output vector containing inst. rate for channel
   double tdres,cf,spont,out;
   double cohc,cihc;
   int species,implnt;
   int ifspike;
   int nstim, nsout,nihcout;
   int nrep;
   cohc = 1.0;
   cihc = 1.0;

//Get Instance Sout vector for synapse data
   nstim = vector_instance_px(vv, &sout);
   nsout = vector_arg_px(1, &stim);

//Get Input arguments
/*   if(ifarg(9)!=1){  //Must be changed if more input arguments added
      hoc_execerror("ERROR: input syntax must be sout.an_zilany_v4( stim, tdres,  cf, fibertype,implnt,cihc,cohc, species,nrep,ihcout[optional])", 0);
      return 0;
   }
*/
 //TDRES  resolution of stim vector
   //Bruce model uses seconds rather than msec
   tdres = (double)(*getarg(2));
   if (tdres > 0.01e-3 || tdres < 0.002e-3){
      //printf("Note: ZilanyBruceV4 resolution should be between 0.01ms (Fs = 100kHz) and 0.002ms (500kHz) for normal usage.\n");
      //tdres = 0.002e-3;
   }
   //CF of fiber
   cf = (double)(*getarg(3));

   if ((cf<80)||(cf>70e3))
   {
      hoc_execerror("ERROR: cf must be between 80 Hz and 70 kHz\n",0);
      return 0;
   }
   //Spontaneous rate of ANF
   spont = (double)(*getarg(4));
   
   //implnt new variable in version 4
   implnt = (long)(*getarg(5));
   if ((implnt<0)||(implnt>1))
   {
      printf("an_zbcatmodel: implnt  must be 0 (actual) or 1 (approx)\n");
      implnt=1;
   }
   
   //Model
   cihc = (double)(*getarg(6));
   if ((cihc<0)||(cihc>1))
   {
      printf("an_zilany_v4: cihc  must be between 0 and 1\n");
      cihc=1;
   }
   cohc = (double)(*getarg(7));
   if ((cohc<0)||(cohc>1))
   {
      printf("an_zilany_v4: cohc  must be between 0 and 1\n");
      cohc=1;
   }

   //Species
   species = (int)(*getarg(8));
   if (species != 1){
     hoc_execerror("an_zilany_v4: species other than cat (1) are not implemented",0);
     return 0;
   }
   //Reps
   nrep = (int)(*getarg(9));
   
if(ifarg(10)) {
  nihcout = vector_arg_px(10, &ihcout);
  if(nihcout != nstim){
    printf("ihcout must be the same size as stim, sout\n"); return 0;}
 } else {
  ihcout = makevector(nstim);
 }

   printf("AN model: Zilany, Carney, Bruce, Nelson and others  (version 4 c2010)\n");
   printf("IHCAN(stim,%.0f,%d,%g,%d,%g,%g,ihcout,%d)\n", cf, nrep, tdres, nstim, cohc, cihc, species);
   IHCAN(stim, cf, nrep, tdres, nstim, cohc, cihc, ihcout,species);
   printf("SingleAN(ihcout,%.0f,%d,%g,%d,%d,%d,sout,%d)\n",cf,nrep,tdres,nstim,spont,implnt,species);
   out= SingleAN_v4_1(ihcout,cf,nrep,tdres,nstim,spont,implnt,sout,species);

   if( !ifarg(10) ) freevector(ihcout);
   return out; 
}

static double ihc_zilany_v4(void *vv)
{
   double *stim;      //Input stimulus vector in pascals
   double *ihcout;   //Output vector containing inst. rate for channel
double tdres,cf;
   double cohc,cihc;
   int species;
   int nstim, nihcout;
   int nrep;
   cohc = 1.0;
   cihc = 1.0;

//Get Instance Sout vector for synapse data
   nihcout = vector_instance_px(vv, &ihcout);
   nstim = vector_arg_px(1, &stim);

//Get Input arguments
   if(ifarg(7)!=1){  //Must be changed if more input arguments added
      hoc_execerror("ERROR: input syntax must be ihcout.ihc_zilany_v4( stim, tdres,  cf,cihc,cohc, species,nrep)", 0);
      return 0;
   }
   //TDRES  resolution of stim vector
   //Bruce model uses seconds rather than msec
   tdres = (double)(*getarg(2));
   if (tdres > 0.01e-3 || tdres < 0.002e-3){
      //printf("Note: ZilanyBruceV4 resolution should be between 0.01ms (Fs = 100kHz) and 0.002ms (500kHz) for normal usage.\n");
      //tdres = 0.002e-3;
   }
   //CF of fiber
   cf = (double)(*getarg(3));

   if ((cf<80)||(cf>70e3))
   {
      hoc_execerror("ERROR: cf must be between 80 Hz and 70 kHz\n",0);
      return 0;
   }
   
   //Model
   cihc = (double)(*getarg(4));
   if ((cihc<0)||(cihc>1))
   {
      printf("an_zilany_v4: cihc  must be between 0 and 1\n");
      cihc=1;
   }
   cohc = (double)(*getarg(5));
   if ((cohc<0)||(cohc>1))
   {
      printf("an_zilany_v4: cohc  must be between 0 and 1\n");
      cohc=1;
   }

   //Species
   species = (int)(*getarg(6));
   if ( species != 1 || species != 9 ){
     printf("an_zilany_v4: species other than cat (1) are not correctly implemented");

   }
   //Reps
   nrep = (int)(*getarg(7));

   printf("AN model: Zilany, Bruce, Nelson and Carney  (version 4 c2010)\n");
   printf("IHCAN(stim,%.0f,%d,%g,%d,%g,%g,ihcout,%d)\n", cf, nrep, tdres, nstim, cohc, cihc,species);
   IHCAN(stim, cf, nrep, tdres, nstim, cohc, cihc, ihcout,species);
   return (double) nihcout; 
}

static double syn_zilany_v4(void *vv)
{

  double *sptime,*sptimes,*ihcout,*sout;   //Input and Output vector containing inst. rate for channel
void** xx;
 void *spks;
   double tdres,cf,out;
   double xspikes;
   int species,fibertype,implnt;
   int ifspike = 0;
   int nihcout,nspikes, nsout;
   int nrep,i;
   out=0;
   ifspike=0;
//Get Instance Sout vector for synapse data
   nsout = vector_instance_px(vv, &sout);
   nihcout = vector_arg_px(1, &ihcout);

//Get Input arguments
   if(ifarg(8)!=1  || ifarg(7)!=1){  //Must be changed if more input arguments added
      hoc_execerror("ERROR: input syntax must be sout.syn_zilany_v4( ihcout, tdres,  cf, fibertype,implnt,cihc,cohc, species,nrep)\n or  sout.syn_zilany_v4( ihcout, tdres,  cf, fibertype,implnt,cihc,cohc, species,nrep,spks)", 0);
      return 0;
   }
   //TDRES  resolution of stim vector
   //Bruce model uses seconds rather than msec
   tdres = (double)(*getarg(2));
   if (tdres > 0.01e-3 || tdres < 0.002e-3){
      //printf("Note: ZilanyBruceV4 resolution should be between 0.01ms (Fs = 100kHz) and 0.002ms (500kHz) for normal usage.\n");
      //tdres = 0.002e-3;
   }
   //CF of fiber
   cf = (double)(*getarg(3));

   if ((cf<80)||(cf>70e3))
   {
      hoc_execerror("ERROR: cf must be between 80 Hz and 70 kHz\n",0);
      return 0;
   }
   //Spontaneous rate of ANF
   fibertype = (int)(*getarg(4));
   if ((fibertype<=0)||(fibertype>3))
   {
      printf("an_zilany_v4: fibertype  must be 1, 2, or 3 \n Default fibertype = 3 (HSR)");
      fibertype=3;
   }
   //New variable in version 4
   implnt = (int)(*getarg(5));
   if ((implnt<0)||(implnt>1))
   {
      printf("an_zbcatmodel: implnt  must be 0 (actual) or 1 (approx)\n");
      implnt=1;
   }
   
   //Species
   species = (int)(*getarg(6));
   if (species != 1){
     hoc_execerror("an_zilany_v4: species other than cat (1) are not implemented",0);
     return 0;
   }
   //Reps
   nrep = (int)(*getarg(7));

   if(ifarg(8)) {
     ifspike = 1;
   
       xx = (void**)(&xspikes);
       *xx = (void*)0;
       if (ifarg(9)) {
	 *xx = vector_arg(8);
       }
       spks = *((void**)(&xspikes));
       if (spks) 
	 vector_resize(spks,0);

   }


   printf("AN model: Zilany, Carney, Bruce, Nelson and others  (version 4 c2010)\n");
   printf("SingleAN(ihcout,%.0f,%d,%g,%d,%d,%d,sout,%d)\n",cf,nrep,tdres,nihcout,fibertype,implnt,sout,species);
    out= SingleAN_v4(ihcout,cf,nrep,tdres,nihcout,fibertype,implnt,sout,species);


    /*======  Spike Generations ======*/
      if(ifspike){
	sptimes  = makevector((long) ceil(out / 0.00075));
	nspikes = SpikeGenerator_v4(sout, tdres, out, nrep, sptimes);
	spks = *((void**)(&xspikes));
	if (spks)  vector_resize(spks, nspikes);
     
	sptime = ((double*) vector_vec(spks));      //Get array ptr to spks Vector	
	for(i = 0; i < nspikes; i++)
	  {
	   sptime[i]  = sptimes[i];
	  };
	out=nspikes;
      	freevector(sptimes);     
      }     
      return out; 
}


/* The spike generator now uses a method coded up by B. Scott Jackson (bsj22@cornell.edu) 
   Scott's original code is available from Laurel Carney's web site at:
   http://www.urmc.rochester.edu/smd/Nanat/faculty-research/lab-pages/LaurelCarney/auditory-models.cfm
*/
/* int SpikeGenerator(double *synouttmp, double tdres, int totalstim, int nrep, double *sptime) 
*/
static double ANFSpikeGenerator3(void *vv)
{  
 
    double  c0,s0,c1,s1,dead;
    int nspikes,k,NoutMax,Nout,deadtimeIndex,randBufIndex;      
    double	deadtimeRnd, endOfLastDeadtime, refracMult0, refracMult1, refracValue0, refracValue1;
    double Xsum, unitRateIntrvl, countTime, DT;    
    double *randNums, *sout;
    double tdres,nsout;
    int totalstim, nrep; 
    double *sptime;
    void      *spks;  //We need to be able to resize this vector
    double   xspikes;

//Get Instance Sout vector for synapse data
   nsout = vector_instance_px(vv, &sout);
   void** xx;
   xx = (void**)(&xspikes);
   *xx = (void*)0;
   if (ifarg(1)) {
      *xx = vector_arg(1);
   }
   spks = *((void**)(&xspikes));
   if (spks) {
      vector_resize(spks, 0);
   }
//Get number of Reps and tdres

   nrep = (int) (*getarg(2));
   tdres = (double)(*getarg(3));


    c0      = 0.5;
    s0      = 0.001;
    c1      = 0.5;
    s1      = 0.0125;
    dead    = 0.00075;

    totalstim = nsout*tdres*1000;
    DT = totalstim * tdres * nrep;  /* Total duration of the rate function */
    Nout = 0;
    NoutMax = (long) ceil(DT/dead);    
   
    spks = *((void**)(&xspikes));
    if(spks){
     vector_resize(spks, NoutMax);
    }
    sptime = ((double*) vector_vec(spks));      //Get array ptr to spks Vector
   
    for (randBufIndex=0; randBufIndex< NoutMax+1;++randBufIndex) 
      randNums[randBufIndex]=scop_random();
    randBufIndex = 0;
    
    /* Calculate useful constants */
    deadtimeIndex = (long) floor(dead/tdres);  /* Integer number of discrete time bins within deadtime */
    deadtimeRnd = deadtimeIndex*tdres;		   /* Deadtime rounded down to length of an integer number of discrete time bins */

    refracMult0 = 1 - tdres/s0;  /* If y0(t) = c0*exp(-t/s0), then y0(t+tdres) = y0(t)*refracMult0 */
    refracMult1 = 1 - tdres/s1;  /* If y1(t) = c1*exp(-t/s1), then y1(t+tdres) = y1(t)*refracMult1 */

    /* Calculate effects of a random spike before t=0 on refractoriness and the time-warping sum at t=0 */
    endOfLastDeadtime = log(randNums[randBufIndex++]) / (sout[0] + dead);  /* End of last deadtime before t=0 */
    if (endOfLastDeadtime< 0) endOfLastDeadtime=0;
    refracValue0 = c0*exp(endOfLastDeadtime/s0);     /* Value of first exponential in refractory function */
    refracValue1 = c1*exp(endOfLastDeadtime/s1);     /* Value of second exponential in refractory function */
    Xsum = sout[0] * (-endOfLastDeadtime + c0*s0*(exp(endOfLastDeadtime/s0)-1) + c1*s1*(exp(endOfLastDeadtime/s1)-1));  
    /* Value of time-warping sum */
    /*  ^^^^ This is the "integral" of the refractory function ^^^^ (normalized by 'tdres') */
    
    /* Calculate first interspike interval in a homogeneous, unit-rate Poisson process (normalized by 'tdres') */
    unitRateIntrvl = -log(randNums[randBufIndex++])/tdres;  
	    /* NOTE: Both 'unitRateInterval' and 'Xsum' are divided (or normalized) by 'tdres' in order to reduce calculation time.  
		This way we only need to divide by 'tdres' once per spike (when calculating 'unitRateInterval'), instead of 
		multiplying by 'tdres' once per time bin (when calculating the new value of 'Xsum').                         */

    countTime = tdres;
    for (k=0; (k<totalstim*nrep) && (countTime<DT); ++k, countTime+=tdres, refracValue0*=refracMult0, refracValue1*=refracMult1)  /* Loop through rate vector */
      {
	if (sout[k]>0)  /* Nothing to do for non-positive rates, i.e. Xsum += 0 for non-positive rates. */
	  {
	    Xsum += sout[k]*(1 - refracValue0 - refracValue1);  /* Add synout*(refractory value) to time-warping sum */
	    
	    if ( Xsum >= unitRateIntrvl )  /* Spike occurs when time-warping sum exceeds interspike "time" in unit-rate process */
	      {
		sptime[Nout] = countTime; Nout = Nout+1;								
		unitRateIntrvl = -log(randNums[randBufIndex++]) /tdres; 
		Xsum = 0;
				
		/* Increase index and time to the last time bin in the deadtime, and reset (relative) refractory function */
		k += deadtimeIndex;
		countTime += deadtimeRnd;
		refracValue0 = c0;
		refracValue1 = c1;
	      }
	  }
      } /* End of rate vector loop */			
    
    vector_resize(spks, Nout);
    nspikes = Nout;  /* Number of spikes that occurred. */
    return(nspikes);
}

static double fast_fGn(void *vv)
{

  void *ptr_ffGn;
   double *vec_ffGn;   //Output vector containing inst. rate for channel
   double tdres,cf,N,Hinput,mu,sigma;
   double xspikes;
   int nout,nsizemax;
   int iarg=1;
   
//Get Instance/out vector 
   nout = vector_instance_px(vv, &vec_ffGn);

   void** xx;
   xx = (void**)(&xspikes);
   *xx = (void*)0;
   if (ifarg(1)) {
      *xx = vector_arg(iarg++);
   }
   ptr_ffGn = *((void**)(&xspikes));
   if (ptr_ffGn) {
      vector_resize(ptr_ffGn, 0);
   }
 

//Get Input arguments
   if( !ifarg(5) && !ifarg(6)){  //Must be changed if more input arguments added
      hoc_execerror("fast_fGn: input syntax must be vec.fast_fGn(N, tdres, Hinput, mu, sigma (optional))", 0);
      return 0;
   }
   if(ifarg(iarg))   N = (double)(*getarg(iarg++));
   nsizemax = (int) N;

   //Bruce model uses seconds rather than msec
      if(ifarg(iarg)) tdres = (double)(*getarg(iarg++));
/*   if (tdres > 0 || tdres < 1){
      printf("Note: ZilanyBruceV4 resolution should be between 0.01ms (Fs = 100kHz) and 0.002ms (500kHz) for normal usage.\n");
      tdres = 0.01e-3;
   }
*/

   if(ifarg(iarg))   Hinput = (double)(*getarg(iarg++));

   if ( ( Hinput <= 0 ) || (Hinput > 2) ){
      hoc_execerror("fast_fGn ERROR: Hinput must be between 0  and 2\n",0);
      return 0;
   }

   if(ifarg(iarg))   mu = (double)(*getarg(iarg++));
   if ( mu < 0 ) {
      printf("mu must be >= 0\n");
      mu=0;
   }

   if(ifarg(iarg)){
     sigma = (double)(*getarg(iarg++));
   }else{
     sigma = -1;
       }

    ptr_ffGn = *((void**)(&xspikes));
    if(ptr_ffGn){
      vector_resize(ptr_ffGn, nsizemax);
      //      vector_resize(ptr_ffGn, 32000);
    }
    vec_ffGn = ((double*) vector_vec(ptr_ffGn));      //Get array ptr to ptr_ffGn Vector
    printf("fast_fGn: calling ffGn(&%x,%d,%g,%g,%g,%g)\n",&vec_ffGn,nsizemax,tdres,Hinput,mu,sigma);
    return ffGn(vec_ffGn,nsizemax,tdres,Hinput,mu,sigma);

}

static double rtresample(void *vv)
{
  void *ptr_resample;
  double *vec1_resample;
   double *vec2_resample;   //Output vector containing inst. rate for channel
   double tdres,cf,N,Hinput,mu,sigma;
   double xresample,resamp;
   int nin,NSizeMax;
   int iarg=1;
   
//Get Instance/out vector 
   nin = vector_instance_px(vv, &vec1_resample);

   void** xx;
   xx = (void**)(&xresample);
   *xx = (void*)0;
   if (ifarg(1)) {
      *xx = vector_arg(iarg++);
   }
   ptr_resample = *((void**)(&xresample));
   if (ptr_resample) {
      vector_resize(ptr_resample, 0);
   }
 
//Get Input arguments
   if(!ifarg(2)){  //Must be changed if more input arguments added
      hoc_execerror("rtresample: input syntax must be vec.rtresample(destvec, resample)", 0);
      return 0;
   }
   if(ifarg(iarg))   resamp = (double)(*getarg(iarg++));
   NSizeMax = (int) (resamp * (double)nin);
   if(NSizeMax <=0 ){  //Must be changed if more input arguments added
      hoc_execerror("rtresample: size not big enough)", 0);
      return 0;
   }

    ptr_resample = *((void**)(&xresample));
    if(ptr_resample){
     vector_resize(ptr_resample, NSizeMax);
    }

    vec2_resample = ((double*) vector_vec(ptr_resample));      //Get array ptr to ptr_resample Vector

    return resample(vec1_resample,vec2_resample,nin,resamp);

}
ENDVERBATIM


PROCEDURE install_an_zilany_v4()
{
VERBATIM
   install_vector_method("an_zilany_v4", an_zilany_v4);
   install_vector_method("ihc_zilany_v4", ihc_zilany_v4);
   install_vector_method("syn_zilany_v4", syn_zilany_v4);
   install_vector_method("ANFSpikeGenerator3", ANFSpikeGenerator3);
   install_vector_method("rtresample", rtresample);
   install_vector_method("fast_fGn", fast_fGn);
ENDVERBATIM
}

   
