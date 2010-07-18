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
#include "../complex.c"
#include "../zilanycarneyv4.c"
#include "ffGn2.c"

//sout.an_zilay_v4(stim, tdres,cf,spont,model, species,ifspike,wavefile)

static double an_zilany_v4(void *vv)
{

   double *stim;      //Input stimulus vector in pascals
   double *ihcout,*sout;   //Output vector containing inst. rate for channel
   double tdres,cf,out;
   double cohc,cihc;
   int species,fibertype,implnt;
   int ifspike;
   int nstim, nsout;
   int nrep;
   cohc = 1.0;
   cihc = 1.0;

//Get Instance Sout vector for synapse data
   nstim = vector_instance_px(vv, &sout);
   nsout = vector_arg_px(1, &stim);

//Get Input arguments
   if(ifarg(8)!=1){  //Must be changed if more input arguments added
      hoc_execerror("ERROR: input syntax must be sout.an_zilany_v4( stim, tdres,  cf, fibertype,implnt,cihc,cohc, species,nrep)", 0);
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

   ihcout = makevector(nstim);
   printf("AN model: Zilany, Carney, Bruce, Nelson and others  (version 4 c2010)\n");
   printf("IHCAN(stim,%.0f,%d,%g,%d,%g,%g,ihcout,%d)\n", cf, nrep, tdres, nstim, cohc, cihc, ihcout,species);
   IHCAN(stim, cf, nrep, tdres, nstim, cohc, cihc, ihcout,species);
   printf("SingleAN(ihcout,%.0f,%d,%g,%d,%d,%d,sout,%d)",cf,nrep,tdres,nstim,fibertype,implnt,sout,species);
   out= SingleAN(ihcout,cf,nrep,tdres,nstim,fibertype,implnt,sout,species);
   freevector(ihcout);
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

  void *spks;
   double *inst;   //Output vector containing inst. rate for channel
   double tdres,cf,N,Hinput,mu,sigma;
   double xspikes;
   int nout,nsizemax;
   int iarg=1;
   
//Get Instance/out vector 
   nout = vector_instance_px(vv, &inst);

   void** xx;
   xx = (void**)(&xspikes);
   *xx = (void*)0;
   if (ifarg(1)) {
      *xx = vector_arg(iarg++);
   }
   spks = *((void**)(&xspikes));
   if (spks) {
      vector_resize(spks, 0);
   }
 

//Get Input arguments
   if( !ifarg(5) && !ifarg(6)){  //Must be changed if more input arguments added
      hoc_execerror("fast_fGn: input syntax must be vec.fast_fGn(N, tdres, Hinput, mu, sigma (optional))", 0);
      return 0;
   }
   N = (double)(*getarg(iarg++));
   nsizemax = (int) N;

   //Bruce model uses seconds rather than msec
   tdres = (double)(*getarg(iarg++));
/*   if (tdres > 0 || tdres < 1){
      printf("Note: ZilanyBruceV4 resolution should be between 0.01ms (Fs = 100kHz) and 0.002ms (500kHz) for normal usage.\n");
      tdres = 0.01e-3;
   }
*/

   Hinput = (double)(*getarg(iarg++));

   if ( ( Hinput <= 0 ) || (Hinput > 2) ){
      hoc_execerror("fast_fGn ERROR: Hinput must be between 0  and 2\n",0);
      return 0;
   }

   mu = (double)(*getarg(iarg++));
   if ( mu < 0 ) {
      printf("mu must be >= 0");
      mu=0;
   }

   if(ifarg(5)){
     sigma = (double)(*getarg(iarg++));
   }else{
     sigma = 0;
       }

    spks = *((void**)(&xspikes));
    if(spks){
     vector_resize(spks, nsizemax);
    }
    inst = ((double*) vector_vec(spks));      //Get array ptr to spks Vector


    return ffGn(inst,nsizemax,tdres,Hinput,mu,sigma);

}

ENDVERBATIM


PROCEDURE install_an_zilany_v4()
{
VERBATIM
   install_vector_method("an_zilany_v4", an_zilany_v4);
   install_vector_method("ANFSpikeGenerator3", ANFSpikeGenerator3);
   install_vector_method("fast_fGn", fast_fGn);
ENDVERBATIM
}

   
