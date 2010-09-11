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

 /* #define DEBUG   */
/* #define _FFGN_ */ 

#include "complex.c"
#include "carneymodel.c"
#include "zilanycarneyv4.c"
#include "resample.c"
#include "ffGn.c"


static double an_zilany_v4(void  *vv)
{

   double *stim;      /*Input stimulus vector in pascals*/
   double *ihcout,*sout;   /*Output vector containing inst. rate for channel*/
   double tdres=0.0,cf=0.0,spont=0.0;
   double cohc=1.0,cihc=1.0;
   int species=0,fibertype=0,implnt=0;
   int ifspike=0,nrep=0,out=0;
   int nstim=0, nsout=0,nihcout=0;


/*Get Instance Sout vector for synapse data*/
   nstim = vector_instance_px(vv, &sout);
   nsout = vector_arg_px(1, &stim);
   if (nstim != nsout){
       printf ("an_zilany_v4: sout must be the same size as stim.\n");
       return 0;
   }
   
/*Get Input arguments*/
   if( ifarg(9)!=1 && ifarg(10)!=1 ){  /*Must be changed if more input arguments added*/
      hoc_execerror("an_zilany_v4: input syntax must be sout.an_zilany_v4( stim, tdres,  cf, fibertype,implnt,cihc,cohc, species,nrep,ihcout[optional])", 0);
      return 0;
   }

 /*TDRES  resolution of stim vector*/
   /*Bruce model uses seconds rather than msec*/
   tdres = (double)(*getarg(2));
   if (tdres > 0.01e-3 || tdres < 0.002e-3){
      printf("\tNote: Zilany Bruce V4 resolution should be between 0.01ms (Fs = 100kHz) for normal usage and 0.002ms (200kHz) for stim above 40kHz.\n");      /*printf("\tNote: ZilanyBruceV4 resolution should be between 0.01ms (Fs = 100kHz) and 0.002ms (500kHz) for normal usage.\n");*/
      /*tdres = 0.002e-3;*/
   }
   /*CF of fiber*/
   cf = (double)(*getarg(3));

   if ((cf<80)||(cf>70e3))
   {
      hoc_execerror("ERROR: cf must be between 80 Hz and 70 kHz\n",0);
      return 0;
   }
   /*Spontaneous rate of ANF*/
   fibertype = (int) round(*getarg(4));
   if ((fibertype<=0)||(fibertype>3))
   {
      printf("an_zilany_v4: fibertype  must be 1, 2, or 3 \n Default fibertype = 3 (HSR)");
      fibertype=3;
   }
   /*New variable in version 4*/
   implnt = (int)round(*getarg(5));
   if ((implnt<0)||(implnt>1))
   {
      printf("an_zbcatmodel: implnt  must be 0 (actual) or 1 (approx)\n");
      implnt=1;
   }
   
   /*Model*/
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

   /*Species*/
   species = (int)round(*getarg(8));
   if (species != 1 && species !=9){
     printf("an_zilany_v4: species other than cat (1 or 9) are not fully implemented");
   }
   /*Reps*/
   nrep = (int)round(*getarg(9));
   if (nrep < 1) {printf("an_zilany_v4: nrep must be 1 or greater"); return 0;}
   
if(ifarg(10)) {
  nihcout = vector_arg_px(10, &ihcout);
  if(nihcout != nstim){
    printf("ihcout must be the same size as stim, sout\n"); return 0;}
 } else {
  ihcout = makevector(nstim);
 }

   printf("AN model: Zilany et al., 2009  (version 4 c2010)\n");
   printf("  IHCAN(stim,%.0f,%d,%g,%d,%g,%g,ihcout,%d)\n", cf, nrep, tdres, nstim, cohc, cihc, species);
   IHCAN(stim, cf, nrep, tdres, nstim, cohc, cihc, ihcout,species);
   printf("  SingleAN(ihcout,%.0f,%d,%g,%d,%d,%d,sout,%d)\n",cf,nrep,tdres,nstim,fibertype,implnt,species);
   out= SingleAN_v4(ihcout,cf,nrep,tdres,nstim,fibertype,implnt,sout,species);

   if( !ifarg(10) ) freevector(ihcout);
   return out; 
}

static double an_zilany_v4_1(void *vv)
{

   double *stim;      /*Input stimulus vector in pascals*/
   double *ihcout,*sout;   /*Output vector containing inst. rate for channel*/
   double tdres,cf,spont,out;
   double cohc,cihc;
   int species,implnt;
   int ifspike;
   int nstim, nsout,nihcout;
   int nrep;
   cohc = 1.0;
   cihc = 1.0;

/*Get Instance Sout vector for synapse data*/
   nstim = vector_instance_px(vv, &sout);
   nsout = vector_arg_px(1, &stim);

/*Get Input arguments*/
/*   if(ifarg(9)!=1){  //Must be changed if more input arguments added
      hoc_execerror("ERROR: input syntax must be sout.an_zilany_v4( stim, tdres,  cf, spont,implnt,cihc,cohc, species,nrep,ihcout[optional])", 0);
      return 0;
   }
*/
 /*TDRES  resolution of stim vector*/
   /*Bruce model uses seconds rather than msec*/
   tdres = (double)(*getarg(2));
   if (tdres > 0.01e-3 || tdres < 0.002e-3){
      printf("Note: Zilany Bruce V4 resolution should be between 0.01ms (Fs = 100kHz) for normal usage and 0.002ms (200kHz) for stim above 40kHz.\n");  
      /*tdres = 0.002e-3;*/
   }
   /*CF of fiber*/
   cf = (double)(*getarg(3));

   if ((cf<80)||(cf>70e3))
   {
      hoc_execerror("ERROR: cf must be between 80 Hz and 70 kHz\n",0);
      return 0;
   }
   /*Spontaneous rate of ANF*/
   spont = (double)(*getarg(4));
   
   /*implnt new variable in version 4*/
   implnt = (long)(*getarg(5));
   if ((implnt<0)||(implnt>1))
   {
      printf("an_zbcatmodel: implnt  must be 0 (actual) or 1 (approx)\n");
      implnt=1;
   }
   
   /*Model*/
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

   /*Species*/
   species = (int)(*getarg(8));
   if (species != 1){
     hoc_execerror("an_zilany_v4: species other than cat (1) are not implemented",0);
     return 0;
   }
   /*Reps*/
   nrep = (int)(*getarg(9));
   if (nrep < 1) {printf("an_zilany_v4: nrep must be 1 or greater"); return 0;}
   
if(ifarg(10)) {
  nihcout = vector_arg_px(10, &ihcout);
  if(nihcout != nstim){
    printf("ihcout must be the same size as stim, sout\n"); return 0;}
 } else {
  ihcout = makevector(nstim);
 }

   printf("\tAN model: Zilany, Carney, Bruce, Nelson and others  (version 4 c2010)\n");
   printf("\tIHCAN(stim,%.0f,%d,%g,%d,%g,%g,ihcout,%d)\n", cf, nrep, tdres, nstim, cohc, cihc, species);
   IHCAN(stim, cf, nrep, tdres, nstim, cohc, cihc, ihcout,species);
   printf("\tSingleAN(ihcout,%.0f,%d,%g,%d,%d,%d,sout,%d)\n",cf,nrep,tdres,nstim,spont,implnt,species);
   out= SingleAN_v4_1(ihcout,cf,nrep,tdres,nstim,spont,implnt,sout,species);

   if( !ifarg(10) ) freevector(ihcout);
   return out; 
}

static double ihc_zilany_v4(void *vv)
{
   double *stim;      /*Input stimulus vector in pascals*/
   double *ihcout;   /*Output vector containing inst. rate for channel*/
double tdres,cf;
   double cohc,cihc;
   int species;
   int nstim, nihcout;
   int nrep;
   cohc = 1.0;
   cihc = 1.0;

/*Get Instance Sout vector for synapse data*/
   nihcout = vector_instance_px(vv, &ihcout);
   nstim = vector_arg_px(1, &stim);

/*Get Input arguments*/
   if(ifarg(7)!=1){  /*Must be changed if more input arguments added*/
      hoc_execerror("ERROR: input syntax must be ihcout.ihc_zilany_v4( stim, tdres,  cf,cihc,cohc, species,nrep)", 0);
      return 0;
   }
   /*TDRES  resolution of stim vector*/
   /*Bruce model uses seconds rather than msec*/
   tdres = (double)(*getarg(2));
   if (tdres > 0.01e-3 || tdres < 0.002e-3){
      printf("Note: Zilany Bruce V4 resolution should be between 0.01ms (Fs = 100kHz) for normal usage and 0.002ms (200kHz) for stim above 40kHz.\n");      /*printf("Note: ZilanyBruceV4 resolution should be between 0.01ms (Fs = 100kHz) and 0.002ms (500kHz) for normal usage.\n");*/
      /*tdres = 0.002e-3;*/
   }
   /*CF of fiber*/
   cf = (double)(*getarg(3));

   if ((cf<80)||(cf>70e3))
   {
      hoc_execerror("ERROR: cf must be between 80 Hz and 70 kHz\n",0);
      return 0;
   }
   
   /*Model*/
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

   /*Species*/
   species = (int)(*getarg(6));
   if ( species != 1 || species != 9 ){
     printf("an_zilany_v4: species other than cat (1) are not correctly implemented");

   }
   /*Reps*/
   nrep = (int)(*getarg(7));
   if (nrep < 1) {printf("an_zilany_v4: nrep must be 1 or greater"); return 0;}

   printf("\tAN model: Zilany, Bruce, Nelson and Carney  (version 4 c2010)\n");
   printf("\tIHCAN(stim,%.0f,%d,%g,%d,%g,%g,ihcout,%d)\n", cf, nrep, tdres, nstim, cohc, cihc,species);
   IHCAN(stim, cf, nrep, tdres, nstim, cohc, cihc, ihcout,species);
   return (double) nihcout; 
}

static double syn_zilany_v4(void *vv)
{

  double *psth,*sptime,*sptimes,*ihcout,*sout;   /*Input and Output vector containing inst. rate for channel*/
  void** xx;
  void *spks,*pst;
   double tdres,cf,out;
   double xspikes,xpsth;
   int species,fibertype,implnt;
   int ifspike,ifpsth;
   int nihcout,nspikes, nsout;
   int nrep,i,ipst;
   out=0;
   ifspike=0;ifpsth=0;
/*Get Instance Sout vector for synapse data*/
   nsout = vector_instance_px(vv, &sout);
   nihcout = vector_arg_px(1, &ihcout);

/*Get Input arguments*/
   if(ifarg(8)!=1  || ifarg(7)!=1 || ifarg(9)!=1 ){  /*Must be changed if more input arguments added*/
      hoc_execerror("ERROR: input syntax must be sout.syn_zilany_v4( ihcout, tdres,  cf, fibertype,implnt, species,nrep)\n or  sout.syn_zilany_v4( ihcout, tdres,  cf, fibertype,implnt, species,nrep,spks [optional],psth [optional])", 0);
      return 0;
   }
   /*TDRES  resolution of stim vector*/
   /*Bruce model uses seconds rather than msec*/
   tdres = (double)(*getarg(2));
   if (tdres > 0.01e-3 || tdres < 0.002e-3){
      printf("Note: Zilany Bruce V4 resolution should be between 0.01ms (Fs = 100kHz) for normal usage and 0.002ms (200kHz) for stim above 40kHz.\n");
      /*printf("Note: ZilanyBruceV4 resolution should be between 0.01ms (Fs = 100kHz) and 0.002ms (500kHz) for normal usage.\n");*/
      /*tdres = 0.002e-3;*/
   }
   /*CF of fiber*/
   cf = (double)(*getarg(3));

   if ((cf<80)||(cf>70e3))
   {
      hoc_execerror("ERROR: cf must be between 80 Hz and 70 kHz\n",0);
      return 0;
   }
   /*Spontaneous rate of ANF*/
   fibertype = (int)(*getarg(4));
   if ((fibertype<=0)||(fibertype>3))
   {
      printf("an_zilany_v4: fibertype  must be 1, 2, or 3 \n Default fibertype = 3 (HSR)");
      fibertype=3;
   }
   /*New variable in version 4*/
   implnt = (int)(*getarg(5));
   if ((implnt<0)||(implnt>1))
   {
      printf("an_zbcatmodel: implnt  must be 0 (actual) or 1 (approx)\n");
      implnt=1;
   }
   
   /*Species*/
   species = (int)(*getarg(6));
   if (species != 1){
     hoc_execerror("an_zilany_v4: species other than cat (1) are not implemented",0);
     return 0;
   }
   /*Reps*/
   nrep = (int)(*getarg(7));
   if (nrep < 1) {printf("an_zilany_v4: nrep must be 1 or greater\n"); return 0;}

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
   if(ifarg(9)) {
     ifpsth = 1;
   
       xx = (void**)(&xpsth);
       *xx = (void*)0;
       if (ifarg(9)) {
	 *xx = vector_arg(9);
       }
       pst = *((void**)(&xpsth));
       if (pst) 
	 vector_resize(pst,0);

   }


   printf("\tAN model: Zilany, Carney, Bruce, Nelson and others  (version 4 c2010)\n");
   printf("\tSingleAN(ihcout,%.0f,%d,%g,%d,%d,%d,sout,%d)\n",cf,nrep,tdres,nihcout,fibertype,implnt,sout,species);
    out= SingleAN_v4(ihcout,cf,nrep,tdres,nihcout,fibertype,implnt,sout,species);


    /*======  Spike Generations ======*/
      if(ifspike){
	  sptimes  = makevector((int) ceil(out / 0.00075)); zero_vector(sptimes,(int) ceil(out / 0.00075));
	  nspikes = SpikeGenerator_v4(sout, tdres, out, nrep, sptimes);
	  spks = *((void**)(&xspikes));
	  if (spks)  vector_resize(spks, nspikes);     
	  sptime = ((double*) vector_vec(spks));      /*Get array ptr to spks Vector	*/
	  for(i = 0; i < nspikes; i++)
	  {
	      sptime[i]  = sptimes[i]*1000.0;
	  };
	  out=nspikes;
      	  freevector(sptimes);     
	  
	  if(ifpsth){
	      pst = *((void**)(&xpsth));
	      if (pst)  vector_resize(pst, nsout);     
	      psth = ((double*) vector_vec(pst));      /*Get array ptr to psth Vector	*/
	      for(i = 0; i < nspikes; i++){        
		  ipst = (int) (fmod(sptime[i],tdres*nsout) / tdres);
		  psth[ipst] = psth[ipst] + 1;       
	      };   
	}
      }     
      return out; 
  }
  
static double psth_zilany_v4(void *vv)
{
   double *stim;      /*Input stimulus vector in pascals*/
   double *ihcout,*sout, *psth;   /*Output vector containing inst. rate for channel*/
   double tdres,cf,spont,out;
   double cohc,cihc;
   int species,fibertype,implnt;
   int ifspike;
   int nstim, nsout,nihcout,npsth;
   int nrep;
   cohc = 1.0;
   cihc = 1.0;

/*Get Instance Sout vector for synapse data*/
   nstim = vector_instance_px(vv, &sout);
   nsout = vector_arg_px(1, &stim);

/*Get Input arguments*/
   if(ifarg(10)!=1){  /*Must be changed if more input arguments added */
      hoc_execerror("ERROR: input syntax must be sout.psth_zilany_v4( stim, tdres,cf,fibertype,implnt,cihc,cohc,species,nrep,psth)", 0);
      return 0;
   }

 /*TDRES  resolution of stim vector*/
   /*Bruce model uses seconds rather than msec*/
   tdres = (double)(*getarg(2));
   if (tdres > 0.01e-3 || tdres < 0.002e-3){
      printf("Note: ZilanyBruceV4 resolution should be between 0.01ms (Fs = 100kHz) for normal usage and 0.002ms (200kHz) for stim above 40kHz.\n");
      /*tdres = 0.002e-3;*/
   }
   /*CF of fiber*/
   cf = (double)(*getarg(3));

   if ((cf<80)||(cf>70e3))
   {
      hoc_execerror("ERROR: cf must be between 80 Hz and 70 kHz\n",0);
      return 0;
   }
   /*type of ANF*/
   fibertype = (long)(*getarg(4));
   if ((fibertype<=0)||(fibertype>3))
   {
      printf("an_zilany_v4: fibertype  must be 1, 2, or 3 \n Default fibertype = 3 (HSR)\n");
      fibertype=3;
   }
   /*Implementation: new variable in version 4*/
   implnt = (long)(*getarg(5));
   if ((implnt<0)||(implnt>1))
   {
      printf("an_zbcatmodel: implnt  must be 0 (actual) or 1 (approx)\n");
      implnt=1;
   }
   
   /*Model*/
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

   /*Species*/
   species = (int)(*getarg(8));
   if (species != 1){
    printf("an_zilany_v4: cat (species=1 or 9) is the default, others are not fully tested\n");
   }
   /*Reps*/
   nrep = (int)(*getarg(9));
   if (nrep < 1) {printf("psth_zilany_v4: nrep must be 1 or greater\n"); return 0;}
   
   if(ifarg(10)) {
     npsth = vector_arg_px(10, &psth);
   }
  ihcout = makevector(nstim);

   printf("\tAN model: Zilany, Carney, Bruce, Nelson and others  (version 4 c2010)\n");
   printf("\tIHCAN(stim,%.0f,%d,%g,%d,%g,%g,ihcout,%d)\n", cf, nrep, tdres, nstim, cohc, cihc, species);
   IHCAN(stim, cf, nrep, tdres, nstim, cohc, cihc, ihcout,species);
   printf("\tPsthAN(ihcout,%.0f,%d,%g,%d,%d,%d,sout,%d)\n",cf,nrep,tdres,nstim,fibertype,implnt,species);
   PsthAN(ihcout,cf,nrep,tdres,nstim,fibertype,implnt,species,sout,psth);

   freevector(ihcout);
   return npsth; 
}


/* The spike generator now uses a method coded up by B. Scott Jackson (bsj22@cornell.edu) Scott's original code is available from Laurel Carney's web site at: http://www.urmc.rochester.edu/smd/Nanat/faculty-research/lab-pages/LaurelCarney/auditory-models.cfm
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
    void      *spks;  /*We need to be able to resize this vector*/
    double   xspikes,stimdur;

/*Get Instance Sout vector for synapse data*/
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
/*Get number of Reps and tdres*/

   if(ifarg(2)) nrep = (int) (*getarg(2));
   if(ifarg(3)) tdres = (double)(*getarg(3));

   c0      = 0.5;
   s0      = 0.001;
   c1       = 0.5;
   s1       = 0.0125;
    dead    = 0.00075;
    stimdur = nsout*tdres*1000;
    totalstim = nsout;
    DT = totalstim * tdres * nrep;  /* Total duration of the rate function */
    Nout = 0;
    NoutMax = (long) ceil(DT/dead);    
    
/*printf("\tSpikeGenerator: nsout %d,\t nreps %d,\ttdres %f\n",nsout,nrep,tdres);*/
    printf("\tANFSpikeGenerator3: resizing spks to %d, nrep %d stimdur %g DT %g\n",NoutMax,nrep,stimdur,DT);
    spks = *((void**)(&xspikes));
    if(spks){
     vector_resize(spks, NoutMax + nrep);
    }
    sptime = ((double*) vector_vec(spks));      /*Get array ptr to spks Vector*/

    Nout= SpikeGenerator_v4(sout, tdres, totalstim, nrep, sptime);
  
    vector_resize(spks, Nout);
    nspikes = Nout;  /* Number of spikes that occurred. */
    return(nspikes);
}

static double fast_fGn(void *vv)
{

  void *ptr_ffGn;
   double *vec_ffGn;   /*Output vector containing inst. rate for channel*/
   double tdres,cf,N,Hinput,mu,sigma;
   double xspikes;
   int nout,nsizemax;
   int iarg=1;
   
/*Get Instance/out vector */
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
 

/*Get Input arguments*/
   if( !ifarg(5) && !ifarg(6)){  /*Must be changed if more input arguments added*/
      hoc_execerror("fast_fGn: input syntax must be vec.fast_fGn(N, tdres, Hinput, mu, sigma (optional))", 0);
      return 0;
   }
   if(ifarg(iarg))   N = (double)(*getarg(iarg++));
   nsizemax = (int) N;

   /*Bruce model uses seconds rather than msec*/
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
      /*      vector_resize(ptr_ffGn, 32000);*/
    }
    vec_ffGn = ((double*) vector_vec(ptr_ffGn));      /*Get array ptr to ptr_ffGn Vector*/
    printf("\tfast_fGn: calling ffGn(&%x,%d,%g,%g,%g,%g)\n",&vec_ffGn,nsizemax,tdres,Hinput,mu,sigma);
    return ffGn(vec_ffGn,nsizemax,tdres,Hinput,mu,sigma);

}

static double rtresample(void *vv)
{
  void *ptr_resample;
  double *vec1_resample;
   double *vec2_resample;   /*Output vector containing inst. rate for channel*/
   double tdres,cf,N,Hinput,mu,sigma;
   double xresample,resamp;
   int nin,NSizeMax;
   int iarg=1;
   
/*Get Instance/out vector */
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
 
/*Get Input arguments*/
   if(!ifarg(2)){  /*Must be changed if more input arguments added*/
      hoc_execerror("rtresample: input syntax must be vec.rtresample(destvec, resample)", 0);
      return 0;
   }
   if(ifarg(iarg))   resamp = (double)(*getarg(iarg++));
   NSizeMax = (int) (resamp * (double)nin);
   if(NSizeMax <=0 ){  /*Must be changed if more input arguments added*/
      hoc_execerror("rtresample: size not big enough)", 0);
      return 0;
   }

    ptr_resample = *((void**)(&xresample));
    if(ptr_resample){
     vector_resize(ptr_resample, NSizeMax);
    }

    vec2_resample = ((double*) vector_vec(ptr_resample));      /*Get array ptr to ptr_resample Vector*/

    return resample(vec1_resample,vec2_resample,nin,resamp);

}
ENDVERBATIM


PROCEDURE install_an_zbcatmodel_v4()
{
VERBATIM
   install_vector_method("an_zilany_v4", an_zilany_v4);
   install_vector_method("an_zilany_v4_1", an_zilany_v4_1);
   install_vector_method("ihc_zilany_v4", ihc_zilany_v4);
   install_vector_method("syn_zilany_v4", syn_zilany_v4);
   install_vector_method("psth_zilany_v4", psth_zilany_v4);
   install_vector_method("ANFSpikeGenerator3", ANFSpikeGenerator3);
   install_vector_method("rtresample", rtresample);
   install_vector_method("fast_fGn", fast_fGn);
ENDVERBATIM
}

   
