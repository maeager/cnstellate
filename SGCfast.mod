:  Spiral Ganglion Cell- NEURON Artificial Cell


NEURON   {
  ARTIFICIAL_CELL SGC_fast

  RANGE spikecount, spont, cf, channel
  RANGE seed, reset, maxsout
  RANGE c0,s0,c1,s1, abs_refr

}
DEFINE MAXSPIKES 100
VERBATIM
extern double* vector_vec();
extern int vector_capacity();
extern void* vector_arg();

extern void vector_resize();

#ifndef __max
#define __max(a,b) (((a) > (b))? (a): (b))
#endif

#define USE_RANDOM_VECTOR 0
int nstim;
ENDVERBATIM


PARAMETER {
   spikecount= 0.0
   stimdur   = 0.0
   nrep      = 1
   spont     = 65.0  <0,150>    : Spontaneous Rate
   channel   = 0
    cf       = 1000    (Hz) <20, 40000>
    stimtdres= 2e-5 (s)        <1e-9,1>
    seed     = 1            : default variable seed for rand num generator
   reset      = 0
   realtime   = 0.0         :seconds
   rsptime    = 0.0
:Poisson random number exponentials and constants
   c0      = 0.5
   s0      = 0.001
   c1       = 0.5
   s1       = 0.0125
   abs_refr    = 0.00075        :absolute refractory period in seconds
   maxsout = 0

}



ASSIGNED {

   spkindex
   eventtime (ms)
   ratein
   spiketimesAddress

}

INITIAL {
   :SGfast(1)   //This must be done after finitialize()
   spkindex = 0
   if (GetSpike()) {
      net_send(eventtime - t, 1)
   }
}


NET_RECEIVE (w){
   if (flag == 0) {:0 progagate synaptic external event
      net_event(t)

   }
   if (flag == 1) {  :1 self event
      net_event(t)

      if (GetSpike())  {
         net_send(eventtime-t, 1)
      }
   }
}


PROCEDURE GetSpike() {
VERBATIM
  { void* vv; int i, size, trigger; double* spksout;
   i = (int)spkindex;
   trigger =0;
   if (i >= 0) {
      vv = *((void**)(&spiketimesAddress));
      if (vv) {
         size = (int)spikecount;
         spksout = vector_vec(vv);
         if (i < size) {
            eventtime = spksout[i];
            spkindex += 1.;
            trigger = 1;

         }else{
             /* spkindex = -1.;  */
         }
      }else{
         /*pkindex = -1.;*/
      }
   }
   return trigger;
  }
ENDVERBATIM
}




:user function to set instantaneous rate function of fiber
PROCEDURE SetFibreRate()
{
VERBATIM
   void **vv;
   double *px;
    int i;
   double max;
   extern void* vector_arg();
   max = 0;
#ifdef DEBUG
   printf("Setting Fibre Rate %g  \n",  &spiketimesAddress);
#endif
   if (ifarg(3)) {
      vv = (void**)(&ratein);
      *vv = (void*)0;
      *vv = vector_arg(1);
      nstim = vector_capacity(*vv);
      
#ifdef DEBUG   
      for (i=0;i<nstim;i++){
         if ((*vv)[i] > max) { max = (*vv)[i];}
         printf("%f ", ratein[i]);
      }
      maxsout = max;
#endif
      vv = (void**)(&spiketimesAddress);
      *vv = (void*)0;
      *vv = vector_arg(2);
      stimtdres = *getarg(3);
      
#ifdef DEBUG
      printf("Instantaneous rate set with %d points, resolution %g and max %g \n" , nstim, stimtdres, maxsout);
#endif
  }
  else{
      printf("Instantaneous rate has NOT been set \n");
  }

ENDVERBATIM
}

PROCEDURE SetSpikes()
{
VERBATIM

   void** vv;
   vv = (void**)(&spiketimesAddress);
   *vv = (void*)0;
   if (ifarg(1)) {
      *vv = vector_arg(1);
      if (ifarg(2)){
         spikecount = *getarg(2);
      } else{
         spikecount = vector_capacity(*vv);
      }
   }
ENDVERBATIM
}

PROCEDURE SpkGen()
{
VERBATIM
   static long seed_time=0;

   double  *spksout, *rate;
   void  *spks, *vrate;
   double rint,prob,randprob;

     int out;
   int  nout, i, size, vecsize;
/*range variables  stimtdres, spontrate, xseed also used */
   realtime   = 0.0;
   rsptime    = 0.0;
   
   
#ifdef DEBUG
   printf("In SpkGen : \n");
#endif
   
   spks  = NULL;
   spks = *((void**)(&spiketimesAddress));     /*void pointer to spiketimesAddress vector */
   vector_resize(spks, 0);         /* delete all previous spikes */
   stimdur  = nstim*stimtdres*1000;
   nout = 0;
   vector_resize(spks, 20);
   vecsize = 20;
   spksout = vector_vec(spks);

   vrate = *((void**)(&ratein));
   rate = vector_vec(vrate);
   
#ifdef DEBUG   
   printf("In SpkGen: nstim %d\n", nstim);
#endif
   spikecount    = 0;
   if (spks) {
       /*
       Initialise random number generator
       On first run reset the random number generator
       Check for user seed to initialise random number generator
       */
       if (reset = 0){
           if(seed != 0)       
           set_seed(seed);  
           reset = 1;
       } 
       
       rsptime =  0.0 - scop_random()*1/spont;
       #ifdef DEBUG
       printf("First rand %f \nSpikeTimes: ",SG.rsptime);
       #endif
       /* Begin spike simulation */
       spikecount = 0.0;
       for(i=0;i<nstim;i++)
       {
           out=0;
           realtime += stimtdres; /* running time */
           /* interval from last spike */
           rint = realtime - rsptime;
           if(rint > abs_refr )
           {
               prob = rate[i]*stimtdres*(1.0-( c0*exp(-(rint - abs_refr)/s0) + c1*exp(-(rint - abs_refr)/s1)));
               randprob = scop_random();
               if( prob > randprob)
               {
		   #ifdef DEBUG
		   printf("Prob %f, Rand %f ", prob, randprob);
		   #endif
		   /*New Spike occurance*/
		   rsptime = realtime;
		   out=1;
               }
               if(out==1){
		   /*Increment spike counter*/
		   spikecount += 1.;
		   size = (int)spikecount;
		   
		   /*Resize vector if needed*/
		   #ifdef DEBUG
		   printf("Vectorsize %d  %d:", size, vecsize);
		   #endif
		   
		   if (size == vecsize){
                       vecsize +=10;
                       vector_resize(spks, vecsize);      /*resize spike vector*/
                       spksout = vector_vec(spks);      /*set spiketimes ptr to resized vector*/
		   }
		   spksout[size-1] = realtime*1000;
                   /*store new spike time in milliseconds */
		   #ifdef DEBUG
		   printf("%.2f ", spksout[size-1]);
		   #endif
		   
               }
           }
	   
	   
       }
       #ifdef DEBUG
       printf("\n");
       printf("\nSize of Spike times %d\n", size);
       #endif
   }else{
      printf("Cannot access ratein and spks");
   }
   /*   resize spike vector*/
   size = (int)spikecount;
   vector_resize(spks, size);
   #ifdef DEBUG
   printf("Returning from SpkGen %d\n\n", (int)spikecount);
   #endif
   
   return spikecount;
ENDVERBATIM
}



PROCEDURE SGfast()
{
VERBATIM
void        *spks, *vrate;
int       NoutMax, i, Incr;
double      *rate,  deadtime;
long      j, k, numRate, nrep;
long       Nout, randBufLen, randBufIndex;
double      meanRate, *spktimes,  *randNums;
double      deadtimeRnd, endOfLastDeadtime, refracMult0, refracMult1, refracValue0, refracValue1;
double      Xsum, unitRateIntrvl, stime;
long      deadtimeIndex;


/*Get "Spks" vector object*/
   spks  = NULL;
   spks = *((void**)(&spiketimesAddress));     /*void pointer to spiketimes vector*/
   vector_resize(spks, 0);         /*delete all previous spikes*/

/* Get the "rate" vector object*/
   vrate = *((void**)(&ratein));       /*Set vrate to ratein Vector*/
   rate = vector_vec(vrate);      /*Get array ptr*/
if (rate == NULL || spks == NULL ) hoc_execerror("SGfast: Unable to locate Rate or Spks",0);


/*Get number of Reps*/
if (ifarg(1)){
   nrep = (int) (*getarg(1));
}else nrep = 1; /*hoc_execerror("SGfast: Must have one argument- nreps",0);*/

/* Get module parameters*/
deadtime = abs_refr;
numRate = nstim;
stimdur = nstim*stimtdres*1000;
   spikecount    = 0;

/* Determine the mean of the rate vector*/
meanRate = 0.0;
for (k=0; k<numRate;k++)
   if (rate[k] > 0.0)  meanRate += rate[k];
meanRate /= numRate;


/* Create output buffer for spike times;*/
   NoutMax = ((int) (meanRate * stimdur * nrep))+1;
   /*printf("NoutMax = %d\n", NoutMax);*/
   if (NoutMax < 10)  NoutMax = 10;
   
   #ifdef DEBUG
   printf("Rate vector mean %g and size  %d\n T dur %g",meanRate,numRate, stimdur);
   #endif
   vector_resize(spks, NoutMax);      /*make space for spikes*/
   spktimes = vector_vec(spks);      /*Get array ptr to spks Vector*/


   Nout = 0;

#if USE_RANDOM_VECTOR
/*Get a vector of pseudo-random numbers.*/
if ( NoutMax < 1 )  hoc_execerror("Random number buffer length must be positive.",0);
randBufLen = (int)(NoutMax+1);  /* Need 1 "extra" pseudo-random number to initialize the process*/
randNums = makevector((int)randBufLen);
for (i=0;i<randBufLen;i++) randNums[i] = scop_random();
#endif
randBufIndex = 0;


/* Calculate useful constants*/
deadtimeIndex = (long) floor(deadtime/stimtdres);   /* Integer number of discrete time bins within deadtime*/
deadtimeRnd = deadtimeIndex*stimtdres;      /* Deadtime rounded down to length of an integer number*/
                     /*of discrete time bins*/

refracMult0 = 1 - stimtdres/s0;  /* If y0(t) = c0*exp(-t/s0), then y0(t+stimtdres) = y0(t)*refracMult0*/
refracMult1 = 1 - stimtdres/s1;  /* If y1(t) = c1*exp(-t/s1), then y1(t+stimtdres) = y1(t)*refracMult1*/

/* Calculate effects of a random spike before t=0 on refractoriness and the time-warping sum at t=0*/
#if USE_RANDOM_VECTOR
endOfLastDeadtime = log( randNums[randBufIndex++] ) / rate[0] + deadtime;  /* End of last deadtime before t=0*/
#else
endOfLastDeadtime = log( scop_random() ) / rate[0] + deadtime;  /* End of last deadtime before t=0*/
#endif
refracValue0 = c0*exp(endOfLastDeadtime/s0);  /* Value of first exponential in refractory function*/
refracValue1 = c1*exp(endOfLastDeadtime/s1);  /* Value of second exponential in refractory function*/
Xsum = rate[0] * ( -endOfLastDeadtime + c0*s0*(exp(endOfLastDeadtime/s0)-1) + c1*s1*(exp(endOfLastDeadtime/s1)-1) );  /* Value of time-warping sum*/
/*  ^^^^ This is the "integral" of the refractory function ^^^^ (normalized by 'stimtdres')*/

/* Calculate first interspike interval in a homogeneous, unit-rate Poisson process (normalized by 'stimtdres')*/
#if USE_RANDOM_VECTOR
unitRateIntrvl = -log( randNums[randBufIndex++] ) / stimtdres;
#else
unitRateIntrvl = -log( scop_random() ) / stimtdres;
#endif
/* NOTE: Both 'unitRateInterval' and 'Xsum' are divided (or normalized) by 'stimtdres' in order to reduce calculation time.
      This way we only need to divide by 'stimtdres' once per spike (when calculating 'unitRateInterval'), instead of
      multiplying by 'stimtdres' once per time bin (when calculating the new value of 'Xsum').
*/
#ifdef DEBUG
printf("\n Stimdur= %g \n \
stimtdres= %g \n \
randBufLen = %d \n \
deadtimeIndex = %d \n \
deadtimeRnd = %g \n \
refracMult0 = %g \n \
refracMult1 = %g \n \
endOfLastDeadtime = %g \n \
refracValue0 = %g \n \
refracValue1 = %g \n \
Xsum = %g \n ", stimdur,stimtdres, \
randBufLen,deadtimeIndex,deadtimeRnd,refracMult0,refracMult1,endOfLastDeadtime, \
refracValue0,refracValue1,Xsum);
#endif

/*Loop through rate vector*/
stime = stimtdres;
k=0;
for (j=0; j<nrep; ++j)
{
   for (; (k<numRate) && (stime<stimdur); \
   ++k, stime+=stimtdres, refracValue0*=refracMult0, refracValue1*=refracMult1)
   {
      if (rate[k] > 0.0)
        {
         Xsum += rate[k] * (1 - refracValue0 - refracValue1);
         if ( Xsum >= unitRateIntrvl )
         {
            spktimes[Nout] = stime*1000;
            if (++Nout >= NoutMax)
            {
               /*Incr = (int)(meanRate * ((stimdur-stime) + (nrep-j-1)*stimdur));*/
               /*if (Incr < 10)  Incr = 10;*/
               Incr = 10;
               NoutMax += Incr;
               vector_resize(spks, NoutMax);      /*resize spike vector*/
               spktimes = vector_vec(spks);      /*set spiketimes ptr to resized vector*/
               #if USE_RANDOM_VECTOR
               freevector(randNums);
               randNums=makevector(Incr);
               for (i=0;i<Incr;++i) randNums[i] = scop_random();
               randBufIndex = 0;
               #endif
            }
         /* Next interspike "stime" in unit-rate process*/
            #if USE_RANDOM_VECTOR
            unitRateIntrvl = -log( randNums[randBufIndex++] ) / stimtdres;
            #else
            unitRateIntrvl = -log( scop_random() ) / stimtdres;
            #endif
            Xsum = 0.0;
         /* Increase index and stime to the last time bin in the deadtime, and reset (relative) refractory function*/
            k += deadtimeIndex;
            stime += deadtimeRnd;
            refracValue0 = c0;
            refracValue1 = c1;
        }
    }
} /* End of rate vector loop*/

/* Reset index and time to begin a new repetion.
Don't just set to zero, since new repetition may start within the
deadtime following the last spike in the previous repetition. */
   stime -= stimdur/1000;
   k -= numRate;
} /* End of nrep loop*/

#ifdef DEBUG
printf("Total spikes before = %d\n", Nout);
#endif
/* 
Delete spike(s) that occur after the last repetition of the rate function ends*/
for (; (Nout>0)&&(spktimes[Nout-1]>stimdur); )  --Nout;

#if USE_RANDOM_VECTOR
freevector(randNums);
#endif

#ifdef DEBUG
printf("Total spikes = %d\n", Nout);
#endif
vector_resize(spks, Nout);
spikecount=Nout;
return Nout;
ENDVERBATIM
}



PROCEDURE SGfast2()
{
VERBATIM
/* 
   Pass the output of Synapse model through the Spike Generator 
   The spike generator now uses a method coded up by B. Scott Jackson
   (bsj22@cornell.edu) Scott's original code is available from Laurel
   Carney's web site at:
   http://www.urmc.rochester.edu/smd/Nanat/faculty-research/lab-pages/LaurelCarney/auditory-models.cfm
*/

 void *spks, *vrate;
 long NoutMax, i;
 double   *rate,  deadtime,tdres,totalstim;
 int      j, k, numRate, nrep;

 long     nspikes, Nout, deadtimeIndex, randBufIndex;
 double deadtimeRnd, endOfLastDeadtime, refracMult0, refracMult1, refracValue0, refracValue1;
 double Xsum, unitRateIntrvl, countTime, DT;

 double *randNums,*spktimes;


/*Get "Spks" vector object*/
   spks  = NULL;
   spks = *((void**)(&spiketimesAddress));     /*void pointer to spike times vector*/
   vector_resize(spks, 0);         /*delete all previous spikes*/

/* Get the "rate" vector object*/
   vrate = *((void**)(&ratein));       /*Set vrate to ratein Vector*/
   rate = vector_vec(vrate);      /*Get array ptr*/
   if (rate == NULL || spks == NULL ) hoc_execerror("SGfast2: Unable to locate Rate or Spks",0);

/*Get number of Reps*/
   if (ifarg(1)){
     nrep = (int) (*getarg(1));
   }else nrep = 1; 

/* Get module parameters*/
   deadtime = abs_refr;
   numRate = nstim;
   stimdur = nstim*stimtdres*1000;
   spikecount    = 0;
   tdres = stimtdres;
   totalstim=nstim;
    
    DT = totalstim * tdres * nrep;  /* Total duration of the rate function */
    Nout = 0;
    NoutMax = (long) ceil(DT / deadtime);
    randNums = makevector(NoutMax + 10);

    vector_resize(spks,NoutMax + nrep);
    spktimes = vector_vec(spks);      /*set spiketimes ptr to resized vector*/

    for (k=0;k<=NoutMax;k++) randNums[k] = scop_random();  /*replace MATLAB call rand(1,N)*/
    randBufIndex = 0;

    /* Calculate useful constants */
    deadtimeIndex = (long) floor(deadtime / tdres);  /* Integer number of discrete time bins within deadtime */
    deadtimeRnd = deadtimeIndex * tdres;   /* Deadtime rounded down to length of an integer number of discrete time bins */

    refracMult0 = 1 - tdres / s0;  /* If y0(t) = c0*exp(-t/s0), then y0(t+tdres) = y0(t)*refracMult0 */
    refracMult1 = 1 - tdres / s1;  /* If y1(t) = c1*exp(-t/s1), then y1(t+tdres) = y1(t)*refracMult1 */

    /* Calculate effects of a random spike before t=0 on refractoriness and the time-warping sum at t=0 */
    endOfLastDeadtime = __max(0, log(randNums[randBufIndex++]) / rate[0] + deadtime); /* End of last deadtime before t=0 */
    refracValue0 = c0 * exp(endOfLastDeadtime / s0); /* Value of first exponential in refractory function */
    refracValue1 = c1 * exp(endOfLastDeadtime / s1); /* Value of second exponential in refractory function */
    Xsum = rate[0] * (-endOfLastDeadtime + c0 * s0 * (exp(endOfLastDeadtime / s0) - 1) + c1 * s1 * (exp(endOfLastDeadtime / s1) - 1));
    /* Value of time-warping sum */
    /*  ^^^^ This is the "integral" of the refractory function ^^^^ (normalized by 'tdres') */

    /* Calculate first interspike interval in a homogeneous, unit-rate Poisson process (normalized by 'tdres') */
    unitRateIntrvl = -log(randNums[randBufIndex++]) / tdres;
    /* NOTE: Both 'unitRateInterval' and 'Xsum' are divided (or normalized) by 'tdres' in order to reduce calculation time.
       This way we only need to divide by 'tdres' once per spike (when calculating 'unitRateInterval'), instead of
       multiplying by 'tdres' once per time bin (when calculating the new value of 'Xsum').                         */
    countTime = tdres;
    
    #ifdef DEBUG
     printf(" c0 %g\ts0 %g\tc1 %g\ts1 %g\tdead %g\n",c0, s0, c1, s1, deadtime);
     printf(" nspikes %ld\tk %ld\tNoutMax %ld\tNout %ld\tdeadtimeIndex %ld\trandBufIndex %ld\n", nspikes, k, NoutMax, Nout, deadtimeIndex, randBufIndex);
     printf("deadtimeRnd %g\tendOfLastDeadtime %g\trefracMult0 %g\trefracMult1 %g\trefracValue0 %g\trefracValue1 %g\n", deadtimeRnd, endOfLastDeadtime, refracMult0, refracMult1, refracValue0, refracValue1);
     printf("Xsum %g\tunitRateIntrvl %g\tcountTime %g\tDT %g\n",Xsum, unitRateIntrvl, countTime, DT);
     #endif
    for (j=0;j<nrep;j++){
      for (k = 0; (k < totalstim) && (countTime < DT); ++k, countTime += tdres, refracValue0 *= refracMult0, refracValue1 *= refracMult1) { /* Loop through rate vector */
        if (rate[k] > 0) { /* Nothing to do for non-positive rates, i.e. Xsum += 0 for non-positive rates. */
	  Xsum += rate[k] * (1 - refracValue0 - refracValue1);  /* Add synout*(refractory value) to time-warping sum */

	  if (Xsum >= unitRateIntrvl) {  /* Spike occurs when time-warping sum exceeds interspike "time" in unit-rate process */
	    spktimes[Nout] = countTime*1000.0;
	    Nout = Nout + 1;
	    unitRateIntrvl = -log(randNums[randBufIndex++]) / tdres;
	    Xsum = 0;

	    /* Increase index and time to the last time bin in the deadtime, and reset (relative) refractory function */
	    k += deadtimeIndex;
	    countTime += deadtimeRnd;
	    refracValue0 = c0;
	    refracValue1 = c1;
	  }
        }
      } /* End of rate vector loop */
      if (nrep ==1) break;
      countTime = tdres;
      spktimes[Nout++] = 0.0;
    }
    freevector(randNums);
    
    
    
    #ifdef DEBUG
    printf("Total spikes before = %d\n", Nout);
    #endif
    
/* Delete spike(s) that occur after the last repetition of the rate function ends*/
/*for (; (Nout>0)&&(spktimes[Nout-1]>stimdur); )  --Nout;*/

#ifdef DEBUG
printf("Total spikes = %d\n", Nout);
#endif

vector_resize(spks, Nout);
spikecount=Nout;
return Nout;
ENDVERBATIM

}
