:  Spiral Ganglion Cell- NEURON Artificial Cell


NEURON   {
  ARTIFICIAL_CELL SGC_fast

  RANGE spikecount,  spont, cf
  RANGE   seed, reset, maxsout
  RANGE  c0,s0,c1,s1, abs_refr

}
DEFINE MAXSPIKES 100
VERBATIM
extern double* vector_vec();
extern int vector_capacity();
extern void* vector_arg();

extern void vector_resize();
#define USE_RANDOM_VECTOR 0
int nstim;
ENDVERBATIM


PARAMETER {
   spikecount   =0.0
   stimdur      =0.0
   nrep      =1
    spont       = 65.0  <0,150>    : Spontaneous Rate
    cf       = 1000    (Hz) <20, 40000>
    stimtdres     = 2e-5 (s)        <1e-9,1>
    seed        = 1            : default variable seed for rand num generator
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
   spiketimes

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
      vv = *((void**)(&spiketimes));
      if (vv) {
         size = (int)spikecount;
         spksout = vector_vec(vv);
         if (i < size) {
            eventtime = spksout[i];
            spkindex += 1.;
            trigger = 1;

         }else{
            //spkindex = -1.;
         }
      }else{
         //spkindex = -1.;
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
//   printf("Setting Fibre Rate %g  \n",  &spiketimes);
   if (ifarg(3)) {

      vv = (void**)(&ratein);
      *vv = (void*)0;
      *vv = vector_arg(1);
      nstim = vector_capacity(*vv);

   /*   for (i=0;i<nstim;i++){
         if ((*vv)[i] > max) { max = (*vv)[i];}
      //   printf("%f ", ratein[i]);
      }
      maxsout = max;

   */
      vv = (void**)(&spiketimes);
      *vv = (void*)0;
      *vv = vector_arg(2);
      stimtdres = *getarg(3);

//   printf("Instantaneous rate set with %d points, resolution %g and max %g \n" , nstim, stimtdres, maxsout);
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
   vv = (void**)(&spiketimes);
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
//range variables  stimtdres, spontrate, xseed also used
   realtime   = 0.0;
   rsptime    = 0.0;


//   printf("In SpkGen : \n");
   spks  = NULL;
   spks = *((void**)(&spiketimes));     //void pointer to spiketimes vector
   vector_resize(spks, 0);         //delete all previous spikes
   stimdur  = nstim*stimtdres*1000;
   nout = 0;
   vector_resize(spks, 20);
   vecsize = 20;
   spksout = vector_vec(spks);

   vrate = *((void**)(&ratein));
   rate = vector_vec(vrate);

//   printf("In SpkGen: nstim %d\n", nstim);
   spikecount    = 0;


   if (spks) {



      //Initialise random number generator
       if (reset = 0){
         if(seed == 0){}//user has initialised random number generator
         else   set_seed(seed);  // user controlled seed

         reset =1;
      } //After set is 1,  do not reset random number generator




      rsptime =  0.0 - scop_random()*1/spont;
   //   printf("First rand %f \nSpikeTimes: ",SG.rsptime);
      //Begin spike simulation
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
            //   printf("Prob %f, Rand %f ", prob, randprob);
               //New Spike occurance
               rsptime = realtime;
               out=1;
            }
            if(out==1){
               //Increment spike counter
               spikecount += 1.;
               size = (int)spikecount;

               //Resize vector if needed
            //   printf("Vectorsize %d  %d:", size, vecsize);

               if (size == vecsize){
                  vecsize +=10;
                  vector_resize(spks, vecsize);      //resize spike vector
                  spksout = vector_vec(spks);      //set spiketimes ptr to resized vector
               }
               spksout[size-1] = realtime*1000;
                  //store new spike time in milliseconds
         //      printf("%.2f ", spksout[size-1]);

            }
                }


      }
//      printf("\n");
   //printf("\nSize of Spike times %d\n", size);
   }else{
      printf("Cannot access ratein and spks");
   }
   //   resize spike vector
   size = (int)spikecount;
   vector_resize(spks, size);
//   printf("Returning from SpkGen %d\n\n", (int)spikecount);
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


//Get "Spks" vector object
   spks  = NULL;
   spks = *((void**)(&spiketimes));     //void pointer to spiketimes vector
   vector_resize(spks, 0);         //delete all previous spikes

// Get the "rate" vector object
   vrate = *((void**)(&ratein));       //Set vrate to ratein Vector
   rate = vector_vec(vrate);      //Get array ptr
if (rate == NULL || spks == NULL ) hoc_execerror("SGfast: Unable to locate Rate or Spks",0);


//Get number of Reps
if (ifarg(1)){
   nrep = (int) (*getarg(1));
}else nrep = 1; //hoc_execerror("SGfast: Must have one argument- nreps",0);

// Get module parameters
deadtime = abs_refr;
numRate = nstim;
stimdur = nstim*stimtdres*1000;
   spikecount    = 0;

// Determine the mean of the rate vector
meanRate = 0.0;
for (k=0; k<numRate;k++)
   if (rate[k] > 0.0)  meanRate += rate[k];
meanRate /= numRate;


// Create output buffer for spike times;
   NoutMax = ((int) (meanRate * stimdur * nrep))+1;
   //printf("NoutMax = %d\n", NoutMax);
   if (NoutMax < 10)  NoutMax = 10;

   //printf("Rate vector mean %g and size  %d\n T dur %g",meanRate,numRate, stimdur);
   vector_resize(spks, NoutMax);      //make space for spikes
   spktimes = vector_vec(spks);      //Get array ptr to spks Vector


   Nout = 0;

#if USE_RANDOM_VECTOR
//Get a vector of pseudo-random numbers.
if ( NoutMax < 1 )  hoc_execerror("Random number buffer length must be positive.",0);
randBufLen = (int)(NoutMax+1);  // Need 1 "extra" pseudo-random number to initialize the process
randNums = makevector((int)randBufLen);
for (i=0;i<randBufLen;i++) randNums[i] = scop_random();
#endif
randBufIndex = 0;


// Calculate useful constants
deadtimeIndex = (long) floor(deadtime/stimtdres);   // Integer number of discrete time bins within deadtime
deadtimeRnd = deadtimeIndex*stimtdres;      // Deadtime rounded down to length of an integer number
                     //of discrete time bins

refracMult0 = 1 - stimtdres/s0;  // If y0(t) = c0*exp(-t/s0), then y0(t+stimtdres) = y0(t)*refracMult0
refracMult1 = 1 - stimtdres/s1;  // If y1(t) = c1*exp(-t/s1), then y1(t+stimtdres) = y1(t)*refracMult1

// Calculate effects of a random spike before t=0 on refractoriness and the time-warping sum at t=0
#if USE_RANDOM_VECTOR
endOfLastDeadtime = log( randNums[randBufIndex++] ) / rate[0] + deadtime;  // End of last deadtime before t=0
#else
endOfLastDeadtime = log( scop_random() ) / rate[0] + deadtime;  // End of last deadtime before t=0
#endif
refracValue0 = c0*exp(endOfLastDeadtime/s0);  // Value of first exponential in refractory function
refracValue1 = c1*exp(endOfLastDeadtime/s1);  // Value of second exponential in refractory function
Xsum = rate[0] * ( -endOfLastDeadtime + c0*s0*(exp(endOfLastDeadtime/s0)-1) + c1*s1*(exp(endOfLastDeadtime/s1)-1) );  // Value of time-warping sum
//  ^^^^ This is the "integral" of the refractory function ^^^^ (normalized by 'stimtdres')

// Calculate first interspike interval in a homogeneous, unit-rate Poisson process (normalized by 'stimtdres')
#if USE_RANDOM_VECTOR
unitRateIntrvl = -log( randNums[randBufIndex++] ) / stimtdres;
#else
unitRateIntrvl = -log( scop_random() ) / stimtdres;
#endif
/* NOTE: Both 'unitRateInterval' and 'Xsum' are divided (or normalized) by 'stimtdres' in order to reduce calculation time.
      This way we only need to divide by 'stimtdres' once per spike (when calculating 'unitRateInterval'), instead of
      multiplying by 'stimtdres' once per time bin (when calculating the new value of 'Xsum').
*/
/*
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
*/
//Loop through rate vector
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
               //Incr = (int)(meanRate * ((stimdur-stime) + (nrep-j-1)*stimdur));
               //if (Incr < 10)  Incr = 10;
               Incr = 10;
               NoutMax += Incr;
               vector_resize(spks, NoutMax);      //resize spike vector
               spktimes = vector_vec(spks);      //set spiketimes ptr to resized vector
               #if USE_RANDOM_VECTOR
               freevector(randNums);
               randNums=makevector(Incr);
               for (i=0;i<Incr;++i) randNums[i] = scop_random();
               randBufIndex = 0;
               #endif
            }
         // Next interspike "stime" in unit-rate process
            #if USE_RANDOM_VECTOR
            unitRateIntrvl = -log( randNums[randBufIndex++] ) / stimtdres;
            #else
            unitRateIntrvl = -log( scop_random() ) / stimtdres;
            #endif
            Xsum = 0.0;
         // Increase index and stime to the last time bin in the deadtime, and reset (relative) refractory function
            k += deadtimeIndex;
            stime += deadtimeRnd;
            refracValue0 = c0;
            refracValue1 = c1;
         }
      }
   } // End of rate vector loop
/* Reset index and time to begin a new repetion.
Don't just set to zero, since new repetition may start within the
deadtime following the last spike in the previous repetition. */
   stime -= stimdur/1000;
   k -= numRate;
} // End of nrep loop

//printf("Total spikes before = %d\n", Nout);
// Delete spike(s) that occur after the last repetition of the rate function ends
for (; (Nout>0)&&(spktimes[Nout-1]>stimdur); )  --Nout;

#if USE_RANDOM_VECTOR
freevector(randNums);
#endif

//printf("Total spikes = %d\n", Nout);
vector_resize(spks, Nout);
spikecount=Nout;
return Nout;
ENDVERBATIM

}
