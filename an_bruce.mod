NEURON {
  SUFFIX nothing
    }
VERBATIM

extern int vector_instance_px();
extern void vector_resize();
extern double* vector_vec();
extern void* vector_arg();
ENDVERBATIM

:  anmodel Vector Method - sout.an4( stim/*wavefile*/, tdres,  cf, spontrate, model, species, ifspike)
     VERBATIM
#include "complex.h"
#include "zbcatmodelv2.c"
#include "carneymodel.h"
   //sout.an_zbcatmodel07(tdres,cf,spont,model, species,ifspike,wavefile)
static double an_zbcatmodel07(void *vv)
{
  double *stim;      //Input stimulus vector in pascals
  double *sout;   //Output vector containing inst. rate for channel
  double tdres,cf,spont;
  double cohc,cihc;
  int species;
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
    hoc_execerror("ERROR: input syntax must be sout.an_zbcatmodel07( stim, tdres,  cf, spontrate, cihc,cohc, species,nrep)", 0);
    return 0;
  }
  //TDRES  resolution of stim vector
  //Bruce model uses seconds rather than msec
  tdres = (double)(*getarg(2));
  if (tdres > 0.01e-3 || tdres < 0.002e-3){
    //printf("Note: ZilanyBruceV2 resolution should be between 0.01ms (Fs = 100kHz) and 0.002ms (500kHz) for normal usage.\n");
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
  if ((spont<0)||(spont>150))
    {
      hoc_execerror("ERROR: spont  must be between 0 and 150 spikes/s\n",0);
      return 0;
    }
  if (spont<0.1)
    {
      //printf("an_zbcatmodel.mod setting spont rate to 0.1 to create effectively low-spont behavior\n");
      spont = 0.1;
    }
  //Model
  cihc = (double)(*getarg(5));
  if ((cihc<0)||(cihc>1))
    {
      printf("an_zbcatmodel: cihc  must be between 0 and 1\n");
      cihc=1;
    }
  cohc = (double)(*getarg(6));
  if ((cohc<0)||(cohc>1))
    {
      printf("an_zbcatmodel: cohc  must be between 0 and 1\n");
      cohc=1;
    }
  //Species
  species = (int)(*getarg(7));
  //Reps
  nrep = (int)(*getarg(8));
  //printf("AN model: Zilany and Bruce 2007()");
  return an_zilanybruce2007(tdres,cf,spont,cihc,cohc,species,nrep,stim,sout,nstim);
}


static double ANFSpikeGenerator(void *vv)
{
  void      *spks;      //We need to be able to resize this vector
  double   *spktimes,*sout;
  double   xspikes;
  double   meanRate,deadtimeRnd, endOfLastDeadtime;
  double   refracMult0, refracMult1, refracValue0, refracValue1;
  int       NspksMax, i, nspks,nsout;
  long      Nspks,j, k,nrep,deadtimeIndex;
  double   Xsum, unitRateIntrvl, stime;
  double   c0,c1,s0,s1,abs_refr,tdres,stimdur;

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

  //printf("SpikeGenerator: nsout %d,\t nreps %d,\ttdres %f\n",nsout,nrep,tdres);
  // Module parameters
  c0      = 0.5;
  s0      = 0.001;
  c1       = 0.5;
  s1       = 0.0125;
  abs_refr    = 0.00075;     //absolute refractory period in seconds

  stimdur = nsout*tdres*1000;
  // Determine the mean of the rate vector
  meanRate = 0.0;
  for (k=0; k<nsout;k++)
    if (sout[k] > 0.0)  meanRate += sout[k];
  meanRate /= nsout;

  // Create output buffer for spike times;
  NspksMax = (int)(meanRate * stimdur * nrep / 1000)+1;
  //printf("NspksMax = %d\n", NspksMax);
  if (NspksMax < 10)  NspksMax = 10;

  //printf("SpikeGenerator: nsout %d,\t nreps %d,\ttdres %f\n",nsout,nrep,tdres);
  //printf("SpikeGenerator: meanrate %g,\t NspksMax %d\t stimdur %g\n",meanRate,NspksMax,stimdur);

  spks = *((void**)(&xspikes));
  if (spks) {
    vector_resize(spks, NspksMax);
  }
  spktimes = ((double*) vector_vec(spks));      //Get array ptr to spks Vector
  printf("Spike Generator: NspksMax %d\n", NspksMax);
  Nspks = 0;

  // Calculate useful constants
  deadtimeIndex = (long) floor(abs_refr/tdres);   // Integer number of discrete time bins within deadtime
  deadtimeRnd = deadtimeIndex*tdres;      // Deadtime rounded down to length of an integer number
  //of discrete time bins

  refracMult0 = 1 - tdres/s0;  // If y0(t) = c0*exp(-t/s0), then y0(t+stimtdres) = y0(t)*refracMult0
  refracMult1 = 1 - tdres/s1;  // If y1(t) = c1*exp(-t/s1), then y1(t+stimtdres) = y1(t)*refracMult1

  // Calculate effects of a random spike before t=0 on refractoriness and the time-warping sum at t=0

  endOfLastDeadtime = log( scop_random() ) / sout[0] + abs_refr;  // End of last deadtime before t=0

  refracValue0 = c0*exp(endOfLastDeadtime/s0);  // Value of first exponential in refractory function
  refracValue1 = c1*exp(endOfLastDeadtime/s1);  // Value of second exponential in refractory function
  Xsum = sout[0] * ( -endOfLastDeadtime + c0*s0*(exp(endOfLastDeadtime/s0)-1) + c1*s1*(exp(endOfLastDeadtime/s1)-1) );  // Value of time-warping sum
  //  ^^^^ This is the "integral" of the refractory function ^^^^ (normalized by 'stimtdres')

  // Calculate first interspike interval in a homogeneous, unit-rate Poisson process (normalized by 'stimtdres')

  unitRateIntrvl = -log( scop_random() ) / tdres;

  //printf("SpikeGenerator: nsout %d,\t nreps %d,\ttdres %f\n",nsout,nrep,tdres);
  //printf("SpikeGenerator: meanrate %g,\t NspksMax %d\t stimdur %g\n",meanRate,NspksMax,stimdur);
  //printf("SpikeGenerator: endOfLastDeadtime %g\t unitRateIntrvl %g\t Xsum %g\n",endOfLastDeadtime,unitRateIntrvl,Xsum);

  //printf("Spike Generator: spktimes[0] %f\t spktimes[NspksMax-1] %f\n",spktimes[0],spktimes[NspksMax-1]);

  //Loop through rate vector
  stime = tdres;
  k=0;
  for (j=0; j<nrep; ++j)
    {
      spktimes[Nspks++] = 0.0;
      for (; (k<nsout) && (stime<stimdur); \
	   ++k, stime+=tdres, refracValue0*=refracMult0, refracValue1*=refracMult1)
	{
	  if (sout[k] > 0.0)
	    {
	      Xsum += sout[k] * (1 - refracValue0 - refracValue1);
	      //printf("SpikeGenerator: unitRateIntrvl %g\t Xsum %g\n",unitRateIntrvl,Xsum);
	      if ( Xsum >= unitRateIntrvl )
		{
		  spktimes[Nspks] = stime*1000;
		  //printf("%f\n",stime*1000);
		  if (++Nspks >= NspksMax)
		    {
		      NspksMax += 10;
		      spks = *((void**)(&xspikes));
		      if (spks) vector_resize(spks, NspksMax);
		      spktimes = ((double*) vector_vec(spks));
		    }
		  // Next interspike "stime" in unit-rate process

		  unitRateIntrvl = -log( scop_random() ) / tdres;
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
      k -= nsout;
    } // End of nrep loop

  // Delete spike(s) that occur after the last repetition of the rate function ends
  for (; (Nspks>0)&&(spktimes[Nspks-1]>stimdur); )  --Nspks;

  vector_resize(spks, Nspks);

  return Nspks;

}
ENDVERBATIM


PROCEDURE install_an_zbcatmodel_v3()
{
  VERBATIM

    install_vector_method("an_zbcatmodel07", an_zbcatmodel07);
  install_vector_method("ANFSpikeGenerator", ANFSpikeGenerator);
  ENDVERBATIM
    }

   
