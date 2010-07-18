/* This is Version 3 of the public distribution of the code for the auditory
   periphery model of:

    Zilany, M.S.A., Bruce, I.C., Nelson, P.C., and Carney, L.H. (2009). "A Phenomenological
        model of the synapse between the inner hair cell and auditory nerve : Long-term adaptation
        with power-law dynamics," Journal of the Acoustical Society of America 126(5): 2390-2412.

   Please cite this paper if you publish any research
   results obtained with this code or any modified versions of this code.

   See the file readme.txt for details of compiling and running the model.

   %%% © Muhammad S.A. Zilany (msazilany@gmail.com), Ian C. Bruce, Paul C. Nelson, and laurel H. Carney October 2008 %%%

*/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>     /* Added for MS Visual C++ compatability, by Ian Bruce, 1999 */
//#include <mex.h>
#include <time.h>


#define MAXSPIKES 1000000
#ifndef TWOPI
#define TWOPI 6.28318530717959
#endif

#ifndef __max
#define __max(a,b) (((a) > (b))? (a): (b))
#endif

#ifndef __min
#define __min(a,b) (((a) < (b))? (a): (b))
#endif

int resample(double *, double *, int , double);
double ffGn(double *yffGn, int N, double tdres, double Hinput, double mu, double sigma);

double SingleAN(double *px, double cf, int nrep, double tdres, int totalstim, double fibertype, double implnt, double *synout, double species)
{

    /*variables for the signal-path, control-path and onward */
    double *synouttmp, *sptime;

    int    i, nspikes, ipst;
    double I, spont;
    double sampFreq = 10e3; /* Sampling frequency used in the synapse */

    /* Declarations of the functions used in the program */
    double Synapse(double *, double, double, int, int, double, double, double, double *, double);
    int    SpikeGenerator(double *, double, int, int, double *);

    /* Allocate dynamic memory for the temporary variables */
    synouttmp  = makevector(totalstim * nrep);
    sptime  = makevector((long) ceil(totalstim * tdres * nrep / 0.00075));

    /* Spontaneous Rate of the fiber corresponding to Fibertype */
    if (fibertype == 1) spont = 0.1;
    if (fibertype == 2) spont = 5.0;
    if (fibertype == 3) spont = 100.0;

    /*====== Run the synapse model ======*/
    I = Synapse(px, tdres, cf, totalstim, nrep, spont, implnt, sampFreq, synouttmp, species);

    /* Wrapping up the unfolded (due to no. of repetitions) Synapse Output */
    for (i = 0; i < I ; i++) {
        ipst = (int)(fmod(i, totalstim));
        synout[ipst] = synout[ipst] + synouttmp[i] / nrep;
    };
    /*======  Spike Generations ======*/
    /*
    nspikes = SpikeGenerator(synouttmp, tdres, totalstim, nrep, sptime);
    for(i = 0; i < nspikes; i++)
    {
    ipst = (int) (fmod(sptime[i],tdres*totalstim) / tdres);
    psth[ipst] = psth[ipst] + 1;
    };
    */

    /* Freeing dynamic memory allocated earlier */
    freevector(sptime); freevector(synouttmp);
    return I;
} /* End of the SingleAN function */

/*-------------------------------------------------------------------------------
  Synapse model: if the time resolution is not small enough, the concentration of
  the immediate pool could be as low as negative, at this time there is an alert
  message print out and the concentration is set at saturated level
 -------------------------------------------------------------------------------*/
double Synapse(double *ihcout, double tdres, double cf, int totalstim, int nrep, double spont, double implnt, double sampFreq, double *synouttmp, double species)
{
    /* Initalize Variables */
    int z, b, n;
    int resamp = (int) ceil(1 / (tdres * sampFreq)); //sampFreq is not the same as 1/tdres
    double incr = 0.0;
    int delaypoint = floor(7500 / (cf / 1e3));   //is this for Cat only

    double alpha1, beta1, I1, alpha2, beta2, I2, binwidth;
    int    k, j, indx, i;
    double synstrength, synslope, CI, CL, PG, CG, VL, PL, VI;
    double cf_factor, PImax, kslope, Ass, Asp, TauR, TauST, Ar_Ast, PTS, Aon, AR, AST, Prest, gamma1, gamma2, k1, k2;
    double VI0, VI1, alpha, beta, theta1, theta2, theta3, vsat, tmpst, tmp, PPI, CIlast, temp;

    double *sout1, *sout2, *synSampOut, *powerLawIn, *exponOut, *TmpSyn;
    double *m1, *m2, *m3, *m4, *m5;
    double *n1, *n2, *n3;

    //    mxArray *randInputArray[4], *randOutputArray[1];
    double *randNums;

    // mxArray *IhcInputArray[3], *IhcOutputArray[1];
    double *sampIHC, *ihcDims;

    exponOut = makevector((long) ceil(totalstim * nrep));
    powerLawIn = makevector((long) ceil(totalstim * nrep + 3 * delaypoint));
    sout1 = makevector((long) ceil((totalstim * nrep + 2 * delaypoint) * tdres * sampFreq));
    sout2 = makevector((long) ceil((totalstim * nrep + 2 * delaypoint) * tdres * sampFreq));
    synSampOut  = makevector((long) ceil((totalstim * nrep + 2 * delaypoint) * tdres * sampFreq));
    TmpSyn  = makevector((long) ceil(totalstim * nrep + 2 * delaypoint));

    m1 = makevector((long) ceil((totalstim * nrep + 2 * delaypoint) * tdres * sampFreq));
    m2 = makevector((long) ceil((totalstim * nrep + 2 * delaypoint) * tdres * sampFreq));
    m3  = makevector((long) ceil((totalstim * nrep + 2 * delaypoint) * tdres * sampFreq));
    m4 = makevector((long) ceil((totalstim * nrep + 2 * delaypoint) * tdres * sampFreq));
    m5  = makevector((long) ceil((totalstim * nrep + 2 * delaypoint) * tdres * sampFreq));

    n1 = makevector((long) ceil((totalstim * nrep + 2 * delaypoint) * tdres * sampFreq));
    n2 = makevector((long) ceil((totalstim * nrep + 2 * delaypoint) * tdres * sampFreq));
    n3 = makevector((long) ceil((totalstim * nrep + 2 * delaypoint) * tdres * sampFreq));

    /*----------------------------------------------------------*/
    /*------- Parameters of the Power-law function -------------*/
    /*----------------------------------------------------------*/
    binwidth = 1 / sampFreq;
    alpha1 = 5e-6 * 100e3; beta1 = 5e-4; I1 = 0;
    alpha2 = 1e-2 * 100e3; beta2 = 1e-1; I2 = 0;

    /*----------------------------------------------------------*/
    /*------- Generating a random sequence ---------------------*/
    /*----------------------------------------------------------
      randInputArray[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
      *mxGetPr(randInputArray[0])= ceil((totalstim*nrep+2*delaypoint)*tdres*sampFreq);
      randInputArray[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
      *mxGetPr(randInputArray[1])= 1/sampFreq;
      randInputArray[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
      *mxGetPr(randInputArray[2])= 0.9; // Hurst index
      randInputArray[3] = mxCreateDoubleMatrix(1, 1, mxREAL);
      *mxGetPr(randInputArray[3])= spont;

      mexCallMATLAB(1, randOutputArray, 4, randInputArray, "ffGn");
      randNums = mxGetPr(randOutputArray[0]);
    */
    randNums = makevector((long) ceil((totalstim * nrep + 2 * delaypoint) * tdres * sampFreq));
    if(!ffGn(randNums, ceil((totalstim*nrep + 2*delaypoint)*tdres*sampFreq), 1000 / sampFreq, 0.9, spont, -1)){
      hoc_execerror("Synapse: error calling ffGn",0);
      return 0;
    }

    /*----------------------------------------------------------*/
    /*----- Double Exponential Adaptation ----------------------*/
    /*----------------------------------------------------------*/
    if (spont == 100) cf_factor = __min(800, pow(10, 0.29 * cf / 1e3 + 0.7));
    if (spont == 5)   cf_factor = __min(50, 2.5e-4 * cf * 4 + 0.2);
    if (spont == 0.1) cf_factor = __min(1.0, 2.5e-4 * cf * 0.1 + 0.15);

    PImax  = 0.6;                /* PI2 : Maximum of the PI(PI at steady state) */
    kslope = (1 + 50.0) / (5 + 50.0) * cf_factor * 20.0 * PImax;

    Ass    = 300 * TWOPI / 2 * (1 + cf / 10e3);    /* Steady State Firing Rate eq.10 */
    if (implnt == 1) Asp = spont * 5;   /* Spontaneous Firing Rate if actual implementation */
    if (implnt == 0) Asp = spont * 4.1; /* Spontaneous Firing Rate if approximate implementation */
    TauR   = 2e-3;               /* Rapid Time Constant eq.10 */
    TauST  = 60e-3;              /* Short Time Constant eq.10 */
    Ar_Ast = 6;                  /* Ratio of Ar/Ast */
    PTS    = 3;                  /* Peak to Steady State Ratio, characteristic of PSTH */

    /* now get the other parameters */
    Aon    = PTS * Ass;                        /* Onset rate = Ass+Ar+Ast eq.10 */
    AR     = (Aon - Ass) * Ar_Ast / (1 + Ar_Ast);      /* Rapid component magnitude: eq.10 */
    AST    = Aon - Ass - AR;                   /* Short time component: eq.10 */
    Prest  = PImax / Aon * Asp;                /* eq.A15 */
    CG  = (Asp * (Aon - Asp)) / (Aon * Prest * (1 - Asp / Ass));    /* eq.A16 */
    gamma1 = CG / Asp;                         /* eq.A19 */
    gamma2 = CG / Ass;                         /* eq.A20 */
    k1     = -1 / TauR;                        /* eq.8 & eq.10 */
    k2     = -1 / TauST;                       /* eq.8 & eq.10 */
    /* eq.A21 & eq.A22 */
    VI0    = (1 - PImax / Prest) / (gamma1 * (AR * (k1 - k2) / CG / PImax + k2 / Prest / gamma1 - k2 / PImax / gamma2));
    VI1    = (1 - PImax / Prest) / (gamma1 * (AST * (k2 - k1) / CG / PImax + k1 / Prest / gamma1 - k1 / PImax / gamma2));
    VI  = (VI0 + VI1) / 2;
    alpha  = gamma2 / k1 / k2;   /* eq.A23,eq.A24 or eq.7 */
    beta   = -(k1 + k2) * alpha; /* eq.A23 or eq.7 */
    theta1 = alpha * PImax / VI;
    theta2 = VI / PImax;
    theta3 = gamma2 - 1 / PImax;

    PL  = ((beta - theta2 * theta3) / theta1 - 1) * PImax;  /* eq.4' */
    PG  = 1 / (theta3 - 1 / PL);                  /* eq.5' */
    VL  = theta1 * PL * PG;                       /* eq.3' */
    CI  = Asp / Prest;                            /* CI at rest, from eq.A3,eq.A12 */
    CL  = CI * (Prest + PL) / PL;                 /* CL at rest, from eq.1 */

    if (kslope >= 0)  vsat = kslope + Prest;
    tmpst  = log(2) * vsat / Prest;
    if (tmpst < 400) synstrength = log(exp(tmpst) - 1);
    else synstrength = tmpst;
    synslope = Prest / log(2) * synstrength;

    k = 0;
    for (indx = 0; indx < totalstim*nrep; ++indx) {
        tmp = synstrength * (ihcout[indx]);
        if (tmp < 400) tmp = log(1 + exp(tmp));
        PPI = synslope / synstrength * tmp;

        CIlast = CI;
        CI = CI + (tdres / VI) * (-PPI * CI + PL * (CL - CI));
        CL = CL + (tdres / VL) * (-PL * (CL - CIlast) + PG * (CG - CL));
        if (CI < 0) {
            temp = 1 / PG + 1 / PL + 1 / PPI;
            CI = CG / (PPI * temp);
            CL = CI * (PPI + PL) / PL;
        };
        exponOut[k] = CI * PPI;
        k = k + 1;
    }
    for (k = 0; k < delaypoint; k++)
        powerLawIn[k] = exponOut[0];
    for (k = delaypoint; k < totalstim*nrep + delaypoint; k++)
        powerLawIn[k] = exponOut[k-delaypoint];
    for (k = totalstim * nrep + delaypoint; k < totalstim*nrep + 3*delaypoint; k++)
        powerLawIn[k] = powerLawIn[k-1];

    /*----------------------------------------------------------*/
    /*------ Downsampling to sampFreq (Low) sampling rate ------*/
    /*----------------------------------------------------------*/
    /*
      IhcInputArray[0] = mxCreateDoubleMatrix(1, k, mxREAL);
      ihcDims = mxGetPr(IhcInputArray[0]);
      for (i=0;i<k;++i)
      ihcDims[i] = powerLawIn[i];
      IhcInputArray[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
      *mxGetPr(IhcInputArray[1])= 1;
      IhcInputArray[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
      *mxGetPr(IhcInputArray[2])= resamp;
      mexCallMATLAB(1, IhcOutputArray, 3, IhcInputArray, "resample");
      sampIHC = mxGetPr(IhcOutputArray[0]);
    */



    /* Simple resampling routine from src/ivoc/ivocvect.cpp */
    int len=0;
    printf("Synapse: calling resample(%d,XX,%d,%g)\n", &powerLawIn[0],  k, 1.0/resamp);
    len=resample(powerLawIn, sampIHC, k, 1/resamp);
    if ( len==0 || ( sampIHC == NULL) ) {
        printf("Synapse: ivoc_resample return error, copying powerLawIn");
	len = k;
	if(sampIHC) freevector(sampIHC); 
	  sampIHC = makevector(k);
        for (indx = 0;indx < k;indx++) sampIHC[indx] = powerLawIn[indx];
    }
    printf("k %d len %d \t old len %g\n",k,len, (totalstim*nrep + 2*delaypoint)*tdres*sampFreq);


    /*** libresample implementation ***/
    //n = (int)(k/resamp);
    //sampIHC = makevector(n);
    //n=resample(powerLawIn,sampIHC,k,n);

    freevector(powerLawIn); freevector(exponOut);
    /*----------------------------------------------------------*/
    /*----- Running Power-law Adaptation -----------------------*/
    /*----------------------------------------------------------*/
    k = 0;
    for (indx = 0; indx < floor((totalstim*nrep + 2*delaypoint)*tdres*sampFreq); indx++) {
        sout1[k]  = __max(0, sampIHC[indx] + randNums[indx] - alpha1 * I1);
        /* sout1[k]  = __max( 0, sampIHC[indx] - alpha1*I1); */   /* No fGn condition */
        sout2[k]  = __max(0, sampIHC[indx] - alpha2 * I2);

        if (implnt == 1) { /* ACTUAL Implementation */
            I1 = 0; I2 = 0;
            for (j = 0; j < k + 1; ++j) {
                I1 += (sout1[j]) * binwidth / ((k - j) * binwidth + beta1);
                I2 += (sout2[j]) * binwidth / ((k - j) * binwidth + beta2);
            }
        } /* end of actual */

        if (implnt == 0) { /* APPROXIMATE Implementation */
            if (k == 0) {
                n1[k] = 1.0e-3 * sout2[k];
                n2[k] = n1[k]; n3[0] = n2[k];
            } else if (k == 1) {
                n1[k] = 1.992127932802320 * n1[k-1] + 1.0e-3 * (sout2[k] - 0.994466986569624 * sout2[k-1]);
                n2[k] = 1.999195329360981 * n2[k-1] + n1[k] - 1.997855276593802 * n1[k-1];
                n3[k] = -0.798261718183851 * n3[k-1] + n2[k] + 0.798261718184977 * n2[k-1];
            } else {
                n1[k] = 1.992127932802320 * n1[k-1] - 0.992140616993846 * n1[k-2] + 1.0e-3 * (sout2[k] - 0.994466986569624 * sout2[k-1] + 0.000000000002347 * sout2[k-2]);
                n2[k] = 1.999195329360981 * n2[k-1] - 0.999195402928777 * n2[k-2] + n1[k] - 1.997855276593802 * n1[k-1] + 0.997855827934345 * n1[k-2];
                n3[k] = -0.798261718183851 * n3[k-1] - 0.199131619873480 * n3[k-2] + n2[k] + 0.798261718184977 * n2[k-1] + 0.199131619874064 * n2[k-2];
            }
            I2 = n3[k];

            if (k == 0) {
                m1[k] = 0.2 * sout1[k];
                m2[k] = m1[k]; m3[k] = m2[k];
                m4[k] = m3[k]; m5[k] = m4[k];
            } else if (k == 1) {
                m1[k] = 0.491115852967412 * m1[k-1] + 0.2 * (sout1[k] - 0.173492003319319 * sout1[k-1]);
                m2[k] = 1.084520302502860 * m2[k-1] + m1[k] - 0.803462163297112 * m1[k-1];
                m3[k] = 1.588427084535629 * m3[k-1] + m2[k] - 1.416084732997016 * m2[k-1];
                m4[k] = 1.886287488516458 * m4[k-1] + m3[k] - 1.830362725074550 * m3[k-1];
                m5[k] = 1.989549282714008 * m5[k-1] + m4[k] - 1.983165053215032 * m4[k-1];
            } else {
                m1[k] = 0.491115852967412 * m1[k-1] - 0.055050209956838 * m1[k-2] + 0.2 * (sout1[k] - 0.173492003319319 * sout1[k-1] + 0.000000172983796 * sout1[k-2]);
                m2[k] = 1.084520302502860 * m2[k-1] - 0.288760329320566 * m2[k-2] + m1[k] - 0.803462163297112 * m1[k-1] + 0.154962026341513 * m1[k-2];
                m3[k] = 1.588427084535629 * m3[k-1] - 0.628138993662508 * m3[k-2] + m2[k] - 1.416084732997016 * m2[k-1] + 0.496615555008723 * m2[k-2];
                m4[k] = 1.886287488516458 * m4[k-1] - 0.888972875389923 * m4[k-2] + m3[k] - 1.830362725074550 * m3[k-1] + 0.836399964176882 * m3[k-2];
                m5[k] = 1.989549282714008 * m5[k-1] - 0.989558985673023 * m5[k-2] + m4[k] - 1.983165053215032 * m4[k-1] + 0.983193027347456 * m4[k-2];
            }
            I1 = m5[k];
        } /* end of approximate implementation */

        synSampOut[k] = sout1[k] + sout2[k];
        k = k + 1;
    }   /* end of all samples */

    freevector(sout1); freevector(sout2);
    freevector(m1); freevector(m2); freevector(m3); freevector(m4); freevector(m5); freevector(n1); freevector(n2); freevector(n3);
    /*----------------------------------------------------------*/
    /*----- Upsampling to original (High 100 kHz) sampling rate --------*/
    /*----------------------------------------------------------*/
    for (z = 0; z < k - 1; ++z) {
        incr = (synSampOut[z+1] - synSampOut[z]) / resamp;
        for (b = 0; b < resamp; ++b) {
            TmpSyn[z*resamp+b] = synSampOut[z] + b * incr;
        }
    }
    for (i = 0;i < totalstim*nrep;++i)
        synouttmp[i] = TmpSyn[i+delaypoint];

    freevector(synSampOut); freevector(TmpSyn);

    return((long) ceil(totalstim*nrep));
}

/* ------------------------------------------------------------------------------------ */
/* Pass the output of Synapse model through the Spike Generator */

/* The spike generator now uses a method coded up by B. Scott Jackson (bsj22@cornell.edu)
   Scott's original code is available from Laurel Carney's web site at:
   http://www.urmc.rochester.edu/smd/Nanat/faculty-research/lab-pages/LaurelCarney/auditory-models.cfm
*/

int SpikeGenerator(double *synouttmp, double tdres, int totalstim, int nrep, double *sptime)
{
    double  c0, s0, c1, s1, dead;
    int     nspikes, k, NoutMax, Nout, deadtimeIndex, randBufIndex;
    double deadtimeRnd, endOfLastDeadtime, refracMult0, refracMult1, refracValue0, refracValue1;
    double Xsum, unitRateIntrvl, countTime, DT;

    double *randNums, *randDims;

    c0      = 0.5;
    s0      = 0.001;
    c1      = 0.5;
    s1      = 0.0125;
    dead    = 0.00075;

    DT = totalstim * tdres * nrep;  /* Total duration of the rate function */
    Nout = 0;
    NoutMax = (long) ceil(totalstim * nrep * tdres / dead);
    randNums = makevector(NoutMax + 1);
    /*
      randInputArray[0] = mxCreateDoubleMatrix(1, 2, mxREAL);
    randDims = mxGetPr(randInputArray[0]);
    randDims[0] = 1;
    randDims[1] = NoutMax+1;
        mexCallMATLAB(1, randOutputArray, 1, randInputArray, "rand");
    randNums = mxGetPr(randOutputArray[0]);
    randBufIndex = 0;
    */
    /* Calculate useful constants */
    deadtimeIndex = (long) floor(dead / tdres);  /* Integer number of discrete time bins within deadtime */
    deadtimeRnd = deadtimeIndex * tdres;   /* Deadtime rounded down to length of an integer number of discrete time bins */

    refracMult0 = 1 - tdres / s0;  /* If y0(t) = c0*exp(-t/s0), then y0(t+tdres) = y0(t)*refracMult0 */
    refracMult1 = 1 - tdres / s1;  /* If y1(t) = c1*exp(-t/s1), then y1(t+tdres) = y1(t)*refracMult1 */

    /* Calculate effects of a random spike before t=0 on refractoriness and the time-warping sum at t=0 */
    endOfLastDeadtime = __max(0, log(randNums[randBufIndex++]) / synouttmp[0] + dead); /* End of last deadtime before t=0 */
    refracValue0 = c0 * exp(endOfLastDeadtime / s0); /* Value of first exponential in refractory function */
    refracValue1 = c1 * exp(endOfLastDeadtime / s1); /* Value of second exponential in refractory function */
    Xsum = synouttmp[0] * (-endOfLastDeadtime + c0 * s0 * (exp(endOfLastDeadtime / s0) - 1) + c1 * s1 * (exp(endOfLastDeadtime / s1) - 1));
    /* Value of time-warping sum */
    /*  ^^^^ This is the "integral" of the refractory function ^^^^ (normalized by 'tdres') */

    /* Calculate first interspike interval in a homogeneous, unit-rate Poisson process (normalized by 'tdres') */
    unitRateIntrvl = -log(randNums[randBufIndex++]) / tdres;
    /* NOTE: Both 'unitRateInterval' and 'Xsum' are divided (or normalized) by 'tdres' in order to reduce calculation time.
       This way we only need to divide by 'tdres' once per spike (when calculating 'unitRateInterval'), instead of
       multiplying by 'tdres' once per time bin (when calculating the new value of 'Xsum').                         */

    countTime = tdres;
    for (k = 0; (k < totalstim*nrep) && (countTime < DT); ++k, countTime += tdres, refracValue0 *= refracMult0, refracValue1 *= refracMult1) { /* Loop through rate vector */
        if (synouttmp[k] > 0) { /* Nothing to do for non-positive rates, i.e. Xsum += 0 for non-positive rates. */
            Xsum += synouttmp[k] * (1 - refracValue0 - refracValue1);  /* Add synout*(refractory value) to time-warping sum */

            if (Xsum >= unitRateIntrvl) {  /* Spike occurs when time-warping sum exceeds interspike "time" in unit-rate process */
                sptime[Nout] = countTime; Nout = Nout + 1;
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


    freevector(randNums);
    nspikes = Nout;  /* Number of spikes that occurred. */
    return(nspikes);
}


/* This is Version 1 of the public distribution of the code for the auditory
   periphery model of:

        Zilany, M. S. A. and Bruce, I. C. (2007). "Representation of the vowel
        /eh/ in normal and impaired auditory nerve fibers: Model predictions of
        responses in cats," Journal of the Acoustical Society of America
        122(1):402-417.

        Zilany, M. S. A. and Bruce, I. C. (2006). "Modeling auditory-nerve
        responses for high sound pressure levels in the normal and impaired
        auditory periphery," Journal of the Acoustical Society of
        America 120(3):1446-1466.

   Please cite these papers if you publish any research
   results obtained with this code or any modified versions of this code.

   See the file readme.txt for details of compiling and running the model.

   %%% (C) Ian C. Bruce (ibruce@ieee.org) and M. S. Arefeen Zilany, June 2006 %%%

*/

void IHCAN(double *px, double cf, int nrep, double tdres, int totalstim, double cohc, double cihc, double *ihcout, double species)
{

    /*variables for middle-ear model */
    double megainmax = 43;
    double *mey1, *mey2, *mey3, meout, c1filterouttmp, c2filterouttmp, c1vihctmp, c2vihctmp;
    double fp, C, m11, m12, m21, m22, m23, m24, m25, m26, m31, m32, m33, m34, m35, m36;

    /*variables for the signal-path, control-path and onward */
    double *ihcouttmp, *tmpgain;
    int    grd;

    double bmplace, centerfreq, gain, taubm, ratiowb, bmTaubm, fcohc, TauWBMax, TauWBMin, tauwb;
    double Taumin[1], Taumax[1], bmTaumin[1], bmTaumax[1], ratiobm[1], lasttmpgain, wbgain, ohcasym, ihcasym, delay;
    int    i, n, delaypoint, grdelay[1], bmorder, wborder;
    double wbout1, wbout, ohcnonlinout, ohcout, tmptauc1, tauc1, rsigma, wb_gain;

    /* Declarations of the functions used in the program */
    double C1ChirpFilt(double, double, double, int, double, double);
    double C2ChirpFilt(double, double, double, int, double, double);
    double WbGammaTone(double, double, double, int, double, double, int);

    double Get_tauwb(double, int, double *, double *);
    double Get_taubm(double, double, double *, double *, double *);
    double gain_groupdelay(double, double, double, double, int *);
    double delay_cat(double cf);

    double OhcLowPass(double, double, double, int, double, int);
    double IhcLowPass(double, double, double, int, double, int);
    double Boltzman(double, double, double, double, double);
    double NLafterohc(double, double, double, double);
    double ControlSignal(double, double, double, double, double);

    double NLogarithm(double, double, double, double);

    /* Allocate dynamic memory for the temporary variables */
    ihcouttmp  = makevector(totalstim * nrep);

    mey1 = makevector(totalstim);
    mey2 = makevector(totalstim);
    mey3 = makevector(totalstim);

    tmpgain = makevector(totalstim);

    /** Calculate the location on basilar membrane from CF */

    bmplace = 11.9 * log10(0.80 + cf / 456.0);

    /** Calculate the center frequency for the control-path wideband filter
    from the location on basilar membrane */

    centerfreq = 456.0 * (pow(10, (bmplace + 1.2) / 11.9) - 0.80); /* shift the center freq */

    /*==================================================================*/
    /*====== Parameters for the gain ===========*/
    gain = 52 / 2 * (tanh(2.2 * log10(cf / 0.6e3) + 0.15) + 1);
    /*gain = 52/2*(tanh(2.2*log10(cf/1e3)+0.15)+1);*/
    if (gain > 60) gain = 60;
    if (gain < 15) gain = 15;
    /*====== Parameters for the control-path wideband filter =======*/
    bmorder = 3;
    Get_tauwb(cf, bmorder, Taumax, Taumin);
    taubm   = cohc * (Taumax[0] - Taumin[0]) + Taumin[0];
    ratiowb = Taumin[0] / Taumax[0];
    /*====== Parameters for the signal-path C1 filter ======*/
    Get_taubm(cf, Taumax[0], bmTaumax, bmTaumin, ratiobm);
    bmTaubm  = cohc * (bmTaumax[0] - bmTaumin[0]) + bmTaumin[0];
    fcohc    = bmTaumax[0] / bmTaubm;
    /*====== Parameters for the control-path wideband filter =======*/
    wborder  = 3;
    TauWBMax = Taumin[0] + 0.2 * (Taumax[0] - Taumin[0]);
    TauWBMin = TauWBMax / Taumax[0] * Taumin[0];
    tauwb    = TauWBMax + (bmTaubm - bmTaumax[0]) * (TauWBMax - TauWBMin) / (bmTaumax[0] - bmTaumin[0]);

    wbgain = gain_groupdelay(tdres, centerfreq, cf, tauwb, grdelay);
    tmpgain[0]   = wbgain;
    lasttmpgain  = wbgain;
    /*===============================================================*/
    /* Nonlinear asymmetry of OHC function and IHC C1 transduction function*/
    ohcasym  = 7.0;
    ihcasym  = 3.0;
    /*===============================================================*/
    /*===============================================================*/
    /* Prewarping and related constants for the middle ear */
    fp = 1e3;  /* prewarping frequency 1 kHz */
    C  = TWOPI * fp / tan(TWOPI / 2 * fp * tdres);
    m11 = C / (C + 693.48);                    m12 = (693.48 - C) / C;
    m21 = 1 / (pow(C, 2) + 11053 * C + 1.163e8);  m22 = -2 * pow(C, 2) + 2.326e8;    m23 = pow(C, 2) - 11053 * C + 1.163e8;
    m24 = pow(C, 2) + 1356.3 * C + 7.4417e8;    m25 = -2 * pow(C, 2) + 14.8834e8;  m26 = pow(C, 2) - 1356.3 * C + 7.4417e8;
    m31 = 1 / (pow(C, 2) + 4620 * C + 909059944); m32 = -2 * pow(C, 2) + 2 * 909059944; m33 = pow(C, 2) - 4620 * C + 909059944;
    m34 = 5.7585e5 * C + 7.1665e7;             m35 = 14.333e7;                 m36 = 7.1665e7 - 5.7585e5 * C;

    for (n = 0;n < totalstim;n++) { /* Start of the loop */
        if (n == 0) { /* Start of the middle-ear filtering section  */
            mey1[0]  = m11 * px[0];
            mey2[0]  = mey1[0] * m24 * m21;
            mey3[0]  = mey2[0] * m34 * m31;
            meout = mey3[0] / megainmax ;
        } else if (n == 1) {
            mey1[1]  = m11 * (-m12 * mey1[0] + px[1]       - px[0]);
            mey2[1]  = m21 * (-m22 * mey2[0] + m24 * mey1[1] + m25 * mey1[0]);
            mey3[1]  = m31 * (-m32 * mey3[0] + m34 * mey2[1] + m35 * mey2[0]);
            meout = mey3[1] / megainmax;
        } else {
            mey1[n]  = m11 * (-m12 * mey1[n-1]  + px[n]         - px[n-1]);
            mey2[n]  = m21 * (-m22 * mey2[n-1] - m23 * mey2[n-2] + m24 * mey1[n] + m25 * mey1[n-1] + m26 * mey1[n-2]);
            mey3[n]  = m31 * (-m32 * mey3[n-1] - m33 * mey3[n-2] + m34 * mey2[n] + m35 * mey2[n-1] + m36 * mey2[n-2]);
            meout = mey3[n] / megainmax;
        }
        ;  /* End of the middle-ear filtering section */

        /* Control-path filter */

        wbout1 = WbGammaTone(meout, tdres, centerfreq, n, tauwb, wbgain, wborder);
        wbout  = pow((tauwb / TauWBMax), wborder) * wbout1 * 10e3 * __max(1, cf / 5e3);

        ohcnonlinout = Boltzman(wbout, ohcasym, 12.0, 5.0, 5.0); /* pass the control signal through OHC Nonlinear Function */
        ohcout = OhcLowPass(ohcnonlinout, tdres, 600, n, 1.0, 2);/* lowpass filtering after the OHC nonlinearity */

        tmptauc1 = NLafterohc(ohcout, bmTaumin[0], bmTaumax[0], ohcasym); /* nonlinear function after OHC low-pass filter */
        tauc1    = cohc * (tmptauc1 - bmTaumin[0]) + bmTaumin[0];  /* time -constant for the signal-path C1 filter */
        rsigma   = 1 / tauc1 - 1 / bmTaumax[0]; /* shift of the location of poles of the C1 filter from the initial positions */

        if (1 / tauc1 < 0.0) printf("IHCAN: The poles are in the right-half plane; system is unstable.\n");

        tauwb = TauWBMax + (tauc1 - bmTaumax[0]) * (TauWBMax - TauWBMin) / (bmTaumax[0] - bmTaumin[0]);

        wb_gain = gain_groupdelay(tdres, centerfreq, cf, tauwb, grdelay);

        grd = grdelay[0];

        if ((grd + n) < totalstim)
            tmpgain[grd+n] = wb_gain;

        if (tmpgain[n] == 0)
            tmpgain[n] = lasttmpgain;

        wbgain      = tmpgain[n];
        lasttmpgain = wbgain;

        /*====== Signal-path C1 filter ======*/

        c1filterouttmp = C1ChirpFilt(meout, tdres, cf, n, bmTaumax[0], rsigma); /* C1 filter output */


        /*====== Parallel-path C2 filter ======*/

        c2filterouttmp  = C2ChirpFilt(meout, tdres, cf, n, bmTaumax[0], 1 / ratiobm[0]); /* parallel-filter output*/

        /*=== Run the inner hair cell (IHC) section: NL function and then lowpass filtering ===*/

        c1vihctmp  = NLogarithm(cihc * c1filterouttmp, 0.1, ihcasym, cf);

        c2vihctmp = -NLogarithm(c2filterouttmp * fabs(c2filterouttmp) * cf / 10 * cf / 2e3, 0.2, 1.0, cf); /* C2 transduction output */

        ihcouttmp[n] = IhcLowPass(c1vihctmp + c2vihctmp, tdres, 3000, n, 1.0, 7);
    }
    ;  /* End of the loop */

    /* Stretched out the IHC output according to nrep (number of repetitions) */

    for (i = 0;i < totalstim*nrep;i++) {
        ihcouttmp[i] = ihcouttmp[(int)(fmod(i,totalstim))];
    };
    /* Adjust total path delay to IHC output signal */
    delay      = delay_cat(cf);
    delaypoint = __max(0, (int) ceil(delay / tdres));

    for (i = delaypoint;i < totalstim*nrep;i++) {
        ihcout[i] = ihcouttmp[i - delaypoint];
    };

    /* Freeing dynamic memory allocated earlier */

    freevector(ihcouttmp);
    freevector(mey1); freevector(mey2); freevector(mey3);
    freevector(tmpgain);

} /* End of the SingleAN function */
/* -------------------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------------------- */
/** Get TauMax, TauMin for the tuning filter. The TauMax is determined by the bandwidth/Q10
    of the tuning filter at low level. The TauMin is determined by the gain change between high
    and low level */

double Get_tauwb(double cf, int order, double *taumax, double *taumin)
{
    double Q10, bw, gain, ratio;

    gain = 52 / 2 * (tanh(2.2 * log10(cf / 0.6e3) + 0.15) + 1);
    /*gain = 52/2*(tanh(2.2*log10(cf/1e3)+0.15)+1);*/

    if (gain > 60) gain = 60;
    if (gain < 15) gain = 15;

    ratio = pow(10, (-gain / (20.0 * order)));  /* ratio of TauMin/TauMax according to the gain, order */

    /*Q10 = pow(10,0.4708*log10(cf/1e3)+0.5469);  /* 75th percentile */
    Q10 = pow(10, 0.4708 * log10(cf / 1e3) + 0.4664); /* 50th percentile */
    /*Q10 = pow(10,0.4708*log10(cf/1e3)+0.3934);  /* 25th percentile */

    bw     = cf / Q10;
    taumax[0] = 2.0 / (TWOPI * bw);

    taumin[0]   = taumax[0] * ratio;

    return 0;
}
/* --------------------------------------------------------------------------- */
double Get_taubm(double cf, double taumax, double *bmTaumax, double *bmTaumin, double *ratio)
{
    double gain, factor, bwfactor;

    gain = 52 / 2 * (tanh(2.2 * log10(cf / 0.6e3) + 0.15) + 1);
    /*gain = 52/2*(tanh(2.2*log10(cf/1e3)+0.15)+1);*/

    if (gain > 60) gain = 60;
    if (gain < 15) gain = 15;

    bwfactor = 0.7;
    factor   = 2.5;

    ratio[0]  = pow(10, (-gain / (20.0 * factor)));

    bmTaumax[0] = taumax / bwfactor;
    bmTaumin[0] = bmTaumax[0] * ratio[0];
    return 0;
}
/* ---------------------------------------------------------------------------- */
/** Pass the signal through the signal-path C1 Tenth Order Nonlinear Chirp-Gammatone Filter */

double C1ChirpFilt(double x, double tdres, double cf, int n, double taumax, double rsigma)
{
    static double C1gain_norm, C1initphase;
    static double C1input[12][4], C1output[12][4];

    double ipw, ipb, rpa, pzero, rzero;
    double sigma0, fs_bilinear, CF, norm_gain, phase, c1filterout;
    int i, r, order_of_pole, half_order_pole, order_of_zero;
    double temp, dy, preal, pimg;

    COMPLEX p[11];

    /* Defining initial locations of the poles and zeros */
    /*======== setup the locations of poles and zeros =======*/
    sigma0 = 1 / taumax;
    ipw    = 1.01 * cf * TWOPI - 50;
    ipb    = 0.2343 * TWOPI * cf - 1104;
    rpa    = pow(10, log10(cf) * 0.9 + 0.55) + 2000;
    pzero  = pow(10, log10(cf) * 0.7 + 1.6) + 500;

    /*===============================================================*/

    order_of_pole    = 10;
    half_order_pole  = order_of_pole / 2;
    order_of_zero    = half_order_pole;

    fs_bilinear = TWOPI * cf / tan(TWOPI * cf * tdres / 2);
    rzero       = -pzero;
    CF          = TWOPI * cf;

    if (n == 0) {
        p[1].x = -sigma0;

        p[1].y = ipw;

        p[5].x = p[1].x - rpa; p[5].y = p[1].y - ipb;

        p[3].x = (p[1].x + p[5].x) * 0.5; p[3].y = (p[1].y + p[5].y) * 0.5;

        p[2]   = compconj(p[1]);    p[4] = compconj(p[3]); p[6] = compconj(p[5]);

        p[7]   = p[1]; p[8] = p[2]; p[9] = p[5]; p[10] = p[6];

        C1initphase = 0.0;
        for (i = 1;i <= half_order_pole;i++) {
            preal     = p[i*2-1].x;
            pimg      = p[i*2-1].y;
            C1initphase = C1initphase + atan(CF / (-rzero)) - atan((CF - pimg) / (-preal)) - atan((CF + pimg) / (-preal));
        };

        /*===================== Initialize C1input & C1output =====================*/

        for (i = 1;i <= (half_order_pole + 1);i++) {
            C1input[i][3] = 0;
            C1input[i][2] = 0;
            C1input[i][1] = 0;
            C1output[i][3] = 0;
            C1output[i][2] = 0;
            C1output[i][1] = 0;
        }

        /*===================== normalize the gain =====================*/

        C1gain_norm = 1.0;
        for (r = 1; r <= order_of_pole; r++)
            C1gain_norm = C1gain_norm * (pow((CF - p[r].y), 2) + p[r].x * p[r].x);

    };

    norm_gain = sqrt(C1gain_norm) / pow(sqrt(CF * CF + rzero * rzero), order_of_zero);

    p[1].x = -sigma0 - rsigma;

    if (p[1].x > 0.0) printf("C1ChirpFilt: The system becomes unstable.\n");

    p[1].y = ipw;

    p[5].x = p[1].x - rpa; p[5].y = p[1].y - ipb;

    p[3].x = (p[1].x + p[5].x) * 0.5; p[3].y = (p[1].y + p[5].y) * 0.5;

    p[2] = compconj(p[1]); p[4] = compconj(p[3]); p[6] = compconj(p[5]);

    p[7] = p[1]; p[8] = p[2]; p[9] = p[5]; p[10] = p[6];

    phase = 0.0;
    for (i = 1;i <= half_order_pole;i++) {
        preal = p[i*2-1].x;
        pimg  = p[i*2-1].y;
        phase = phase - atan((CF - pimg) / (-preal)) - atan((CF + pimg) / (-preal));
    };

    rzero = -CF / tan((C1initphase - phase) / order_of_zero);

    if (rzero > 0.0) printf("C1ChirpFilt: The zeros are in the right-half plane.\n");

    /*%==================================================  */
    /*each loop below is for a pair of poles and one zero */
    /*%      time loop begins here                         */
    /*%==================================================  */

    C1input[1][3] = C1input[1][2];
    C1input[1][2] = C1input[1][1];
    C1input[1][1] = x;

    for (i = 1;i <= half_order_pole;i++) {
        preal = p[i*2-1].x;
        pimg  = p[i*2-1].y;

        temp  = pow((fs_bilinear - preal), 2) + pow(pimg, 2);


        /*dy = (input[i][1] + (1-(fs_bilinear+rzero)/(fs_bilinear-rzero))*input[i][2]
                              - (fs_bilinear+rzero)/(fs_bilinear-rzero)*input[i][3] );
        dy = dy+2*output[i][1]*(fs_bilinear*fs_bilinear-preal*preal-pimg*pimg);

        dy = dy-output[i][2]*((fs_bilinear+preal)*(fs_bilinear+preal)+pimg*pimg);*/

        dy = C1input[i][1] * (fs_bilinear - rzero) - 2 * rzero * C1input[i][2] - (fs_bilinear + rzero) * C1input[i][3]
             + 2 * C1output[i][1] * (fs_bilinear * fs_bilinear - preal * preal - pimg * pimg)
             - C1output[i][2] * ((fs_bilinear + preal) * (fs_bilinear + preal) + pimg * pimg);

        dy = dy / temp;

        C1input[i+1][3] = C1output[i][2];
        C1input[i+1][2] = C1output[i][1];
        C1input[i+1][1] = dy;

        C1output[i][2] = C1output[i][1];
        C1output[i][1] = dy;
    }

    dy = C1output[half_order_pole][1] * norm_gain;  /* don't forget the gain term */
    c1filterout = dy / 4.0;   /* signal path output is divided by 4 to give correct C1 filter gain */

    return (c1filterout);
}

/* ---------------------------------------------------------------------------- */
/** Parallelpath C2 filter: same as the signal-path C1 filter with the OHC completely impaired */

double C2ChirpFilt(double xx, double tdres, double cf, int n, double taumax, double fcohc)
{
    static double C2gain_norm, C2initphase;
    static double C2input[12][4];  static double C2output[12][4];

    double ipw, ipb, rpa, pzero, rzero;

    double sigma0, fs_bilinear, CF, norm_gain, phase, c2filterout;
    int    i, r, order_of_pole, half_order_pole, order_of_zero;
    double temp, dy, preal, pimg;

    COMPLEX p[11];

    /*================ setup the locations of poles and zeros =======*/

    sigma0 = 1 / taumax;
    ipw    = 1.01 * cf * TWOPI - 50;
    ipb    = 0.2343 * TWOPI * cf - 1104;
    rpa    = pow(10, log10(cf) * 0.9 + 0.55) + 2000;
    pzero  = pow(10, log10(cf) * 0.7 + 1.6) + 500;
    /*===============================================================*/

    order_of_pole    = 10;
    half_order_pole  = order_of_pole / 2;
    order_of_zero    = half_order_pole;

    fs_bilinear = TWOPI * cf / tan(TWOPI * cf * tdres / 2);
    rzero       = -pzero;
    CF          = TWOPI * cf;

    if (n == 0) {
        p[1].x = -sigma0;

        p[1].y = ipw;

        p[5].x = p[1].x - rpa; p[5].y = p[1].y - ipb;

        p[3].x = (p[1].x + p[5].x) * 0.5; p[3].y = (p[1].y + p[5].y) * 0.5;

        p[2] = compconj(p[1]); p[4] = compconj(p[3]); p[6] = compconj(p[5]);

        p[7] = p[1]; p[8] = p[2]; p[9] = p[5]; p[10] = p[6];

        C2initphase = 0.0;
        for (i = 1;i <= half_order_pole;i++) {
            preal     = p[i*2-1].x;
            pimg      = p[i*2-1].y;
            C2initphase = C2initphase + atan(CF / (-rzero)) - atan((CF - pimg) / (-preal)) - atan((CF + pimg) / (-preal));
        };

        /*===================== Initialize C2input & C2output =====================*/

        for (i = 1;i <= (half_order_pole + 1);i++) {
            C2input[i][3] = 0;
            C2input[i][2] = 0;
            C2input[i][1] = 0;
            C2output[i][3] = 0;
            C2output[i][2] = 0;
            C2output[i][1] = 0;
        }

        /*===================== normalize the gain =====================*/

        C2gain_norm = 1.0;
        for (r = 1; r <= order_of_pole; r++)
            C2gain_norm = C2gain_norm * (pow((CF - p[r].y), 2) + p[r].x * p[r].x);
    };

    norm_gain = sqrt(C2gain_norm) / pow(sqrt(CF * CF + rzero * rzero), order_of_zero);

    p[1].x = -sigma0 * fcohc;

    if (p[1].x > 0.0) printf("C2ChirpFilt: The system becomes unstable.\n");

    p[1].y = ipw;

    p[5].x = p[1].x - rpa; p[5].y = p[1].y - ipb;

    p[3].x = (p[1].x + p[5].x) * 0.5; p[3].y = (p[1].y + p[5].y) * 0.5;

    p[2] = compconj(p[1]); p[4] = compconj(p[3]); p[6] = compconj(p[5]);

    p[7] = p[1]; p[8] = p[2]; p[9] = p[5]; p[10] = p[6];

    phase = 0.0;
    for (i = 1;i <= half_order_pole;i++) {
        preal = p[i*2-1].x;
        pimg  = p[i*2-1].y;
        phase = phase - atan((CF - pimg) / (-preal)) - atan((CF + pimg) / (-preal));
    };

    rzero = -CF / tan((C2initphase - phase) / order_of_zero);
    if (rzero > 0.0) printf("C2ChirpFilt: The zeros are in the right-hand plane.\n");
    /*%==================================================  */
    /*%      time loop begins here                         */
    /*%==================================================  */

    C2input[1][3] = C2input[1][2];
    C2input[1][2] = C2input[1][1];
    C2input[1][1] = xx;

    for (i = 1;i <= half_order_pole;i++) {
        preal = p[i*2-1].x;
        pimg  = p[i*2-1].y;

        temp  = pow((fs_bilinear - preal), 2) + pow(pimg, 2);

        /*dy = (input[i][1] + (1-(fs_bilinear+rzero)/(fs_bilinear-rzero))*input[i][2]
                              - (fs_bilinear+rzero)/(fs_bilinear-rzero)*input[i][3] );
        dy = dy+2*output[i][1]*(fs_bilinear*fs_bilinear-preal*preal-pimg*pimg);

        dy = dy-output[i][2]*((fs_bilinear+preal)*(fs_bilinear+preal)+pimg*pimg);*/

        dy = C2input[i][1] * (fs_bilinear - rzero) - 2 * rzero * C2input[i][2] - (fs_bilinear + rzero) * C2input[i][3]
             + 2 * C2output[i][1] * (fs_bilinear * fs_bilinear - preal * preal - pimg * pimg)
             - C2output[i][2] * ((fs_bilinear + preal) * (fs_bilinear + preal) + pimg * pimg);

        dy = dy / temp;

        C2input[i+1][3] = C2output[i][2];
        C2input[i+1][2] = C2output[i][1];
        C2input[i+1][1] = dy;

        C2output[i][2] = C2output[i][1];
        C2output[i][1] = dy;

    };

    dy = C2output[half_order_pole][1] * norm_gain;
    c2filterout = dy / 4.0;

    return (c2filterout);
}

/* -------------------------------------------------------------------------------------------- */
/** Pass the signal through the Control path Third Order Nonlinear Gammatone Filter */

double WbGammaTone(double x, double tdres, double centerfreq, int n, double tau, double gain, int order)
{
    static double wbphase;
    static COMPLEX wbgtf[4], wbgtfl[4];

    double delta_phase, dtmp, c1LP, c2LP, out;
    int i, j;

    if (n == 0) {
        wbphase = 0;
        for (i = 0; i <= order;i++) {
            wbgtfl[i] = compmult(0, compexp(0));
            wbgtf[i]  = compmult(0, compexp(0));
        }
    }

    delta_phase = -TWOPI * centerfreq * tdres;
    wbphase += delta_phase;

    dtmp = tau * 2.0 / tdres;
    c1LP = (dtmp - 1) / (dtmp + 1);
    c2LP = 1.0 / (dtmp + 1);
    wbgtf[0] = compmult(x, compexp(wbphase));                /* FREQUENCY SHIFT */

    for (j = 1; j <= order; j++)                             /* IIR Bilinear transformation LPF */
        wbgtf[j] = comp2sum(compmult(c2LP * gain, comp2sum(wbgtf[j-1], wbgtfl[j-1])),
                            compmult(c1LP, wbgtfl[j]));
    out = REAL(compprod(compexp(-wbphase), wbgtf[order])); /* FREQ SHIFT BACK UP */

    for (i = 0; i <= order;i++) wbgtfl[i] = wbgtf[i];
    return(out);
}

/* ----------------------------------------------------------------------------- */
/** Calculate the gain and group delay for the Control path Filter */

double gain_groupdelay(double tdres, double centerfreq, double cf, double tau, int *grdelay)
{
    double tmpcos, dtmp2, c1LP, c2LP, tmp1, tmp2, wb_gain;

    tmpcos = cos(TWOPI * (centerfreq - cf) * tdres);
    dtmp2 = tau * 2.0 / tdres;
    c1LP = (dtmp2 - 1) / (dtmp2 + 1);
    c2LP = 1.0 / (dtmp2 + 1);
    tmp1 = 1 + c1LP * c1LP - 2 * c1LP * tmpcos;
    tmp2 = 2 * c2LP * c2LP * (1 + tmpcos);

    wb_gain = pow(tmp1 / tmp2, 1.0 / 2.0);

    grdelay[0] = (int)floor((0.5 - (c1LP * c1LP - c1LP * tmpcos) / (1 + c1LP * c1LP - 2 * c1LP * tmpcos)));

    return(wb_gain);
}
/* -------------------------------------------------------------------------------------------- */
/** Calculate the delay (basilar membrane, synapse, etc. for cat) */

double delay_cat(double cf)
{
    double A0, A1, x, delay;

    A0    = 3.0;
    A1    = 12.5;
    x     = 11.9 * log10(0.80 + cf / 456.0);      /* cat mapping */
    delay = A0 * exp(-x / A1) * 1e-3;

    return(delay);
}
/* -------------------------------------------------------------------------------------------- */
/* Get the output of the OHC Nonlinear Function (Boltzman Function) */

double Boltzman(double x, double asym, double s0, double s1, double x1)
{
    double shift, x0, out1, out;

    shift = 1.0 / (1.0 + asym);  /* asym is the ratio of positive Max to negative Max*/
    x0    = s0 * log((1.0 / shift - 1) / (1 + exp(x1 / s1)));

    out1 = 1.0 / (1.0 + exp(-(x - x0) / s0) * (1.0 + exp(-(x - x1) / s1))) - shift;
    out = out1 / (1 - shift);

    return(out);
}  /* output of the nonlinear function, the output is normalized with maximum value of 1 */

/* -------------------------------------------------------------------------------------------- */
/* Get the output of the OHC Low Pass Filter in the Control path */

double OhcLowPass(double x, double tdres, double Fc, int n, double gain, int order)
{
    static double ohc[4], ohcl[4];

    double c, c1LP, c2LP;
    int i, j;

    if (n == 0) {
        for (i = 0; i < (order + 1);i++) {
            ohc[i] = 0;
            ohcl[i] = 0;
        }
    }

    c = 2.0 / tdres;
    c1LP = (c - TWOPI * Fc) / (c + TWOPI * Fc);
    c2LP = TWOPI * Fc / (TWOPI * Fc + c);

    ohc[0] = x * gain;
    for (i = 0; i < order;i++)
        ohc[i+1] = c1LP * ohcl[i+1] + c2LP * (ohc[i] + ohcl[i]);
    for (j = 0; j <= order;j++) ohcl[j] = ohc[j];
    return(ohc[order]);
}
/* -------------------------------------------------------------------------------------------- */
/* Get the output of the IHC Low Pass Filter  */

double IhcLowPass(double x, double tdres, double Fc, int n, double gain, int order)
{
    static double ihc[8], ihcl[8];

    double C, c1LP, c2LP;
    int i, j;

    if (n == 0) {
        for (i = 0; i < (order + 1);i++) {
            ihc[i] = 0;
            ihcl[i] = 0;
        }
    }

    C = 2.0 / tdres;
    c1LP = (C - TWOPI * Fc) / (C + TWOPI * Fc);
    c2LP = TWOPI * Fc / (TWOPI * Fc + C);

    ihc[0] = x * gain;
    for (i = 0; i < order;i++)
        ihc[i+1] = c1LP * ihcl[i+1] + c2LP * (ihc[i] + ihcl[i]);
    for (j = 0; j <= order;j++) ihcl[j] = ihc[j];
    return(ihc[order]);
}
/* -------------------------------------------------------------------------------------------- */
/* Get the output of the Control path using Nonlinear Function after OHC */

double NLafterohc(double x, double taumin, double taumax, double asym)
{
    double R, dc, R1, s0, x1, out, minR;

    minR = 0.05;
    R  = taumin / taumax;

    if (R < minR) minR = 0.5 * R;
    else       minR = minR;

    dc = (asym - 1) / (asym + 1.0) / 2.0 - minR;
    R1 = R - minR;

    /* This is for new nonlinearity */
    s0 = -dc / log(R1 / (1 - minR));

    x1  = fabs(x);
    out = taumax * (minR + (1.0 - minR) * exp(-x1 / s0));
    if (out < taumin) out = taumin;
    if (out > taumax) out = taumax;
    return(out);
}
/* -------------------------------------------------------------------------------------------- */
/* Get the output of the IHC Nonlinear Function (Logarithmic Transduction Functions) */

double NLogarithm(double x, double slope, double asym, double cf)
{
    double corner, strength, xx, splx, asym_t;

    corner    = 80;
    strength  = 20.0e6 / pow(10, corner / 20);

    xx = log(1.0 + strength * fabs(x)) * slope;

    if (x < 0) {
        splx   = 20 * log10(-x / 20e-6);
        asym_t = asym - (asym - 1) / (1 + exp(splx / 5.0));
        xx = -1 / asym_t * xx;
    };
    return(xx);
}
/* -------------------------------------------------------------------------------------------- */


/**********************************************************************

  resample function ammended from testresample.c

  Real-time library interface by Dominic Mazzoni

  Based on resample-1.7:
    http://www-ccrma.stanford.edu/~jos/resample/

  License: LGPL - see the file LICENSE.txt for more information

**********************************************************************/
/*#include "libresample.h"



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
   printf("-- srclen: %d  factor: %g srcblk: %d dstblk: %d expected %d\n",srclen,  factor, srcblocksize, dstblocksize,expectedlen);


    src = (float*) malloc((unsigned) (srclen*sizeof(float)));
   for (i=0;i<srclen;i++) src[i] = (float) source[i];
   dst = (float*) malloc((unsigned) ((dstlen+100)*sizeof(float)));//makevector(dstlen+100);

   printf(" source  %x\t src %x\t dst %x\n",&source[0],&src[0],&dst[0]);

   handle = resample_open(1, factor, factor);
   fwidth = resample_get_filter_width(handle);
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
      hoc_execerror(" resample_process returned an error ", 0);
      return 0;
   }

   if (out <= 0) {
     hoc_execerror(" resample_process returned less than zero samples",0);
     free((char*)src);
     free((char*)dst);
     return 0;
   }

   lendiff = abs(out - expectedlen);
   if (lendiff > (int)(2*factor + 1.0)) {
      printf("resample:   Expected ~%d, got %d samples out\n", expectedlen, out);
   }
   printf("resample_close 1: dst[0] %g\t dst[end] %g\t dstlen %d\n",dst[0],dst[dstlen+99], dstlen);
   int len =(int)(srclen * factor);
   dest = makevector(len);
zero_vector(dest,len);
 printf("resample_close 2:dstlen %d, len %d,  dest[0] %g\t dest[end] %g\n",dstlen, len,dest[0],dest[len-1]);
   
 for (i=0;i<len;i++) {
   if ( isnan(dst[i]) ) {len = i-1; break;}
   dest[i] = (double) dst[i];

   //   printf("dst[%d]\t %g\n",i,dst[i]);

 }   

   printf("resample_close 3:  len %d \tdest[0] %g\t dest[len] %g\n",len,dest[0],dest[len-1]);   
   free((char*) src);
   free((char*) dst);
   printf("resample_close\n");
   return len;
}
*/
/*  resample apapted from src/ivoc/ivocvect.cpp */

 
int resample(double* v1, double * ans, int capacity, double f)
{
    //  double * temp;
    int i;

    printf("resample: &v1 %x\t capacity %d  factor %g destlen %d\n", &v1[0], capacity, f, (int)(capacity * f));

    if (f < 0 || f > capacity / 2) {
        printf("ivoc_resample: nyquist error with factor");
        //return 0;
    }
    int n = (int)(capacity * f);
    ans = makevector(n);
    for (i = 0; i < n; i++) ans[i] = v1[(int)(i/f)];
    //*ans = temp;
    return n;
}

