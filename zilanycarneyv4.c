/* This is Version 4 of the public distribution of the code for the auditory
   periphery model of:

   Zilany, M.S.A., Bruce, I.C., Nelson, P.C., and Carney, L.H. (2009). "A Phenomenological
   model of the synapse between the inner hair cell and auditory nerve : Long-term adaptation
   with power-law dynamics," Journal of the Acoustical Society of America 126(5): 2390-2412.

   Please cite this paper if you publish any research
   results obtained with this code or any modified versions of this code.

   See the file readme.txt for details of compiling and running the model.

   %%% © Muhammad S.A. Zilany (msazilany@gmail.com), Ian C. Bruce, Paul C. Nelson, and laurel H. Carney October 2008 %%%

*/
//#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>
//#include <math.h>     /* Added for MS Visual C++ compatability, by Ian Bruce, 1999 */
//#include <mex.h>
//#include <time.h>


/* Declarations of the functions used in the program */
/*
  extern double C1ChirpFilt(double, double, double, int, double, double);
  extern double C2ChirpFilt(double, double, double, int, double, double);
  extern double WbGammaTone(double, double, double, int, double, double, int);
  extern double gain_groupdelay(double, double, double, double, int *);
  extern double Get_tauwb(double,double, int, double *, double *); // Calc gain in IHCAN then pass to tauwb
  extern double Get_taubm(double,double, double, double *, double *, double *);

  extern double OhcLowPass(double, double, double, int, double, int);
  extern double IhcLowPass(double, double, double, int, double, int);
  extern double Boltzman(double, double, double, double, double);
  extern double NLafterohc(double, double, double, double);
  extern double ControlSignal(double, double, double, double, double);
  extern double NLogarithm(double, double, double); // extra double does nothing in this version

  extern double cochlea_f2x(int , double);
  extern double cochlea_x2f(int , double);
  extern double delay_cat(double , int);
*/

/* Declarations of the functions used in the SingleAN_v4 program */
double Synapse_v4(double *, double, double, int, int, double, double, double, double *, double);
int    SpikeGenerator_v4(double *, double, int, int, double *);
double dbl_exp_adaptation(double cf, double spont);

#ifdef _FFGN_
 double ffGn(double *yffGn, int N, double tdres, double Hinput, double mu, double sigma);
#endif


void IHCAN(double *px, double cf, int nrep, double tdres, int totalstim, 
	   double cohc, double cihc, double *ihcout, double species)
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

  nrep = 1; /* restrict stim should be what is presented in px, not multiples */

  /* Allocate dynamic memory for the temporary variables */
  ihcouttmp  = makevector(totalstim*nrep);

  mey1 = makevector(totalstim);
  mey2 = makevector(totalstim);
  mey3 = makevector(totalstim);

  tmpgain = makevector(totalstim);

  /** Calculate the location on basilar membrane from CF */

  bmplace = cochlea_f2x(species, cf);// 11.9 * log10(0.80 + cf / 456.0);

  /** Calculate the center frequency for the control-path wideband filter
      from the location on basilar membrane */

  centerfreq = cochlea_x2f(species, bmplace + 1.2); /* shift the center freq */

  /*==================================================================*/
  /*====== Parameters for the gain ===========*/
  gain = 52 / 2 * (tanh(2.2 * log10(cf / 0.6e3) + 0.15) + 1);
  /*gain = 52/2*(tanh(2.2*log10(cf/1e3)+0.15)+1);*/
  if (gain > 60) gain = 60;
  if (gain < 15) gain = 15;

  /*====== Parameters for the control-path wideband filter =======*/
  bmorder = 3;
  Get_tauwb(cf, gain, bmorder, Taumax, Taumin);
  taubm   = cohc * (Taumax[0] - Taumin[0]) + Taumin[0];
  ratiowb = Taumin[0] / Taumax[0];
  /*====== Parameters for the signal-path C1 filter ======*/
  Get_taubm(cf, gain, Taumax[0], bmTaumax, bmTaumin, ratiobm);
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

    if (1 / tauc1 < 0.0) printf("\t\t\tIHCAN: The poles are in the right-half plane; system is unstable.\n");

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

    c1vihctmp  = NLogarithm(cihc * c1filterouttmp, 0.1, ihcasym); // no need for cf

    c2vihctmp = -NLogarithm(c2filterouttmp * fabs(c2filterouttmp) * cf / 10 * cf / 2e3, 0.2, 1.0); /* C2 transduction output EDIT removed cf in NLog function*/

    ihcouttmp[n] = IhcLowPass(c1vihctmp + c2vihctmp, tdres, 3000, n, 1.0, 7);

    if (isnan(ihcouttmp[i])) {
      printf("\t\t\t\tIHCAN: nan at %d\n\tError in an_zilany V4\n",i);
      // return ;
    }

  }
  ;  /* End of the loop */
#ifdef DEBUG
  printf("\t\t\t\tIHCAN: End of the loop.\n");
#endif
  /* Stretched out the IHC output according to nrep (number of repetitions) */

  for (i = 0;i < totalstim*nrep;i++) {
    ihcouttmp[i] = ihcouttmp[(int)(fmod(i,totalstim))];
  };

  /* Adjust total path delay to IHC output signal */

  delay      = delay_species(cf,species);
  delaypoint = __max(0, (int) ceil(delay / tdres));

  for (i = delaypoint;i < totalstim*nrep;i++) {
    ihcout[i] = ihcouttmp[i - delaypoint];
  };

  /* Freeing dynamic memory allocated earlier */

  freevector(ihcouttmp);
  freevector(mey1);
  freevector(mey2);
  freevector(mey3);
  freevector(tmpgain);
#ifdef DEBUG
  printf("\t\t\t\tIHCAN: done.\n");
#endif
} /* End of the IHCAN function */


double SingleAN_v4_1(double *px, double cf, int nrep, double tdres, int totalstim, double spont, double implnt, double *synout, double species)
{

  /*variables for the signal-path, control-path and onward */
  double *synouttmp;
  int    i, nspikes, ipst;
  double nsout;
  double sampFreq = 10e3; /* Sampling frequency used in the synapse */

  /* Allocate dynamic memory for the temporary variables */
  synouttmp  = makevector(totalstim * nrep);

  /*====== Run the synapse model ======*/
  nsout = Synapse_v4(px, tdres, cf, totalstim, nrep, spont, implnt, sampFreq, synouttmp, species);

  /* Wrapping up the unfolded (due to no. of repetitions) Synapse Output */
  for (i = 0; i < nsout ; i++) {
    ipst = (int)(fmod(i, totalstim));
    synout[ipst] = synout[ipst] + synouttmp[i] / nrep;
  };

  /* Freeing dynamic memory allocated earlier */
  freevector(synouttmp);
  return nsout;
} /* End of the SingleAN function */


/*-------------------------------------------------------------------------------
  Synapse model: if the time resolution is not small enough, the concentration of
  the immediate pool could be as low as negative, at this time there is an alert
  message print out and the concentration is set at saturated level
  -------------------------------------------------------------------------------*/
double Synapse_v4(double *ihcout, double tdres, double cf, int totalstim, int nrep, double spont, double implnt, double sampFreq, double *synouttmp, double species)
{
  /* Initalize Variables */
  double rmean, rstd;
  int z, b, n;
  int resamp = (int) ceil(1 / (tdres * sampFreq)); //sampFreq is not the same as 1/tdres
  double incr = 0.0;    // int delaypoint = floor(7500 / (cf / 1e3));   
    
  double delay      = delay_species(cf, species);
  int delaypoint = __max(0, (int) ceil(delay / tdres));  // from version 2
#ifdef DEBUG    
  printf("\t\t\tSynapse_v4: resamp %d delaypoint %d\n",resamp,delaypoint);
#endif
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

  nrep=1; /* do overlapping in the input stimulus vector */

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
  /*----------------------------------------------------------*/
#ifdef DEBUG
  printf("\t\t\tSynapse: Generating a random sequence\n");
#endif
  int Nrand = ((int) ceil((totalstim * nrep + 2 * delaypoint) * tdres * sampFreq));

  randNums = makevector(Nrand); zero_vector(randNums,Nrand);
    
#ifdef _FFGN_ 
  if (!(ffGn(randNums, Nrand, 1 / sampFreq, 0.9, spont, -1.0))) {
    hoc_execerror("\t\t\tSynapse: error calling ffGn", 0);
    return 0;
  }
  printf("\t\t\tSynapse: Completed ffGn\n");

  for (indx = 0;indx < Nrand;indx++)  {
    if (isnan(randNums[indx])){
      printf("\t\t\tSynapse: found NaN  %d\n",indx);
      return 0;
    }
    rmean += randNums[indx];

  }
  for (indx = 0;indx < Nrand;indx++) rstd += pow((randNums[indx] - rmean), 2);
  printf("\t\t\tSynapse: Completed ffGn: mean  %g\t stdev %g\n", rmean, sqrt(rstd / Nrand));
#else
  //     for (indx = 0;indx < Nrand;indx++)  randNums[indx]=-spont;
  //  for (indx = 0;indx < Nrand;indx++)  randNums[indx]=0.0;//spont+exp(-indx);
#endif

  /*----------------------------------------------------------*/
  /*----- Double Exponential Adaptation -----(Replacement)----*/
  /*----------------------------------------------------------*/
  cf_factor = dbl_exp_adaptation(cf,spont);

#ifdef DEBUG
  printf("\t\t\tSynapse: Generating parameters\n");
#endif

  /*----------------------------------------------------------*/
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

#ifdef DEBUG
  printf("\t\t\tSynapse: Generating exponOut\n");
#endif

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

#ifdef DEBUG
  printf("\t\t\tSynapse: Generating powerLawIn\n");
#endif

  for (k = 0; k < delaypoint; k++)
    powerLawIn[k] = exponOut[0];
  for (k = delaypoint; k < totalstim*nrep + delaypoint; k++)
    powerLawIn[k] = exponOut[k-delaypoint];
  for (k = totalstim * nrep + delaypoint; k < totalstim*nrep + 3*delaypoint; k++)
    powerLawIn[k] = powerLawIn[k-1];

  /*----------------------------------------------------------*/
  /*------ Downsampling to sampFreq (Low) sampling rate ------*/
  /*----------------------------------------------------------*/
#ifdef DEBUG
  printf("\t\t\tSynapse: Downsampling to sampFreq (Low) sampling rate\n");
#endif

  /* Resampling routine from libresample examples ***/
  int len = 0;
#ifdef DEBUG
  printf("\t\t\tSynapse: calling resample(%d,NULL,%d,%g)\n", &powerLawIn[0],  k, 1.0 / resamp);
#endif
  sampIHC = makevector( (int)((double)k / resamp) );
  len = resample(powerLawIn, sampIHC, k, 1.0/resamp);
  if (len == 0 || (sampIHC == NULL)) {
    printf("\t\t\tSynapse: resample return error, copying powerLawIn\n");
    len = k;
    if (sampIHC) {
      /* printf("\t\t\tsampIHC address %x", &sampIHC[0]); */ freevector(sampIHC);
    }
    sampIHC = makevector( (int)(k / resamp) );
    for (indx = 0;indx < (int)(k / resamp);indx++) sampIHC[indx] = powerLawIn[(int)round(indx /resamp)];
  }
#ifdef DEBUG
  printf("\t\t\tSynapse: resample done, k %d len %d \t old len %d sampIHC[0] %g\n", k, len, floor((totalstim*nrep + 2*delaypoint)*tdres*sampFreq), sampIHC[0]);

#endif
  freevector(powerLawIn); freevector(exponOut);
  /*----------------------------------------------------------*/
  /*----- Running Power-law Adaptation -----------------------*/
  /*----------------------------------------------------------*/
#ifdef DEBUG
  printf("\t\t\tSynapse: Running Power-law Adaptation\n");
#endif

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

  /*----------------------------------------------------------*/
  /*----- Upsampling to original (High 100 kHz) sampling rate */
  /*----------------------------------------------------------*/
#ifdef DEBUG
  printf("\t\t\tSynapse: Upsampling to original (High 100 kHz) sampling rate\n");
#endif

  for (z = 0; z < k - 1; ++z) {
    incr = (synSampOut[z+1] - synSampOut[z]) / resamp;
    for (b = 0; b < resamp; ++b) {
      TmpSyn[z*resamp+b] = synSampOut[z] + b * incr;
    }
  }

  for (i = 0;i < totalstim*nrep;++i)
    synouttmp[i] = TmpSyn[i+delaypoint];

#ifdef DEBUG
  printf("\t\t\tSynapse_v4: done!\n");
  printf("\t\t\t alpha1 %g beta1 %g I1 %g alpha2 %g beta2 %g I2 %g binwidth %g\n",alpha1, beta1, I1, alpha2, beta2, I2, binwidth);
  printf("\t\t\t synstrength %g synslope %g CI %g CL %g PG %g CG %g VL %g PL %g VI %g\n", synstrength, synslope, CI, CL, PG, CG, VL, PL, VI);
  printf("\t\t\t cf_factor %g PImax %g kslope %g Ass %g Asp %g TauR %g TauST %g Ar_Ast %g PTS %g Aon %g AR %g AST %g Prest %g gamma1 %g gamma2 %g k1 %g k2 %g\n",      cf_factor, PImax, kslope, Ass, Asp, TauR, TauST, Ar_Ast, PTS, Aon, AR, AST, Prest, gamma1, gamma2, k1, k2);
  printf("\t\t\t VI0 %g VI1 %g alpha %g beta %g theta1 %g theta2 %g theta3 %g vsat %g tmpst %g tmp %g PPI %g CIlast %g temp %g\n",     VI0, VI1, alpha, beta, theta1, theta2, theta3, vsat, tmpst, tmp, PPI, CIlast, temp);

  printf("\t\t\tsout1[0] %g sout2[0] %g \
    m1[0] %g m2[0] %g m3[0] %g m4[0] %g m5[0] %g \
    n1[0] %g n2[0] %g n3[0] %g synSampOut[0] %g TmpSyn[0] %g\n",
	 sout1[0],sout2[0],
	 m1[0],m2[0],m3[0],m4[0],m5[0],
	 n1[0],n2[0],n3[0],synSampOut[0],TmpSyn[0]) ;

#endif
  freevector(sout1); freevector(sout2);
  freevector(m1); freevector(m2); freevector(m3); freevector(m4); freevector(m5); 
  freevector(n1); freevector(n2); freevector(n3);

  freevector(synSampOut); freevector(TmpSyn);
  freevector(randNums);
  return((long) ceil(totalstim*nrep));
}



/* ------------------------------------------------------------------------------------ */
/* Pass the output of Synapse model through the Spike Generator */

/* The spike generator now uses a method coded up by B. Scott Jackson (bsj22@cornell.edu)
   Scott's original code is available from Laurel Carney's web site at:
   http://www.urmc.rochester.edu/smd/Nanat/faculty-research/lab-pages/LaurelCarney/auditory-models.cfm
*/

int SpikeGenerator_v4(double *synouttmp, double tdres, int totalstim, int nrep, double *sptime)
{
  double  c0, s0, c1, s1, dead;
  int j;
  int     nspikes, k, NoutMax, Nout, deadtimeIndex, randBufIndex;
  double deadtimeRnd, endOfLastDeadtime, refracMult0, refracMult1, refracValue0, refracValue1;
  double Xsum, unitRateIntrvl, countTime, DT;

  double *randNums;

  c0      = 0.5;
  s0      = 0.001;
  c1      = 0.5;
  s1      = 0.0125;
  dead    = 0.00075;

  DT = totalstim * tdres * nrep;  /* Total duration of the rate function */
  Nout = 0;
  NoutMax = (long) ceil(totalstim * nrep * tdres / dead);
  randNums = makevector(NoutMax + 1);

  for (k=0;k<=NoutMax;k++) randNums[k] = scop_random();  //replace mexCallMATLAB to rand(1,NoutMax)
  randBufIndex = 0;

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

  printf("\t\t\t c0 %g\ts0 %g\tc1 %g\ts1 %g\tdead %g\n",c0, s0, c1, s1, dead);
  printf("\t\t\t nspikes %ld\tk %ld\tNoutMax %ld\tNout %ld\tdeadtimeIndex %ld\trandBufIndex %ld\n", nspikes, k, NoutMax, Nout, deadtimeIndex, randBufIndex);
  printf("\t\t\tdeadtimeRnd %g\tendOfLastDeadtime %g\trefracMult0 %g\trefracMult1 %g\trefracValue0 %g\trefracValue1 %g\n", deadtimeRnd, endOfLastDeadtime, refracMult0, refracMult1, refracValue0, refracValue1);
  printf("\t\t\tXsum %g\tunitRateIntrvl %g\tcountTime %g\tDT %g\n",Xsum, unitRateIntrvl, countTime, DT);



  for (j=0;j<nrep;j++){
    countTime = tdres;
    for (k = 0; (k < totalstim) && (countTime < DT); ++k, countTime += tdres, refracValue0 *= refracMult0, refracValue1 *= refracMult1) { /* Loop through rate vector */
      if (synouttmp[k] > 0) { /* Nothing to do for non-positive rates, i.e. Xsum += 0 for non-positive rates. */
	Xsum += synouttmp[k] * (1 - refracValue0 - refracValue1);  /* Add synout*(refractory value) to time-warping sum */

	if (Xsum >= unitRateIntrvl) {  /* Spike occurs when time-warping sum exceeds interspike "time" in unit-rate process */
	  sptime[Nout] = countTime;
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
    sptime[Nout++] = 0.0;
  }


  freevector(randNums);
  nspikes = Nout;  /* Number of spikes that occurred. */
  return(nspikes);
}


double SingleAN_v4(double *px, double cf, int nrep, double tdres, int totalstim, double fibertype, double implnt, double *synout, double species)
{

  /*variables for the signal-path, control-path and onward */
  double *synouttmp, *sptime;

  int    i, nspikes, ipst;
  double I, spont;
  double sampFreq = 10e3; /* Sampling frequency used in the synapse */

  /* Allocate dynamic memory for the temporary variables */
  synouttmp  = makevector(totalstim * nrep);
  //sptime  = makevector((long) ceil(totalstim * tdres * nrep / 0.00075));

  /* Spontaneous Rate of the fiber corresponding to Fibertype */
    
  if(fibertype < 0){ spont = fabs(fibertype);}
  else if (fibertype == 1){ spont = 0.1;}
  else if (fibertype == 2){ spont = 5.0;}
  else if (fibertype == 3){ spont = 100.0;}

  /*====== Run the synapse model ======*/
  I = Synapse_v4(px, tdres, cf, totalstim, nrep, spont, implnt, sampFreq, synouttmp, species);

  /* Wrapping up the unfolded (due to no. of repetitions) Synapse Output */
  for (i = 0; i < I ; i++) {
    ipst = (int)(fmod(i, totalstim));
    synout[ipst] = synout[ipst] + synouttmp[i] / nrep;
  };

  /* Freeing dynamic memory allocated earlier */
  //freevector(sptime);
  freevector(synouttmp);
  return I;
} /* End of the SingleAN function */

void PsthAN(double *px, double cf, int nrep, double tdres, int totalstim, double fibertype, double implnt, double species,double *synout, double *psth)
{	
  
  /*variables for the signal-path, control-path and onward */
  double *synouttmp,*sptime;
  int    i,nspikes,ipst;
  double I,spont;
  double sampFreq = 10e3; /* Sampling frequency used in the synapse */
        
    
  /* Allocate dynamic memory for the temporary variables */
  synouttmp  = makevector(totalstim*nrep);
  sptime  = makevector((long) ceil(totalstim*tdres*nrep/0.00075));  	
	   
  /* Spontaneous Rate of the fiber corresponding to Fibertype */    
  if (fibertype==1) spont = 0.1;
  if (fibertype==2) spont = 5.0;
  if (fibertype==3) spont = 100.0;
    
  /*====== Run the synapse model ======*/    
  I = Synapse_v4(px, tdres, cf, totalstim, nrep, spont, implnt, sampFreq, synouttmp,species);
            
  /* Wrapping up the unfolded (due to no. of repetitions) Synapse Output */
  for(i = 0; i <I ; i++)
    {
      // ipst = (int) (fmod(i,totalstim));
      // synout[ipst] = synout[ipst] + synouttmp[i]/nrep;
      synout[i] = synouttmp[i];
    };    
  /*======  Spike Generations ======*/
  printf("\t\t\tPsthAN: calling SpikeGenerator_v4");
  nspikes = SpikeGenerator_v4(synouttmp, tdres, totalstim, nrep, sptime);
  for(i = 0; i < nspikes; i++)
    {        
      //	ipst = (int) (fmod(sptime[i],tdres*totalstim) / tdres);
      ipst = (int) (fmod(sptime[i],tdres*totalstim) * sampFreq);
      psth[ipst] = psth[ipst] + 1;       
    };

  /* Freeing dynamic memory allocated earlier */

  freevector(sptime); freevector(synouttmp); 

} /* End of the PSTH version of SingleAN function */




double dbl_exp_adaptation(double cf, double spont)
{
  /*----------------------------------------------------------*/
  /*----- Double Exponential Adaptation ----------------------*/
  /*----------------------------------------------------------*/
  if (spont >= 18){ return  __min(800.0, pow(10, 0.29 * cf / 1e3 + 0.7));
  } else if (spont > 2) { return __min(50.0, 2.5e-4 * cf * 4 + 0.2);
  } else return  __min(1.0, 2.5e-4 * cf * 0.1 + 0.15);

}


