/* This is Version 2 of the public distribution of the code for the auditory
   periphery model of:

   Zilany, M. S. A. and Bruce, I. C. (2006). "Modeling auditory-nerve
   responses for high sound pressure levels in the normal and impaired
   auditory periphery," Journal of the Acoustical Society of
   America 120(3):1446-1466,

   Zilany, M. S. A. and Bruce, I. C. (2007). "Representation of the vowel
   /eh/ in normal and impaired auditory nerve fibers: Model predictions of
   responses in cats," Journal of the Acoustical Society of America
   122(1):402-417.

   Please cite the Zilany and Bruce (2006, 2007) papers if you publish any research
   results obtained with this code or any modified versions of this code.

   Note that the only difference between the models presented in the two papers is
   the function for cochlear amplifier (CA) gain versus characteristic frequency (CF),
   given by Eq. (6) in Zilany and Bruce (2006) and by Eq. (1) in Zilany and Bruce (2007).
   Since version 1.1 of the code, the new CA gain versus CF function of Zilany and Bruce
   (2007) is used by default.  If you wish to use the old function, then you will need to
   uncomment line 288 of this code and comment out line 289, and then recompile the code.

   See the file readme.txt for details of compiling and running the model.

   %%% (C) Ian C. Bruce (ibruce@ieee.org), M. S. Arefeen Zilany, Rasha Ibrahim, June 2006 - December 2007 %%%

*/

//#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>
//#include <math.h>
//#include <mex.h>
//#include <time.h>

#include "complex.h"

//#include "spikegen.h"
//#include "gsl_rng.h"

#ifndef TWOPI
#define TWOPI 6.28318530717959
#endif

#ifndef __max
#define __max(a,b) (((a) > (b))? (a): (b))
#endif

#ifndef __min
#define __min(a,b) (((a) < (b))? (a): (b))
#endif


double cochlea_f2x(int , double);
double cochlea_x2f(int , double);
double delay_cat(double , int);

double Synapse(double, double, double, double, int);
//int    SpikeGenerator(double *, double, int, double, double *);


int an_zilanybruce2007(double binwidth, double cf, double spont, double cihc, double cohc, int species, int nrep, double *px, double *synout, int totalstim)
{
  double *timeout, *meout, *c1filterout, *c2filterout, *c1vihc,  *c2vihc,  *ihcout;

  /*variables for middle-ear model */
  double megainmax = 41.1405;
  double *mey1, *mey2, *mey3;
  double fp, C, m11, m12, m21, m22, m23, m24, m25, m26, m31, m32, m33, m34, m35, m36;

  /*variables for the signal-path, control-path and onward */
  double *c1filterouttmp, *c2filterouttmp, *c1vihctmp, *c2vihctmp, *ihcouttmp, *synouttmp, *tmpgain, *sptime;
  double    *grd; //makevector only allows double* vectors

  double bmplace, centerfreq, CAgain, taubm, ratiowb, bmTaubm, fcohc, TauWBMax, TauWBMin, tauwb;
  double Taumin[1], Taumax[1], bmTaumin[1], bmTaumax[1], ratiobm[1], lasttmpgain, wbgain, ohcasym, ihcasym, delay;
  int    i, n, delaypoint, grdelay[1], bmorder, wborder, nspikes, ipst;
  double wbout1, wbout, ohcnonlinout, ohcout, tmptauc1, tauc1, rsigma, wb_gain;

  /* Declarations of the functions used in the program */
  double C1ChirpFilt(double, double, double, int, double, double);
  double C2ChirpFilt(double, double, double, int, double, double);
  double WbGammaTone(double, double, double, int, double, double, int);
  double gain_groupdelay(double, double, double, double, int *);
  double Get_tauwb(double, double, int, double *, double *);
  double Get_taubm(double, double, double, double *, double *, double *);

  double OhcLowPass(double, double, double, int, double, int);
  double IhcLowPass(double, double, double, int, double, int);
  double Boltzman(double, double, double, double, double);
  double NLafterohc(double, double, double, double);
  double ControlSignal(double, double, double, double, double);
  double NLogarithm(double, double, double);


  /* Allocate dynamic memory for the temporary variables */

  c1filterouttmp   = makevector(totalstim);
  c2filterouttmp   = makevector(totalstim);
  c1vihctmp = makevector(totalstim);
  c2vihctmp = makevector(totalstim);
  ihcouttmp = makevector(totalstim);
  synouttmp = makevector(totalstim);

  mey1 = makevector(totalstim);
  mey2 = makevector(totalstim);
  mey3 = makevector(totalstim);

  grd     = makevector(totalstim);
  tmpgain = makevector(totalstim);
  meout = makevector(totalstim);
  // sptime  = (double*)mxCalloc((long) ceil(totalstim*binwidth*nrep/0.00075),sizeof(double));
  /** Calculate the location on basilar membrane from CF */
  // cochlea-freq mapping function using data from Zhang et al. 2001 and Heinz et al. 2001
  bmplace = cochlea_f2x(species, cf);

  /** Calculate the center frequency for the control-path wideband filter
      from the location on basilar membrane */

  centerfreq = cochlea_x2f(species, bmplace + 1.2); /* shift the center freq */

  /*==================================================================*/
  /*====== Parameters for the CAgain ===========*/
  /* CAgain = 52/2*(tanh(2.2*log10(cf/1e3)+0.15)+1); */ /* CA gain function used in Zilany and Bruce (2006) */
  CAgain = 52 / 2 * (tanh(2.2 * log10(cf / 600) + 0.15) + 1);    /* CA gain function used in Zilany and Bruce (2007) */
  if (CAgain < 15) CAgain = 15;

  /*====== Parameters for the control-path wideband filter =======*/
  bmorder = 3;
  Get_tauwb(cf, CAgain, bmorder, Taumax, Taumin);
  taubm   = cohc * (Taumax[0] - Taumin[0]) + Taumin[0];
  ratiowb = Taumin[0] / Taumax[0];
  /*====== Parameters for the signal-path C1 filter ======*/
  Get_taubm(cf, CAgain, Taumax[0], bmTaumax, bmTaumin, ratiobm);
  bmTaubm  = cohc * (bmTaumax[0] - bmTaumin[0]) + bmTaumin[0];
  fcohc    = bmTaumax[0] / bmTaubm;
  /*====== Parameters for the control-path wideband filter =======*/
  wborder  = 3;
  TauWBMax = Taumin[0] + 0.2 * (Taumax[0] - Taumin[0]);
  TauWBMin = TauWBMax / Taumax[0] * Taumin[0];
  tauwb    = TauWBMax + (bmTaubm - bmTaumax[0]) * (TauWBMax - TauWBMin) / (bmTaumax[0] - bmTaumin[0]);

  wbgain = gain_groupdelay(binwidth, centerfreq, cf, tauwb, grdelay);
  tmpgain[0]   = wbgain;
  lasttmpgain  = wbgain;
  /*===============================================================*/
  /* Nonlinear asymmetry of OHC function and IHC C1 transduction function*/
  ohcasym  = 7.0;
  ihcasym  = 3.0;
  /*===============================================================*/
  /* Prewarping and related constants for the middle ear */
  fp = 1e3;  /* prewarping frequency 1 kHz */
  C  = TWOPI * fp / tan(TWOPI / 2 * fp * binwidth);
  m11 = C / (C + 693.48);                    m12 = (693.48 - C) / C;
  m21 = 1 / (pow(C, 2) + 11053 * C + 1.163e8);  m22 = -2 * pow(C, 2) + 2.326e8;     m23 = pow(C, 2) - 11053 * C + 1.163e8;
  m24 = pow(C, 2) + 1356.3 * C + 7.4417e8;    m25 = -2 * pow(C, 2) + 14.8834e8;   m26 = pow(C, 2) - 1356.3 * C + 7.4417e8;
  m31 = 1 / (pow(C, 2) + 4620 * C + 909059944); m32 = -2 * pow(C, 2) + 2 * 909059944; m33 = pow(C, 2) - 4620 * C + 909059944;
  m34 = 5.7585e5 * C + 7.1665e7;             m35 = 14.333e7;                  m36 = 7.1665e7 - 5.7585e5 * C;

  for (n = 0;n < totalstim;n++) { /* Start of the loop */
    if (n == 0) { /* Start of the middle-ear filtering section  */
      mey1[0]  = m11 * px[0];
      mey2[0]  = mey1[0] * m24 * m21;
      mey3[0]  = mey2[0] * m34 * m31;
      meout[0] = mey3[0] / megainmax ;
    }

    else if (n == 1) {
      mey1[1]  = m11 * (-m12 * mey1[0] + px[1]       - px[0]);
      mey2[1]  = m21 * (-m22 * mey2[0] + m24 * mey1[1] + m25 * mey1[0]);
      mey3[1]  = m31 * (-m32 * mey3[0] + m34 * mey2[1] + m35 * mey2[0]);
      meout[1] = mey3[1] / megainmax;
    } else {
      mey1[n]  = m11 * (-m12 * mey1[n-1]  + px[n]         - px[n-1]);
      mey2[n]  = m21 * (-m22 * mey2[n-1] - m23 * mey2[n-2] + m24 * mey1[n] + m25 * mey1[n-1] + m26 * mey1[n-2]);
      mey3[n]  = m31 * (-m32 * mey3[n-1] - m33 * mey3[n-2] + m34 * mey2[n] + m35 * mey2[n-1] + m36 * mey2[n-2]);
      meout[n] = mey3[n] / megainmax;
    }
    ;  /* End of the middle-ear filtering section */

    //     timeout[n] = n*binwidth;

    /* Control-path filter */

    wbout1 = WbGammaTone(meout[n], binwidth, centerfreq, n, tauwb, wbgain, wborder);
    wbout  = pow((tauwb / TauWBMax), wborder) * wbout1 * 10e3 * __max(1, cf / 5e3);

    ohcnonlinout = Boltzman(wbout, ohcasym, 12.0, 5.0, 5.0); /* pass the control signal through OHC Nonlinear Function */
    ohcout = OhcLowPass(ohcnonlinout, binwidth, 600, n, 1.0, 2);/* lowpass filtering after the OHC nonlinearity */

    tmptauc1 = NLafterohc(ohcout, bmTaumin[0], bmTaumax[0], ohcasym); /* nonlinear function after OHC low-pass filter */
    tauc1    = cohc * (tmptauc1 - bmTaumin[0]) + bmTaumin[0];  /* time -constant for the signal-path C1 filter */
    rsigma   = 1 / tauc1 - 1 / bmTaumax[0]; /* shift of the location of poles of the C1 filter from the initial positions */

    if (1 / tauc1 < 0.0) printf("The poles are in the right-half plane; system is unstable.\n");

    tauwb = TauWBMax + (tauc1 - bmTaumax[0]) * (TauWBMax - TauWBMin) / (bmTaumax[0] - bmTaumin[0]);

    wb_gain = gain_groupdelay(binwidth, centerfreq, cf, tauwb, grdelay);

    grd[n] = grdelay[0];

    if (floor(grd[n] + n) < totalstim)
      tmpgain[(int)floor(grd[n] + n)] = wb_gain;

    if (tmpgain[n] == 0)
      tmpgain[n] = lasttmpgain;

    wbgain      = tmpgain[n];
    lasttmpgain = wbgain;

    /*====== Signal-path C1 filter ======*/

    c1filterouttmp[n] = C1ChirpFilt(meout[n], binwidth, cf, n, bmTaumax[0], rsigma); /* C1 filter output */


    /*====== Parallel-path C2 filter ======*/

    c2filterouttmp[n]  = C2ChirpFilt(meout[n], binwidth, cf, n, bmTaumax[0], 1 / ratiobm[0]); /* parallel-filter output*/

    /*=== Run the inner hair cell (IHC) section: NL function and then lowpass filtering ===*/

    c1vihctmp[n]  = NLogarithm(cihc * c1filterouttmp[n], 0.1, ihcasym);

    c2vihctmp[n] = -NLogarithm(c2filterouttmp[n] * fabs(c2filterouttmp[n]) * cf / 10 * cf / 2e3, 0.2, 1.0); /* C2 transduction output */

    ihcouttmp[n] = IhcLowPass(c1vihctmp[n] + c2vihctmp[n], binwidth, 3800, n, 1.0, 7);

    /*====== Run the synapse model ======*/

    synouttmp[n] = Synapse(ihcouttmp[n], binwidth, cf, spont, n);

  }
  ;  /* End of the loop */

  /* Adjust total path delay to all signals after BM */

  delay      = delay_cat(cf, species);
  delaypoint = __max(0, (int) ceil(delay / binwidth));

  for (i = 0;i < delaypoint;i++) {
    synout[i] = spont;
  };
  for (i = delaypoint;i < totalstim;i++) {
    //   c1filterout[i]  = c1filterouttmp[i-delaypoint];
    //   c2filterout[i]  = c2filterouttmp[i-delaypoint];
    //   c1vihc[i] = c1vihctmp[i-delaypoint];
    //   c2vihc[i] = c2vihctmp[i-delaypoint];
    //   ihcout[i] = ihcouttmp[i-delaypoint];
    synout[i] = synouttmp[i-delaypoint];
  };

  /* Spike Generation */

  //  nspikes = SpikeGenerator(synout, binwidth, totalstim, nrep, sptime);
  //
  //  for(i = 0; i < nspikes; i++)
  //  {
  //   ipst = (int) (fmod(sptime[i],binwidth*totalstim) / binwidth);
  //         psth[ipst] = psth[ipst] + 1;
  //  };


  /* Freeing dynamic memory allocated earlier */

  freevector(c1filterouttmp);
  freevector(c2filterouttmp);
  freevector(c1vihctmp);
  freevector(c2vihctmp);
  freevector(ihcouttmp);
  freevector(synouttmp);

  freevector(mey1);
  freevector(mey2);
  freevector(mey3);

  freevector(meout);
  freevector(grd);
  freevector(tmpgain);

} /* End of the SingleAN function */


/*---------------------------------------------------------------------------------
  Synapse model: if the time resolution is not small enough, the concentration of
  the immediate pool could be as low as negative, at this time there is an alert 
  message print out and the concentration is set at saturated level  
  --------------------------------------------------------------------------------*/
double Synapse(double x, double binwidth, double cf, double spont, int n)
{
  static double synstrength, synslope, CI, CL, PG, CG, VL, PL, VI;

  double cf_factor, PImax, kslope, Ass, Asp, TauR, TauST, Ar_Ast, PTS, Aon, AR, AST, Prest, gamma1, gamma2, k1, k2;
  double VI0, VI1, alpha, beta, theta1, theta2, theta3, vsat, tmpst, tmp, PPI, CIlast, temp, out;
  double cfsat, cfslope, cfconst;

  if (n == 0) {
    /*        cf_factor = __min(1e3,pow(10,0.29*cf/1e3 + 0.4));
     */
    if (spont >= 50)
      cf_factor = __min(1e3, pow(10, 0.29 * cf / 1e3 + 0.4));
    else {
      cfslope = pow(spont, 0.19) * pow(10, -0.87);
      cfconst = 0.1 * pow(log10(spont), 2) + 0.56 * log10(spont) - 0.84;
      cfsat = pow(10, cfslope * 8965.5 / 1e3 + cfconst);    /*find saturation at saturation freq: 8965.5 Hz*/
      cf_factor = __min(cfsat, pow(10, cfslope * cf / 1e3 + cfconst));
    }
    ;                                       /*added by Tim Zeyl June 14 2006*/

    PImax  = 0.6;                /* PI2 : Maximum of the PI(PI at steady state) */
    kslope = (1 + 50.0) / (5 + 50.0) * cf_factor * 20.0 * PImax;

    Ass    = 350;                /* Steady State Firing Rate eq.10 */
    Asp    = spont;              /* Spontaneous Firing Rate  eq.10 */
    TauR   = 2e-3;               /* Rapid Time Constant eq.10 */
    TauST  = 60e-3;              /* Short Time Constant eq.10 */
    Ar_Ast = 6;                  /* Ratio of Ar/Ast */
    PTS    = 1.0 + 9.0 * 50.0 / (9.0 + 50.0);    /* Peak to Steady State Ratio, characteristic of PSTH */

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

    if (kslope >= 0)  vsat    = kslope + Prest;
    tmpst  = log(2) * vsat / Prest;
    if (tmpst < 400) synstrength = log(exp(tmpst) - 1);
    else synstrength = tmpst;
    synslope = Prest / log(2) * synstrength;

    if (spont < 0) spont = 50;
  };

  tmp = synstrength * x;
  if (tmp < 400) tmp = log(1 + exp(tmp));
  PPI = synslope / synstrength * tmp;

  CIlast = CI;
  CI = CI + (binwidth / VI) * (-PPI * CI + PL * (CL - CI));
  CL = CL + (binwidth / VL) * (-PL * (CL - CIlast) + PG * (CG - CL));

  if (CI < 0) {
    temp = 1 / PG + 1 / PL + 1 / PPI;
    CI = CG / (PPI * temp);
    CL = CI * (PPI + PL) / PL;
  };

  out = CI * PPI;
  return(out);
}


