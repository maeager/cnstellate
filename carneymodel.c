

#include "carneymodel.h"

/* --------------------------------------------------------------------------
** Calculate the location on Basilar Membrane from best frequency
*/
double cochlea_f2x(int species, double f)
{
  double x;
  switch (species) {
  case 0: /*/human */
    x = (1.0 / 0.06) * log10((f / 165.4) + 0.88);
    break;
  case 2: /* Rat data from earlab.bu.edu */
    x = (8.03 / 0.928) * log10(1.0 + f / 7613.3);
    break;
  default:
  case 1: /*/cat */
    x = 11.9 * log10(0.80 + f / 456.0);
    break;
  };
  return(x);
}
/* --------------------------------------------------------------------------------
** Calculate the best frequency from the location on basilar membrane
*/
double cochlea_x2f(int species, double x)
{
  double f;
  switch (species) {
  case 0: /* human */
    if((x>35)||(x<0)) hoc_execerror("BM distance out of human range, [in cochlea_x2f(...)]",0); 
    f = 165.4 * (pow(10, (0.06 * x)) - 0.88);
    break;
  case 2: //rat
    f = 7613.3 * (pow(10, 0.924 * x / 8.03) - 1.0);
    break;
  default:
  case 1: /*/cat */
    f = 456.0 * (pow(10, x / 11.9) - 0.80);
    break;
  };
  return(f);
}


/* -------------------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------------------- */
/** Get TauMax, TauMin for the tuning filter. The TauMax is determined by the bandwidth/Q10
    of the tuning filter at low level. The TauMin is determined by the gain change between high
    and low level */

double Get_tauwb(double cf, double CAgain, int order, double *taumax, double *taumin)
{

  double Q10, bw, ratio;

  ratio = pow(10, (-CAgain / (20.0 * order)));  /* ratio of TauMin/TauMax according to the gain, order */

  /* Q10 = pow(10,0.4708*log10(cf/1e3)+0.5469); */ /* 75th percentile */
  Q10 = pow(10, 0.4708 * log10(cf / 1e3) + 0.4664); /* 50th percentile */
  /* Q10 = pow(10,0.4708*log10(cf/1e3)+0.3934); */ /* 25th percentile */

  bw     = cf / Q10;
  taumax[0] = 2.0 / (TWOPI * bw);

  taumin[0]   = taumax[0] * ratio;

  return 0;
}
/* -------------------------------------------------------------------------------------------- */
double Get_taubm(double cf, double CAgain, double taumax, double *bmTaumax, double *bmTaumin, double *ratio)
{
  double factor, bwfactor;

  bwfactor = 0.7;
  factor   = 2.5;

  ratio[0]  = pow(10, (-CAgain / (20.0 * factor)));

  bmTaumax[0] = taumax / bwfactor;
  bmTaumin[0] = bmTaumax[0] * ratio[0];
  return 0;
}
/* -------------------------------------------------------------------------------------------- */
/** Pass the signal through the signal-path C1 Tenth Order Nonlinear Chirp-Gammatone Filter */

double C1ChirpFilt(double x, double binwidth, double cf, int n, double taumax, double rsigma)
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
  ipw    = 1.01 * cf * TWOPI-50;
  ipb    = 0.2343 * TWOPI * cf-1104;
  rpa    = pow(10, log10(cf) * 0.9 + 0.55) + 2000;
  pzero  = pow(10, log10(cf) * 0.7 + 1.6) + 500;

  /*===============================================================*/

  order_of_pole    = 10;
  half_order_pole  = order_of_pole / 2;
  order_of_zero    = half_order_pole;

  fs_bilinear = TWOPI * cf / tan(TWOPI * cf * binwidth / 2);
  rzero       = -pzero;
  CF          = TWOPI * cf;

  if (n == 0) 
    {
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

  if (p[1].x > 0.0) printf("The system becomes unstable.\n");

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

  if (rzero > 0.0) printf("The zeros are in the right-half plane.\n");

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

/* -------------------------------------------------------------------------------------------- */
/** Parallelpath C2 filter: same as the signal-path C1 filter with the OHC completely impaired */

double C2ChirpFilt(double xx, double binwidth, double cf, int n, double taumax, double fcohc)
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

  fs_bilinear = TWOPI * cf / tan(TWOPI * cf * binwidth / 2);
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

  if (p[1].x > 0.0) printf("The system becomes unstable.\n");

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
  if (rzero > 0.0) printf("The zeros are in the right-hand plane.\n");
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

double WbGammaTone(double x, double binwidth, double centerfreq, int n, double tau, double gain, int order)
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

  delta_phase = -TWOPI * centerfreq * binwidth;
  wbphase += delta_phase;

  dtmp = tau * 2.0 / binwidth;
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

/* -------------------------------------------------------------------------------------------- */
/** Calculate the gain and group delay for the Control path Filter */

double gain_groupdelay(double binwidth, double centerfreq, double cf, double tau, int *grdelay)
{
  double tmpcos, dtmp2, c1LP, c2LP, tmp1, tmp2, wb_gain;

  tmpcos = cos(TWOPI * (centerfreq - cf) * binwidth);
  dtmp2 = tau * 2.0 / binwidth;
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

double delay_cat(double cf, int species)
{
  /* DELAY THE WAVEFORM (delay buf1, tauf, ihc for display purposes)  */
  /* Note: Latency vs. CF for click responses is available for Cat only (not human) */
  /* Use original fit for Tl (latency vs. CF in msec) from Carney & Yin '88
     and then correct by .75 cycles to go from PEAK delay to ONSET delay */
  /* from Carney and Yin '88 */
  double A0, A1, x, delay;

  switch (species){
  case 0: /*Human*/
    x = cochlea_f2x(species, cf); /* human mapping */
    if ((x<5)||(x>35)) printf("delay_cat: x out of range.\n");
    delay = 4.915 + 0.3631*exp(0.11324*x);   /* 5<x<35 m */
    break;

  case 9: /* Cat */
  case 1:
  default:
    /*  
	delay (msec) = A0 * exp( -x /A1 ) * 1e-3 - 1.0/CF 
	where A0 = 8.13(ms), A1 = 6.49 (nm)
    */
    A0    = 3.0;   
    A1    = 12.5;
    x = cochlea_f2x(species, cf); /* cat mapping */
    delay = A0 * exp(-x / A1) * 1e-3;
  }
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

double OhcLowPass(double x, double binwidth, double Fc, int n, double gain, int order)
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

  c = 2.0 / binwidth;
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

double IhcLowPass(double x, double binwidth, double Fc, int n, double gain, int order)
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


  C = 2.0 / binwidth;
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
  if (out < taumin) out = taumin; if (out > taumax) out = taumax;
  return(out);
}
/* -------------------------------------------------------------------------------------------- */
/* Get the output of the IHC Nonlinear Function (Logarithmic Transduction Functions) */

double NLogarithm(double x, double slope, double asym)
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
