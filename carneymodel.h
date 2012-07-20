#ifndef _carney_model_
#define _carney_model_

#include "complex.h"

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



double cochlea_f2x(int species, double f);
double cochlea_x2f(int species, double x);
double Get_tauwb(double cf, double CAgain, int order, double *taumax, double *taumin);
double Get_taubm(double cf, double CAgain, double taumax, double *bmTaumax, double *bmTaumin, double *ratio);
double C1ChirpFilt(double x, double binwidth, double cf, int n, double taumax, double rsigma);
double C2ChirpFilt(double xx, double binwidth, double cf, int n, double taumax, double fcohc);
double WbGammaTone(double x, double binwidth, double centerfreq, int n, double tau, double gain, int order);
double gain_groupdelay(double binwidth, double centerfreq, double cf, double tau, int *grdelay);
double delay_cat(double cf, int species);
double Boltzman(double x, double asym, double s0, double s1, double x1);
double OhcLowPass(double x, double binwidth, double Fc, int n, double gain, int order);
double IhcLowPass(double x, double binwidth, double Fc, int n, double gain, int order);
double NLafterohc(double x, double taumin, double taumax, double asym);
double NLogarithm(double x, double slope, double asym);

#endif
