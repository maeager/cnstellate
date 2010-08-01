/*
 * Run Time Utility
 *
 */
#ifndef _runmodel_h
#define _runmodel_h

//#include <stdlib.h>
#include <stdio.h>
#include <math.h>

typedef struct _Tstim T_stim;

struct _Tstim {
    /*/General information about the model */
    double tdres, cf, spont;
    int species, model;
    /*/About filter banks */
    int banks; /*/for filterbank */
    double delx;  /*/for filter bank */
    double cfhi, cflo, cfinc;
    /*/About stimulus */
    int stim; /*/stimulus type */
    char wavefile[100]; /*/get wave from the file */
    double reptim; /*/repetition time */
    int trials; /*/for spike generator */

    int nstim, nrep; /*/runtime usage */
    double *buf;
};


#endif









#ifndef _COMPLEX_H
#define _COMPLEX_H

/* COMPLEX.H header file
 * use for complex arithmetic in C
 (part of them are from "C Tools for Scientists and Engineers" by L. Baker)
*/
extern double cos(double);
extern double sin(double);
struct __COMPLEX {
    double x, y;
};
typedef struct __COMPLEX COMPLEX;

/* for below, X, Y are complex structures, and one is returned */

/*//real part of the complex multiplication */
#define CMULTR(X,Y) ((X).x*(Y).x-(X).y*(Y).y)
/*//image part of the complex multiplication */
#define CMULTI(X,Y) ((X).y*(Y).x+(X).x*(Y).y)
/*// used in the Division : real part of the division */
#define CDRN(X,Y) ((X).x*(Y).x+(Y).y*(X).y)
/*// used in the Division : image part of the division */
#define CDIN(X,Y) ((X).y*(Y).x-(X).x*(Y).y)
/*// used in the Division : denumerator of the division */
#define CNORM(X) ((X).x*(X).x+(X).y*(X).y)
/*//real part of the complex */
#define CREAL(X) ((double)((X).x))
/*//conjunction value */
#define CONJG(z,X) {(z).x=(X).x;(z).y= -(X).y;}
/*//conjunction value */
#define CONJ(X) {(X).y= -(X).y;}
/*//muliply : z could not be same variable as X or Y, same rule for other Macro */
#define CMULT(z,X,Y) {(z).x=CMULTR((X),(Y));(z).y=CMULTI((X),(Y));}
/*//division */
#define CDIV(z,X,Y){double d=CNORM(Y); (z).x=CDRN(X,Y)/d; (z).y=CDIN(X,Y)/d;}
/*//addition */
#define CADD(z,X,Y) {(z).x=(X).x+(Y).x;(z).y=(X).y+(Y).y;}
/*//subtraction */
#define CSUB(z,X,Y) {(z).x=(X).x-(Y).x;(z).y=(X).y-(Y).y;}
/*//assign */
#define CLET(to,from) {(to).x=(from).x;(to).y=(from).y;}
/*//abstract value(magnitude) */
#define cabs(X) sqrt((X).y*(X).y+(X).x*(X).x)
/*//real to complex */
#define CMPLX(X,real,imag) {(X).x=(real);(X).y=(imag);}
/*//multiply with real */
#define CTREAL(z,X,real) {(z).x=(X).x*(real);(z).y=(X).y*(real);}

#define CEXP(z,phase) {(z).x = cos(phase); (z).y = sin(phase); }
/* implementation using function : for compatibility */
COMPLEX compdiv(COMPLEX ne, COMPLEX de);
COMPLEX compexp(double theta);
COMPLEX compmult(double scalar, COMPLEX compnum);
COMPLEX compprod(COMPLEX compnum1, COMPLEX compnum2);
COMPLEX comp2sum(COMPLEX summand1, COMPLEX summand2);
COMPLEX comp3sum(COMPLEX summand1, COMPLEX summand2, COMPLEX summand3);
COMPLEX compsubtract(COMPLEX complexA, COMPLEX complexB);
double  REAL(COMPLEX compnum);

#endif
/**
   This file defines several filters that is used in the program
 */
#ifndef _FILTER_H
#define _FILTER_H
//#include <math.h>
//#include <stdlib.h>
//#include "complex.h"

#ifndef TWOPI
#define TWOPI 6.2831853
#endif
/*// Maximum order of the filter */
#define MAX_ORDER 10

/*/######################################################################## */
/* interface the users */
typedef struct __LowPass TLowPass;
typedef struct __HighPass THighPass;
typedef struct __GammaTone TGammaTone;
/*// Construct the filter with time_resolustion(1/F_s), cutoff frequency, gain, filter order */
TLowPass* getLowPass(double _tdres, double _Fc, double _gain, int _LPorder);
void initLowPass(TLowPass* p, double _tdres, double _Fc, double _gain, int _LPorder);
/*// */
TGammaTone* getGammaTone(double _tdres, double _Fshift, double _tau, double _gain, int _order);
void initGammaTone(TGammaTone* p, double _tdres, double _Fshift, double _tau, double _gain, int _order);
/*// Init the filter with time_resolution(1/F_s), cut-off freq., gain,order */
THighPass* getHighPass(double _tdres, double _Fc, double _gain, int _HPorder);
void initHighPass(THighPass* p, double _tdres, double _Fc, double _gain, int _HPorder);
/* /########################################################################################## */
/*/ ----------------------------------------------------------------------------
 **
   Lowpass Filter
*/

struct __LowPass {
    /*//input one signal to the filter and get the output */
    double(*run)(TLowPass *p, double x);
    /*//input multiple signal to the filter and get the multiple output */
    void (*run2)(TLowPass *p, const double *in, double *out, const int length);

    /*/time-domain resolution,cut-off frequency,gain, filter order */
    double tdres, Fc, gain;
    int Order;
    /*/ parameters used in calculation */
    double c1LP, c2LP, hc[MAX_ORDER], hcl[MAX_ORDER];
};

/*/ ----------------------------------------------------------------------------
 **
  Highpass Filter
*/
struct __HighPass {

    /*//input one signal to the filter and get the output */
    double(*run)(THighPass *p, double x);
    /*//input multiple signal to the filter and get the multiple output */
    void (*run2)(THighPass *p, const double *in, double *out, const int length);
    /*time-domain resolution, cut-off frequency, gain of the filter */
    double tdres, Fc, gain;
    int Order;       /*/order of the hi-pass */
    double c1HP, c2HP, hc[MAX_ORDER], hcl[MAX_ORDER];
};

/*/ ----------------------------------------------------------------------------
 **
   Gammatone filter
 */
struct __GammaTone {
    /*//input one signal to the filter and get the output */
    double(*run)(TGammaTone *p, double x);
    /*//input multiple signal to the filter and get the multiple output */
    void (*run2)(TGammaTone *p, const double *in, double *out, const int length);

    double phase;
    /* Cutoff Freq(tau), Shift Freq, ... */
    double tdres, tau;
    double F_shift, delta_phase;
    double gain, c1LP, c2LP;
    COMPLEX gtf[MAX_ORDER], gtfl[MAX_ORDER];
    int Order;

    /*// Set the tau of the gammatone filter, this is useful for time-varying filter */
    void (*settau)(TGammaTone *p, double _tau);
};
#endif
#ifndef _HC_H
#define _HC_H
//#include <math.h>
//#include <stdlib.h>
//#include "filters.h"

typedef struct __NonLinear TNonLinear;
typedef struct __HairCell THairCell;
THairCell* getHairCell(double tdres, double cutoff, int order);
TNonLinear* getAfterOhcNL(double taumin, double taumax, double dc, double minR);

/*/AfterOHC NonLinear */
double runAfterOhcNL(TNonLinear* p, double x);
void runAfterOhcNL2(TNonLinear* p, const double *in, double *out, const int length);
/*/OHC NonLinear function */
void init_boltzman(TNonLinear* p, double _corner, double _slope, double _strength,
                   double _x0, double _s0, double _x1, double _s1, double _asym);
double runBoltzman(TNonLinear *p, double x);
void runBoltzman2(TNonLinear *p, const double *in, double *out, const int length);
/*/IHC NonLinear function */
double runIHCNL(TNonLinear * p, double x);
void runIHCNL2(TNonLinear* p, const double *in, double *out, const int length);
/*/IHCPPI Nonlinear */
double runIHCPPI(TNonLinear* p, double x);
void runIHCPPI2(TNonLinear* p, const double *in, double *out, const int length);
/*/HairCell run */
double runHairCell(THairCell* p, double x);
void runHairCell2(THairCell* p, const double *in, double *out, const int length);

/*
 *
 * ################# Structure Implementation #######################
 *
 */
struct __NonLinear {
    double(*run)(TNonLinear* p, double x);
    void (*run2)(TNonLinear *p, const double *in, double *out, const int length);
    /*/For OHC Boltzman */
    double p_corner, p_slope, p_strength, x0, s0, x1, s1, shift;
    double Acp, Bcp, Ccp;
    /*/For AfterOHCNL */
    double dc, minR, A, TauMax, TauMin; /*/s0 also used by AfterOHCNL */
    /*/For IHC nonlinear function */
    double A0, B, C, D;
    /*/For IHCPPI nonlinear */
    double psl, pst, p1, p2;
};

struct __HairCell {
    double(*run)(THairCell *p, double x);
    void (*run2)(THairCell* p, const double *in, double *out, const int length);
    TLowPass hclp;
    TNonLinear hcnl;
    /*/Boltzman Like NonLinear */
    void (*setOHC_NL)(THairCell*, double, double, double, double, double, double, double, double);
    void (*setIHC_NL)();
};

#endif

#ifndef _CMPA_H
#define _CMPA_H
//#include "filters.h"
//#include "hc.h"
//#include "synapse.h"
//#include "complex.h"

#ifndef TWOPI
#define TWOPI 6.2831853
#endif
#define Broad_ALL  0x80
#define Linear_ALL  0x40
#define NonLinear_ALL  0x20

#define FeedBack_NL  (NonLinear_ALL|0x01)
#define FeedForward_NL  (NonLinear_ALL|0x02)
#define Broad_Linear  (Broad_ALL|Linear_ALL|0x01)
#define Sharp_Linear  (Linear_ALL|0x02)
#define Broad_Linear_High (Broad_ALL|Linear_ALL|0x03)

/*/############################################################################## */
typedef struct __AuditoryNerve TAuditoryNerve;
typedef struct __BasilarMembrane TBasilarMembrane;
double cochlea_f2x(int species, double f);
double cochlea_x2f(int species, double x);

double runAN(TAuditoryNerve* p, double x);
void runAN2(TAuditoryNerve *p, const double *in, double *out, const int length);

void initAuditoryNerve(TAuditoryNerve *p, int model, int species, double tdres, double cf, double spont);
void initBasilarMembrane(TBasilarMembrane* bm, int model, int species, double tdres, double cf);

/*/############################################################################## */

/** The class the define the basic structure of the time-varing filter

    the class consists of following components: \n
    1. tuning filter(bmfilter), the tau is controlled by the control path\n
    2. wideband pass filter(wbfilter)\n
    3. outer hair cell model(ohc)\n
    4. nonlinear function after the outer haircell(afterohc)
 */
struct __BasilarMembrane { /* class of basilar membrane */

    double(*run)(TBasilarMembrane *p, double x);
    void (*run2)(TBasilarMembrane *p, const double *in, double *out, const int length);

    int bmmodel; /* determine if the bm is broad_linear, sharp_linear or other */
    double tdres;
    int bmorder, wborder;

    double tau, TauMax, TauMin;
    double TauWB, TauWBMin;
    double A, B;
    /* --------Model -------------- */
    TGammaTone bmfilter; /*/NonLinear Filter */
    TGammaTone gfagain; /*/Linear Filter */
    TGammaTone wbfilter; /*/Control Path filter */
    THairCell ohc;
    TNonLinear afterohc;
};

/** Class of the auditory nerve fiber, this is a complete model of the fiber

    The class consists of all the parts of the auditory nerve fiber\n
    1. time-varying tuning filter with control path: TBasilarMembrane(bm, 3rd order)\n
    2. the gammatone filter after the 3rd-order nonlinear filter(gffilter)\n
    3. inner hair cell model(ihc)\n
    4. synapse model, from Westman,1986(syn)\n
    5. spike generator, from Carney,1993(sg)
 */
struct __AuditoryNerve {
    /*/ Run Function */
    double(*run)(TAuditoryNerve *p, double x);
    void (*run2)(TAuditoryNerve *p, const double *in,  double *out,  const int length);
    /*/ Model Structor */
    TBasilarMembrane bm;
    THairCell ihc;
    TNonLinear ihcppi; /*/From ihc->ppi */
    Tsynapse syn;
    /*/  TSpikeGenerator sg; */
    /*/ Model Parameters */
    double tdres, cf, spont;
    int species, model;
    /* This parameter indicates if we are using sout only or spikes */
    int ifspike;
};

#endif





#ifndef _synapse_h
#define _synapse_h
typedef struct __tsynapse Tsynapse;
struct __tsynapse {
    double(*run)(Tsynapse *p, double x);
    void (*run2)(Tsynapse *p, const double *in, double *out, const int length);
    double cf, tdres;
    double spont;
    double PTS, Ass, Ar_over_Ast, Pimax, tauR, tauST;
    double Prest, PPIlast, PL, PG, CIrest, CIlast, CLrest, CLlast, CG, VI, VL;

    double Vsat;
};
int runSynapse(Tsynapse *pthis, const double *in, double *out, const int length);
void runsyn_dynamic(Tsynapse *pthis, const double *in, double *out, const int length);
double run1syn_dynamic(Tsynapse *pthis, double x);
int initSynapse(Tsynapse *pthis);
#endif


/**
complex.cpp includes all of the COMPLEX math functions needed for model programs
*/

//#include <stdlib.h>
//#include <math.h>
//#include "complex.h"

/*/divide */
COMPLEX compdiv(COMPLEX ne, COMPLEX de)
{
    double d;
    COMPLEX z;
    d = de.x * de.x + de.y * de.y;
    z.x = (ne.x * de.x + ne.y * de.y) / d;
    z.y = (ne.y * de.x - ne.x * de.y) / d;
    return(z);
}
/*/ this returns a complex number equal to exp(i*theta) */
COMPLEX compexp(double theta)
{
    COMPLEX dummy;
    dummy.x = cos(theta);
    dummy.y = sin(theta);
    return dummy;
}
/*/ Multiply a complex number by a scalar */
COMPLEX compmult(double scalar, COMPLEX compnum)
{
    COMPLEX answer;
    answer.x = scalar * compnum.x;
    answer.y = scalar * compnum.y;
    return answer;
}
/*/ Find the product of 2 complex numbers */
COMPLEX compprod(COMPLEX compnum1, COMPLEX compnum2)
{
    COMPLEX answer;
    answer.x = (compnum1.x * compnum2.x) - (compnum1.y * compnum2.y);
    answer.y = (compnum1.x * compnum2.y) + (compnum1.y * compnum2.x);
    return answer;
}
/*/ add 2 complex numbers */
COMPLEX comp2sum(COMPLEX summand1, COMPLEX summand2)
{
    COMPLEX answer;
    answer.x = summand1.x + summand2.x;
    answer.y = summand1.y + summand2.y;
    return answer;
}
/*/ add three complex numbers */
COMPLEX comp3sum(COMPLEX summand1, COMPLEX summand2, COMPLEX summand3)
{
    COMPLEX answer;
    answer.x = summand1.x + summand2.x + summand3.x;
    answer.y = summand1.y + summand2.y + summand3.y;
    return answer;
}

/*  subtraction: complexA - complexB */
COMPLEX compsubtract(COMPLEX complexA, COMPLEX complexB)
{
    COMPLEX answer;
    answer.x = complexA.x - complexB.x;
    answer.y = complexA.y - complexB.y;
    return answer;
}
/* Get the real part of the complex */
double REAL(COMPLEX compnum)
{
    return compnum.x;
};

/* Get the imaginary part of the complex */
double IMAG(COMPLEX compnum)
{
    return(compnum.y);
}

/* Get the conjugate of the complex signal */
COMPLEX compconj(COMPLEX complexA)
{
    COMPLEX answer;
    answer.x = complexA.x;
    answer.y = -complexA.y;
    return (answer);
}


//#include <stdlib.h>
//#include "filters.h"
//#include "complex.h"
/**
    The low-pass filter is in the form of cascade of first order lowpass filter\\
    Each of them is implemented as y(i) = c1*y(i-1)+c2*(x(i)+x(i-1))

    This function initializes the filter:\\
    1. initialize/reset the filter\\
    2. calculation c1 c2 thru bilinear transformation\\
    3. set the gain and order of the filter\\

    @author Xuedong Zhang
 */
/**
 * Internal function used by TLowPass
 */
double runLowPass(TLowPass *p, double x);
void runLowPass2(TLowPass *p, const double *in, double *out, const int length);
void initLowPass(TLowPass* res, double _tdres, double _Fc, double _gain, int _LPorder);
TLowPass* getLowPass(double _tdres, double _Fc, double _gain, int _LPorder)
{
    TLowPass *res = (TLowPass*) hoc_Emalloc(sizeof(TLowPass));
    initLowPass(res, _tdres, _Fc, _gain, _LPorder);
    return(res);
};

void initLowPass(TLowPass* res, double _tdres, double _Fc, double _gain, int _LPorder)
{
    double c;
    int i;
    res->tdres = _tdres;
    c = 2.0 / _tdres;
    res->Fc = _Fc;
    res->Order = _LPorder;
    res->c1LP = (c - TWOPI * _Fc) / (c + TWOPI * _Fc);
    res->c2LP = TWOPI * _Fc / (TWOPI * _Fc + c);
    for (i = 0; i <= res->Order; i++) res->hc[i] = res->hcl[i] = 0.0;
    res->gain = _gain;

    res->run = runLowPass;
    res->run2 = runLowPass2;
    return;
};
/**
  This function runs the low-pass filter
   @author Xuedong Zhang
 */
double runLowPass(TLowPass *p, double x)
{
    register int i;
    register int pOrder = p->Order;
    p->hc[0] = x * p->gain;

    for (i = 0; i < pOrder; i++)
        p->hc[i+1] = p->c1LP * p->hcl[i+1] + p->c2LP * (p->hc[i] + p->hcl[i]);

    for (i = 0; i <= pOrder;i++) p->hcl[i] = p->hc[i];

    return(p->hc[pOrder]);
};
/**
  This function runs the low-pass filter
   @author Xuedong Zhang
 */
void runLowPass2(TLowPass *p, const double *in, double *out, const int length)
{
    register int loopSig, loopLP;
    int pOrder = p->Order;
    double *hc, *hcl, c1LP, c2LP;
    double gain;
    gain = p->gain;
    c1LP = p->c1LP;
    c2LP = p->c2LP;
    hc = p->hc;
    hcl = p->hcl;
    for (loopSig = 0; loopSig < length; loopSig++) {
        hc[0] = in[loopSig] * gain;

        for (loopLP = 0; loopLP < pOrder; loopLP++)
            hc[loopLP+1] = c1LP * hcl[loopLP+1] + c2LP * (hc[loopLP] + hcl[loopLP]);

        for (loopLP = 0; loopLP <= pOrder;loopLP++) hcl[loopLP] = hc[loopLP];
        out[loopSig] = hc[pOrder];
    };
    return;
};

/**

    The high-pass filter is in the form of a cascade of first order highpass filter\\
    Each of them is implemented as y(i) = c1*y(i-1)+c2*(x(i)-x(i-1))

    This function does several things:\\
    1. initialize/reset the filter\\
    2. calculation c1 c2 thru bilinear transformation\\
    3. set the gain and order of the filter\\

    @author Xuedong Zhang
*/

/**
 * Internal function used by THighPass
 */
double runHighPass(THighPass *p, double x);
void runHighPass2(THighPass *p, const double *in, double *out, const int length);
void initHighPass(THighPass* p, double _tdres, double _Fc, double _gain, int _HPorder);

THighPass* getHighPass(double _tdres, double _Fc, double _gain, int _HPorder)
{
    THighPass* res = (THighPass*) hoc_Emalloc(sizeof(THighPass));
    initHighPass(res, _tdres, _Fc, _gain, _HPorder);
    return(res);
};

void initHighPass(THighPass* p, double _tdres, double _Fc, double _gain, int _HPorder)
{
    double c;
    int i;
    THighPass *res = p;
    res->tdres = _tdres;
    c = 2.0 / _tdres;
    res->Fc = _Fc;
    res->Order = _HPorder;
    res->c1HP = (c - TWOPI * _Fc) / (c + TWOPI * _Fc);
    res->c2HP =  c / (TWOPI * _Fc + c);
    for (i = 0; i <= res->Order; i++) res->hc[i] = res->hcl[i] = 0.0;
    res->gain = _gain;

    res->run = runHighPass;
    res->run2 = runHighPass2;
    return;
};
/**

   This function get the filtering output of the high-pass filter
   @author Xuedong Zhang
 */
double runHighPass(THighPass *p, double x)
{
    register int i;
    int pOrder = p->Order;
    p->hc[0] = x * p->gain;

    for (i = 0; i < pOrder;i++)
        p->hc[i+1] = p->c1HP * p->hcl[i+1] + p->c2HP * (p->hc[i] - p->hcl[i]);

    for (i = 0; i <= pOrder;i++) p->hcl[i] = p->hc[i];

    return(p->hc[pOrder]);
};
/**

   This function get the filtering output of the high-pass filter
   @author Xuedong Zhang
 */
void runHighPass2(THighPass *p, const double *in, double *out, const int length)
{
    register int loopSig, loopHP;
    int pOrder = p->Order;
    double *hc = p->hc;
    double *hcl = p->hcl;
    double c1HP = p->c1HP;
    double c2HP = p->c2HP;
    double gain = p->gain;

    for (loopSig = 0; loopSig < length; loopSig++) {
        hc[0] = in[loopSig] * gain;

        for (loopHP = 0; loopHP < pOrder;loopHP++)
            hc[loopHP+1] = c1HP * hcl[loopHP+1] + c2HP * (hc[loopHP] - hcl[loopHP]);

        for (loopHP = 0; loopHP <= pOrder;loopHP++) hcl[loopHP] = hc[loopHP];

        out[loopSig] = hc[pOrder];
    };
    return;
};

/**
    The gammatone filter is in the form of cascade of first order lowpass filter
    shifted by the center frequency (from Carney,1993)\\
    Each of the low pass filter is implemented as y(i) = c1*y(i-1)+c2*(x(i)-x(i-1))

    This function does several things:\\
    1. initialize/reset the filter\\
    2. calculation c1 c2 thru bilinear transformation from tau\\
    3. set the gain and order of the filter\\

    @author Xuedong Zhang
*/
void setGammaToneTau(TGammaTone *p, double tau);
double runGammaTone(TGammaTone *p, double x);
void runGammaTone2(TGammaTone *p, const double *in, double *out, const int length);
void initTGammaTone(TGammaTone* res, double _tdres, double _Fshift, double _tau, double _gain, int _order);
TGammaTone* getGammaTone(double _tdres, double _Fshift, double _tau, double _gain, int _order)
{
    TGammaTone *res = (TGammaTone*) hoc_Emalloc(sizeof(TGammaTone));
    initGammaTone(res, _tdres, _Fshift, _tau, _gain, _order);
    return(res);
};

void initGammaTone(TGammaTone* res, double _tdres, double _Fshift, double _tau, double _gain, int _order)
{
    double c;
    int i;
    res->tdres = _tdres;
    res->F_shift = _Fshift;
    res->delta_phase = -TWOPI * _Fshift * _tdres;
    res->phase = 0;
    res->tau = _tau;

    c = 2.0 / _tdres; /* for bilinear transformation */
    res->c1LP = (_tau * c - 1) / (_tau * c + 1);
    res->c2LP = 1 / (_tau * c + 1);
    res->gain = _gain;
    res->Order = _order;
    for (i = 0; i <= res->Order; i++)
        res->gtf[i].x = res->gtfl[i].x = res->gtf[i].y = res->gtfl[i].y = 0.0;

    res->run = runGammaTone;
    res->run2 = runGammaTone2;
    res->settau = setGammaToneTau;
    return;
};
/**

   Reset the tau of the gammatone filter\\
   it recalculate the c1 c2 used by the filtering function
 */
void setGammaToneTau(TGammaTone *p, double tau)
{
    double dtmp;
    p->tau = tau;
    dtmp = tau * 2.0 / p->tdres;
    p->c1LP = (dtmp - 1) / (dtmp + 1);
    p->c2LP = 1.0 / (dtmp + 1);
};
/**

   Pass the signal through the gammatone filter\\
   1. shift the signal by centeral frequency of the gamma tone filter\\
   2. low pass the shifted signal \\
   3. shift back the signal \\
   4. take the real part of the signal as output
   @author Xuedong Zhang
 */
double runGammaTone(TGammaTone *p, double x)
{
    int i, j;
    double out;
    COMPLEX c1, c2, c_phase;
    x *= p->gain;
    p->phase += p->delta_phase;

    CEXP(c_phase, p->phase); /*/ FREQUENCY SHIFT */
    CTREAL(p->gtf[0], c_phase, x);
    for (j = 1; j <= p->Order; j++) {    /*/ IIR Bilinear transformation LPF */
        CADD(c1, p->gtf[j-1], p->gtfl[j-1]);
        CTREAL(c2, c1, p->c2LP);
        CTREAL(c1, p->gtfl[j], p->c1LP);
        CADD(p->gtf[j], c1, c2);
    };
    CONJ(c_phase); /*/ FREQ SHIFT BACK UP */
    CMULT(c1, c_phase, p->gtf[p->Order]);
    out = CREAL(c1);
    for (i = 0; i <= p->Order;i++) p->gtfl[i] = p->gtf[i];
    return(out);
};

void runGammaTone2(TGammaTone *p, const double *in, double *out, const int length)
{
    int register loopSig, loopGT;
    double x;
    COMPLEX c1, c2, c_phase;

    for (loopSig = 0; loopSig < length; loopSig++) {
        x = p->gain * in[loopSig];
        p->phase += p->delta_phase;

        CEXP(c_phase, p->phase); /*/ FREQUENCY SHIFT */
        CTREAL(p->gtf[0], c_phase, x);
        for (loopGT = 1; loopGT <= p->Order; loopGT++) {    /*/ IIR Bilinear transformation LPF */
            CADD(c1, p->gtf[loopGT-1], p->gtfl[loopGT-1]);
            CTREAL(c2, c1, p->c2LP);
            CTREAL(c1, p->gtfl[loopGT], p->c1LP);
            CADD(p->gtf[loopGT], c1, c2);
        };
        CONJ(c_phase); /* FREQ SHIFT BACK UP */
        CMULT(c1, c_phase, p->gtf[p->Order]);
        for (loopGT = 0; loopGT <= p->Order;loopGT++) p->gtfl[loopGT] = p->gtf[loopGT];

        out[loopSig] = CREAL(c1);
    };
    return;
};

//#include <math.h>
//#include <stdlib.h>
//#include "filters.h"
//#include "hc.h"

/*
 *
 * ############## Implemenation of OHC_NL ######################
 *
 */
double runBoltzman(TNonLinear *p, double x);
void runBoltzman2(TNonLinear *p, const double *in, double *out, const int length);
void init_boltzman(TNonLinear* p, double _corner, double _slope, double _strength,
                   double _x0, double _s0, double _x1, double _s1, double _asym)
{
    /* asym is the ratio of positive Max to negative Max*/

    p->p_corner = _corner;
    p->p_slope = _slope;
    p->p_strength = _strength;
    p->x0 = _x0;
    p->s0 = _s0;
    p->x1 = _x1;
    p->s1 = _s1;
    if (_asym < 0) p->shift = 1.0 / (1.0 + exp(_x0 / _s0) * (1 + exp(_x1 / _s1)));
    else {
        p->shift = 1.0 / (1.0 + _asym);
        p->x0 = _s0 * log((1.0 / p->shift - 1) / (1 + exp(_x1 / _s1)));
        p->Bcp = p->p_slope / p->p_strength;
        p->Acp = exp((-p->p_corner - 20 * log10(20e-6)) * p->p_strength);
        p->Ccp = 20 * p->p_strength / log(10);
    };
    p->run = runBoltzman;
    p->run2 = runBoltzman2;
    return;
};

double runBoltzman(TNonLinear *p, double x)
{
    /*// get the output of the first nonlinear function */
    double xx, out;
    xx = fabs(x);
    if (x > 0)
        xx = p->Bcp * log(1 + p->Acp * pow(xx, p->Ccp));
    else if (x < 0)
        xx = -p->Bcp * log(1 + p->Acp * pow(xx, p->Ccp));
    else xx = 0;

    /*// get the output of the second nonlinear function(Boltzman Function) */
    out = 1.0 / (1.0 + exp(-(xx - p->x0) / p->s0) * (1.0 + exp(-(xx - p->x1) / p->s1))) - p->shift;
    return(out / (1 - p->shift));
};

void runBoltzman2(TNonLinear *p, const double *in, double *out, const int length)
{
    /*// get the output of the first nonlinear function */
    double x, xx, out1;
    int register i;
    for (i = 0; i < length; i++) {
        x = in[i];
        xx = fabs(x);
        if (x > 0)
            xx = p->Bcp * log(1 + p->Acp * pow(xx, p->Ccp));
        else if (x < 0)
            xx = -p->Bcp * log(1 + p->Acp * pow(xx, p->Ccp));
        else xx = 0;
        /*// get the output of the second nonlinear function(Boltzman Function) */
        out1 = 1.0 / (1.0 + exp(-(xx - p->x0) / p->s0) * (1.0 + exp(-(xx - p->x1) / p->s1))) - p->shift;
        out[i] = out1 / (1 - p->shift);
    };
    return;
};
/*
 *
 * ################## Implementation of THairCell ###################
 *
 */
double runHairCell(THairCell *p, double x);
void runHairCell2(THairCell *p, const double *in, double *out, const int length);

double runHairCell(THairCell *p, double x)
{
    double y;
    y = p->hcnl.run(&(p->hcnl), x);
    return(p->hclp.run(&(p->hclp), y));
};

void runHairCell2(THairCell *p, const double *in, double *out, const int length)
{
    p->hcnl.run2(&(p->hcnl), in, out, length);
    p->hclp.run2(&(p->hclp), out, out, length);
};

/*
 *
 * ############## Implementation InterFace about AfterOHC_NL ###############
 *
 */
double runAfterOhcNL(TNonLinear* p, double x);
void runAfterOhcNL2(TNonLinear* p, const double *in, double *out, const int length);
void initAfterOhcNL(TNonLinear* p, double taumin, double taumax, double dc, double minR);
TNonLinear* getAfterOhcNL(double taumin, double taumax, double dc, double minR)
{
    TNonLinear* p = (TNonLinear*)hoc_Emalloc(sizeof(TNonLinear));
    initAfterOhcNL(p, taumin, taumax, dc, minR);
    return p;
};

void initAfterOhcNL(TNonLinear* p, double taumin, double taumax, double dc, double minR)
{
    double R;

    p->TauMax = taumax;
    p->TauMin = taumin;

    R = taumin / taumax;
    if (R < minR) p->minR = 0.5 * R;
    else p->minR   = minR;
    p->A = p->minR / (1 - p->minR); /* makes x = 0; output = 1; */
    p->dc = dc;
    R = R - p->minR;
    /*/ This is for new nonlinearity */
    p->s0 = -dc / log(R / (1 - p->minR));

    p->run = runAfterOhcNL;
    p->run2 = runAfterOhcNL2;

    return;
};

double runAfterOhcNL(TNonLinear* p, double x)
{
    /** output of the nonlinearity
        out = TauMax*(minR+(1.0-minR)*exp(-x1/s0));\\
        if the input is zero, the output is TauMax,\\
        if the input is dc, the output is TauMin,\\
        if the input is too large, the output is pushed to TauMax*minR
        if the input is negative, the output is the 2*TauMax-out (keep
        the control signal continuous)
     */
    double out;
    double x1 = fabs(x);
    out = p->TauMax * (p->minR + (1.0 - p->minR) * exp(-x1 / p->s0));
    return(out);
};

void runAfterOhcNL2(TNonLinear* p, const double *in, double *out, const int length)
{
    int register i;
    double x1;

    for (i = 0;i < length;i++) {
        x1 = fabs(in[i]);
        out[i] = p->TauMax * (p->minR + (1.0 - p->minR) * exp(-x1 / p->s0));
    };
    return;
};

/*
 *
 * ############## Implementation InterFace about IHCPPI ###############
 *
 */

/*/IHCPPI Nonlinear */
double runIHCPPI(TNonLinear* p, double x)
{
    double PPI;
    double temp;
    temp = p->p2 * x;
    if (temp > 400)
        PPI = p->p1 * temp;
    else
        PPI = p->p1 * log(1. + exp(temp)); /*/ soft-rectifier */
    return(PPI);
};

void runIHCPPI2(TNonLinear* p, const double *in, double *out, const int length)
{
    int register i;
    double PPI;
    double temp;
    for (i = 0;i < length;i++) {
        temp = p->p2 * in[i];
        if (temp < 400) {
            PPI = p->p1 * log(1. + exp(temp)); /*/ soft-rectifier */
        } else {
            PPI = p->p1 * temp;
        }
        out[i] = PPI;
    };
    return;
};

/*
 *
 * ############## Implementation InterFace about IHC NL ###############
 *
 */

double runIHCNL(TNonLinear * p, double x)
{
    double temp, dtemp, tempA;

    if (x >= 0) {
        tempA = p->A0;
    } else {
        dtemp = pow(-x, p->C);
        tempA = -p->A0 * (dtemp + p->D) / (3 * dtemp + p->D);
    };
    temp = tempA * log(fabs(x) * p->B + 1.0);

    return(temp);
};

void runIHCNL2(TNonLinear* p, const double *in, double *out, const int length)
{
    int register i;
    double temp, dtemp, tempA;
    for (i = 0; i < length; i++) {
        /*/begin Vsp -> Vihc */
        temp = in[i];
        if (temp >= 0) {
            tempA = p->A0;
        } else {
            dtemp = pow(-temp, p->C);
            tempA = -p->A0 * (dtemp + p->D) / (3 * dtemp + p->D);
        };
        out[i] = tempA * log(fabs(temp) * p->B + 1.0);
    };
    return;
};

// To do : 1. change the option that about # of model
//         2. add the option that can generate the synapse output / spike rate
//         3. check the parameter for CAT/HUMAN, and maybe more
//         4. Species can be specified as others as the same in the paper

/* Model numbers as used in ARLO Heinz et al., Fig. 4 */
//  model == 1 bmmodel = FeedForward_NL;
//  model == 2 bmmodel = FeedBack_NL;
//  model == 3 bmmodel = Sharp_Linear;
//  model == 4 bmmodel = Broad_Linear;
//  model == 5 bmmodel = Broad_Linear_High;

//#include "cmpa.h"
//#include "hc.h"
//#include "synapse.h"
//#include "complex.h"
//#include "filters.h"
//#include <stdlib.h>
double Get_tau(int species, double cf, int order, double* taumax, double* taumin, double* taurange);
double erbGM(double);
double cochlea_f2x(int species, double f);
double cochlea_x2f(int species, double x);

void runAN2(TAuditoryNerve *p, const double *in, double *out, const int length);
double runBasilarMembrane(TBasilarMembrane *bm, double x);
void run2BasilarMembrane(TBasilarMembrane *bm, const double *in, double *out, const int length);
/*
 * different species have different parameters !!!!
 * 1. Get_Tau(...)
 * 2. Synapse setup in this function
 * 3. IHC low pass filter cutoff freq
 *
 * */
void initAuditoryNerve(TAuditoryNerve *p, int model, int species, double tdres, double cf, double spont)
{ /* default */
    Tsynapse *syn;

    double Kcf, temp;
    double p1, p2, pst, psl;

    p->run = runAN;
    p->run2 = runAN2;

    syn = &(p->syn);

    p->model = model;
    p->species = species;
    p->cf = cf;
    p->spont = spont;
    p->tdres = tdres;

    p->run = runAN;
    p->run2 = runAN2;
    /*
     * Init the basilar membrane
     * */
    initBasilarMembrane(&(p->bm), model, species, tdres, cf);

    /* Now add support for different species */

    /* This line outputs the appropriate sout rate based on whether or not spikes are being generated (with refractoriness) */
    if (p->ifspike == 1) syn->Ass = 350;
    else syn->Ass = 130; /* borrow this for Mike G. Heinz */
    switch (species) {
    case 1 : /* Cat - low CFs only (Carney '93) */
    case 9 : /* Cat - all CFs "Universal"  (Zhang et al 2001) */
    case 2 : /* Rat Model - Michael Eager */
        /*
           * Init the synapse
           * */
        syn->tdres = tdres;
        syn->cf = cf;
        syn->spont = spont;
        syn->Pimax = 0.6;
        syn->PTS = 1 + 9 * spont / (9 + spont);
        syn->Ar_over_Ast = 6.0;
        syn->tauR = 0.002;
        syn->tauST = 0.060;
        initSynapse(syn);  /*This function calcuate all other parameters as described in the appendix*/

        /*
         * Init the ihcppi rectify function and parameters
         * !!!!Use the result from initSynapse(syn)
         *
         * */
        Kcf = 2.0 + 3.0 * log10(cf / 1000);
        if (Kcf < 1.5) Kcf = 1.5;
        syn->Vsat = syn->Pimax * Kcf * 20.0 * (1 + spont) / (5 + spont);
        /* Important !!! The original code seems use the following equation!!!
           This is not shown in the paper but Prest is very small compare to the other part < 0.01*Vsat
           and the results is the same
         */
        /* syn->Vsat = syn->Pimax*Kcf*20.0*(1+spont)/(5+spont)+syn->Prest; */

        break;
    case 0 :
    default: /* Human */
        /*
         * Init the synapse
         * */
        syn->tdres = tdres;
        syn->cf = cf;
        syn->spont = spont;
        syn->Pimax = 0.6;
        syn->PTS = 1 + 9 * spont / (9 + spont);   /* Peak to Steady State Ratio, characteristic of PSTH */
        syn->Ar_over_Ast = 6.0;
        syn->tauR = 0.002;
        syn->tauST = 0.060;
        initSynapse(syn);  /*This function calcuate all other parameters as described in the appendix*/

        /*
         * Init the ihcppi rectify function and parameters
         * !!!!Use the result from initSynapse(syn)
         *
         * */
        Kcf = 2.0 + 1.3 * log10(cf / 1000); /*/ From M. G. Heinz */
        if (Kcf < 1.5) Kcf = 1.5;
        syn->Vsat = syn->Pimax * Kcf * 60.0 * (1 + spont) / (6 + spont);
        /* Important !!! The original code seems use the following equation!!!
           This is not shown in the paper but Prest is very small compare to the other part < 0.01*Vsat
           and the results is the same, so we use the upper equation
         */
        /* syn->Vsat = syn->Pimax*Kcf*60.0*(1+spont)/(6+spont)+syn->Prest; */

        break;
    };
    temp = log(2) * syn->Vsat / syn->Prest;
    if (temp < 400) pst = log(exp(temp) - 1);
    else pst = temp;
    psl = syn->Prest * pst / log(2);
    p2 = pst;
    p1 = psl / pst;

#ifdef DEBUG
    printf("\ninitSyanpse : p1=%f; p2=%f", p1, p2);
#endif

    p->ihcppi.psl = psl;
    p->ihcppi.pst = pst;
    p->ihcppi.p1 = p1;
    p->ihcppi.p2 = p2;
    p->ihcppi.run = runIHCPPI; /*These are functions used to get ihcppi output */
    p->ihcppi.run2 = runIHCPPI2;  /*Defined in hc.c, hc.h */

    /*
     * Init Inner Hair Cell Model
     */
    /* For Human Only !!!!!!! Hair Cell low pass filter (tdres,cutoff,gain,order) */
    switch (species) {
    case 0:
    default: /* Human */
        initLowPass(&(p->ihc.hclp), tdres, 4500.0, 1.0, 7);
        break;
    case 1:  /* Cat, this version use the Laurel's revcor data to determine taumin, taumax - only good for low freqs*/
    case 9:  /* Also for universal CAT; as in JASA 2001 Zhang et al. - good for all CFs */
    case 2:
        initLowPass(&(p->ihc.hclp), tdres, 3800.0, 1.0, 7);
        break;
    };
    p->ihc.run = runHairCell;
    p->ihc.run2 = runHairCell2;
    p->ihc.hcnl.A0 = 0.1;  /* Inner Hair Cell Nonlinear Function */
    p->ihc.hcnl.B = 2000.;
    p->ihc.hcnl.C = 1.74;
    p->ihc.hcnl.D = 6.87e-9;
    p->ihc.hcnl.run = runIHCNL;
    p->ihc.hcnl.run2 = runIHCNL2;
    return;
}

double runAN(TAuditoryNerve *p, double x)
{
    double bmout, ihcout, ppiout, sout;
    bmout = p->bm.run(&(p->bm), x);
    ihcout = p->ihc.run(&(p->ihc), bmout);
    ppiout = p->ihcppi.run(&(p->ihcppi), ihcout);
    sout = p->syn.run(&(p->syn), ppiout);
    return(sout);
}

void runAN2(TAuditoryNerve *p, const double *in, double *out, const int length)
{
    p->bm.run2(&(p->bm), in, out, length);
    p->ihc.run2(&(p->ihc), out, out, length);
    p->ihcppi.run2(&(p->ihcppi), out, out, length);
    p->syn.run2(&(p->syn), out, out, length);
    return;
}

/* User Get_Tau, initGammaTone, initBoltzman*/
void initBasilarMembrane(TBasilarMembrane* bm, int model, int species, double tdres, double cf)
{ /*
   *
   * ##### Get Basilar Membrane ########
   * 1. Get a structure of BasilarMembrane
   * 2. Specify wide band filter in the BasilarMembrane Model
   *    //WB filter not used for all model versions, but it's computed here anyway
   * 3. Specify the OHC model in the control path
   *    3.2 Specify the NL function used in the OHC
   * 4. Specify the AfterOHC NL in the control path
   * */

    int bmmodel;
    double taumax, taumin, taurange; /* general */
    double x, centerfreq, tauwb, tauwbmin, dtmp, wb_gain; /* for wb */
    double absdb, ohcasym; /* for ohc */
    double dc, R, minR; /* for afterohc */
    TNonLinear *p;

    /* Model numbers as used in ARLO Heinz et al., Fig. 4 */
    if (model == 1) bmmodel = FeedForward_NL;
    else if (model == 2) bmmodel = FeedBack_NL;
    else if (model == 3) bmmodel = Sharp_Linear;
    else if (model == 4) bmmodel = Broad_Linear;
    else if (model == 5) bmmodel = Broad_Linear_High;
    bm->bmmodel = bmmodel;
    bm->tdres = tdres;

    /*
     *  Determine taumax,taumin,order here
     */
    bm->run = runBasilarMembrane;
    bm->run2 = run2BasilarMembrane;

    bm->bmorder = 3;
    Get_tau(species, cf, bm->bmorder, &taumax, &taumin, &taurange);

    bm->TauMax = taumax;
    bm->TauMin = taumin;
    if (bm->bmmodel&Broad_ALL)
        bm->tau = taumin;
    else
        bm->tau = taumax;

    initGammaTone(&(bm->bmfilter), tdres, cf, bm->tau, 1.0, bm->bmorder);
    initGammaTone(&(bm->gfagain), tdres, cf, taumin, 1.0, 1);
    /*
     * Get Wbfilter parameters
     */
    x = cochlea_f2x(species, cf);
    centerfreq = cochlea_x2f(species, x + 1.2); /* shift the center freq Qing use 1.1 shift */
    bm->wborder = 3;
    tauwb = taumin + 0.2 * (taumax - taumin);
    tauwbmin = tauwb / taumax * taumin;
    dtmp = tauwb * TWOPI * (centerfreq - cf);
    wb_gain = pow((1 + dtmp * dtmp), bm->wborder / 2.0);

    bm->TauWB = tauwb;
    bm->TauWBMin = tauwbmin;
    initGammaTone(&(bm->wbfilter), tdres, centerfreq, tauwb, wb_gain, bm->wborder);
    bm->A = (tauwb / taumax - tauwbmin / taumin) / (taumax - taumin);
    bm->B = (taumin * taumin * tauwb - taumax * taumax * tauwbmin) / (taumax * taumin * (taumax - taumin));

    /*
     * Init OHC model
     */
    absdb = 20; /* The value that the BM starts compression */
    initLowPass(&(bm->ohc.hclp), tdres, 800.0, 1.0, 3); /* This is now use in both Human & Cat MODEL */
    /*/ parameter into boltzman is corner,slope,strength,x0,s0,x1,s1,asym
    // The corner determines the level that BM have compressive nonlinearity */
    ohcasym = 7.0;
    /*/set OHC Nonlinearity as boltzman function combination */
    init_boltzman(&(bm->ohc.hcnl), absdb - 12, 0.22, 0.08, 5, 12, 5, 5, ohcasym);
    bm->ohc.run = runHairCell;
    bm->ohc.run2 = runHairCell2;

    /*
     * Init AfterOHC model
     */
    p = &(bm->afterohc);
    dc = (ohcasym - 1) / (ohcasym + 1.0) / 2.0;
    dc -= 0.05;
    minR = 0.05;
    p->TauMax = taumax;
    p->TauMin = taumin;
    R = taumin / taumax;
    if (R < minR) p->minR = 0.5 * R;
    else p->minR   = minR;
    p->A = p->minR / (1 - p->minR); /* makes x = 0; output = 1; */
    p->dc = dc;
    R = R - p->minR;
    /* This is for new nonlinearity */
    p->s0 = -dc / log(R / (1 - p->minR));
    p->run = runAfterOhcNL;
    p->run2 = runAfterOhcNL2;
    return;
}

/** pass the signal through the tuning filter.
    using different model
    Sharp_Linear | Broad_Linear | Broad_Linear_High | FeedBack_NL | FeedForward_NL
*/
double runBasilarMembrane(TBasilarMembrane *bm, double x)
{
    double out;
    const int length = 1;
    run2BasilarMembrane(bm, &x, &out, length);
    return(out);
}

void run2BasilarMembrane(TBasilarMembrane *bm, const double *in, double *out, const int length)
{
    register int i;
    double wb_gain, dtmp, taunow;
    double wbout, ohcout;
    double x, x1, out1;

    for (i = 0; i < length; i++) {
        x = in[i];
        x1 = bm->bmfilter.run(&(bm->bmfilter), x); /*/ pass the signal through the tuning filter */
        switch (bm->bmmodel) {
        default:
        case Sharp_Linear: /* / if linear model, needn't get the control signal */
        case Broad_Linear:
            out1 = x1; break;
        case Broad_Linear_High: /* /adjust the gain of the tuning filter */
            out1 = x1 * pow((bm->TauMin / bm->TauMax), bm->bmfilter.Order);
            break;
        case FeedBack_NL: /* /get the output of the tuning filter as the control signal */
            wbout = x1;

            ohcout = bm->ohc.run(&(bm->ohc), wbout); /*/ pass the control signal through OHC model */
            /*/ pass the control signal through nonliearity after OHC */
            bm->tau = bm->afterohc.run(&(bm->afterohc), ohcout);
            /*/ set the tau of the tuning filter */
            bm->bmfilter.settau(&(bm->bmfilter), bm->tau);
            /*  Gain Control of the tuning filter */
            out1 = pow((bm->tau / bm->TauMax), bm->bmfilter.Order) * x1;
            break;
        case FeedForward_NL:
            /*/get the output of the wide-band pass as the control signal */
            wbout = bm->wbfilter.run(&(bm->wbfilter), x);
            /*/ scale the tau for wide band filter in control path */
            taunow = bm->A * bm->tau * bm->tau - bm->B * bm->tau;

            bm->wbfilter.settau(&(bm->wbfilter), taunow); /*/set the tau of the wide-band filter*/
            /*/ normalize the gain of the wideband pass filter as 0dB at CF */
            dtmp = taunow * TWOPI * (bm->wbfilter.F_shift - bm->bmfilter.F_shift);
            wb_gain = pow((1 + dtmp * dtmp), bm->wbfilter.Order / 2.0);
            bm->wbfilter.gain = wb_gain;

            ohcout = bm->ohc.run(&(bm->ohc), wbout); /*/ pass the control signal through OHC model*/
            /*/ pass the control signal through nonliearity after OHC */
            bm->tau = bm->afterohc.run(&(bm->afterohc), ohcout);
            /*/ set the tau of the tuning filter */
            bm->bmfilter.settau(&(bm->bmfilter), bm->tau);
            /*/ Gain Control of the tuning filter */
            out1 = pow((bm->tau / bm->TauMax), bm->bmfilter.Order) * x1;
            break;
        };
        out[i] = out1;
    };
    bm->gfagain.run2(&(bm->gfagain), out, out, length);

    return;
}


/* ********************************************************************************************
 * * Get TauMax, TauMin for the tuning filter. The TauMax is determined by the bandwidth/Q10
    of the tuning filter at low level. The TauMin is determined by the gain change between high
    and low level
    Also the calculation is different for different species
 */
double Get_tau(int species, double cf, int order, double* _taumax, double* _taumin, double* _taurange)
{
    double taumin, taumax, taurange;
    double ss0, cc0, ss1, cc1;
    double Q10, bw, gain, ratio;
    double xcf, x1000;
    gain = 20 + 42 * log10(cf / 1e3);          /*/ estimate compression gain of the filter */
    if (gain > 70) gain = 70;
    if (gain < 15) gain = 15;
    ratio = pow(10, (-gain / (20.0 * order))); /*/ ratio of TauMin/TauMax according to the gain, order */

    /*/ Calculate the TauMax according to different species */
    switch (species) {
    case 1:
        /* Cat parameters: Tau0 vs. CF (from Carney & Yin 88) - only good for low CFs*/
        /* Note: Tau0 is gammatone time constant for 'high' levels (80 dB rms) */
        /* parameters for cat Tau0 vs. CF */
        xcf = cochlea_f2x(species, cf);    /* position of cf unit; from Liberman's map */
        x1000 = cochlea_f2x(species, 1000); /* position for 1 Khz */
        ss0 = 6.; cc0 = 1.1; ss1 = 2.2; cc1 = 1.1;
        taumin = (cc0 * exp(-xcf / ss0) + cc1 * exp(-xcf / ss1)) * 1e-3;     /* in sec */
        taurange = taumin * xcf / x1000;
        taumax = taumin + taurange;
        break;
    case 0:
        /* Human Parameters: From Mike Heinz: */
        /* Bandwidths now are based on Glasberg and Moore's (1990) ERB=f(CF,level) equations  */
        taumax =  1. / (TWOPI * 1.019 * erbGM(cf));
        break;
    case 9:
    case 2:
        /* Universal species from data fitting : From Xuedong Zhang,Ian (JASA 2001) */
        /* the Q10 determine the taumax(bandwidths at low level) Based on Cat*/
        Q10 = pow(10, 0.4708 * log10(cf / 1e3) + 0.4664);
        bw = cf / Q10;
        taumax = 2.0 / (TWOPI * bw);
    }
    ; /*/end of switch */
    taumin =  taumax * ratio;
    taurange = taumax - taumin;
    *_taumin = taumin;
    *_taumax = taumax;
    *_taurange = taurange;
    return 0;
}

/*/ --------------------------------------------------------------------------------
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
        x  = (8.03 / 0.928) * log10(1.0 + f / 7613.3);
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
    case 0: /*/human
      //      if((x>35)||(x<0)) error("BM distance out of human range, [in cochlea_x2f(...)]"); */
        f = 165.4 * (pow(10, (0.06 * x)) - 0.88);
        break;
    case 2:
        f = 7613.3 * (10 ^(0.924 * x / 8.03) - 1.0);
        break;
    default:
    case 1: /*/cat */
        f = 456.0 * (pow(10, x / 11.9) - 0.80);
        break;
    };
    return(f);
}
/*/ --------------------------------------------------------------------------------
 ** Calculate the erb at 65dB */
double erb51_65(double CF)
{
    double erb1000, erbCf, erb;
    erb1000 = 24.7 * (4.37 * 1000 / 1000 + 1);
    erbCf = 24.7 * (4.37 * CF / 1000 + 1);
    erb = (0.5 * erbCf) / (1 - (0.38 / 4000) * erb1000 * (65 + 10 * log10(erbCf) - 51)) + .5 * erbCf;
    return(erb);
}
/*/ --------------------------------------------------------------------------------
 ** Calculate the erb using GM's method */
double erbGM(double CF)
{
    double erbCf;
    erbCf = (1 / 1.2) * 24.7 * (4.37 * CF / 1000 + 1);
    return(erbCf);
}

/*/ ---------------------------------------------------------------------------
 ** Calculate the delay(basilar membrane, synapse for cat*/
double delay_cat(double cf)
{
    /* DELAY THE WAVEFORM (delay buf1, tauf, ihc for display purposes)  */
    /* Note: Latency vs. CF for click responses is available for Cat only (not human) */
    /* Use original fit for Tl (latency vs. CF in msec) from Carney & Yin '88
       and then correct by .75 cycles to go from PEAK delay to ONSET delay */
    double A0 = 8.13; /* from Carney and Yin '88 */
    double A1 = 6.49;
    double x = cochlea_f2x(1, cf); /*/cat mapping */
    double delay = A0 * exp(-x / A1) * 1e-3 - 1.0 / cf;
    return(delay);
}

//#include <stdio.h>
//#include <stdlib.h>
//#include <math.h>
//#include "synapse.h"
//#include "filters.h"

/*/interface to all the code */
int runSynapse(Tsynapse *pthis, const double *in, double *out, int length)
{
    int error = 0;
    initSynapse(pthis);
    ihc_nl(pthis, in, out, length);
    runsyn_dynamic(pthis, in, out, length);
    return(error);
}
/*/function [V,PPI] = ihcmodel(p,Fs)
//% Just the IHC part of ANmod2 - use this after gammatone filter.
//!!!!! This function should be called after initSynapse(); */
int ihc_nl(Tsynapse *pthis, const double *in, double *out, const int length)
{
    int error = 0;
    static double A0 = 0.1;
    static double B = 2000.;
    static double C = 1.74;
    static double D = 6.87e-9;
    register int i;
    double temp, tempA, dtemp;
    double p1, p2, Vihc, PPI, pst, psl;
    double spont, Kcf, Pimax;

    TLowPass ihclowpass;

    double cf = pthis->cf;
    for (i = 0; i < length; i++) {
        /*/begin Vsp -> Vihc */
        temp = in[i];
        if (temp >= 0) {
            tempA = A0;
        } else {
            dtemp = pow(-temp, C);
            tempA = -A0 * (dtemp + D) / (3 * dtemp + D);
        };
        out[i] = tempA * log(fabs(temp) * B + 1.0);
    };

    /*/ Low Pass Filter Goes here !!!!!!!!!!!!!!!!!!!!! */
    initLowPass(&ihclowpass, pthis->tdres, 4500.0, 1.0, 7);
    ihclowpass.run2(&ihclowpass, out, out, length);

    /*/ begin the Vinc -> PPI */
    Kcf = 2.0 + 1.3 * log10(cf / 1000); /*/ From M. G. Heinz */
    if (Kcf < 1.5) Kcf = 1.5;
    /*/ Important !!! This is not shown in the paper
    // For Prest is so small compare to the other part */
    pthis->Vsat = pthis->Pimax * Kcf * 60.0 * (1 + pthis->spont) / (6 + pthis->spont) + pthis->Prest;

    temp = log(2) * pthis->Vsat / pthis->Prest;
    if (temp < 400) pst = log(exp(temp) - 1);
    else pst = temp;
    psl = pthis->Prest * pst / log(2);
    p2 = pst;
    p1 = psl / pst;
    for (i = 0;i < length;i++) {
        Vihc = out[i];
        PPI = p1 * log(1. + exp(p2 * Vihc)); /*/ soft-rectifier */
        out[i] = PPI;
    };

    return(error);
};

int initSynapse(Tsynapse *pthis)
{
    int error = 0;
    double PTS, Ass, Aon, Ar_over_Ast, Ar, Ast;
    double Pimax, spont;
    double Prest;
    double CG, gamma1, gamma2, tauR, tauST, kappa1, kappa2, VI0, VI1, VI;
    double alpha, beta, theta1, theta2, theta3;
    double PL, PG, VL, Cirest, CLrest;

    double CLlast, CIlast, PPIlast;
    double tdres, CInow, CLnow;
    register int i;

    pthis->run = run1syn_dynamic;
    pthis->run2 = runsyn_dynamic;
    tdres = pthis->tdres;
    /*/begin the Synapse dynamic */
    PTS = pthis->PTS;
    Ass = pthis->Ass; /* For Human, from M. G. Heinz */
    Aon = PTS * Ass;
    Ar_over_Ast = pthis->Ar_over_Ast;
    Ar = (Aon - Ass) * (Ar_over_Ast) / (1. + Ar_over_Ast);  /*/%???? Ar on both sides */
    Ast = Aon - Ass - Ar;
    Pimax = pthis->Pimax;
    spont = pthis->spont; /*/50 in default */
    Prest = Pimax * spont / Aon;
    CG = spont * (Aon - spont) / (Aon * Prest * (1. - spont / Ass));
    gamma1 = CG / spont;
    gamma2 = CG / Ass;
    tauR = pthis->tauR;
    tauST = pthis->tauST;
    kappa1 = -(1. / tauR);
    kappa2 = -(1. / tauST);

    VI0 = (1. - Pimax / Prest) / (gamma1 * (Ar * (kappa1 - kappa2) / CG / Pimax + kappa2 / Prest / gamma1 - kappa2 / Pimax / gamma2));
    VI1 = (1. - Pimax / Prest) / (gamma1 * (Ast * (kappa2 - kappa1) / CG / Pimax + kappa1 / Prest / gamma1 - kappa1 / Pimax / gamma2));
    VI = (VI0 + VI1) / 2.;

    alpha = gamma2 / (kappa1 * kappa2);
    beta = -(kappa1 + kappa2) * alpha;
    theta1 = alpha * Pimax / VI;
    theta2 = VI / Pimax;
    theta3 = gamma2 - 1. / Pimax;
    PL = (((beta - theta2 * theta3) / theta1) - 1.) * Pimax;
    PG = 1. / (theta3 - 1. / PL);
    VL = theta1 * PL * PG;
    Cirest = spont / Prest;
    CLrest = Cirest * (Prest + PL) / PL;
    pthis->Prest = Prest;
    pthis->CIrest = Cirest;
    pthis->CLrest = CLrest;
    pthis->VI = VI;
    pthis->VL = VL;
    pthis->PL = PL;
    pthis->PG = PG;
    pthis->CG = CG;

    pthis->PPIlast = Prest;
    pthis->CLlast = CLrest;
    pthis->CIlast = Cirest;
    return(error);
};

void runsyn_dynamic(Tsynapse *pthis, const double *in, double *out, int length)
{
    /* !!!!!!!!!!
     * The double *in and double *out could be the same pointer,
     * SO be careful about this
     */
    int error;
    int register i;
    double PPIlast, PL, PG, CIlast, CLlast, CG, VI, VL;
    double tdres;
    double CInow, CLnow;

    error = 0;
    tdres = pthis->tdres;
    PL = pthis->PL;
    PG = pthis->PG;
    CG = pthis->CG;
    VI = pthis->VI;
    VL = pthis->VL;
    CIlast = pthis->CIlast;
    CLlast = pthis->CLlast;
    PPIlast = pthis->PPIlast;

    for (i = 0; i < length;i++) {
        CInow = CIlast + (tdres / VI) * ((-PPIlast * CIlast) + PL * (CLlast - CIlast));
        CLnow = CLlast + (tdres / VL) * (-PL * (CLlast - CIlast) + PG * (CG - CLlast));
        PPIlast = in[i];
        CIlast = CInow;
        CLlast = CLnow;
        out[i] = CInow * PPIlast;
    };

    pthis->CIlast = CIlast;
    pthis->CLlast = CLlast;
    pthis->PPIlast = pthis->PPIlast;
    return;
};

double run1syn_dynamic(Tsynapse *pthis, double x)
{
    double out;
    runsyn_dynamic(pthis, &x, &out, 1);
    return out;
};

//#include "runmodel.h"
//#include "cmpa.h"
//#include "complex.h"
//#include "filters.h"

int an_arlo(double tdres, double cf, double spont, int model, int species, int ifspike,
            const double *in, double *out, int length)
{
    TAuditoryNerve anf;
    anf.ifspike = ifspike;
    initAuditoryNerve(&anf, model, species, tdres, cf, spont);
    runAN2(&anf, in, out, length);
    return(0);
};

/*/ ---------------------------------------------------------------- */

int parsecommandline(T_stim *ptm, int argc, char* argv[])
{
    /*int i;
    int needhelp = 0;
    char *para;
    i = 1;
    if(argc ==1 ) needhelp = 1;
    while(i<argc)
    { para = argv[i];
      if((*para)=='-')
        { para++;
    if(strcmp(para,"tdres")==0)
       ptm->tdres = (double)(atof(argv[++i]));
    else if(strcmp(para,"cf")==0)
       ptm->cf = (double)(atof(argv[++i]));
    else if(strcmp(para,"spont")==0)
       ptm->spont = (double)(atof(argv[++i]));
    else if(strcmp(para,"species")==0)
       ptm->species = _atoi(argv[++i]);
    else if(strcmp(para,"model")==0)
       ptm->model = atoi(argv[++i]);
    else if(strcmp(para,"fibers")==0)
       ptm->banks = atoi(argv[++i]);
    else if(strcmp(para,"delx")==0)
     ptm->delx = (double)(atof(argv[++i]));
    else if(strcmp(para,"cfhi")==0)
     ptm->cfhi = (double)(atof(argv[++i]));
    else if(strcmp(para,"cflo")==0)
     ptm->cflo = (double)(atof(argv[++i]));
    else if(strcmp(para,"reptim")==0)
       ptm->reptim = (double)(atof(argv[++i]));
    else if(strcmp(para,"trials")==0)
     ptm->nrep = atoi(argv[++i]);
    else if(strcmp(para,"wavefile")==0)
     {
       ptm->stim = 11;
       strcpy(ptm->wavefile,argv[++i]);
     }
    else if(strcmp(para,"help")==0)
     needhelp = 1;
    else
     { printf("\nUnkown parameters --> %s",para); needhelp = 1; break; };
        }
      else { printf("\nUnkown parameters --> %s",para); needhelp = 1; break; };
      i++;
    };
    if(needhelp==1)
      {
        printf("\n This program accept following parameters:\n");

        printf("\n -species #(0) --> input the species(0=human,1=cat(LF only,JASA '93),9=cat (all CFs, JASA 2001))");
        printf("\n -model #(1)   --> anmodel(1:Nonlinear_w/comp&supp,");
        printf("\n                           2:Nonlinear_w/comp, w/o supp,");
        printf("\n                           3:linear sharp,");
        printf("\n                           4:linear broad, low threshold)");
        printf("\n                           5:linear broad, high threshold");
        printf("\n -cf #(1000)   --> character frequency of the an tested(center cf for filter banks)");
        printf("\n -spont #(50)  --> spontaneous rate of the fier");
        printf("\n -tdres #(2e-6)--> time domain resolution(second)");

        printf("\n\nFor filter banks >>>>>>>>");
        printf("\n -fibers #(1) --># of filters, use this option with cf,cflo,cfhi,delx");
        printf("\n -cfhi #(-1)   --> highest cf to go(specify cfhi,cflo will recalculate cf&delx");
        printf("\n -cflo #(-1)   --> lowest cf to go");
        printf("\n -delx #(0.05mm) --> distance between filters along basilar membrane");
        printf("\n                     !!!!!!!!!");

        printf("\n\nAbout stimulus>>>>>>>>");
        printf("\n -wavefile filename(click) --> specify the stimulus wavefile name(click)");
        printf("\n -reptim #(0.02) --> the time you want to run the model(20msec)");
        printf("\n -trials #(0)    --> spike generation trials");
        printf("\n");
        //exit(1);
      };*/
    return(0);
};


