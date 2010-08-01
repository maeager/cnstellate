/* anmod3m.c */
#define usemiddleear 1         /* 1-> include middle ear model; 0-> don't include middle ear */



typedef struct Tanstruct Tanmodel;

struct Tanstruct {

//General information about the model
    double tdres, cf, spont, fp1;
    int species, model;
    int nstim;
//  runtime usage
    int    ifspike;
    int    control_type;
// 0 ->without suppression
// 1 ->with suppression
    double PI;
//      double *soundin;
// double *meout;
// double *soundout;
// double *control_signal;
// double *ihc_out;
// double *sout;

    // locations of poles
    double ta;
//Pa in Fig. 3 of the paper
    double tb;
//Pb in Fig. 3 of the paper
    double rgain;
// location of the pole closest to imaginary axis
    double nlgain;
// gain for the control signal
    double zero_r;
// Location of zeros
    int delayn;
// forced delay for AN model

};


typedef struct {
    double realpart;
    double imgpart;
}  mycomplex;

mycomplex complexplus(mycomplex complex01, mycomplex complex02)
{
    mycomplex complex_result;
    complex_result.realpart = complex01.realpart + complex02.realpart;
    complex_result.imgpart = complex01.imgpart + complex02.imgpart;
    return (complex_result);
};

mycomplex complexmulti(mycomplex complex01, mycomplex complex02)
{
    mycomplex complex_result;
    complex_result.realpart = complex01.realpart * complex02.realpart
                              - complex01.imgpart * complex02.imgpart;
    complex_result.imgpart = complex01.realpart * complex02.imgpart
                             + complex01.imgpart * complex02.realpart;
    return (complex_result);
};

mycomplex myconj(mycomplex complex01)
{
    mycomplex complex_result;
    complex_result.realpart = complex01.realpart;
    complex_result.imgpart = -complex01.imgpart;
    return (complex_result);
};

int middleear(Tanmodel *Tan, double *stim, double *meout)
{
    int error_number = 1;

    int pole_order = 4;         // two pairs of poles
    int half_pole_order = 2;
    int zero_order = 2;         // two zeros

    double gain_norm;

    double fs_bilinear; //bilinear transformation frequency
    double kkd;
    double bpinput[4][4];
    double bpoutput[4][4];
    double zerobp;

    double preal;
    double pimg;

    double a1;
    double a2;
    double b1[3];
    double b2[3];
    double dy;

    mycomplex p_11, p_12, p_21, p_22; // poles in s plain
    mycomplex p[5];

    int n, i, j;

    //%========== locations of poles ===================%

    p_11.realpart = -250 * 2 * Tan->PI;
    p_11.imgpart = 400 * 2 * Tan->PI;

    p_21.realpart = -2000 * 2 * Tan->PI;
    p_21.imgpart = 6000 * 2 * Tan->PI;

    p_12 = myconj(p_11);
    p_22 = myconj(p_21);

    p[1] = p_11;
    p[2] = p_12;
    p[3] = p_21;
    p[4] = p_22;

    fs_bilinear = 2.0 / Tan->tdres;

    kkd = 1.0;
    for (i = 1; i <= half_pole_order; i = i + 1) {
        kkd = kkd / ((fs_bilinear - p[i*2].realpart) * (fs_bilinear - p[i*2].realpart)
                     + p[i*2].imgpart * p[i*2].imgpart);
    }

    /* ========= setup zeros =================*/
    zerobp = -200;
    for (i = 1; i <= zero_order; i = i + 1) {
        kkd = kkd * (fs_bilinear - zerobp);
    }

    /*========= initialize some variables======*/
    for (i = 1; i <= 3; i = i + 1) {
        for (j = 1; j <= 3; j = j + 1) {
            bpinput[i][j] = 0.0;
            bpoutput[i][j] = 0.0;
        }
    }

    /*======= normalize gain at 1000Hz========== */
    gain_norm = 1.0;
    for (i = 1; i <= pole_order; i = i + 1) {
        gain_norm = gain_norm * sqrt(p[i].realpart * p[i].realpart
                                     + (p[i].imgpart - 1000 * 2 * Tan->PI) * (p[i].imgpart - 1000 * 2 * Tan->PI));
    }
    for (i = 1; i <= zero_order; i = i + 1) {
        gain_norm = gain_norm / sqrt(2 * Tan->PI * zerobp * 2 * Tan->PI * zerobp + 1000.0 * 2 * Tan->PI * 1000.0 * 2 * Tan->PI);
    }

    /* filter coefficients */

    a1 = (1 - (fs_bilinear + zerobp) / (fs_bilinear - zerobp));
    a2 = - (fs_bilinear + zerobp) / (fs_bilinear - zerobp);
    for (i = 1; i <= half_pole_order; i = i + 1) {
        b1[i] = 2 * (fs_bilinear * fs_bilinear
                     - p[i*2].realpart * p[i*2].realpart - p[i*2].imgpart * p[i*2].imgpart)
                / ((fs_bilinear - p[i*2].realpart) * (fs_bilinear - p[i*2].realpart)
                   + p[i*2].imgpart * p[i*2].imgpart);
        b2[i] = -((fs_bilinear + p[i*2].realpart) * (fs_bilinear + p[i*2].realpart)
                  + p[i*2].imgpart * p[i*2].imgpart)
                / ((fs_bilinear - p[i*2].realpart) * (fs_bilinear - p[i*2].realpart)
                   + p[i*2].imgpart * p[i*2].imgpart);
    }


    for (n = 1; n <= Tan->nstim ; n = n + 1) {
        bpinput[1][3] = bpinput[1][2];
        bpinput[1][2] = bpinput[1][1];
        bpinput[1][1] = stim[n-1];

        for (i = 1; i <= half_pole_order; i = i + 1) {
            dy = bpinput[i][1] + a1 * bpinput[i][2] + a2 * bpinput[i][3];
            dy = dy + b1[i] * bpoutput[i][1] + b2[i] * bpoutput[i][2];

            bpinput[i+1][3] = bpoutput[i][2];
            bpinput[i+1][2] = bpoutput[i][1];
            bpinput[i+1][1] = dy;

            bpoutput[i][2] = bpoutput[i][1];
            bpoutput[i][1] = dy;
        }

        meout[n-1] = dy * kkd * gain_norm;

    }
    return (error_number);
}

#define z_order  1

void controlpath(Tanmodel *T, double *meout, double *control_signal)
{


    int i, n;

    double fs_bilinear;  /* bilinear frequency =2.0/tdres */

    double x_cf;           /* location of CF on frequency map   */
    double f_shift;        /* frequency shift corresponding to 1.2mm distance along the BM */
    double wbw;            /* bandwidth of the wide-band filter */
    double gain_norm_bp;   /* normalization factor for Bandpass filter   */
    double temp01_bp, temp02_bp;  /* temporary variable for normalization*/
    mycomplex p[7];      /* poles in control space   */
    mycomplex pd[7];     /* poles in discrete domain */
    double z[7];         /* zeros in control space   */
    double zd[7];        /* zeros in discrete domain */
    double kkd;          /* a variable used in gain control */
    double pda[6], pdb[6];
    double goutput[10][4];
    double ginput[10][4];
    double dy;           /* output of digital filter */
    double tempdouble01; /* temporary variables      */
    double tempdouble02;
    double wbw2pi;       /* this is equal to 2*PI*wbw */

    double preal;
    double pimg;

    double control_nl = 0.0;

    double p_corner;     /* parameters for the first nonlinearity */
    double p_strength;
    double p_slope;
    double splx;
    double xx;           /* output of the first nonlinearity */
    double acp;          /* parameters for Zhang's(2001) first NL */
    double bcp;
    double ccp;


    double potential;    /* parameters for the second nonlinearity */
    double s0, s1, x0, x1;
    double shift;
    double asym;

    double conlp[5][2];
    double bw_conlp;
    double bw_conlp_2pi;

    /* parameters for the feedback LP filter */
    double fblp[2][2];   /* 1st order LP */
    double bw_fblp;      /* bandwidth of the feedback LP */
    double bw_fblp_2pi;  /* bandwidth multiplied by 2pi  */


    /* bilinear transformation frequency */
    fs_bilinear = 2.0 / T->tdres;

    //%=========================================================
    //% band-pass filter
    //%=========================================================

    //frequency map
    x_cf = 11.9 * log10(0.8 + T->cf / 456);
    //%the peak of the suppression curve is shifted by 1.2mm for all CF
    f_shift = (pow(10, ((x_cf + 1.2) / 11.9)) - 0.8) * 456 - T->cf;

    // wbw controls the bandwidth of the wide-bandpass filter
    wbw = T->cf / 4.0;

    /* initial locations of poles and zeros */
    p[1].realpart = -2 * T->PI * wbw;
    p[1].imgpart = 2 * T->PI * (T->cf + f_shift);
    p[2].realpart = -2 * T->PI * wbw;
    p[2].imgpart = -2 * T->PI * (T->cf + f_shift);

    p[3] = p[1];
    p[4] = p[2];

    p[5] = p[1];
    p[6] = p[2];

    z[1] = 0;       /* only one zero is used */
    z[2] = 0;
    z[3] = 0;
    z[4] = 0;
    z[5] = 0;

    for (i = 1; i <= 8; i = i + 1) {
        ginput[i][1] = 0.0;
        ginput[i][2] = 0.0;
        ginput[i][3] = 0.0;
        goutput[i][1] = 0.0;
        goutput[i][2] = 0.0;
        goutput[i][3] = 0.0;
    }

    /* setup for the first NL */
    acp = 100;
    bcp = 2.5;
    ccp = 0.60;

    /* setup for the control LP */
    bw_conlp = 800;
    bw_conlp_2pi = bw_conlp * 2 * T->PI;

    conlp[0][0] = 0.0;
    conlp[0][1] = 0.0;
    conlp[1][0] = 0.0;
    conlp[1][1] = 0.0;
    conlp[2][0] = 0.0;
    conlp[2][1] = 0.0;
    conlp[3][0] = 0.0;
    conlp[3][1] = 0.0;

    /* setup for the feedback LP */
    bw_fblp = 500;
    bw_fblp_2pi = bw_fblp * 2 * T->PI;
    fblp[0][0] = 0.0;
    fblp[0][1] = 0.0;
    fblp[1][0] = 0.0;
    fblp[1][1] = 0.0;

    for (n = 1; n <= T->nstim; n = n + 1) {
        /*=========================================================*/
        /* band-pass filter                                        */
        /*=========================================================*/

        ginput[1][3] = ginput[1][2];
        ginput[1][2] = ginput[1][1];
        ginput[1][1] = meout[n-1];     // meout is X0 in paper Fig. 1.


        /* this part needs to be modified to be a loop,         */
        /* if you want different poles at different locations. */
        /* wbw and wbw2pi are always positive                  */
        wbw2pi = -(p[1].realpart - control_nl);

        /* normalize the gain at cf */
        wbw = wbw2pi / 2.0 / T->PI;
        temp01_bp = sqrt(wbw * wbw + f_shift * f_shift);
        temp02_bp = sqrt((2 * T->cf + f_shift) * (2 * T->cf + f_shift) + (wbw * wbw));
        gain_norm_bp = 2.0 * T->PI * temp01_bp * 2.0 * T->PI * (temp02_bp);
        gain_norm_bp = gain_norm_bp * gain_norm_bp * gain_norm_bp;

        /* normalization factor related to zero(s)     */

        gain_norm_bp = gain_norm_bp /
                       pow(sqrt(2 * T->PI * z[1] * 2 * T->PI * z[1] + 2 * T->PI * T->cf * 2 * T->PI * T->cf), z_order);

        kkd = 1.0;

        for (i = 1; i <= 3; i = i + 1) { /* 3*2 poles */
            preal = p[i*2].realpart - control_nl;
            pimg = p[i*2].imgpart;

            tempdouble01 = (fs_bilinear - preal);
            tempdouble01 = tempdouble01 * tempdouble01 + pimg * pimg;

            dy = (ginput[i][1] + 2 * ginput[i][2] + ginput[i][3]);
            dy = dy + 2 * goutput[i][1] * (fs_bilinear * fs_bilinear - preal * preal - pimg * pimg);

            dy = dy - goutput[i][2] * ((fs_bilinear + preal) * (fs_bilinear + preal) + pimg * pimg);

            dy = dy / tempdouble01;

            ginput[i+1][3] = goutput[i][2];
            ginput[i+1][2] = goutput[i][1];
            ginput[i+1][1] = dy;

            goutput[i][2] = goutput[i][1];
            goutput[i][1] = dy;

        }


        for (i = 4; i <= 3 + z_order; i = i + 1) { /* zeros */
            goutput[i][1] = ginput[i][1] * (fs_bilinear - z[i-3])
                            - ginput[i][2] * (fs_bilinear + z[i-3])
                            - goutput[i][1];
            ginput[i+1][2] = ginput[i+1][1];
            ginput[i+1][1] = goutput[i][1];

        }

        dy = goutput[3+z_order][1];

        dy = dy * kkd * gain_norm_bp;   /* don't forget the gain term */

        /* this dy is X1 in paper fig. 1 */

        /*=========================================================*/
        /*  the first nonlinearity                                 */
        /*=========================================================*/

        /* the first NL of Zhang(2001) */
        if (dy >= 0) {
            xx = bcp * log(1.0 + acp * pow(dy, ccp));
        } else {
            xx = -bcp * log(1.0 + acp * pow(-dy, ccp));
        }

        /* this xx is X2 in paper fig. 1 */

        /*==========================================================*/
        /*  the second nonlinearity                                 */
        /*==========================================================*/

        asym = 7.0;
        s0 = 8.0;
        x1 = 5.0;
        s1 = 3.0;
        shift = 1.0 / (1.0 + asym);

        x0 = s0 * log((1.0 / shift - 1) / (1 + exp(x1 / s1)));

        potential = 1.0 / (1.0 + exp(-(xx - x0) / s0) * (1.0 + exp(-(xx - x1) / s1))) - shift;

        /* this potential is Y in paper figure 1 */

        potential = potential * T->nlgain;

        //%==========================================================
        //%   low pass again w/ cut-off freq = 800 Hz
        //%==========================================================

        conlp[0][0] = conlp[0][1];
        conlp[0][1] = potential;

        for (i = 1; i <= 3; i = i + 1) {
            conlp[i][1] = (conlp[i-1][1] + conlp[i-1][0]
                           + conlp[i][0] * (fs_bilinear - bw_conlp_2pi))
                          / (fs_bilinear + bw_conlp_2pi);
        }
        for (i = 1; i <= 3; i = i + 1) {
            conlp[i][0] = conlp[i][1];
        }
        control_signal[n-1] = conlp[3][1] * 800 * 2 * T->PI * 800 * 2 * T->PI * 800 * 2 * T->PI * 1.5;

        /*===================================================================*/
        /* feedback low-pass filter                                          */
        /* parameters for the feedback LP filter                             */

        fblp[0][0] = fblp[0][1];
        fblp[0][1] = control_signal[n-1];

        for (i = 1; i <= 1; i = i + 1) {
            fblp[i][1] = (fblp[i-1][1] + fblp[i-1][0]
                          + fblp[i][0] * (fs_bilinear - bw_fblp_2pi))
                         / (fs_bilinear + bw_fblp_2pi);
        }
        for (i = 1; i <= 1; i = i + 1) {
            fblp[i][0] = fblp[i][1];
        }
        control_nl = fblp[1][1] * 500 * 2 * T->PI * 10;

    }    // =================================================end of time loop

}


void signalpath(Tanmodel *Tan, double* control_signal, double* meout, double* soundout)
{
    int n, i;

    int order_of_pole;
    int half_order_pole;
    int order_of_zero;

    double aa;
    double zeroa; /* location of the zeros */
    double fs_bilinear; /* bilinear transformation frequency */


    double bfp; /* best frequency location on imaginary axis */

    mycomplex poles[21];   /* can use up to 10 pairs of poles */

    double gain_norm;

    double dy;

    double bpinput[22][4];
    double bpoutput[22][4];
    double preal;
    double pimg;

    double tempdouble01;


    fs_bilinear = 2.0 / Tan->tdres;

    /*================ setup the locations of poles and zeros =======*/

    order_of_pole = 20;           /* 5 pairs of poles */
    half_order_pole = 10;
    order_of_zero = 10;

    poles[1].realpart = -(Tan->rgain);
    poles[1].imgpart = Tan->fp1 * 2 * Tan->PI;

    poles[5].realpart = poles[1].realpart - Tan->ta;
    poles[5].imgpart = poles[1].imgpart - Tan->tb;

    poles[3].realpart = (poles[1].realpart + poles[5].realpart) * 0.5;
    poles[3].imgpart = (poles[1].imgpart +  poles[5].imgpart) * 0.5;

    poles[2] = myconj(poles[1]);
    poles[4] = myconj(poles[3]);
    poles[6] = myconj(poles[5]);

    poles[7] = poles[1];
    poles[8] = poles[2];
    poles[9] = poles[5];
    poles[10] = poles[6];

    for (i = 1; i <= 10; i = i + 1) {
        poles[i+10] = poles[i];
    }

    zeroa = -(Tan->zero_r);

    /*===================== normalize the gain =====================*/

    bfp = 2 * Tan->PI * Tan->cf;

    gain_norm = 1.0;

    for (n = 1; n <= order_of_pole; n = n + 1) {
        tempdouble01 = bfp - poles[n].imgpart;
        gain_norm = gain_norm *
                    (tempdouble01 * tempdouble01 + poles[n].realpart * poles[n].realpart);
    }


    gain_norm = sqrt(gain_norm);

    gain_norm = gain_norm / pow((sqrt(bfp * bfp + zeroa * zeroa)), order_of_zero);


    for (i = 1; i <= half_order_pole + order_of_zero + 1; i = i + 1) {
        bpinput[i][1] = 0.0;
        bpinput[i][2] = 0.0;
        bpinput[i][3] = 0.0;
        bpoutput[i][1] = 0.0;
        bpoutput[i][2] = 0.0;
        bpoutput[i][3] = 0.0;
    }

    /*%==================================================  */
    /*%      time loop begins here                         */
    /*%==================================================  */
    for (n = 1; n <= Tan->nstim ; n = n + 1) {
        aa = -Tan->rgain - control_signal[n-1];
        if (aa >= 0) {
        stayhere:  aa = 100;       // you can choose to print an error msg here
            goto stayhere;
        }

        //%========== update pole locations ===================%
        poles[1].realpart = aa;
        poles[1].imgpart = Tan->fp1 * 2 * Tan->PI;

        poles[5].realpart = poles[1].realpart - Tan->ta;
        poles[5].imgpart = poles[1].imgpart - Tan->tb;

        poles[3].realpart = (poles[1].realpart + poles[5].realpart) * 0.5;
        poles[3].imgpart = (poles[1].imgpart +  poles[5].imgpart) * 0.5;

        poles[2] = myconj(poles[1]);
        poles[4] = myconj(poles[3]);
        poles[6] = myconj(poles[5]);

        poles[7] = poles[1];
        poles[8] = poles[2];
        poles[9] = poles[5];
        poles[10] = poles[6];

        for (i = 1; i <= 10; i = i + 1) {
            poles[i+10] = poles[i];
        }

        bpinput[1][3] = bpinput[1][2];
        bpinput[1][2] = bpinput[1][1];
        bpinput[1][1] = meout[n-1];

        for (i = 1; i <= half_order_pole; i = i + 1) { // 10*2 poles
            preal = poles[i*2].realpart;
            pimg = poles[i*2].imgpart;

            tempdouble01 = (fs_bilinear - preal);
            tempdouble01 = tempdouble01 * tempdouble01 + pimg * pimg;

            dy = (bpinput[i][1] +
                  (1 - (fs_bilinear + zeroa) / (fs_bilinear - zeroa)) * bpinput[i][2]
                  - (fs_bilinear + zeroa) / (fs_bilinear - zeroa) * bpinput[i][3]);

            dy = dy + 2 * bpoutput[i][1] * (fs_bilinear * fs_bilinear - preal * preal - pimg * pimg);

            dy = dy - bpoutput[i][2] * ((fs_bilinear + preal) * (fs_bilinear + preal) + pimg * pimg);

            dy = dy / tempdouble01;

            bpinput[i+1][3] = bpoutput[i][2];
            bpinput[i+1][2] = bpoutput[i][1];
            bpinput[i+1][1] = dy;

            bpoutput[i][2] = bpoutput[i][1];
            bpoutput[i][1] = dy;

        }

        dy = bpoutput[half_order_pole][1] * pow((fs_bilinear - zeroa), 10);

        dy = dy * gain_norm;     /* don't forget the gain term */
        soundout[n-1] = dy;

        /*each loop below is for a pair of poles and one zero */
        /*and for each loop, a zero at -1 (-infinite in control space) is added*/
        /* so that the order of zeros is same as the order of poles */

        soundout[n-1] = soundout[n-1] / 3;
        /* signal path output is divided by 3 to increase the threshold at CF */

    } /*end of time loop */

}


void ihczxd2001(Tanmodel *Tan, double* soundout, double* ihc_out, double* sout)
{


    //extern double sout_average;

    //double sout_sum;   //used for the calculation of sout_average;
    //double tscale;     //used for the calculation of sout_average;

    int i, j;           // loop counters
    double tempdouble;
    double tempA;


    int ihc_lp_order;   // order of the low-pass filter in IHC
    double ihc_lp_gain; // gain of the low-pass filter in IHC
    double ihc_lp_Fc;   // cut-off freq of the low-pass filter in IHC
    double ihc_lp[15];  // tmp for ihc lowpass
    double ihc_c1LP;
    double ihc_c2LP;
    double ihc_lp_c;
    double hcl[15];

    double A0 = 0.1;
    double B = 2000.;
    double C = 1.74;
    double D = 6.87e-9;
    double dtemp;

    double Ass;
    double PTS = 8.627;
    double Aon;
    double Ar_over_Ast = 6.;
    double Ar;
    double Ast;
    double Pimax = 0.6;
    double spont = 50.0;
    double Kcf;
    double Vsat;
    double pst;
    double psl;
    double p1;
    double p2;

    double Prest;
    double CG;
    double gamma1;
    double gamma2;
    double kappa1;
    double kappa2;

    double tauR = 0.002; //human
    double tauST = 0.060;

    double VI0;
    double VI1;
    double VI;

    double alpha;
    double beta;
    double theta1;
    double theta2;
    double theta3;

    double PL;
    double PG;
    double VL;

    double Cirest;
    double CLrest;

    double Vihc;
    double PPI;
    double PPIlast;

    double CInow;
    double CLnow;
    double CIlast;
    double CLlast;

    FILE *psave02;

    if (Tan->ifspike == 1) Ass = 350;
    else Ass = 130; /* read Frank's cmpa.c */

    for (i = 1; i <= Tan->nstim; i++) {
        /*/begin Vsp -> Vihc */
        tempdouble = soundout[i-1];
        if (tempdouble >= 0) {
            tempA = A0;
        } else {
            dtemp = pow(-tempdouble, C);
            tempA = -A0 * (dtemp + D) / (3 * dtemp + D);
        };
        ihc_out[i-1] = tempA * log(fabs(tempdouble) * B + 1.0);

    };

    //-----------------------------------------------------------------------
    //      low pass filter
    //-----------------------------------------------------------------------
    // initialize the parameters of the low-pass
    // from zxd's synapse.c and filters.c
    ihc_lp_order = 7;
    ihc_lp_gain = 1.0;
    ihc_lp_Fc = 3800.0;       //this would be 3800Hz if for cat.  4500 for human

    ihc_lp_c = 2.0 / Tan->tdres;
    ihc_c1LP = (ihc_lp_c - 2 * Tan->PI * ihc_lp_Fc) / (ihc_lp_c + 2 * Tan->PI * ihc_lp_Fc);
    ihc_c2LP = 2 * Tan->PI * ihc_lp_Fc / (2 * Tan->PI * ihc_lp_Fc + ihc_lp_c);

    for (j = 1; j <= ihc_lp_order + 1; j = j + 1)
        hcl[j] = 0.0;

    for (i = 1; i <= Tan->nstim; i = i + 1) {
        /*hc[0] = in[loopSig]*gain;

        for(loopLP=0; loopLP<pOrder; loopLP++)
        hc[loopLP+1] = c1LP * hcl[loopLP+1] + c2LP*(hc[loopLP]+hcl[loopLP]);

        for(loopLP=0; loopLP<=pOrder;loopLP++) hcl[loopLP] = hc[loopLP];
        out[loopSig] = hc[pOrder];   */

        // the above part is copied from zxd
        ihc_lp[1] = ihc_out[i-1] * ihc_lp_gain;

        for (j = 1; j <= ihc_lp_order; j = j + 1)
            ihc_lp[j+1] = ihc_c1LP * hcl[j+1] + ihc_c2LP * (ihc_lp[j] + hcl[j]);

        for (j = 1; j <= ihc_lp_order + 1; j = j + 1)
            hcl[j] = ihc_lp[j];

        ihc_out[i-1] = ihc_lp[ihc_lp_order+1];

    };


    //=======================================================================

    /*/begin the Synapse dynamic */
    Aon = PTS * Ass;
    Ar = (Aon - Ass) * (Ar_over_Ast) / (1. + Ar_over_Ast);  /*/%???? Ar on both sides */
    Ast = Aon - Ass - Ar;
    Prest = Pimax * spont / Aon;
    CG = spont * (Aon - spont) / (Aon * Prest * (1. - spont / Ass));
    gamma1 = CG / spont;
    gamma2 = CG / Ass;

    kappa1 = -(1. / tauR);
    kappa2 = -(1. / tauST);

    VI0 = (1. - Pimax / Prest) /
          (gamma1 * (Ar * (kappa1 - kappa2) / CG / Pimax + kappa2 / Prest / gamma1 - kappa2 / Pimax / gamma2));
    VI1 = (1. - Pimax / Prest) /
          (gamma1 * (Ast * (kappa2 - kappa1) / CG / Pimax + kappa1 / Prest / gamma1 - kappa1 / Pimax / gamma2));
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

    PPIlast = Prest;
    CLlast = CLrest;
    CIlast = Cirest;


    /* from Frank's cmpa.c for species==9 universal */
    Kcf = 2.0 + 3.0 * log10(Tan->cf / 1000);
    if (Kcf < 1.5) Kcf = 1.5;
    Vsat = Pimax * Kcf * 20.0 * (1 + spont) / (5 + spont);

    tempdouble = log(2) * Vsat / Prest;
    //if(tempdouble<400)
    if (tempdouble < 100) {    //zxd used 400, I think 100 is enough and safer.
        pst = log(exp(tempdouble) - 1);
    } else {
        pst = tempdouble;
    }

    psl = Prest * pst / log(2);
    p2 = pst;
    p1 = psl / pst;

    for (i = 1; i <= Tan->nstim; i++) {
        Vihc = ihc_out[i-1];
        tempdouble = p2 * Vihc;
        if (tempdouble < 100) {
            PPI = p1 * log(1. + exp(p2 * Vihc)); /*/ soft-rectifier */
        } else {
            PPI = p1 * tempdouble;
        }

        ihc_out[i-1] = PPI;

    };


    for (i = 1; i <= Tan->nstim; i++) {
        CInow = CIlast + (Tan->tdres / VI) * ((-PPIlast * CIlast) + PL * (CLlast - CIlast));
        CLnow = CLlast + (Tan->tdres / VL) * (-PL * (CLlast - CIlast) + PG * (CG - CLlast));
        PPIlast = ihc_out[i-1];
        CIlast = CInow;
        CLlast = CLnow;
        sout[i-1] = CInow * PPIlast;
    };


};

int initialise(Tanmodel *T,  double  cf, double tdres, double spont, double NLsuppression, double species, double ifspike, double length)
{
    int error_number = 1;

    double rgain80;
    double average_control = 0.3357;
    rgain80 =  pow(10, log10(cf) * 0.5732 + 1.5220);

    T->PI = 3.14159265358979;
    T->fp1 = 1.0854 * cf - 106.0034;
    T->ta =     pow(10, log10(cf) * 1.0230 + 0.1607);
    T->tb =     pow(10, log10(cf) * 1.4292 - 1.1550) - 1000;
    T->rgain =  pow(10, log10(cf) * 0.4 + 1.9);   /* 75% */
    T->nlgain = (rgain80 - T->rgain) / average_control;
    T->zero_r = pow(10, log10(cf) * 1.5 - 0.9);

    T->cf = cf;
    T->tdres = 2e-5;                /* sampling time (20ns, by default) */
    T->ifspike = ifspike;
    T->nstim = length;

    T->model = NLsuppression;

    T->control_type = 1;    /* 1 ->with suppression    */

    T->species = species;
//T->sout = sout;
//T->soundin = stim;
//T->meout = makevector(length);
//T->soundout= makevector(length);
//T->control_signal= makevector(length);
//T->ihc_out= makevector(length);
//*********************************
    T->delayn = 25; /* 0.5ms delay for all CFs */

    return(error_number);

}






double an_tan(double tdres, double cf, double spont, double NLsuppression, double species, double ifspike, double *stim, double *sout, int length)
{
    int error;
    int n;

    error = 0;
    Tanmodel Tan;
    double *meout;
    double *soundout;
    double *control_signal;
    double *ihc_out;

    meout = makevector(length);
    soundout = makevector(length);
    control_signal = makevector(length);
    ihc_out = makevector(length);


    /* parameter setup   */
    initialise(&Tan,  cf, tdres, spont, NLsuppression, species, ifspike, length);

    /*Middle Ear Model*/
    middleear(&Tan, stim, meout);
    controlpath(&Tan, meout, control_signal);
    signalpath(&Tan , control_signal, meout, soundout);
    ihczxd2001(&Tan, soundout, ihc_out, sout);

    printf("nlgain = %f\n", Tan.nlgain); /* This gain, due to the nonlinear control path, varies with CF */
    freevector(meout);
    freevector(soundout);
    freevector(control_signal);
    freevector(ihc_out);
    return error;
}


