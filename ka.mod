TITLE klt.mod  The low threshold conductance of cochlear nucleus neurons

COMMENT

NEURON implementation of Jason Rothman's measurements of VCN conductances.

This file implements the transient potassium current found in ventral cochlear
nucleus "Type I" cells, which are largely "stellate" or "multipolar" cells  (Manis and
Marx, 1991; Rothman and Manis, 2003a,b; Manis et al, 1996). The current is likely
 mediated by Kv4.2 potassium channel subunits, but this has not been directly
demonstrated. The specific implementation is described in Rothman and Manis, J.
Neurophysiol. 2003, in the appendix. Measurements were made from isolated
neurons from adult guinea pig, under reasonably stringent voltage clamp conditions.
 The measured current is sensitive to 4-aminopyridine.
Original implementation by Paul B. Manis, April (JHU) and Sept, (UNC)1999.

File split implementaiton, April 1, 2004.

Contact: pmanis@med.unc.edu

Modified by Michael Eager March 31, 2005 to use with CVode

ENDCOMMENT

NEURON {
        SUFFIX ka
        USEION k READ ek WRITE ik
        RANGE gkabar, gka, ik
        GLOBAL ainf, binf, cinf, atau, btau, ctau
}
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
        (nA) = (nanoamp)
}



:INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
        v (mV)
        celsius  (degC)  : model is defined on measurements made at room temp (22C) in Baltimore
        dt (ms)
        ek = -77 (mV)
        gkabar = 0.00477 (mho/cm2) <0,1e9>
   qt    (/ms)
   q10 = 3
   q = 0.19245 (ms)
}

STATE {
        a b c
}

ASSIGNED {
    ik (mA/cm2)
    gka (mho/cm2)
    ainf binf cinf
    atau (ms) btau (ms) ctau (ms)
}

:LOCAL aexp, bexp, cexp

BREAKPOINT {
   SOLVE states METHOD derivimplicit

   gka = gkabar*(a^4)*b*c
       ik = gka*(v - ek)

}



INITIAL {

   qt = q10^((celsius - 22)/10)
   q=1/qt
    rates(v)

    a = ainf
    b = binf
    c = cinf
}

DERIVATIVE states {
   rates(v)      :
       a' = (ainf - a) / atau
       b' = (binf - b) / btau
       c' = (cinf - c) / ctau

}

:PROCEDURE states() {  :Computes state variables m, h, and n
:   trates(v)      :             at the current v and dt.
:   a = a + aexp*(ainf-a)
:   b = b + bexp*(binf-b)
:   c = c + cexp*(cinf-c)
:VERBATIM
:   return 0;
:ENDVERBATIM
:}



PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.

   TABLE ainf, binf, cinf, atau, btau,ctau
   DEPEND dt, celsius FROM -150 TO 150 WITH 300

    ainf = (1 / (1 + exp(-1*(v + 31) / 6)))^0.25
    binf = 1 / (1 + exp((v + 66) / 7))^0.5
    cinf = 1 / (1 + exp((v + 66) / 7))^0.5

    atau =  (q)*((100 / (7*exp((v+60) / 14) + 29*exp(-(v+60) / 24))) + 0.1)
    btau =  (q)*((1000 / (14*exp((v+60) / 27) + 29*exp(-(v+60) / 24))) + 1)
    ctau = (q)*((90 / (1 + exp((-66-v) / 17))) + 10)
}
