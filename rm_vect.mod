TITLE rm_vect.mod Conductances of cochlear nucleus neurons

COMMENT

NEURON implementation of Jason Rothman's measurements of VCN
conductances.

This file implements the average brain sodium current used in the
Rothman model. In the absence of direct measurements in the VCN,
this is a fair assumption. The model differs from the one used in
Rothman et al, (1993) in that the steep voltage dependence of
recovery from inactivation in that model is missing. This may
affect the refractory period. To use the other model, use
najsr.mod instead.

Original implementation by Paul B. Manis, April (JHU) and Sept,
(UNC)1999.

File split implementaiton, April 1, 2004. Contact:
pmanis@med.unc.edu


Modified by Michael Eager March 31, 2005 to use with CVode


ENDCOMMENT

UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
        (nA) = (nanoamp)
}

NEURON {
        SUFFIX rm2
        USEION na READ ena WRITE ina
        USEION k READ ek WRITE ik
    NONSPECIFIC_CURRENT ih
    NONSPECIFIC_CURRENT il

    RANGE gnabar,gkhtbar,ghbar,gleak, erev,eh


}


PARAMETER {
        v (mV)
        celsius =   22      (degC)  : model is defined on measurements made at room temp in Baltimore
        dt (ms)
        ena     =   55      (mV)
        gnabar  =   0.23    (mho/cm2) <0,1e9>
        ek      =  -70      (mV)
        gkhtbar =   0.0159  (mho/cm2) <0,1e9>
        nf      =   0.85                   <0,1> :proportion of n vs p kinetics
        gleak   =   .001        (mho/cm2)
        erev    =   -65  (mV)
        ghbar   =   0.00318 (mho/cm2) <0,1e9>
        eh      = -43 (mV)
        q10=3
        qt  (/ms)
        q = 0.19245 (ms)
}

STATE {
        m h  n p   r
}

ASSIGNED {
    il (mA/cm2)
    ina (mA/cm2)
 ik (mA/cm2)
    ih (mA/cm2)
    gna (mho/cm2)
    gkht (mho/cm2)
    gh (mho/cm2)

    minf hinf pinf ninf rinf
    mexp hexp nexp pexp rexp

}



BREAKPOINT {
    SOLVE states

    ina = gna*(v - ena)
    ik = gkht*(v - ek)
    il = gleak*(v - erev)
    ih = gh*(v - eh)

}


UNITSOFF
INITIAL {
    rates(v)
    m = minf
    h = hinf
    p = pinf
    n = ninf
    r = rinf
}

PROCEDURE states() {  :Computes state variables m, h, and n
        rates(v)      :             at the current v and dt.
        m = m + mexp*(minf-m)
        h = h + hexp*(hinf-h)
        n = n + nexp*(ninf-n)
        p = p + pexp*(pinf-m)
        r = r + rexp*(rinf-h)

        gna = gnabar*m*m*m*h
        gkht = gkhtbar*(nf*(n^2) + (1-nf)*p)
        gh = ghbar*r
        }

PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
    LOCAL  q10, tinc, tau,q
    TABLE minf, hinf, mexp,hexp,pinf,ninf, pexp,nexp,rexp,rinf DEPEND dt,celsius FROM -100 TO 80 WITH 200
    :for celsius=37 where q10 =3, tested at 22 C
    q10 = 3^((celsius - 22)/10)
      tinc = -dt * q10
: average sodium channel
    minf = 1 / (1+exp(-(v + 38) / 7))
    hinf = 1 / (1+exp((v + 65) / 6))
    tau =  ((10 / (5*exp((v+60) / 18) + 36*exp(-(v+60) / 25))) + 0.04)
    mexp   = 1 - exp(tinc/tau)
    tau =  ((100 / (7*exp((v+60) / 11) + 10*exp(-(v+60) / 25))) + 0.6)
    hexp  = 1 - exp(tinc/tau)
:rectifying potassium channel
    ninf =   (1 + exp(-(v + 15) / 5))^-0.5
    pinf =  1 / (1 + exp(-(v + 23) / 6))
    tau = ((100 / (11*exp((v+60) / 24) + 21*exp(-(v+60) / 23))) + 0.7)
    nexp  = 1 - exp(tinc/tau)
    tau = ((100 / (4*exp((v+60) / 32) + 5*exp(-(v+60) / 22))) + 5)
    pexp  = 1 - exp(tinc/tau)
:hyperpolarising mixed-cation current
    rinf = 1 / (1+exp((v + 76) / 7))
    tau = ((100000 / (237*exp((v+60) / 12) + 17*exp(-(v+60) / 14))) + 25)
    rexp  = 1 - exp(tinc/tau)
}
UNITSON
