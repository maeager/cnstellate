TITLE klt.mod  The low threshold conductance of cochlear nucleus neurons

COMMENT

NEURON implementation of Jason Rothman's measurements of VCN conductances.

This file implements the average brain sodium current used in the Rothman model.
In the absence of direct measurements in the VCN, this is a fair assumption.
The model differs from the one used in Rothman et al, (1993) in that the steep
voltage dependence of recovery from inactivation in that model is missing. This
may affect the refractory period. To use the other model, use najsr.mod instead.

Original implementation by Paul B. Manis, April (JHU) and Sept, (UNC)1999.

File split implementaiton, April 1, 2004.
Contact: pmanis@med.unc.edu


Modified by Michael Eager March 31, 2005 to use with CVode


ENDCOMMENT

UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
        (nA) = (nanoamp)
}

NEURON {
    SUFFIX rm
    USEION na READ ena WRITE ina
    USEION k READ ek WRITE ik
    NONSPECIFIC_CURRENT ih
    NONSPECIFIC_CURRENT il
    
    RANGE gleak, erev
    RANGE gnabar, gna, ina,gkhtbar, gkht, ik,ghbar, gh, ih
    GLOBAL hinf, minf, htau, mtau, ninf, pinf, ntau, ptau, rinf, rtau
}


PARAMETER {
    v (mV)
    celsius  (degC)  : model is defined on measurements made at room temp in Baltimore
    dt (ms)
    ena (mV)
    gnabar =  0.07958 (mho/cm2) <0,1e9>
    ek  (mV)
    gkhtbar = 0.01592 (mho/cm2) <0,1e9>
    nf = 0.85 <0,1> :proportion of n vs p kinetics
    gleak = .001   (mho/cm2)
    erev = -65   (mV)
    ghbar = 0.00318 (mho/cm2) <0,1e9>
    eh = -43 (mV)
    q10 = 3
    qt  (/ms)
    q = 0.19245 (ms)
}

STATE {
    m h  n p   r
}

ASSIGNED {
   il (mA/cm2)
    ina (mA/cm2)
    gna (mho/cm2)
    minf hinf
    mtau (ms) 
    htau (ms)
    ik (mA/cm2)
    gkht (mho/cm2)
    pinf ninf
    ptau (ms) 
    ntau (ms)
    gh (mho/cm2)
    ih (mA/cm2)
    rinf
    rtau (ms)

}



BREAKPOINT {
    SOLVE states METHOD cnexp
    gna = gnabar*(m^3)*h
    ina = gna*(v - ena)
    gkht = gkhtbar*(nf*(n^2) + (1-nf)*p)
    ik = gkht*(v - ek)
    il = gleak*(v - erev)
    gh = ghbar*r
    ih = gh*(v - eh)

}



INITIAL {
    qt = q10^((celsius - 22)/10) : if you don't like room temp, it can be changed!
    q=1/qt
    rates(v)
    m = minf
    h = hinf
    p = pinf
    n = ninf
    r = rinf
}



DERIVATIVE states {
    rates(v)
    m' = (minf - m) / mtau
    h' = (hinf - h) / htau
    n' = (ninf - n) / ntau
    p' = (pinf - p) / ptau
    r' = (rinf - r) / rtau
}



PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.

   TABLE minf, hinf, mtau,htau,pinf,ninf, ptau,ntau,rtau,rinf
   DEPEND dt FROM -150 TO 80 WITH 300
   :for celsius=37 where q10 =3, tested at 22 C

: average sodium channel
    minf = 1 / (1+exp(-(v + 38) / 7))
    hinf = 1 / (1+exp((v + 65) / 6))
    mtau =  (q)*((10 / (5*exp((v+60) / 18) + 36*exp(-(v+60) / 25))) + 0.04)
    htau =  (q)*((100 / (7*exp((v+60) / 11) + 10*exp(-(v+60) / 25))) + 0.6)
:rectifying potassium channel
    ninf =   (1 + exp(-(v + 15) / 5))^-0.5
    pinf =  1 / (1 + exp(-(v + 23) / 6))
    ntau = (q)*((100 / (11*exp((v+60) / 24) + 21*exp(-(v+60) / 23))) + 0.7)
    ptau = (q)*((100 / (4*exp((v+60) / 32) + 5*exp(-(v+60) / 22))) + 5)
:hyperpolarising mixed-cation current
    rinf = 1 / (1+exp((v + 76) / 7))
    rtau = q*((100000 / (237*exp((v+60) / 12) + 17*exp(-(v+60) / 14))) + 25)
}
