TITLE klt.mod  The low threshold conductance of cochlear nucleus neurons

COMMENT

NEURON implementation of Jason Rothman's measurements of VCN conductances.

This file implements the low threshold potassium current found in several brainstem
 nuclei of the auditory system, including the spherical and globular bushy cells
  (Manis and Marx, 1991; Rothman and Manis, 2003a,b) and octopus cells (Bal and
  Oertel, 2000) of the ventral cochlear nucleus, principal cells of the medial
  nucleus of the trapzoid body (Brew and Forsythe, 1995, Wang and Kaczmarek,
  1997) and neurons of the medial superior olive. The current is likely mediated by
  heteromultimers of Kv1.1 and Kv1.2 potassium channel subunits. The specific
  implementation is described in Rothman and Manis, J. Neurophysiol. 2003, in the
  appendix. Measurements were made from isolated neurons from adult guinea pig,
  under reasonably stringent voltage clamp conditions. The measured current is
  sensitive to the mamba snake toxin dendrotoxin-I.


Similar conductrances are found in the homologous neurons of the avian auditory
system (Reyes and Rubel; Zhang and Trussell; Rathouz and Trussell), and the
conductance described here, in the absence of more detailed kinetic measurements
, is probably suitable for use in modeling that system.


Original implementation by Paul B. Manis, April (JHU) and Sept, (UNC)1999.

File split implementation, February 28, 2004.

Contact: pmanis@med.unc.edu

Modified by Michael Eager March 31, 2005 to use with CVode

ENDCOMMENT

UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
        (nA) = (nanoamp)
}

NEURON {
        SUFFIX klt
        USEION k READ ek WRITE ik
        RANGE gkltbar, gklt, ik
        GLOBAL winf, zinf, wtau, ztau
}

:INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
        v (mV)
        celsius = 22 (degC)  : model is defined on measurements made at room temp in Baltimore
        dt (ms)
        ek = -77 (mV)
        gkltbar = 0.01592 (mho/cm2) <0,1e9>
        zss = 0.5   <0,1>   : steady state inactivation of glt
   q10=3
   qt   (/ms)
   q = 0.19245 (ms)
}

STATE {
        w z
}

ASSIGNED {
    ik (mA/cm2)
    gklt (mho/cm2)
    winf zinf
    wtau (ms) ztau (ms)
    }

:LOCAL wexp, zexp

BREAKPOINT {
   SOLVE states METHOD derivimplicit

   gklt = gkltbar*(w^4)*z
       ik = gklt*(v - ek)

}

DERIVATIVE states {
   rates(v)      :
       w' = (winf - w) / wtau
       z' = (zinf - z) / ztau
}

:UNITSOFF

INITIAL {
   qt = q10^((celsius - 22)/10) : if you don't like room temp, it can be changed!
   q=1/qt
    rates(v)
    w = winf
    z = zinf
}

:PROCEDURE states() {  :Computes state variables m, h, and n
:   trates(v)      :             at the current v and dt.
:   w = w + wexp*(winf-w)
:   z = z + zexp*(zinf-z)
:VERBATIM
:   return 0;
:ENDVERBATIM
:}



PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.

   TABLE winf,  zinf, wtau,ztau
   DEPEND dt, celsius FROM -150 TO 150 WITH 300
   

    winf = (1 / (1 + exp(-(v + 48) / 6)))^0.25
    zinf = zss + ((1-zss) / (1 + exp((v + 71) / 10)))

    wtau =  (q)*((100 / (6*exp((v+60) / 6) + 16*exp(-(v+60) / 45))) + 1.5)
    ztau =  (q)*((1000 / (exp((v+60) / 20) + exp(-(v+60) / 8))) + 50)
}


:UNITSON
