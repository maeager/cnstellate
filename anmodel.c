/* Created by Language version: 6.0.2 */
/* NOT VECTORIZED */
#include <stdio.h>
#include <math.h>
#include "scoplib.h"
#undef PI

#include "md1redef.h"
#include "section.h"
#include "nrnoc_ml.h"
#include "md2redef.h"

#if METHOD3
extern int _method3;
#endif

#undef exp
#define exp hoc_Exp
extern double hoc_Exp();
/*SUPPRESS 761*/
/*SUPPRESS 762*/
/*SUPPRESS 763*/
/*SUPPRESS 765*/
extern double *getarg();
static double *_p; static Datum *_ppvar;

#define delta_t dt

#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
static int hoc_nrnpointerindex =  -1;
/* external NEURON variables */
extern double dt;
extern double t;
/* declaration of user functions */
static int _hoc_install_anmodel();
static int _mechtype;
extern int nrn_get_mechtype();
static _hoc_setdata()
{
    Prop *_prop, *hoc_getdata_range();
    _prop = hoc_getdata_range("anmodel");
    _p = _prop->param; _ppvar = _prop->dparam;
    ret(1.);
}
/* connect user functions to hoc names */
static IntFunc hoc_intfunc[] = {
    "setdata_anmodel", _hoc_setdata,
    "install_anmodel", _hoc_install_anmodel,
    0, 0
};
/* declare global and static user variables */
/* some parameters have upper and lower limits */
static HocParmLimits _hoc_parm_limits[] = {
    0, 0, 0
};
static HocParmUnits _hoc_parm_units[] = {
    0, 0
};
static double v = 0;
/* connect global user variables to hoc */
static DoubScal hoc_scdoub[] = {
    0, 0
};
static DoubVec hoc_vdoub[] = {
    0, 0, 0
};
static double _sav_indep;
static nrn_alloc(), nrn_init(), nrn_state();
/* connect range variables in _p that hoc is supposed to know about */
static char *_mechanism[] = {
    "6.0.2",
    "anmodel",
    0,
    0,
    0,
    0
};

static nrn_alloc(_prop)
Prop *_prop;
{
    Prop *prop_ion, *need_memb();
    double *_p; Datum *_ppvar;
    _p = nrn_prop_data_alloc(_mechtype, 0);
    /*initialize range parameters*/
    _prop->param = _p;
    _prop->param_size = 0;

}
static _initlists();
_anmodel_reg()
{
    int _vectorized = 0;
    _initlists();
    hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
    ivoc_help("help ?1 anmodel /cygdrive/c/nrn57/work/modfiles/anmodel.mod\n");
}
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse = 1;
static _modl_cleanup()
{
    _match_recurse = 1;
}
static install_anmodel();

/*VERBATIM*/
extern double* hoc_pgetarg();
extern int vector_instance_px();
extern void* vector_arg();
extern char *gargstr();


/*VERBATIM*/
#include "./HeinzARLO/ARLO.c"


//sout.an4(tdres,cf,spont,model, species,ifspike,wavefile)
#define NSTIMMAX 100000 // enough points for durations up to 250 msec
#define MAXCHARBUF 100
static double an4(void *vv)
{
//stimulus file  FILE *fpstim;
    double *stim;  //Input stimulus vector

    double tdres, cf, spont;
    int model;
    int species;
    int ifspike;
    int nowstim, totalstim;
    double  x;
    int     length, mrows, ncols;

    double *sout;   //Output vector containing all data for each channel (#banks X datapoints)
    int nsout, nstim, i;

    char *wavefile;
    double  cfhi, cflo;
    int banks;

// stim = makevector(NSTIMMAX);


//Get Instance Sout matrix for synapse data across filterbanks
    nsout = vector_instance_px(vv, &sout);
    nstim = vector_arg_px(1, &stim);

//Get Input arguments
    if (ifarg(7) != 1) {  //Must be changed if more input arguments added
        hoc_execerror("ERROR: input syntax must be sout.an4( stim, tdres,  cf, spontrate, model, species, ifspike)", 0);
        return 0;
    }
    //TDRES  resolution of wavefile
    tdres = (double)(*getarg(2));
    //cf of fiber
    cf = (int)(*getarg(3));
    //Spontaneous rate of ANF
    spont = (double)(*getarg(4));
    //Model
    model = (double)(*getarg(5));
    //Species
    species = (double)(*getarg(6));
    //Spike generator on =1
    ifspike = (double)(*getarg(7));
    /* //Wavefile
         wavefile=gargstr(1);
     totalstim=0;
     for(nowstim=0; nowstim<NSTIMMAX; nowstim++) stim[nowstim] = 0.0;
      printf("\nReading in stimulus \n");
       fpstim = fopen(wavefile,"r");
       nowstim=0;
       while((fscanf(fpstim,"%le",&stim[nowstim]) != EOF) & (nowstim<NSTIMMAX)) nowstim++;
       printf("%d points read in\n",nowstim);
       if(totalstim < nowstim) totalstim=nowstim;
       fclose(fpstim);
       if(nowstim>=NSTIMMAX) hoc_execerror("Entire stim not read in: increase NSTIMMAX");
       if(totalstim>=NSTIMMAX) hoc_execerror("Entire reptim too long: increase NSTTMMAX");
     length = totalstim;
    */
    length = nstim;

    //Number of filter banks
//   banks = (int)(*getarg(3));

    //lowest cf
//   cflo = (double)(*getarg(4));
    //highest cf
//   cfhi = (double)(*getarg(5));
// if(cfhi<cflo) hoc_execerror("cfhi should be greater than cflo", 0);


    return an_arlo(tdres, cf, spont, model, species, ifspike, stim, sout, length);

}


static int  install_anmodel()
{

    /*VERBATIM*/

    install_vector_method("an4", an4);

// install_vector_method("SpkGen", SpkGen);

    return 0;
}
static int _hoc_install_anmodel()
{
    double _r;
    _r = 1.;
    install_anmodel() ;
    ret(_r);
}

static initmodel()
{
    int _i; double _save;_ninits++;
    {

    }
}

static nrn_init(_ml, _type) _Memb_list* _ml; int _type;
{
    Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
    _cntml = _ml->_nodecount;
    for (_iml = 0; _iml < _cntml; ++_iml) {
        _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
        if (use_cachevec) {
            _v = VEC_V(_ni[_iml]);
        } else
#endif
        {
            _nd = _ml->_nodelist[_iml];
            _v = NODEV(_nd);
        }
        v = _v;
        initmodel();
    }
}

static double _nrn_current(_v) double _v;
{
    double _current = 0.;v = _v;{
    } return _current;
}

static nrn_state(_ml, _type) _Memb_list* _ml; int _type;
{
    double _break, _save;
    Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
    _cntml = _ml->_nodecount;
    for (_iml = 0; _iml < _cntml; ++_iml) {
        _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
        _nd = _ml->_nodelist[_iml];
#if CACHEVEC
        if (use_cachevec) {
            _v = VEC_V(_ni[_iml]);
        } else
#endif
        {
            _nd = _ml->_nodelist[_iml];
            _v = NODEV(_nd);
        }
        _break = t + .5 * dt; _save = t; delta_t = dt;
        v = _v;
        {
        }}

}

static terminal() {}

static _initlists()
{
    int _i; static int _first = 1;
    if (!_first) return;
    _first = 0;
}
