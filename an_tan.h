#ifndef AN_TAN_H
#define AN_TAN_H

#define usemiddleear 1         /* 1-> include middle ear model; 0-> don't include middle ear */



typedef struct Tanstruct Tanmodel;

struct Tanstruct{

//General information about the model 
  	double tdres, cf, spont,fp1;
  	int species, model;
  	int nstim; 
//  runtime usage 
  	int    ifspike;
  	int    control_type;              
// 0 ->without suppression 
// 1 ->with suppression    
	double PI;
//      double *soundin;
//	double *meout;
//	double *soundout;
//	double *control_signal;
//	double *ihc_out;
//	double *sout;

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

//double an_tan(double tdres,double cf,double spont, double NLsuppression, double species, double ifspike,double *stim,double *sout,int length);

//The following files from Tan and Carney 2003 have been combined
// into  anmod3m.c
//#include "mycomplex.h"
//#include "middleear.c"
//#include "controlpath.c"
//#include "signalpath.c"
//#include "ihczxd2001.c"
#include "anmod3m.c"
#endif
