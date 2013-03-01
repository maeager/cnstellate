:  Spiral Ganglion Cell- NEURON Artificial Cell
:  VecStim version

NEURON   {
  ARTIFICIAL_CELL SGC_VecStim
  RANGE spikecount,  spont, cf


}

VERBATIM
extern double* vector_vec();
extern int vector_capacity();
extern void* vector_arg();
extern void vector_resize();

ENDVERBATIM


PARAMETER {
   spikecount   =0.0
    spont       = 50.0  <0,150>    : Spontaneous Rate
    cf       = 1000    (Hz) <20, 40000>

}



ASSIGNED {
   spkindex
   eventtime (ms)
   spiketimes

}

INITIAL {
   spkindex = 0
   if (GetSpike()) {
      net_send(eventtime - t, 1)
   }
}


NET_RECEIVE (w){
   if (flag == 0) {:0 progagate synaptic external event
      net_event(t)

   }
   if (flag == 1) {  :1 self event
      net_event(t)

      if (GetSpike())  {
         net_send(eventtime-t, 1)
      }
   }
}


FUNCTION GetSpike() {
VERBATIM
  { void* vv; int i, size, trigger; double* spksout;
   i = (int)spkindex;
   trigger =0;
   if (i >= 0) {
      vv = *((void**)(&spiketimes));
      if (vv) {
         size = (int)spikecount;
         spksout = vector_vec(vv);
         if (i < size) {
            eventtime = spksout[i];
            spkindex += 1.;
            trigger = 1;

         }else{
            spkindex = -1.;
         }
      }else{
         spkindex = -1.;
      }
   }
   return trigger;
  }
ENDVERBATIM
}

PROCEDURE SetSpikes() {
VERBATIM

   void** vv;
   vv = (void**)(&spiketimes);
   *vv = (void*)0;
   if (ifarg(1)) {
      *vv = vector_arg(1);
   }
ENDVERBATIM
}
