#!/bin/sh
rm -f i686/an*
rm -f i686/.libs/an*
nrnivmodl
(cd i686; gcc -shared .libs/flushf.o .libs/SGC_VecStim.o .libs/SGCfast.o \
.libs/an_bruce.o \
.libs/an_zilany_v4.o \
.libs/ka.o .libs/klt.o .libs/rm.o .libs/rm_vect.o .libs/mod_func.o \
  -Wl,-rpath -Wl,/home/eagerm/Neuron/nrn-7.2/nrngui/i686/lib -Wl,-rpath -Wl,/home/eagerm/Neuron/nrn-7.2/nrngui/i686/lib -L/usr/local/lib -L/home/eagerm/Neuron/nrn-7.2/nrngui/i686/lib /home/eagerm/Neuron/nrn-7.2/nrngui/i686/lib/libnrnoc.so /home/eagerm/Neuron/nrn-7.2/nrngui/i686/lib/liboc.so /home/eagerm/Neuron/nrn-7.2/nrngui/i686/lib/libmemacs.so /home/eagerm/Neuron/nrn-7.2/nrngui/i686/lib/libnrnmpi.so /home/eagerm/Neuron/nrn-7.2/nrngui/i686/lib/libscopmath.so /home/eagerm/Neuron/nrn-7.2/nrngui/i686/lib/libsparse13.so -lreadline -lncurses /home/eagerm/Neuron/nrn-7.2/nrngui/i686/lib/libnrniv.so /home/eagerm/Neuron/nrn-7.2/nrngui/i686/lib/libivoc.so /home/eagerm/Neuron/nrn-7.2/nrngui/i686/lib/libneuron_gnu.so /home/eagerm/Neuron/nrn-7.2/nrngui/i686/lib/libmeschach.so /home/eagerm/Neuron/nrn-7.2/nrngui/i686/lib/libsundials.so -lm -ldl \
../libresample.a \
-march=atom -mtune=atom -mssse3 -mfpmath=sse  -pthread -Wl,-soname -Wl,libnrnmech.so.0 -o .libs/libnrnmech.so.0.0.0)
