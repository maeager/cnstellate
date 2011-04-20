#!/bin/sh

# Force rebuilding of AN model files
rm -f i686/an*
rm -f i686/.libs/an*

#Check libresample
(cd libresample-0.1.3; make)

NRNIVPATH=${HOME}/src/neuron/nrnmpi/i686

${NRNIVPATH}/bin/nrnivmodl -loadflags `pwd`/libresample-0.1.3/libresample.a


# previous version recompiled the libnrnmech library with the libresample library included
# ${NRNIVPATH}/bin/nrnivmodl &&  \
# (cd i686; gcc -shared .libs/flushf.o .libs/SGC_VecStim.o .libs/SGCfast.o \
#   .libs/an_bruce.o \
#   .libs/an_zilany_v4.o \
#   .libs/ka.o .libs/klt.o .libs/rm.o .libs/rm_vect.o .libs/mod_func.o \
#   -Wl,-rpath -Wl,${NRNIVPATH}/lib -Wl,-rpath -Wl,${NRNIVPATH}/lib -L${NRNIVPATH}/lib \
#   ${NRNIVPATH}/lib/libnrnoc.so ${NRNIVPATH}/lib/liboc.so ${NRNIVPATH}/lib/libmemacs.so \
#   ${NRNIVPATH}/lib/libnrnmpi.so ${NRNIVPATH}/lib/libscopmath.so ${NRNIVPATH}/lib/libsparse13.so \
#   -lreadline -lncurses ${NRNIVPATH}/lib/libnrniv.so ${NRNIVPATH}/lib/libivoc.so \
#   ${NRNIVPATH}/lib/libneuron_gnu.so ${NRNIVPATH}/lib/libmeschach.so ${NRNIVPATH}/lib/libsundials.so \
#   -lm -ldl ../libresample-0.1.3/libresample.a -march=atom -mtune=atom -mssse3 -mfpmath=sse  -pthread \
#   -Wl,-soname -Wl,libnrnmech.so.0 -o .libs/libnrnmech.so.0.0.0)

