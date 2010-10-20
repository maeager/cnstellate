#!/bin/sh
rm -f i686/an*
rm -f i686/.libs/an*
/usr/local/nrn/i686/bin/nrnivmodl
(cd i686; gcc -shared .libs/flushf.o .libs/SGC_VecStim.o .libs/SGCfast.o \
.libs/an_bruce.o \
.libs/an_zilany_v4.o \
.libs/ka.o .libs/klt.o .libs/rm.o .libs/rm_vect.o .libs/mod_func.o \
  -Wl,-rpath -Wl,/usr/local/nrn/i686/lib -Wl,-rpath -Wl,/usr/local/nrn/i686/lib -L/usr/local/lib -L/usr/local/nrn/i686/lib /usr/local/nrn/i686/lib/libnrnoc.so /usr/local/nrn/i686/lib/liboc.so /usr/local/nrn/i686/lib/libmemacs.so /usr/local/nrn/i686/lib/libnrnmpi.so /usr/local/nrn/i686/lib/libscopmath.so /usr/local/nrn/i686/lib/libsparse13.so -lreadline -lncurses /usr/local/nrn/i686/lib/libnrniv.so /usr/local/nrn/i686/lib/libivoc.so /usr/local/nrn/i686/lib/libneuron_gnu.so /usr/local/nrn/i686/lib/libmeschach.so /usr/local/nrn/i686/lib/libsundials.so -lm -ldl \
../libresample.a \
-march=atom -mtune=atom -mssse3 -mfpmath=sse  -pthread -Wl,-soname -Wl,libnrnmech.so.0 -o .libs/libnrnmech.so.0.0.0)
