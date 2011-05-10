ARCH=$(shell arch)
NRNMODL=$(ARCH)/bin/nrnivmodl
MPIMODL=$(HOME)/src/neuron/nrnmpi/$(NRNMODL)
GUIMODL=$(HOME/src/neuron/nrngui/$(NRNMODL)
MODLFLAGS=-loadflags $(shell pwd)/libresample-0.1.3/libresample.a

all: clean-all libresample gui mpi

clean-mpi:
	-rm -rf *.[0-9]*
	-rm -f PI*

clean: clean-mpi
	-rm -f *.ras *.hist *.connect *.curvs
	-rm -f out.dat *.[0-9]* *~ *.bak *.sav 

clean-all: clean
	-rm -rf mpi gui
	-(cd libresample-0.1.3; make clean)

remove:
	-rm -f i686/an*
	-rm -f i686/.libs/an*

libresample:
	(cd libresample-0.1.3; make)

rebuild: remove libresample

gui:
	$(GUIMODL)
	mv $(ARCH) gui
	cd gui
	sed -i 's/cnstellate\/i686/cnstellate\/gui/g' *

mpi:
	$(MPIMODL)
	mv $(ARCH) mpi
	cd mpi
	sed -i 's/cnstellate\/i686/cnstellate\/mpi/g' *
