

NRNHOME="/usr/local/nrn"
ARCH=$(shell arch)
NRNMODL=$(ARCH)/bin/nrnivmodl
MPIMODL=$(HOME)/src/neuron/nrnmpi/$(NRNMODL)
GUIMODL=/usr/local/nrn/$(NRNMODL)
MODLFLAGS=-loadflags "$(shell pwd)/libresample-0.1.3/libresample.a"

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

remove-an:
	-rm -f $(ARCH)/an*
	-rm -f $(ARCH)/.libs/an*

libresample:
	(cd libresample-0.1.3; make)

rebuild-an: remove-an libresample

gui: libresample
	$(GUIMODL) $(MODLFLAGS)
	[ -d gui ] && rm -rf gui
	mv $(ARCH) gui
	(cd gui; sed -i 's_cnstellate/$(ARCH)_cnstellate/gui_g' special)

mpi: libresample
	$(MPIMODL) $(MODLFLAGS)
	[ -d mpi ] && rm -rf mpi
	mv $(ARCH) mpi
	(cd mpi; find -type f -exec sed -i 's_cnstellate/$(ARCH)_cnstellate/mpi_g' {} \; )
