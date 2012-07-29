

NRNHOME="/usr/local/nrn"
ARCH=$(shell arch)
NRNMODL=$(ARCH)/bin/nrnivmodl
MPIMODL=$(HOME)/src/neuron/nrnmpi/$(NRNMODL)
GUIMODL=/usr/local/nrn/$(NRNMODL)
MODLFLAGS=-loadflags "$(shell pwd)/libresample-0.1.3/libresample.a"

all: clean-all libresample gui mpi

clean-mpi:
	-rm -rf mpi 
	-rm -rf $(ARCH)

clean: clean-mpi
	-rm -f *.ras *.hist *.connect *.curvs
	-rm -f out.dat *.[0-9]* *~ *.bak *.sav 

clean-all: clean
	-rm -rf mpi gui
	-(cd libresample-0.1.3; make clean)


gui:  # libresample-0.1.3/libresample.a
	$(GUIMODL) $(MODLFLAGS)
	[ -d gui ] && rm -rf gui
	mv $(ARCH) gui
	(cd gui; sed -i 's_cnstellate/$(ARCH)_cnstellate/gui_g' special)

mpi: libresample-0.1.3/libresample.a
	$(MPIMODL) $(MODLFLAGS)
	[ -d mpi ] && rm -rf mpi
	mv $(ARCH) mpi
	(cd mpi; find -type f -exec sed -i 's_cnstellate/$(ARCH)_cnstellate/mpi_g' {} \; )

libresample-0.1.3/libresample.a: libresample_0.1.3.tar.gz
	(cd libresample-0.1.3; make CFLAGS="-fPIC $CFLAGS -fPIC" CXXFLAGS="-fPIC $CXXFLAGS -fPIC")

libresample_0.1.3.tar.gz:
	wget -c http://ftp.debian.org/debian/pool/main/libr/libresample/libresample_0.1.3.orig.tar.gz
	tar zxvf libresample_0.1.3.orig.tar.gz



