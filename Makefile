

NRNHOME="/home/eagerm/src/neuron/"
ARCH=$(shell arch)
NRNMODL=$(ARCH)/bin/nrnivmodl
MPIMODL=$(NRNHOME)/nrnmpi/$(NRNMODL)
GUIMODL=$(NRNHOME)/nrngui/$(NRNMODL)

MODLFLAGS=-loadflags "$(shell pwd)/libresample-0.1.3/libresample.a"
#MODLFLAGS=''


all: gui mpi

clean-mpi:
	-[ -d mpi ] && rm -rf mpi 

clean-gui:
	-[ -d gui ] && rm -rf gui

clean: clean-mpi clean-gui 
	-[ -d $(ARCH) ] && rm -rf $(ARCH)
	-rm -f *.ras *.hist *.connect *.curvs
	-rm -f out.dat *.[0-9]* *~ *.bak *.sav 

clean-libresample:
	make -C libresample-0.1.3 clean


clean-all: clean clean-libresample


gui:  clean-gui libresample-0.1.3/libresample.a
	$(GUIMODL) $(MODLFLAGS)
	mv $(ARCH) gui
	sed -i 's#cnstellate/$(ARCH)#cnstellate/gui#g' gui/special

mpi: clean-mpi libresample-0.1.3/libresample.a
	$(MPIMODL) $(MODLFLAGS)
	mv $(ARCH) mpi
	sed -i 's#cnstellate/$(ARCH)#cnstellate/mpi#g' mpi/special


libresample-0.1.3/libresample.a:
	(cd libresample-0.1.3; make CFLAGS="-fPIC ${CFLAGS} -fPIC" CXXFLAGS="-fPIC ${CXXFLAGS} -fPIC" libresample.a)


libresample_0.1.3.tar.gz:
	wget -c http://ftp.debian.org/debian/pool/main/libr/libresample/libresample_0.1.3.orig.tar.gz
	tar zxvf libresample_0.1.3.orig.tar.gz



