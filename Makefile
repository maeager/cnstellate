ARCH=$(shell arch)
NRNMODL = $HOME/Neuron/nrn/nrnmpi/$(ARCH)/bin/nrnivmodl
GUIMODL=$HOME/Neuron/nrn/nrngui/$(ARCH)/bin/nrnivmodl

all: clean-all gui mpi

clean-mpi:
	-rm -rf *.[0-9]*
	-rm -f PI*
	
clean: clean-mpi
	-rm -f *.ras *.hist *.connect *.curvs
	-rm -f out.dat *.[0-9]* *~ *.bak *.sav 

clean-all: clean
	-rm -rf mpi gui

gui:
	$(GUIMODL)
	mv $(ARCH) gui
	cd gui
	sed -i 's/cnstellate\/i686/cnstellate\/gui/g' *
	
mpi:
	$(NRNMODL)
	mv $(ARCH) mpi
	cd mpi
	sed -i 's/cnstellate\/i686/cnstellate\/mpi/g' *
