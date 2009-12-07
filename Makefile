ARCH=$(shell arch)
NRNMODL = ~/neuron/nrn-7.1/nrnmpi/$(ARCH)/bin/nrnivmodl
IVMODL=~/neuron/nrngui/$(ARCH)/bin/nrnivmodl

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
	$(IVMODL)
	mv $(ARCH) gui
	
mpi:
	$(NRNMODL)
	mv $(ARCH) mpi
