
NRNMODL = /home/meager/nrn-6.1/nrnmpi/ia64/bin/nrnivmodl
IVMODL=/home/meager/nrn-6.1/nrniv/ia64/bin/nrnivmodl

all: gui mpi

clean-mpi:
	-rm -rf *.[0-9]*
	-rm -f PI*
	
clean: clean-mpi
	-rm -f *.ras *.hist *.connect *.curvs
	-rm -f out.dat *.[0-9]* *~ *.bak *.sav 

clean-all: clean
	-rm -rf ia64 gui

gui:
	$(IVMODL)
	mv ia64 gui
	
mpi:
	$(NRNMODL)
