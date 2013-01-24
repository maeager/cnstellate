load '../Responses/default.gnu'
set xlabel "time (msec)"
set ylabel "Membrane Voltage (mV)"
#set title "CT2 optimisation results"
set terminal postscript eps enhanced defaultplex \
   leveldefault color \
   solid dashlength 1.0 linewidth 2.0 butt noclip \
   palfuncparam 2000,0.003 \
   "Helvetica" 18 

set output "plotCT2.eps"
set border 3
unset x2tics
unset y2tics

whsrts = `awk '/hsr/  {print $3}' TStellate_CT2.Fit.dat` 
gleak = `awk '/gleak/  {print $3}' TStellate_CT2.Fit.dat` 
erev = `awk '/erev/  {print $3}' TStellate_CT2.Fit.dat` 
error = `awk -F'=' '/final error/  {print $2}' TStellate_CT2.Fit.dat` 

set print "plotCT2.params"
print "w_{HSR}          ",whsrts
print ""
print "{g}_{leak}          ",gleak
print ""
print "E_{rev}            ",erev
print ""
print "Error           ",error
unset print

set label 1 system("cat plotCT2.params") at 60,-30

plot 'TStellate_CT2.Fit.dat' u 1:2 t "Test" w l, \
     '' u 1:3 t "Reference" w l
