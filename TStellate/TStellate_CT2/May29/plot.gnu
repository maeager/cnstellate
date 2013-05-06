load '../Responses/default.gnu'
set xlabel "time (msec)"
set ylabel "Membrane Voltage (mV)"
#set title "CT1 optimisation results"
set terminal postscript eps enhanced defaultplex \
   leveldefault color \
   solid dashlength 1.0 linewidth 2.0 butt noclip \
   palfuncparam 2000,0.003 \
   "Helvetica" 18 

set output "plot.eps"
set border 3
unset x2tics
unset y2tics

wlsrts = `awk '/lsr/  {print $3}' TStellate.Fit.dat` 
whsrts = `awk '/hsr/  {print $3}' TStellate.Fit.dat` 
wdsrts = `awk '/ds/  {print $3}' TStellate.Fit.dat` 
wtvrts = `awk '/tv/  {print $3}' TStellate.Fit.dat` 
wglgts = `awk '/glg/  {print $3}' TStellate.Fit.dat` 
#gkht = `awk '/gkhtbar/  {print $3}' TStellate.Fit.dat` 
gleak = `awk '/gleak/  {print $3}' TStellate.Fit.dat` 
#erev = `awk '/erev/  {print $3}' TStellate.Fit.dat` 
error = `awk -F'=' '/final error/  {print $2}' TStellate.Fit.dat` 

set print "plot.params"
print "w_{LSR}          ",wlsrts
print ""
print "w_{HSR}          ",whsrts
print ""
print "w_{DS}           ",wlsrts
print ""
print "w_{TV}           ",whsrts
print ""
print "w_{GLG}          ",wlsrts
print ""
print "{g}_{leak}          ",gleak
print ""
print "Error           ",error
unset print

set label 1 system("cat plot.params") at 60,-30

plot 'TStellate.Fit.dat' u 1:2 t "Test" w l, \
     '' u 1:3 t "Reference" w l
