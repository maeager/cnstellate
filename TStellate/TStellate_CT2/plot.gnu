load '../Responses/default.gnu'
set xlabel "time (msec)"
set ylabel "Membrane Voltage (mV)"
#set title "CT1 optimisation results"
#set term post enh mono
set terminal postscript eps enhanced defaultplex \
   leveldefault color \
   solid dashlength 1.0 linewidth 2.0 butt noclip \
   palfuncparam 2000,0.003 \
   "Helvetica" 18 

set output "plot.eps"
set border 3
unset x2tics
unset y2tics

whsrts = `awk '/param.*hsr/  {print $3}' TStellate.Fit.dat` 
wlsrts = `awk '/param.*lsr/  {print $3}' TStellate.Fit.dat` 
wdsts = `awk '/param.*ds/  {print $3}' TStellate.Fit.dat` 
wtvts = `awk '/param.*tv/  {print $3}' TStellate.Fit.dat` 
wglgts = `awk '/param.*glg/  {print $3}' TStellate.Fit.dat` 
gleak = `awk '/leak/  {print $3}' TStellate.Fit.dat`
error = `awk -F'=' '/final error/  {print $2}' TStellate.Fit.dat` 

set print "plot.params"
print "w_{HSR}          ",whsrts
print ""
print "w_{LSR}          ",wlsrts
print ""
print "w_{DS}          ",wdsts
print ""
print "w_{TV}          ",wtvts
print ""
print "w_{GLG}          ",wglgts
print ""
print "g_{leak}          ",gleak
print ""
print "Error           ",error
unset print
# Position a red square with lower left at 0,0 and upper right at 2,3

set style rect fc lt -1 fs solid 0.15 noborder

set obj rect from -10, graph 0 to 0 , graph 1
set obj rect from 60 , graph 0 to 62, graph 1
set obj rect from 5  , graph 0 to 10, graph 1
set obj rect from 20 , graph 0 to 30, graph 1

set obj rect from 20 , graph 0 to 22, graph 1 fc lt -1 fs solid 0.25 noborder
set obj rect from 50 , graph 0 to 52, graph 1 fc lt -1 fs solid 0.25 noborder

set label 1 system("cat plot.params") at 60,-30
plot 'TStellate.Fit.dat' u 1:3 t "Test" w l lc "black", \
     '' u 1:2 t "Reference" w l lc "black
