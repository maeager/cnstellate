set terminal postscript eps enhanced defaultplex \
   leveldefault mono \
   solid dashlength 1.0 linewidth 2.0 butt noclip \
   palfuncparam 2000,0.003 \
   "Helvetica" 18 
set border 3
unset x2tics
unset y2tics
set xlabel "time (ms)" font "Helvetica,28" 
set ylabel "spike count" font "Helvetica,28" offset +2,0
set key off
set output "psthsingle50.0.eps"
plot  "<grep '^50' ./50/psth.0.dat" using 2:3 w boxes
set output "psthsingle50.1.eps"
plot "<grep '^50' ./50/psth.1.dat" using 2:3 w boxes
set output "psthsingle50.2.eps"
plot  "<grep '^50' ./50/psth.2.dat" using 2:3 w boxes
set output "psthsingle50.3.eps"
plot "<grep '^50' ./50/psth.3.dat" using 2:3 w boxes
set output "psthsingle90.0.eps"
plot  "<grep '^50' ./90/psth.0.dat" using 2:3 w boxes
set output "psthsingle90.1.eps"
plot  "<grep '^50' ./90/psth.1.dat" using 2:3 w boxes
set output "psthsingle90.2.eps"
plot  "<grep '^50' ./90/psth.2.dat" using 2:3 w boxes
set output "psthsingle90.3.eps"
plot "<grep '^50' ./90/psth.3.dat" using 2:3 w boxes
