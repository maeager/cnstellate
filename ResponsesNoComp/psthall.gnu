#!/usr/local/bin/gnuplot
set xrange [0:80]
set terminal postscript eps enhanced defaultplex \
   leveldefault color colortext \
   solid dashlength 1.0 linewidth 2.0 butt noclip \
   palfuncparam 2000,0.003 \
   "Helvetica" 18
#set ticslevel 0
set pm3d map
set xlabel "Time (ms)" font "Helvetica,24" 
set key off
#
#set yrange [0:100]
set ylabel "Channel" font "Helvetica,24" offset +1 

set output "psthallVAR-0.eps"
splot "./VAR/psth.0.dat" using 2:1:3
set output "psthallVAR-1.eps"
splot "./VAR/psth.1.dat" using 2:1:3
set output "psthallVAR-2.eps"
splot "./VAR/psth.2.dat" using 2:1:3
set output "psthallVAR-3.eps"
splot "./VAR/psth.3.dat" using 2:1:3
