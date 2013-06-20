#!/usr/local/bin/gnuplot
set yrange [0:90]
set xrange [0:80]
set terminal postscript eps enhanced defaultplex \
   leveldefault color colortext \
   solid dashlength 1.0 linewidth 2.0 butt noclip \
   palfuncparam 2000,0.003 \
   "Helvetica" 18
#set ticslevel 0
set pm3d map
#set palette rgbformulae 22,13,-31
set xlabel "time (ms)" font "Helvetica,28" 
set ylabel "Level (dB SPL)" font "Helvetica,28" 
set key off
set output "psthVlevel.0.eps"
splot "psth.0.dat" using 1:2:3
set output "psthVlevel.1.eps"
splot "psth.1.dat" using 1:2:3
set output "psthVlevel.2.eps"
splot "psth.2.dat" using 1:2:3
set output "psthVlevel.3.eps"
splot "psth.3.dat" using 1:2:3

