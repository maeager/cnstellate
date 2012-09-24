#!/usr/local/bin/gnuplot

set terminal postscript eps enhanced defaultplex \
   leveldefault color colortext \
   solid dashlength 1.0 linewidth 2.0 butt noclip \
   palfuncparam 2000,0.003 \
   "Helvetica" 18 
set yrange [0:90]
set ylabel "Level (dB SPL)" font "Helvetica,28"
set xlabel "Frequency (Hz)" font "Helvetica,28"
#set ticslevel 0
set pm3d map
set palette rgbformulae 22,13,-31
set logscale x 10
set xrange [200:64000]
set key off
#default response area based on CN network  separation
set output "response_area.0.eps"
splot "response_area.0.dat" using 2:1:3


#
# f1/f2
# set xrange [-2:1]
# set output "response_area_div.eps"
# splot "response_area.0.dat" using ((5180)/($2)):1:3
#splot "response_area.0.dat" using ($2<5810?(5180-$2)/$2:-(($2-5180)/$2)):1:3

#log2(f1/f2)
unset logscale
set xrange [-4:2]
set output "response_area_log2.0.eps"
set xlabel "Frequency (Octave re CF)" font "Helvetica,28"
splot "response_area.0.dat" using (log10(5180/$2)/log10(2)):1:3
