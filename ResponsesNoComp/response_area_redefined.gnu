#!/usr/local/bin/gnuplot
set yrange [0:90]
set xrange [-2:1]
set term postscript eps enhanced color
set output "response_area_redefined.eps"
#set ticslevel 0
set pm3d map
set palette rgbformulae 22,13,-31
#set logscale x 10
#splot "response_area.0.dat" using (log10($2/5180)/log10(2)):1:3
splot "response_area.0.dat" using (($2-5180)/($2)):1:3
#splot "response_area.0.dat" using ($2<5810?(5180-$2)/$2:($2-5180)/$2):1:3
#