#!/usr/local/bin/gnuplot
set yrange [50:800]
set xrange [200:64000]
set term postscript eps enhanced color
set output "response_area.0.eps"
#set ticslevel 0
set pm3d map
set palette rgbformulae 22,13,-31
set logscale x 10
splot "response_area.0.dat" using 2:1:3
#
