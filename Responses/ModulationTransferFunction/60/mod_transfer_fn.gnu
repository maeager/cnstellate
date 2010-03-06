#!/usr/local/bin/gnuplot
set term postscript eps enhanced color
set output "mtf_rate_raster.0.eps"
set nokey
set noxtics
set noborder
set noytics
#set tmargin 0
#set rmargin 0
plot "raster.0.dat" using 1:2 with dots linetype 1
