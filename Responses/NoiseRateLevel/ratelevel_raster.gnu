#!/usr/local/bin/gnuplot
set term postscript eps enhanced color
set output "ratelevel_raster.0.eps"
set nokey
set noxtics
set noborder
set noytics
#set tmargin 0
#set rmargin 0
plot [0:80][0:90] "raster.0.dat" using 1:2 with dots linetype 1
