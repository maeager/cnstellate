#! /usr/bin/gnuplot 
# tMTF and rMTF  plot for AN
# gnuplot is available from http://www.gnuplot.info

set term postscript portrait enhanced mono solid "Helvetica" 8
set output "ANF_MTF.eps"


set size 0.48,0.25
unset key
unset x2tics
set border 3
set ytics nomirror
set xtics nomirror
unset y2tics
set tmargin 0
set rmargin 0

#  ANF rMTF
set label 1 '{/Helvetica=14 A}' at screen 0,0.49
set origin 0.01,0.25
plot "ANF_MTF.dat" using 1:2 with lines linewidth 2

set label 2 '{/Helvetica=14 B}' at screen 0,0.49
set origin 0.01,0.0
plot "ANF_MTF.dat" using 1:3 with lines linewidth 2


