#!/usr/local/bin/gnuplot
set yrange [0:90]
set xrange [0:325]
set term postscript eps enhanced color
set output "psthVlevel.0.eps"
#set ticslevel 0
set pm3d map
#set palette rgbformulae 22,13,-31
splot "psth.0.dat" using 1:2:3
#
unset yrange
unset xrange
set yrange [0:100]
set xrange [0:325]
set output "psth250.0.eps"
splot "./250/psth.0.dat" using 2:1:3
set output "psth250.1.eps"
splot "./250/psth.1.dat" using 2:1:3
set output "psth250.2.eps"
splot "./250/psth.2.dat" using 2:1:3
set output "psth250.3.eps"
splot "./250/psth.3.dat" using 2:1:3

set output "psth750.0.eps"
splot "./750/psth.0.dat" using 2:1:3
set output "psth750.1.eps"
splot "./750/psth.1.dat" using 2:1:3
set output "psth750.2.eps"
splot "./750/psth.2.dat" using 2:1:3
set output "psth750.3.eps"
splot "./750/psth.3.dat" using 2:1:3
