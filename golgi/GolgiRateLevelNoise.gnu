# Golgi Rate level- noise plot
# gnuplot is available from http://www.gnuplot.info

set term postscript portrait enhanced mono solid "Helvetica" 8
set output "GolgiRateLevelNoise.eps"


set size 0.48,0.2
unset key
unset x2tics
set border 3
set ytics nomirror
set xtics nomirror
unset y2tics
set tmargin 0
set rmargin 0

#  
set label 1 '{/Helvetica=14 A}' at screen 0,0.49
set origin 0.01,0.4
plot "GolgiRateLevelNoise.dat" using 1:2 with lines linewidth 2




