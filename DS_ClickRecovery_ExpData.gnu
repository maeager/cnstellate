# DS Click Recovery
#
# gnuplot is available from http://www.gnuplot.info

set term postscript enhanced mono solid "Helvetica" 12

set border 3
set ytics nomirror
set xtics nomirror
unset y2tics
unset x2tics
set tmargin 0
set rmargin 0
set size noratio
set xlabel "Delay (ms)"
set ylabel " Normalised Rate"
set xrange [0:16]
set yrange [0:1.1]

set output "DS_ClickRecoveryExpData.eps"
plot "DS_ClickRecovery_ExpData.dat" using 1:2 with linespoints pointsize 3 linewidth 5, "DS_ClickRecovery_ExpData.dat" using 1:3 with linespoints pointsize 3 linewidth 5,"DS_ClickRecovery_ExpData.dat" using 1:4 with linespoints pointsize 3 linewidth 5
