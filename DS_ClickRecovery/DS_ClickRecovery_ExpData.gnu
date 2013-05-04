# DS Click Recovery
#
# gnuplot is available from http://www.gnuplot.info

set term postscript enhanced color solid "Helvetica" 16

set border 3
set ytics nomirror
set xtics nomirror
unset y2tics
unset x2tics
set tmargin 0
set rmargin 0
set size noratio
set xlabel "Delay (ms)" font "Helvetica,22"
set ylabel " Normalised Rate" "Helvetica,22"
set xrange [0:20]
set yrange [0:1.2]

#set output "DSClickRecoveryExpData.eps"
#plot "DS_ClickRecovery_ExpData.dat" using 1:2 with linespoints pointsize 3 linewidth 5, "DS_ClickRecovery_ExpData.dat" using 1:3 with #linespoints pointsize 3 linewidth 5,"DS_ClickRecovery_ExpData.dat" using 1:4 with linespoints pointsize 3 linewidth 5


set output "DS_ClickRecovery_result.eps"
plot "DS_ClickRecovery.Fit.dat" using 1:2 with linespoints title "No GABA" pointsize 3 linewidth 5, "DS_ClickRecovery.Fit.dat" using 1:3 title "Target" with linespoints pointsize 3 linewidth 5,"DS_ClickRecovery.Fit.dat" using 1:5 with linespoints title "Result" pointsize 3 linewidth 5 
