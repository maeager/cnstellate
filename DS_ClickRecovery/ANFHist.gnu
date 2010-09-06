

set term postscript enhanced color solid "Helvetica" 12
set output "ANF_histogram.eps"

# set border 3
# set ytics nomirror
# set xtics nomirror
# unset y2tics
# unset x2tics
# set tmargin 0
# set rmargin 0
# set size noratio
set xlabel "time (ms)" font "Helvetica,22"
set ylabel "Freq. Channel" font "Helvetica,22"
set xrange [0:100]
set yrange [0:1500]

set pm3d map
set palette rgbformulae 22,13,-31
splot "ANFHist.dat" matrix
