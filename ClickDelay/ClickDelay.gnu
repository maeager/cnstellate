# gnuplot file to produce a figure similar to Figure 2
# in Davison A.P., Feng J. and Brown D. (2003) J. Neurophysiol.

# Note that the figure produced will not match exactly the
# published figure due to differences in the sequence of
# random numbers used to set up the network.

# gnuplot is available from http://www.gnuplot.info

set term postscript portrait enhanced mono solid "Helvetica" 8
set output "ClickDelay.eps"

set size 0.49,0.5
set multiplot
set size 0.48,0.1
unset key
unset xtics
unset border
unset ytics
set tmargin 0
set rmargin 0

#  stellate cell
set label 1 '{/Helvetica=14 A}' at screen 0,0.49
#set arrow 1 from 1200,-30 to 1000,-30 nohead linewidth 2
#set arrow 2 from 1200,-30 to 1200,20 nohead linewidth 2
#set label 2 "200 ms" at 1000,-45
#set label 3 " 50 mV" at 1200,-10 left
set origin 0.01,0.4
plot [0:0.8][-75:0] "ClickDelay.curvs" using 1:3 with lines linewidth 0.5

# granule cell
#unset arrow
#unset label 2
#unset label 3
set label 1 '{/Helvetica=14 B}' at screen 0,0.39
set origin 0.01,0.3
plot [0:0.8][-75:0] "ClickDelay.curvs" using 1:4 with lines linewidth 0.5

# dstellate
#unset arrow
#unset label 2
#unset label 3
set label 1 '{/Helvetica=14 C}' at screen 0,0.29
set origin 0.01,0.2
plot [0:0.8][-75:0] "ClickDelay.curvs" using 1:2 with lines linewidth 0.5

# Raster
#unset arrow
#unset label 2
unset border
set label 1 '{/Helvetica=14 D}' at screen 0,0.19
set origin 0.01,0.1
plot [0:80][0:10000] "ClickDelay.an.HSR.ras" using 4:3 with dots, "ClickDelay.an.HSR.ras" using 4:($3) with dots linetype 1

# Histogram
#set border 1 linewidth 0.5
set label 1 '{/Helvetica=14 E}' at screen 0,0.09
set origin 0.01,0.0
set arrow 1 from 10,1000 to 10,1100 nohead linewidth 2
set label 2 "100 spikes per second" at 15,1010 right
plot [0:80][0:4000] "ClickDelay.an.HSR.hist" using 0:1 with steps linewidth 0.5




