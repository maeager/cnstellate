#!/usr/bin/gnuplot -persist
# set terminal postscript eps enhanced defaultplex \
#    leveldefault mono solid\
#    dashlength 1.0 linewidth 2.0 butt noclip \
#    palfuncparam 2000,0.003 \
#    "Helvetica" 18
# #color solid
# set output "multiplotCV.eps"

#set size 0.49,0.5
set multiplot layout 2,2
set style fill solid 1.0
set ylabel "Rate (sp/s)" font "Helvetica,18"
set xlabel "Time (ms)" font "Helvetica,18"

plot [-10:90] 'psth_0.25.dat' i 0 u 1:2 w boxes lc "black"
set ylabel "Rate (sp/s)" font "Helvetica,18"
set xlabel "Interval (ms)" font "Helvetica,18"
#plot [0:15] 'psth_0.25.dat' i 2 u 1:2 w boxes lc "black"

#plot [0:15] 'psth_0.2.dat' i 2 u 1:3 w boxes lc "black"


set xtics nomirror out ("0-10" 0,"10-20" 1,"20-30" 2,"30-40" 3,"40-50" 4)
set xrange [-0.5:4.5]
set yrange [0.1:0.5]
set border 3
#set ytics nomirror out 0.1,0.1,0.5
set ylabel "ISI (ms)" font "Helvetica,18"
set xlabel "Time Window" font "Helvetica,18"
set style line 1 lc rgb "#0060ad" lt 1 lw 2 pt 7 ps 1.5 # --beautiful blue line and filled circle
set style line 2 lc rgb '#dd181f' lt 1 lw 2 pt 5 ps 1.5 # --beautiful red line and square
set style line 3 lc rgb '#00ad06' lt 1 lw 2 pt 9 ps 1.5 # --green line and triangle
set yrange [0:10]

plot  'psth_10.dat' i 1 u 1:3 w lp ls 3, \
 'psth_10.dat' i 1 u 1:4 w lp ls 2
set ylabel "CV " font "Helvetica,18"
set yrange [0.1:0.5]
set ytics nomirror out 0.1,0.1,0.5
#plot 'psth_10.dat' i 1 u 1:2 w lp ls 1

#load '../Responses/default.gnu'
set xlabel "time (msec)"
set ylabel "Membrane Voltage (mV)"
set border 3
unset tics
unset x2tics
unset y2tics
set style rect fc lt -1 fs solid 0.15 noborder
set obj rect from -10, graph 0 to 0 , graph 1
set obj rect from 60 , graph 0 to 62, graph 1
set obj rect from 5  , graph 0 to 10, graph 1
set obj rect from 20 , graph 0 to 30, graph 1
set obj rect from 20 , graph 0 to 22, graph 1 fc lt -1 fs solid 0.25 noborder
set obj rect from 50 , graph 0 to 52, graph 1 fc lt -1 fs solid 0.25 noborder
set tics nomirror out
set autoscale  y
set autoscale x

plot 'TStellate.Fit.dat' u 1:3 t "Test" w l lc "black", \
     '' u 1:2 t "Reference" w l lc "black

#load '../Responses/default.gnu'
set xtics nomirror out ("0-10" 0,"10-20" 1,"20-30" 2,"30-40" 3,"40-50" 4)
 set xrange [-0.5:4.5]
set yrange [0.1:0.5]
set border 3
#set ytics nomirror out 0.1,0.1,0.5
set ylabel "ISI (ms)" font "Helvetica,18"
set xlabel "Time Window" font "Helvetica,18"
set yrange [0:5]
set ylabel "CV " font "Helvetica,18"
set yrange [0.1:0.5]
set ytics nomirror out 0.1,0.1,0.5
plot 'psth_10.dat' i 1 u 1:2 w lp ls 1, \
     '../TStellate/PaoliniBalancedInh-Fig2.png.dat' i 0 u 1:2:3 w yerr lc rgb '#0060ad', \
     '' i 0 u 1:2 w lp lc rgb '#0060ad' lt 1 lw 2 pt 5 ps 1.5

unset multiplot

# set output "psthcv.eps"

# set dummy jw
# set grid x y2
# set key center top title " "
# set xrange [0 : 60]
# set x2range [0 : 60]
# set ytics nomirror
# set y2tics
# set tics out
# set autoscale  y
# set autoscale y2
# set style fill solid 1.0
# set ylabel "Rate (sp/s)" font "Helvetica,18"
# set xlabel "Time (ms)" font "Helvetica,18"
# set y2label "CV " font "Helvetica,18"
# #set y2range [0.1:0.5]
# plot  'psth_1.dat'  i 0 u 1:2 axes x1y1 w boxes lc "black" , 'psth_1e+01.dat' i 1 u (($1*10)+5):2 axes x2y2 w lp ls 1



# #          Stacked Plot Demo
# #
# # Set top and bottom margins to 0 so that there is no space between plots.
# # Fix left and right margins to make sure that the alignment is perfect.
# # Turn off xtics for all plots except the bottom one.
# # In order to leave room for axis and tic labels underneath, we ask for
# # a 4-plot layout but only use the top 3 slots.
# #
# set tmargin 0
# set bmargin 0
# set lmargin 3
# set rmargin 3
# unset xtics
# unset ytics

# set multiplot layout 4,1 title "Auto-layout of stacked plots\n"

# set key autotitle column nobox samplen 1
# unset title
# set style data boxes
# set yrange [0 : 800000]

# set yrange [0:5]
# set ylabel "ISI (ms)" font "Helvetica,18"
# plot  'psth_1e+01.dat' i 1 u (($1*10)+5):3 w lp ls 3, \
#  'psth_1e+01.dat' i 1 u (($1*10)+5):4 w lp ls 2
# set ylabel "CV " font "Helvetica,18"
# set yrange [0.1:0.5]
# #set ytics nomirror out 0.1,0.1,0.5
# plot 'psth_1e+01.dat' i 1 u (($1*10)+5):2 w lp ls 1

# set xtics nomirror
# set tics scale 0
# set xlabel "Time (ms)"
# set tmargin 0
# set bmargin 0
# set lmargin 3
# set rmargin 3
# unset xtics
# unset ytics

# set multiplot layout 4,1 title "Auto-layout of stacked plots\n"

# set key autotitle column nobox samplen 1
# unset title
# set style data boxes
# set yrange [0 : 800000]

# set yrange [0:5]
# set ylabel "ISI (ms)" font "Helvetica,18"
# plot  'psth_1e+01.dat' i 1 u (($1*10)+5):3 w lp ls 3, \
#  'psth_1e+01.dat' i 1 u (($1*10)+5):4 w lp ls 2
# set ylabel "CV " font "Helvetica,18"
# set yrange [0.1:0.5]
# #set ytics nomirror out 0.1,0.1,0.5
# plot 'psth_1e+01.dat' i 1 u (($1*10)+5):2 w lp ls 1

# set xtics nomirror
# set tics scale 0
# set xlabel "Time (ms)"

# set style fill solid 1.0
# plot  'psth_1.dat'  i 0 u 1:2 w boxes lc "black"
# unset yrange
# set autoscale  y
# set style fill solid 1.0
# plot  'psth_1.dat'  i 0 u 1:2 w boxes lc "black"
# unset multiplot
