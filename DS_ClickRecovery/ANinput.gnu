set terminal postscript eps enhanced defaultplex \
   leveldefault mono solid\
   dashlength 1.0 linewidth 2.0 butt noclip \
   palfuncparam 2000,0.003 \
   "Helvetica" 18 
#color solid 
set output "ANinput.eps"
unset key
unset xtics
unset border
unset ytics
set tmargin 0
set rmargin 0

set multiplot layout 4,1
set label 1 '{/Helvetica=22 Stimulus }' at screen -0.1,0.9
plot 'Stim.dat' u 1:2 w l
#set label 1 '{/Helvetica=22 AN Model }' at screen -0.1,0.75
#plot 'ANFrate.dat' u 1:2 w l
set label 1 '{/Helvetica=22 HSR }' at screen -0.1,0.65
plot 'ANFhist.dat' u 1:($2-245) w l
#plot "< awk '{print $50}' ANFHist.dat" u 0:1 w l
set label 1 '{/Helvetica=22 GLG }' at screen -0.1,0.4
plot 'GolgiPSTH.dat' u 1:2 w l

set xtics out nomirror
#set ytics in nomirror
#set border 3
set xlabel "Time (ms)" font "Helvetica,18"

set label 1 '{/Helvetica=22 DS }' at screen -0.1,0.2
plot [0:300] 'DSPSTH.dat' u 1:2 w l

unset multiplot

!fixbb ANinput.eps