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

set multiplot layout 5,1

plot 'Stim.dat' u 1:2 w l
plot 'ANFrate.dat' u 1:2 w l
plot 'ANFhist.dat' u 1:($2-245) w l
#plot "< awk '{print $50}' ANFHist.dat" u 0:1 w l
plot 'GolgiPSTH.dat' u 1:2 w l

set xtics out nomirror
#set ytics in nomirror
#set border 3
plot [0:300] 'DSPSTH.dat' u 1:2 w l
