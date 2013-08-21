
#set terminal postscript eps enhanced defaultplex \
#   leveldefault color colortext \
##   solid dashlength 1.0 linewidth 2.0 butt noclip \
#   palfuncparam 2000,0.003 \
#   "Helvetica" 18
#set output "psth-250.0.eps"

set tmargin 0
set bmargin 0
set lmargin 3
set rmargin 3
unset xtics
unset ytics

#set tic scale 0
#unset xtics
#unset yticsls 
#unset x2tics
#unset y2tics
set pm3d map
unset title
unset ylabel
unset xlabel
unset key
unset border
set palette gray
#unset cbtics
#set colorbox user size .03, .6 noborder
unset colorbox
 splot  '250/psth.0.dat' u 2:1:3
#!fixbb psth-250.0.eps
