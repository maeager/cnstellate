# 
# gnuplot is available from http://www.gnuplot.info
load '../Responses/default.gnu'

set terminal postscript eps enhanced defaultplex \
   leveldefault mono solid\
   dashlength 1.0 linewidth 2.0 butt noclip \
   palfuncparam 2000,0.003 \
   "Helvetica" 18 

set output "baseline.eps"

#set size 0.49,0.5
set multiplot layout 4,2

unset key
unset xtics
unset border
unset ytics
set tmargin 0
set rmargin 0

#  D stellate cell
#set label 1 '{/Helvetica=14 A}' at screen 0,0.49
#set origin 0.01,0.4
set arrow 1 from 5,-65 to 5,-55 nohead lw 1
plot [10:110][-70:-40] "lsr_single_noNa.dat" w l
unset arrow 1
plot [10:110][-70:-40] "lsr_average_noNa.dat" u ($1+20):2 w l

set arrow 1 from 5,-65 to 5,-55 nohead lw 1
plot [10:110][-70:-40] "hsr_single_noNa.dat" w l
unset arrow 1
plot [10:110][-70:-40] "hsr_average_noNa.dat" u ($1+20):2 w l

set arrow 1 from 5,-65 to 5,-55 nohead lw 1
plot [10:110][-70:-40] "ds_single.dat" w l
unset arrow 1
plot [10:110][-70:-40] "ds_average.dat" u ($1+20):2 w l

set arrow 1 from 5,-65 to 5,-55 nohead lw 1
set arrow 2 from 20,-70 to 70,-70 nohead lw 3
plot [10:110][-70:-40] "golgi_single.dat" w l
unset arrow 1

plot [10:110][-70:-40] "golgi_average.dat" u ($1+20):2 w l

unset multiplot

set output "baseline_overlap.eps"
set multiplot layout 2,2
unset key
unset xtics
unset border
unset ytics
set tmargin 0
set rmargin 0
unset arrow 2
#  D stellate cell
#set label 1 '{/Helvetica=14 A}' at screen 0,0.49
#set origin 0.01,0.4
set arrow 1 from 5,-65 to 5,-55 nohead lw 1
plot [10:110][-70:-40] "lsr_single_noNa.dat" w l lw 1, \
 "lsr_average_noNa.dat" u ($1+20):2 w l lw 2

unset arrow 1
plot [10:110][-70:-40] "hsr_single_noNa.dat" w l lw 1, \
 "hsr_average_noNa.dat" u ($1+20):2 w l lw 2

set arrow 2 from 20,-72 to 70,-72 nohead lw 3
set arrow 1 from 5,-70 to 5,-65 nohead lw 1

plot [10:110][-70:-55] "ds_single.dat" w l lw 1,\
 "ds_average.dat" u ($1+20):2 w l lw 2
unset arrow 1
plot [10:110][-70:-55] "golgi_single.dat" w l lw 1, \
 "golgi_average.dat" u ($1+20):2 w l lw 2



unset multiplot

set output "baseline_exc.eps"
set multiplot layout 1,2 scale 1,0.3
unset key
unset xtics
unset border
unset ytics
set tmargin 0
set rmargin 0
unset arrow 2

set label 1 '{/Helvetica=20 A}' at 10,-50
set label 3 '{/Helvetica=14 10 mV}' at 0,-65 right rotate 
set arrow 1 from 5,-65 to 5,-55 nohead lw 1
plot [10:110][-70:-50] "lsr_single_noNa.dat" w l lw 1, \
 "lsr_average_noNa.dat" u ($1+20):2 w l lw 2
set label 1 '{/Helvetica=20 B}' at 10,-50
unset arrow 1
unset label 3
plot [10:110][-70:-50] "hsr_single_noNa.dat" w l lw 1, \
 "hsr_average_noNa.dat" u ($1+20):2 w l lw 2


unset multiplot

set output "baseline_inh.eps"
set multiplot layout 1,2 scale 1,0.3
unset key
unset xtics
unset border
unset ytics
set tmargin 0
set rmargin 0

set arrow 2 from 20,-72 to 70,-72 nohead lw 3
set arrow 1 from 5,-70 to 5,-65 nohead lw 1
set label 1 '{/Helvetica=20 C}' at 10,-60
set label 3 '{/Helvetica=14 5 mV}' at 0,-65 right rotate 
set label 4 '{/Helvetica=14 50 msec}' at 20,-73

plot [10:110][-75:-60] "ds_single.dat" w l lw 1,\
 "ds_average.dat" u ($1+20):2 w l lw 2
unset arrow 1
set label 1 '{/Helvetica=20 D}' at 10,-60
unset label 3
unset label 4
plot [10:110][-75:-60] "golgi_single.dat" w l lw 1, \
 "golgi_average.dat" u ($1+20):2 w l lw 2


