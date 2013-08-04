set tmargin 0
set bmargin 0
set lmargin 3
set rmargin 3
unset xtics
unset ytics

#set tic scale 0
#unset xtics
#unset ytics
#unset x2tics
#unset y2tics
set pm3d map
unset title
unset ylabel
unset xlabel
unset key
unset border
#unset cbtics
#set colorbox user size .03, .6 noborder
unset colorbox
set multiplot layout 5,3
# splot [0:275][20:80] '150/psth.0.dat' u 2:1:3
 splot  '50/periodhist.0.dat' u 2:1:3
 splot  '100/periodhist.0.dat' u 2:1:3
 splot  '150/periodhist.0.dat' u 2:1:3
 splot  '200/periodhist.0.dat' u 2:1:3
 splot  '250/periodhist.0.dat' u 2:1:3
 splot  '300/periodhist.0.dat' u 2:1:3
 splot  '350/periodhist.0.dat' u 2:1:3
 splot  '400/periodhist.0.dat' u 2:1:3
 splot  '450/periodhist.0.dat' u 2:1:3
 splot  '500/periodhist.0.dat' u 2:1:3
 splot  '550/periodhist.0.dat' u 2:1:3
 splot  '600/periodhist.0.dat' u 2:1:3
 splot  '650/periodhist.0.dat' u 2:1:3
 splot  '700/periodhist.0.dat' u 2:1:3
 splot  '750/periodhist.0.dat' u 2:1:3
 splot  '800/periodhist.0.dat' u 2:1:3
 unset multiplot
