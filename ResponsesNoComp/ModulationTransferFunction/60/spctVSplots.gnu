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
 splot  '50/spctVS.0.dat' u 2:1:3
 splot  '100/spctVS.0.dat' u 2:1:3
 splot  '150/spctVS.0.dat' u 2:1:3
 splot  '200/spctVS.0.dat' u 2:1:3
 splot  '250/spctVS.0.dat' u 2:1:3
 splot  '300/spctVS.0.dat' u 2:1:3
 splot  '350/spctVS.0.dat' u 2:1:3
 splot  '400/spctVS.0.dat' u 2:1:3
 splot  '450/spctVS.0.dat' u 2:1:3
 splot  '500/spctVS.0.dat' u 2:1:3
 splot  '550/spctVS.0.dat' u 2:1:3
 splot  '600/spctVS.0.dat' u 2:1:3
 splot  '650/spctVS.0.dat' u 2:1:3
 splot  '700/spctVS.0.dat' u 2:1:3
 splot  '750/spctVS.0.dat' u 2:1:3
 splot  '800/spctVS.0.dat' u 2:1:3
 unset multiplot
