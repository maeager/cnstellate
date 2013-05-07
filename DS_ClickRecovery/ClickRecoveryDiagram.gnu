reset
set terminal postscript eps enhanced size 7,4 defaultplex \
    leveldefault color solid\
    dashlength 1.0 linewidth 2.0 butt noclip \
    palfuncparam 2000,0.003 \
    "Helvetica" 22 

set output "ClickRecoveryDiagram.eps"
unset key
unset xtics
unset border
unset ytics
set tmargin 0
set rmargin 0
set xtics out nomirror
set multiplot layout 2,1
set label 1 '{/Helvetica=32 Stimulus }' at screen -0.05,0.9


plot [0:100] 'Stim.dat' u ($1*1000):2 w l lc rgb "black"

#set ytics in nomirror
#set border 3
set xlabel "Time (ms)" font "Helvetica,30"
set label 1 '{/Helvetica=32 DS Cell}' at screen -0.05,0.35
set label 2 '{/Helvetica=32   PSTH}' at screen -0.05,0.3
set object 1 rect from  13, -10 to  15,40 fs empty border rgb "blue"
set object 2 rect from  29, -10 to  31,40 fs empty border rgb "red"
set object 3 rect from  83, -10 to  85,40 fs empty border rgb "blue"
set object 4 rect from  85, -10 to  87,40 fs empty border rgb "red"

set style fill solid border lc rgb "black"
set arrow from 0,0 to 100,0 nohead 
plot [0:100] 'DSPSTH.dat' u 1:2 w fillsteps lc rgb "black"

unset multiplot

!fixbb ClickRecoveryDiagram.eps
!epstopdf ClickRecoveryDiagram.eps

