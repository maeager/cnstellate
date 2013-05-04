
load '../default.gnu'

set terminal postscript eps enhanced defaultplex \
   leveldefault color \
   solid dashlength 1.0 linewidth 2.0 butt noclip \
   palfuncparam 2000,0.003 \
   "Helvetica" 18 
#set border 3

unset x2tics
unset y2tics
set xlabel "time (ms)" font "Helvetica,28" 
set ylabel "spike count" font "Helvetica,28" offset +2,0
REPS=20

set label 1 "VAR dB" at graph 0.85,0.9 font "Helvetica,32"
plot [0:*][0:10] "<grep '^50' ./VAR/psth.3.dat" using 2:3 \
     w boxes fs solid 1 rgb "black"




