
set terminal postscript eps enhanced defaultplex \
   leveldefault mono solid\
   dashlength 1.0 linewidth 2.0 butt noclip \
   palfuncparam 2000,0.003 \
   "Helvetica" 18 
#color solid 
set output "PaoliniCV.eps"

set xtics nomirror out ("0-10" 0,"10-20" 1,"20-30" 2,"30-40" 3,"40-50" 4)
set xrange [-0.5:4.5]
set yrange [0.1:0.5]
set border 3
set xlabel "Time Window" font "Helvetica,18"
set ytics nomirror out 0.1,0.1,0.5
set ylabel "Average CV (+-S.E.)" font "Helvetica,18"
set style line 1 lc rgb "#0060ad" lt 1 lw 2 pt 7 ps 1.5 # --beautiful blue line and filled circle
set style line 2 lc rgb '#dd181f' lt 1 lw 2 pt 5 ps 1.5 # --beautiful red line and square
set style line 3 lc rgb '#00ad06' lt 1 lw 2 pt 9 ps 1.5 # --green line and triangle

plot 'PaoliniBalancedInh-Fig2.png.dat' i 0 u 1:2:3 w yerr lc rgb '#0060ad', '' i 0 u 1:2 w lp ls 1, \
     '' i 1 u 1:2:3 w yerr lc rgb '#dd181f', '' i 1 u 1:2 w lp ls 2, \
     '' i 2 u 1:2:3 w yerr lc rgb '#00ad06' ps 0, '' i 2 u 1:2 w lp ls 3
     
