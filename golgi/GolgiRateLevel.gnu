# Golgi Rate level
#
# gnuplot is available from http://www.gnuplot.info

set term postscript enhanced mono solid "Helvetica" 12



# set size 0.48,0.2
# set nokey

set border 3
set ytics nomirror
set xtics nomirror
unset y2tics
unset x2tics
set tmargin 0
set rmargin 0
set size noratio
set xlabel "Tone Loudness (SPL dB)"
set ylabel " Discharge Rate (sp/s)"
set xrange [-20:100]
set yrange [0:120]
#  
# set label 1 '{/Helvetica=14 A}' at screen 0,0.49
# set origin 0.01,0.4
set output "GolgiRateLevel.eps"
plot "GolgiRateLevelFit.dat" using 1:2 title "Ghoshal and Kim (1996) S03-07" smooth csplines with lines linewidth 1.5, "GolgiRateLevelFit.dat" using 1:2 title "" with points linewidth 3, "GolgiRateLevel.dat" using 1:2 title "Fitted Output" smooth csplines with lines,  "GolgiRateLevel.dat" using 1:2 title "" with points  linewidth 3

# plot "GolgiRateLevelFit.dat"  using 1:2 title "" smooth csplines with lines linewidth 1, \
#  "GolgiRateLevelFit.dat" using 1:2 title "" with points pointtype 1 linewidth 3, \
#  "GolgiRateLevel.dat" using 1:2 title "" smooth csplines with lines linewidth 1, \
#  "GolgiRateLevel.dat" using 1:2 title "" with points pointtype 7, \
#  "GolgiRateLevelFitOutput.dat" using 1:2 title "" with points pointtype 7 pointsize 3, \
#  "GolgiRateLevelFitOutput.dat" using 1:3 title "" with points pointtype 1 pointsize 3 linewidth 5

# "GolgiRateLevelFitOutput.dat" using 1:2 title "Fitted Output" smooth csplines with lines linewidth 2, \
# "GolgiRateLevelFitOutput.dat" using 1:3 title "Ghoshal and Kim (1996) S03-07" smooth csplines with lines linewidth 2,  \

set output "GolgiRateLevelActualFit.eps"
# plot "GolgiRateLevelFitOutput.dat" using 1:2 title "Fitted Output" smooth csplines with lines linewidth 2, "GolgiRateLevelFitOutput.dat" using 1:2 title "" with points linewidth 5, "GolgiRateLevelFitOutput.dat" using 1:3 title "Ghoshal and Kim (1996) S03-07" smooth csplines with lines linewidth 2,  "GolgiRateLevelFitOutput.dat" using 1:3 title "" with points  linewidth 5
plot  "GolgiRateLevelFitOutput.dat" using 1:3 title "Ghoshal and Kim (1996) S03-07" with linespoints pointtype 7 linewidth 2, "GolgiRateLevelFitOutput.dat" using 1:2 title "Fitted Output" with linespoints pointtype 5 linewidth 2

# set output "GolgiRateLevelActualFit2.eps"
# plot "GolgiRateLevelFitOutput.dat" using 1:2 title "Fitted Output" smooth csplines with lines linewidth 2, "GolgiRateLevelFitOutput.dat" using 1:2 title "" with points pointtype 7 pointsize 3, "GolgiRateLevelFitOutput.dat" using 1:3 title "Ghoshal and Kim (1996) S03-07" smooth csplines with lines linewidth 2,  "GolgiRateLevelFitOutput.dat" using 1:3 title "" with points  pointtype 1 pointsize 5 linewidth 3

set output "GolgiRateLevel2.eps"
plot "GolgiRateLevel_Opt.dat" using 1:2 title "Exp Data" with linespoints lw 2, "GolgiRateLevel_Opt.dat" using 1:3 title "Fitted Data" with linespoints lw 2, "GolgiRateLevel_Opt.dat" using 1:4 title "LSR" with linespoints, "GolgiRateLevel2.dat" using 1:2 title "Exp Data with spikes" with linespoints lw 2
