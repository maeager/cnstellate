# 
# gnuplot is available from http://www.gnuplot.info

set term postscript portrait enhanced mono solid "Helvetica" 8
set output "vowel_baseline.ras.eps"

set size 0.49,0.5
set multiplot
set size 0.48,0.1
set nokey
set noxtics
set noborder
set noytics
set tmargin 0
set rmargin 0

#  D stellate cell
set label 1 '{/Helvetica=14 A}' at screen 0,0.49
set origin 0.01,0.4
plot [0:600][0:100] "vowel_baseline.ds.ras" using 4:3 with dots #, "vowel_baseline.ds.ras" using 4:($3-7000) with dots linetype 1


# T stellate
#set noarrow
#set nolabel 2
#set nolabel 3
set label 1 '{/Helvetica=14 B}' at screen 0,0.39
set origin 0.01,0.3
plot [0:600][0:200] "vowel_baseline.ras" using 4:3 with dots #, "vowel_baseline.ras" using 4:($3-7000) with dots linetype 1

# Raster
#set noarrow
#set nolabel 2
#set noborder
set label 1 '{/Helvetica=14 C}' at screen 0,0.29
set origin 0.01,0.2
plot [0:600][0:2000] "vowel_baseline.an.HSR.ras" using 4:3 with dots #, "vowel_baseline.an.HSR.ras" using 4:($3-500) with dots linetype 1

