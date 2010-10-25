

set term postscript eps enhanced mono
set output "rasters.eps"

set size 0.49,0.5
set multiplot
set size 0.48,0.1
unset key
unset xtics
unset border
unset ytics
set tmargin 0
set rmargin 0


set origin 0.01,0.4
set label 1 "HSR" at screen 0.05,0.44
plot "./5810/anHSR_raster.dat" using 4:($3) with dots linetype 1

set origin 0.01,0.3
set label 1 "LSR" at screen 0.05,0.34
plot "./5810/anLSR_raster.dat" using 4:($3) with dots linetype 1

set origin 0.01,0.2
set label 1 "DS" at screen 0.05,0.24
plot "./5810/ds_raster.dat" using 4:($3) with dots linetype 1

set origin 0.01,0.1
set label 1 "TV" at screen 0.05,0.14
plot "./5810/tv_raster.dat" using 4:($3) with dots linetype 1

set origin 0.01,0.0
set label 1 "TS" at screen 0.05,0.04
plot "./5810/ts_raster.dat" using 4:($3) with dots linetype 1
