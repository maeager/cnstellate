

set term postscript eps enhanced mono
set output "rasters.eps"

set size 0.49,0.5
set multiplot
set size 0.48,0.1
set nokey
set noxtics
set noborder
set noytics
set tmargin 0
set rmargin 0


set origin 0.01,0.4
plot "./90/anHSR_raster.dat" using 4:($3) with dots linetype 1

set origin 0.01,0.3
plot "./90/anLSR_raster.dat" using 4:($3) with dots linetype 1

set origin 0.01,0.2
plot "./90/ds_raster.dat" using 4:($3) with dots linetype 1

set origin 0.01,0.1
plot "./90/tv_raster.dat" using 4:($3) with dots linetype 1

set origin 0.01,0.0
plot  "./90/ts_raster.dat" using 4:($3) with dots linetype 1
