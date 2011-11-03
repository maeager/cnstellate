#!/usr/bin/env gnuplot -p
#
# Plotting auditory nerve fibre rates in pseudo 3D
#
# AUTHOR: Michael Eager

reset

# wxt
set terminal wxt size 350,262 enhanced font 'Verdana,10' persist
# png
#set terminal pngcairo size 350,262 enhanced font 'Verdana,10'
#set output 'pseudo3d.png'
# svg
#set terminal svg size 350,262 fname 'Verdana, Helvetica, Arial, sans-serif' \
#fsize '10'
#set output 'pseudo3d.svg'

# color definitions
set style line 1 lc rgb '#ffffff' lt 1 lw 1.5 # --- white
# remove border on top and right and set color to gray
set style line 11 lc rgb '#808080' lt 1
set border 3 front ls 11
set tics nomirror out scale 0.75

unset key
set xrange [0:1300]
set yrange [-2:49]
set xtics 300
set ytics 0,5,30

set xlabel 'Time (ms)'
set ylabel 'Centre Frequency (ERB)'
set style fill solid 1.0 border rgb 'black'

erb(f) = 11.17268 * log(1+(46.06538*f)/(f+14678.49))
f(x) = x>0 ? erb(x) : -erb(-x)
plot for [ii=:1:-1] 'bmm.txt' u (f(column(ii))+ii) w filledcu y1=-2 ls 1


set pm3d map
splot '<  sed -ne "9~6p" /media/data/sounds/zbmodelv4_0/Cat/Click_Recovery_1000_2600_8000_8200_15000_15800_21000_21400_28000_28300.dat' matrix u ($1/50000):2:3 w pm3d

