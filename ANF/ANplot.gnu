#!/usr/bin/env gnuplot -p
#
# Plotting auditory nerve fibre rates in pseudo 3D
#
# AUTHOR: Michael Eager

reset
load 'cf.gpi'
# eps
set terminal postscript eps size 7,5 enhanced color \
    font 'Helvetica,20'
set output 'ANplot.eps'

# wxt
# set terminal wxt size 350,262 enhanced font 'Verdana,10' persist
# png
#set terminal pngcairo size 350,262 enhanced font 'Verdana,10'
#set output 'pseudo3d.png'
# svg
#set terminal svg size 350,262 fname 'Verdana, Helvetica, Arial, sans-serif' \
#fsize '10'
#set output 'pseudo3d.svg'

# # color definitions
# set style line 1 lc rgb '#ffffff' lt 1 lw 1.5 # --- white
# # remove border on top and right and set color to gray
# set style line 11 lc rgb '#808080' lt 1
# set border 3 front ls 11
# set tics nomirror out scale 0.75

# unset key
# set xrange [0:1300]
# set yrange [-2:49]
# set xtics 300
# set ytics 0,5,30

# set xlabel 'Time (ms)'
# set ylabel 'Centre Frequency (ERB)'
# set style fill solid 1.0 border rgb 'black'

# erb(f) = 11.17268 * log(1+(46.06538*f)/(f+14678.49))
# f(x) = x>0 ? erb(x) : -erb(-x)
# plot for [ii=25:1:-1] 'bmm.txt' u (f(column(ii))+ii) w filledcu y1=-2 ls 1


# set pm3d map
# splot '<  sed -ne "9~6p" /media/data/sounds/bruce/Cat/FilteredFile_Notch-sb10-1oct_L50_Del20_Dur49.dat' matrix u ($1/50000):2:3 w pm3d


# #Click_Recovery_1000_2600_8000_8200_15000_15800_21000_21400_28000_28300.dat

# reset


set multiplot # title "HSR auditory nerve fibres\nResponse to Notch Noise"

#
# First plot  (large)
#
set lmargin at screen 0.35
set rmargin at screen 0.85
set bmargin at screen 0.25
set tmargin at screen 0.90

set pm3d
set palette rgbformulae 7, 5, 15
set view map

set xrange [ 0.02 : 0.08 ]
set yrange [ 0 : 99 ]
set zrange [ 0 : 4000 ]

set tics nomirror out

unset key
set xlabel 'Time (ms)'



set palette defined (0 "white",50 "#00ffff", 100 "blue",500 "yellow",2000 "red",5000 "#990000")
set cbrange [0:5000]
set cblabel "Firing Rate (spikes/s)"
# set logscale cb
unset ytics

set parametric

#splot '<  sed -ne "11~6p" /media/data/sounds/bruce/Cat/FilteredFile_Notch-sb6-1oct_L50_Del10_Dur95.dat' matrix u ($1/50000):2:3 w pm3d
splot '<  sed -ne "11~6p" /media/data/sounds/bruce/FilteredFile_Notch-sb10-1oct_L50_Del20_Dur49.dat' matrix u ($1/50000):2:3 w pm3d


unset pm3d
unset key


# Average the rows 
#
# second plot  (tall and narrow; at left of main plot)
#
set lmargin at screen 0.10
set rmargin at screen 0.30
set border 3
set ytics nomirror in 0,10,99
set xtics nomirror out 0,100,400

set xrange [*:*]
set parametric

set ylabel 'Centre Frequency (kHz)'

set ytics axis nomirror norotate ("0.2  0" 0, sprintf("%.2f  25",cf(25)/1000.0) 25,sprintf("%.2f  50",cf(50)/1000.0) 50,sprintf("%.2f  75",cf(75)/1000.0) 75,"48  99" 99)


#plot '<sed -ne "11~6p" /media/data/sounds/bruce/Cat/FilteredFile_Notch-sb6-1oct_L50_Del10_Dur95.dat| awk "{sum=0; for(i=1; i<=NF; i++){sum+=\$i}; sum/=NF; print sum}"' u 1:0 w l 
plot '<sed -ne "11~6p" /media/data/sounds/bruce/FilteredFile_Notch-sb10-1oct_L50_Del20_Dur49.dat| awk "{sum=0; for(i=1; i<=NF; i++){sum+=\$i}; sum/=NF; print sum}"'  u 1:0 t "HSR" w l, \
     '<sed -ne "9~6p" /media/data/sounds/bruce/FilteredFile_Notch-sb10-1oct_L50_Del20_Dur49.dat| awk "{sum=0; for(i=1; i<=NF; i++){sum+=\$i}; sum/=NF; print sum}"' u 1:0 t "LSR" w l lc "black"



unset parametric

#
# third plot  (short and wide; at bottom of main plot)
#
# set lmargin at screen 0.20
# set rmargin at screen 0.85
# set bmargin at screen 0.10
# set tmargin at screen 0.25

# set xrange [ -15.00 : 15.00 ]
# set yrange [ * : * ]
# set xtics
# unset ytics

# y = 0
# plot sin(sqrt(x**2+y**2))/sqrt(x**2+y**2)

unset multiplot
