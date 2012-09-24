# gnuplot file to produce a figure similar to Figure 2
# in Davison A.P., Feng J. and Brown D. (2003) J. Neurophysiol.

# Note that the figure produced will not match exactly the
# published figure due to differences in the sequence of
# random numbers used to set up the network.

# gnuplot is available from http://www.gnuplot.info

set term postscript eps enhanced color 
set output "vowel_image.eps"

set border 1 linewidth 0.5

set pm3d map

plot "vowel_b6t_HSR.sout" using 1:2:3 




