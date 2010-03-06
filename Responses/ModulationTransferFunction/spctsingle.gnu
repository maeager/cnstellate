set terminal postscript eps enhanced defaultplex \
   leveldefault color colortext\
   solid dashlength 1.0 linewidth 2.0 butt noclip \
   palfuncparam 2000,0.003 \
   "Helvetica" 18 
set pm3d map
set output ''spct.MOD.CELL.eps"
splot "MOD/spct.CELL.dat" using 2:(log10($1)):(log10($3))