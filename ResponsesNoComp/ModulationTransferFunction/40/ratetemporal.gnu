set terminal postscript eps enhanced defaultplex \
   leveldefault color colortext \
   solid dashlength 1.0 linewidth 2.0 butt noclip \
   palfuncparam 2000,0.003 \
   "Helvetica" 18 size 5,3
#   fontfile "uhvr8a.pfb" 18 

   set output "ratetemporal.2.eps"
#set yrange [0:90]
set xlabel "f_m (Hz)" font "Helvetica,24"
set ylabel "Channel Position" font "Helvetica,24"

set pm3d map
#set logscale x 10
set colorbox noborder
set multiplot layout 1,2
set xtics out ( "100" 100, "" 200, "300" 300, "" 400, "500" 500, "" 600, "700" 700, "" 800)

#set logscale y 10
set cbrange [0:*]
#set palette model RGB
#set palette defined
#set palette defined (0 "blue", 150 "white", 300 "red")
set palette rgbformulae 22,13,-31
splot [50:800][0:99] 'response_area.2.dat' u 1:2:($4*5)
#unset palette
unset ylabel
unset logscale y
set cbrange [0:1]
#set palette model HSV rgbformulae 3,2,2
#set palette model XYZ rgbformulae 7,5,15
#set palette defined ( 0 0 0 0, 1 1 1 1 )
set palette rgbformulae 7,5,15
splot [50:800][0:99] 'vsSPIKES.2.dat' matrix u ($1*50+50):2:3
unset multiplot
!fixbb ratetemporal.2.eps

