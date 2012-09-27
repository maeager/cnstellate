#!/bin/bash
# temporal.sh: Reorganise temporal information and plot using gnuplot
#
# Michael Eager

set -eu

ls [0-9]*| grep :|sed 's/://'| sort -n> freq.dat

> TS_vs.dat
{
for i in `ls  */vs.0.dat|sort -n`
do 
	grep '^50' $i | awk '{print $2}'
	awk '{print $2}' $i > /tmp/vs.dat
	paste TS_vs.dat /tmp/vs.dat > TS_vs.dat
done
} > tstellate.dat
paste freq.dat tstellate.dat > TS_vs_5810.dat
rm -f tstellate.dat
gnuplot TS_temporal.gnu


> DS_vs.dat
{
for i in `ls */vs.2.dat|sort -n`
do 
	grep '^50' $i | awk '{print $2}'
	awk '{print $2}' $i > /tmp/vs.dat
	paste DS_vs.dat /tmp/vs.dat > DS_vs.dat
done
} > dstellate.dat
paste freq.dat dstellate.dat > DS_vs_5810.dat
rm -f dstellate.dat
sed 's/TS/DS/' TS_temporal.gnu | gnuplot

> TV_vs.dat
{
for i in `ls */vs.1.dat|sort -n`; 
do 
    grep '^50' $i | awk '{print $2}'; done
	awk '{print $2}'$i > /tmp/vs.dat
	paste TV_vs.dat /tmp/vs.dat > TV_vs.dat
}> tuberculo.dat
paste freq.dat tuberculo.dat > TV_vs_5810.dat
rm -f tuberculo.dat
sed 's/TS/TV/' TS_temporal.gnu | gnuplot

> G_vs.dat
{
for i in `ls */vs.3.dat|sort -n`
do 
	grep '^50' $i | awk '{print $2}'
	awk '{print $2}' $i > /tmp/vs.dat
	paste G_vs.dat /tmp/vs.dat > G_vs.dat
done
} > golgi.dat
paste freq.dat golgi.dat > G_vs_5810.dat
rm -f golgi.dat
sed 's/TS/G/' TS_temporal.gnu | gnuplot



display TS_tMTF.eps &
display DS_tMTF.eps &
display TV_tMTF.eps &
display G_tMTF.eps &


#rm -f *.dat
