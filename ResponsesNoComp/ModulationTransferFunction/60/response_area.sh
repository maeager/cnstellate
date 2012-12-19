#!/bin/bash
# response_area.sh: Reorganise data into 2D files and plot using gnuplot
#
# Michael Eager

set -eu


rm -f response_area.0.dat
{
for i in $(ls */rateplace.0.dat|sort -n)
do 
        modfreq=$(echo $i| sed 's/\(.*\)\/rateplace.0.dat/\1/')
#	echo $modfreq
	grep '^[0-9]' $i | awk '{print '$modfreq',$1,$2,$3}'
	echo ""
done
}>> response_area.0.dat
gnuplot response_area.gnu &

rm -f response_area.1.dat
{
for i in $(ls */rateplace.1.dat|sort -n)
do 
        modfreq=$(echo $i| sed 's/\(.*\)\/rateplace.1.dat/\1/')
#	echo $modfreq
	grep '^[0-9]' $i | awk '{print '$modfreq',$1,$2,$3}' 
	echo ""
done
}>> response_area.1.dat
sed 's/response_area.0/response_area.1/' response_area.gnu | gnuplot  &

rm -f response_area.2.dat
{
for i in $(ls */rateplace.2.dat|sort -n)
do 
        modfreq=$(echo $i| sed 's/\(.*\)\/rateplace.2.dat/\1/')
#	echo $modfreq
	grep '^[0-9]' $i | awk '{print '$modfreq',$1,$2,$3}' 
	echo ""
done
}>> response_area.2.dat
sed 's/response_area.0/response_area.2/' response_area.gnu | gnuplot  &

rm -f response_area.3.dat
{
for i in $(ls */rateplace.3.dat|sort -n)
do 
        modfreq=$(echo $i| sed 's/\(.*\)\/rateplace.3.dat/\1/')
#	echo $modfreq
	grep '^[0-9]' $i | awk '{print '$modfreq',$1,$2,$3}' 
	echo ""
done
}>> response_area.3.dat
sed 's/response_area.0/response_area.3/' response_area.gnu | gnuplot &



#display response_area.0.eps &
#display  response_area.1.eps &
#display  response_area.2.eps &
#display  response_area.3.eps &

#gnuplot response_area_redefined.gnu
#display  response_area_redefined.eps &



#  Rate and temporal plots


sed 's/\.2\./.0./' ratetemporal.gnu |gnuplot -p
sed 's/\.2\./.1./' ratetemporal.gnu |gnuplot -p
sed 's/\.2\./.2./' ratetemporal.gnu |gnuplot -p
sed 's/\.2\./.3./' ratetemporal.gnu |gnuplot -p
sed 's/\.2\./.4./' ratetemporal.gnu |gnuplot -p

sed 's/\.0\./.0./' periodhist-multiplot.gnu |gnuplot -p
sed 's/\.0\./.1./' periodhist-multiplot.gnu |gnuplot -p
sed 's/\.0\./.2./' periodhist-multiplot.gnu |gnuplot -p
sed 's/\.0\./.3./' periodhist-multiplot.gnu |gnuplot -p
