#!/bin/bash
# rMTF.sh: Reorganise rate MTF information and plot using gnuplot
#
# Michael Eager

set -eu

ls [0-9]*| grep :|sed 's/://'| sort -n> freq.dat

touch TS_rate.dat
{
for i in $(ls  */rateplace.0.dat|sort -n)
do 
	grep 5810 $i | awk '{print $3}' 
	awk '{print $3}' $i > /tmp/rate.dat
	paste TS_rate.dat /tmp/rate.dat > TS_rate.dat
done
} > tstellate.dat
paste freq.dat tstellate.dat > TS_rate_5810.dat
rm -f tstellate.dat
gnuplot TS_rate.gnu


touch DS_rate.dat
{
for i in $(ls  */rateplace.2.dat|sort -n)
do 
	grep 5810 $i | awk '{print $3}' 
	awk '{print $3}' $i > /tmp/rate.dat
	paste DS_rate.dat /tmp/rate.dat > DS_rate.dat
done
} > dstellate.dat
paste freq.dat dstellate.dat > DS_rate_5810.dat
rm -f dstellate.dat
sed 's/TS/DS/' TS_rate.gnu | gnuplot 

touch TV_rate.dat
{
for i in $(ls */rateplace.1.dat|sort -n); 
do 
    grep 5810 $i | awk '{print $3}' 
	awk '{print $3}' $i > /tmp/rate.dat
	paste TV_rate.dat /tmp/rate.dat > TV_rate.dat
done
}> tuberculo.dat
paste freq.dat tuberculo.dat > TV_rate_5810.dat
rm -f tuberculo.dat
sed 's/TS/TV/' TS_rate.gnu | gnuplot

{
for i in $(ls */rateplace.3.dat|sort -n)
do 
	grep 5810 $i | awk '{print $3}' 
	awk '{print $3}' $i > /tmp/rate.dat
	paste G_rate.dat /tmp/rate.dat > G_rate.dat
done
} > golgi.dat
paste freq.dat golgi.dat > G_rate_5810.dat
rm -f golgi.dat
sed 's/TS/G/' TS_rate.gnu | gnuplot


display TS_rateMTF.eps &
display DS_rateMTF.eps &
display TV_rateMTF.eps &
display G_rateMTF.eps &


## RASTER
rm -f raster.0.dat
{
for i in $(ls */ts_raster.dat|sort -n)
do 
        spl=$(echo $i| sed 's/\(.*\)\/ts_raster.dat/\1/')
#	echo $spl
	grep '^50' $i | awk '{print $4,'$spl'}'
	echo ""
done
}>> raster.0.dat
#gnuplot rate_raster.gnu
#display rate_raster.0.eps &

##PSTHs
rm -f psth.0.dat
{
for i in $(ls */psth.0.dat|sort -n)
do 
        spl=$(echo $i| sed 's/\(.*\)\/psth.0.dat/\1/')
#	echo $spl
	grep '^50' $i | awk '{print $2,'$spl',$3}'
	echo ""
done
}>> psth.0.dat
gnuplot psth.gnu
display psth.0.eps &



{
for i in $(ls  */periodhist.0.dat|sort -n)
do 
        fm=$(echo $i| sed 's/\(.*\)\/periodhist.0.dat/\1/')
	grep '^50' $i | awk '{print $2,'$fm',$3}'
	echo ""
done
} > periodhist5810.0.dat

gnuplot ../periodhist5810.gnu
display periodhist5810.0.eps &



#gnuplot rasters.gnu
#display rasters.eps &

