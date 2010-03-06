#!/bin/bash

ls [0-9]*| grep :|sed 's/://'| sort -n> freq.dat

touch TS_rate.dat
{
for i in `ls  */rateplace.0.dat|sort -n`
do 
	grep 5810 $i | awk '{print $3}'
	awk '{print $3}'$i > /tmp/rate.dat
	paste TS_rate.dat /tmp/rate.dat > TS_rate.dat
done
} > tstellate.dat
paste freq.dat tstellate.dat > TS_rate_5810.dat
rm -f tstellate.dat
gnuplot TS_rate.gnu


touch DS_rate.dat
{
for i in `ls  */rateplace.2.dat|sort -n`
do 
	grep 5810 $i | awk '{print $3}'
	awk '{print $3}'$i > /tmp/rate.dat
	paste DS_rate.dat /tmp/rate.dat > DS_rate.dat
done
} > dstellate.dat
paste freq.dat dstellate.dat > DS_rate_5810.dat
rm -f dstellate.dat
gnuplot DS_rate.gnu

touch TV_rate.dat
{
for i in `ls */rateplace.1.dat|sort -n`; 
do 
    grep 5810 $i | awk '{print $3}'; done
	awk '{print $3}'$i > /tmp/rate.dat
	paste TV_rate.dat /tmp/rate.dat > TV_rate.dat
}> tuberculo.dat
paste freq.dat tuberculo.dat > TV_rate_5810.dat
rm -f tuberculo.dat
gnuplot TV_rate.gnu

{
for i in `ls */rateplace.3.dat|sort -n`
do 
	grep 5810 $i | awk '{print $3}'
	awk '{print $3}'$i > /tmp/rate.dat
	paste G_rate.dat /tmp/rate.dat > G_rate.dat
done
} > golgi.dat
paste freq.dat golgi.dat > G_rate_5810.dat
rm -f golgi.dat
gnuplot G_rate.gnu



{
for i in `ls  */periodhist.0.dat|sort -n`
do 
        fm=`echo $i| sed 's/\(.*\)\/periodhist.dat/\1/'`
	grep '^50' $i | awk '{print $2,'$fm',$3}'
done
} > periodhist5810.0.dat












# gnuplot TS_response_area.gnu

display TS_rate.eps &
display DS_rate.eps &
display TV_rate.eps &
display G_rate.eps &

rm -f raster.0.dat
{
for i in `ls */ts_raster.dat|sort -n`
do 
        spl=`echo $i| sed 's/\(.*\)\/ts_raster.dat/\1/'`
#	echo $spl
	grep '^50' $i | awk '{print $4,'$spl'}'
	echo ""
done
}>> raster.0.dat
gnuplot rate_raster.gnu
display rate_raster.0.eps &

rm -f psth.0.dat
{
for i in `ls */psth.0.dat|sort -n`
do 
        spl=`echo $i| sed 's/\(.*\)\/psth.0.dat/\1/'`
#	echo $spl
	grep '^50' $i | awk '{print $2,'$spl',$3}'
	echo ""
done
}>> psth.0.dat
#gnuplot psth.gnu
#display psth.0.eps &

gnuplot rasters.gnu
display rasters.eps &

