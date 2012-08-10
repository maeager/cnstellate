#!/bin/bash

ls [0-9]*| grep :|sed 's/://'| sort -n> freq.dat

{
for i in `ls -rt */rateplace.0.dat`
do 
	grep 5810 $i | awk '{print $3}'
done
} > tstellate.dat
paste freq.dat tstellate.dat > TS_masked_5810.dat
rm -f tstellate.dat
gnuplot TS_masked.gnu

{
for i in `ls -rt */rateplace.2.dat`
do 
	grep 5810 $i | awk '{print $3}'
done
} > dstellate.dat
paste freq.dat dstellate.dat > DS_masked_5810.dat
rm -f dstellate.dat
gnuplot DS_masked.gnu

{
for i in `ls -rt */rateplace.1.dat`; do grep 5810 $i | awk '{print $3}'; done
}> tuberculo.dat
paste freq.dat tuberculo.dat > TV_masked_5810.dat
rm -f tuberculo.dat
gnuplot TV_masked.gnu

{
for i in `ls -rt */rateplace.3.dat`
do 
	grep 5810 $i | awk '{print $3}'
done
} > golgi.dat
paste freq.dat golgi.dat > G_masked_5810.dat
rm -f golgi.dat
gnuplot G_masked.gnu

# gnuplot TS_response_area.gnu

display TS_masked.eps &
display DS_masked.eps &
display TV_masked.eps &
display G_masked.eps &

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
gnuplot masked_raster.gnu
display masked_raster.0.eps &

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


