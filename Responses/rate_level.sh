#!/bin/bash

ls [0-9]*| grep :|sed 's/://'| sort -n> level.dat

{
for i in $(ls -rt */rateplace.0.dat)
do 
	grep 5810 $i | awk '{print $3}'
done
} > tstellate.dat
paste level.dat tstellate.dat > TS_ratelevel_5810.dat
rm -f tstellate.dat
gnuplot TS_ratelevel.gnu

{
for i in $(ls -rt */rateplace.2.dat)
do 
	grep 5810 $i | awk '{print $3}'
done
} > dstellate.dat
paste level.dat dstellate.dat > DS_ratelevel_5810.dat
rm -f dstellate.dat
gnuplot DS_ratelevel.gnu

{
for i in $(ls -rt */rateplace.1.dat); do 
    grep 5810 $i | awk '{print $3}'; 
done
}> tuberculo.dat
paste level.dat tuberculo.dat > TV_ratelevel_5810.dat
rm -f tuberculo.dat
gnuplot TV_ratelevel.gnu

{
for i in $(ls -rt */rateplace.3.dat)
do 
	grep 5810 $i | awk '{print $3}'
done;}> golgi.dat
paste level.dat golgi.dat > G_ratelevel_5810.dat
rm -f golgi.dat
gnuplot G_ratelevel.gnu

# gnuplot TS_response_area.gnu

display TS_ratelevel.eps &
display DS_ratelevel.eps &
display TV_ratelevel.eps &
display G_ratelevel.eps &

# rm -f raster.0.dat
# {
# for i in $(ls */ts_raster.dat|sort -n)
# do 
#         spl=$(echo $i| sed 's/\(.*\)\/ts_raster.dat/\1/')
# #	echo $spl
# 	grep '^50' $i | awk '{print $4,'$spl'}'
# 	echo ""
# done
# }>> raster.0.dat
# gnuplot ratelevel_raster.gnu


rm -f psth.0.dat
{ for i in $(ls */psth.0.dat|sort -n); do 
        spl=$(echo $i| sed 's/\(.*\)\/psth.0.dat/\1/')
#	echo $spl
	grep '^50' $i | awk '{print $2,'$spl',$3}'
	echo ""
done; }>> psth.0.dat

rm -f psth.1.dat
{ for i in $(ls */psth.1.dat|sort -n); do 
        spl=$(echo $i| sed 's/\(.*\)\/psth.1.dat/\1/')
#	echo $spl
	grep '^50' $i | awk '{print $2,'$spl',$3}'
	echo ""
done; }>> psth.1.dat

rm -f psth.2.dat
{
for i in $(ls */psth.2.dat|sort -n); do 
        spl=$(echo $i| sed 's/\(.*\)\/psth.2.dat/\1/')
#	echo $spl
	grep '^50' $i | awk '{print $2,'$spl',$3}'
	echo ""
done; }>> psth.2.dat

rm -f psth.3.dat
{for i in $(ls */psth.3.dat|sort -n); do 
        spl=$(echo $i| sed 's/\(.*\)\/psth.3.dat/\1/')
#	echo $spl
	grep '^50' $i | awk '{print $2,'$spl',$3}'
	echo ""
done; }>> psth.3.dat
gnuplot psth.gnu

display psthVlevel.0.eps &
display psthVlevel.1.eps &
display psthVlevel.2.eps &
display psthVlevel.3.eps &

gnuplot rasters.gnu
display rasters.eps &

