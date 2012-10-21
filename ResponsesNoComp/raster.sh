#!/bin/bash

set -eu

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
gnuplot ratelevel_raster.gnu

rm -f raster.1.dat
{
for i in `ls */tv_raster.dat|sort -n`
do 
        spl=`echo $i| sed 's/\(.*\)\/tv_raster.dat/\1/'`
#	echo $spl
	grep '^50' $i | awk '{print $4,'$spl'}'
	echo ""
done
}>> raster.1.dat
gnuplot ratelevel_raster.gnu

rm -f raster.2.dat
{
for i in `ls */ds_raster.dat|sort -n`
do 
        spl=`echo $i| sed 's/\(.*\)\/ds_raster.dat/\1/'`
#	echo $spl
	grep '^50' $i | awk '{print $4,'$spl'}'
	echo ""
done
}>> raster.2.dat
gnuplot ratelevel_raster.gnu

rm -f raster.3.dat
{
for i in `ls */glg_raster.dat|sort -n`
do 
        spl=`echo $i| sed 's/\(.*\)\/glg_raster.dat/\1/'`
#	echo $spl
	grep '^50' $i | awk '{print $4,'$spl'}'
	echo ""
done
}>> raster.3.dat
gnuplot ratelevel_raster.gnu

rm -f raster.4.dat
{
for i in `ls */anHSR_raster.dat|sort -n`
do 
        spl=`echo $i| sed 's/\(.*\)\/anHSR_raster.dat/\1/'`
#	echo $spl
	grep '^50' $i | awk '{print $4,'$spl'}'
	echo ""
done
}>> raster.4.dat
gnuplot ratelevel_raster.gnu


rm -f raster.5.dat
{
for i in `ls */anLSR_raster.dat|sort -n`
do 
        spl=`echo $i| sed 's/\(.*\)\/anLSR_raster.dat/\1/'`
#	echo $spl
	grep '^50' $i | awk '{print $4,'$spl'}'
	echo ""
done
}>> raster.5.dat
gnuplot ratelevel_raster.gnu


gnuplot rasters.gnu
display rasters.eps &
