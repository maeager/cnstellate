#!/bin/bash

{
for i in `ls -rt */rateplace.0.dat`
do 
	grep 5810 $i | awk '{print $3}'
done
} > tstellate.dat
paste level.dat tstellate.dat > TS_ratelevel_5810.dat
rm -f tstellate.dat
gnuplot TS_ratelevel.gnu

{
for i in `ls -rt */rateplace.2.dat`
do 
	grep 5810 $i | awk '{print $3}'
done
} > dstellate.dat
paste level.dat dstellate.dat > DS_ratelevel_5810.dat
rm -f dstellate.dat
gnuplot DS_ratelevel.gnu

{
for i in `ls -rt */rateplace.1.dat`; do grep 5810 $i | awk '{print $3}'; done
}> tuberculo.dat
paste level.dat tuberculo.dat > TV_ratelevel_5810.dat
rm -f tuberculo.dat
gnuplot TV_ratelevel.gnu

{
for i in `ls -rt */rateplace.3.dat`
do 
	grep 5810 $i | awk '{print $3}'
done
} > golgi.dat
paste level.dat golgi.dat > G_ratelevel_5810.dat
rm -f golgi.dat
gnuplot G_ratelevel.gnu


