#!/bin/bash
# Organise rate data into format for pm3d/splot in gnuplot
#
#  - (C) Michael Eager (mick.eager@gmail.com)

set -eu


rm -f response_area.[0-3].dat

for cell in $(seq 0 1 3); do
{
for i in `ls -- */rateplace.$cell.dat|sort -n`
do 
        spl=`echo $i| sed 's/\(.*\)\/.*/\1/'`
#	echo $spl
	grep '^[0-9]' -- $i | awk '{print '$spl',$1,$2,$3}'
	echo ""
done
}>> response_area.$cell.dat
done

