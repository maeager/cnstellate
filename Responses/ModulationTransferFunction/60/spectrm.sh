#!/bin/bash


# period histograms of each stimulus and its specrogram
for cell in $(seq 0 1 3);
do
    {
	for fm in $(seq 50 50 800);do 
	    grep '^50' $fm/spct.$cell.dat | awk '{print $2,'$fm',$3}'
	    echo ""
	done
    } > spct5810.$cell.dat
   sed -e 's/CELL/'$cell'/' ../spct5810.gnu | gnuplot

    for fm in $(seq 50 50 800);
    do
	sed -e 's/MOD/'$fm'/' -e 's/CELL/'$cell'/' ../spctsingle.gnu | gnuplot
    done
done
