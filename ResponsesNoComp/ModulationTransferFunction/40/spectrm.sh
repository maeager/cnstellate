#!/bin/bash
# temporal.sh: Reorganise period histograms and spectrum information and plot using gnuplot
#
# Michael Eager

set -eu


# period histograms of each stimulus and its specrogram
for cell in $(seq 0 1 3);
do
    {
	for fm in $(seq 50 50 800);do 
	    grep '^50' $fm/spct.$cell.dat | awk '{print $2,'$fm',$3}'
	    echo ""
	done
    } > spct5810.$cell.dat
   sed -e 's/CELL/'$cell'/' ../../spct5810.gnu | gnuplot

    for fm in $(seq 50 50 800);
    do
	#sed -e 's/MOD/'$fm'/' -e 's/CELL/'$cell'/' ../../spctsingle.gnu | gnuplot
	sed -e 's/MOD/'$fm'/' -e 's/CELL/'$cell'/' ../../periodhistsingle.gnu | gnuplot
	sed -e 's/MOD/'$fm'/' -e 's/CELL/'$cell'/' ../../psthsingle.gnu | gnuplot
    done
done


#  specrogram VS
for cell in $(seq 0 1 3);
do
    {
	for fm in $(seq 50 50 800);do 
	    awk '{if($2 != 0) {print $1, '$fm', $3/$2 } else {print $1, '$fm', "0"}}' $fm/spctVS.$cell.dat
	    echo ""
	done
    } > spctVS.$cell.dat
#   sed -e 's/CELL/'$cell'/' ../../spctVS.gnu | gnuplot

    {
	for fm in $(seq 50 50 800);do 
	    awk '{if($2 == 0){MAX0=$3}; if($2 == 1){print $1, '$fm', $3/MAX0} }' $fm/spctFULL.$cell.dat
	    echo ""
	done
    } > spctFULL.$cell.dat


done

