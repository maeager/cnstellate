#!/bin/bash
# siMTF.sh: Reorganise synchronisation index information and plot using gnuplot
#
# Michael Eager

set -eu


ls [0-9]*| grep :|sed 's/://'| sort -n> freq.dat


for cell in `seq 0 1 3`; do
#    awk '! /#/ {print $1}' ./50/vsSPIKES.$cell.dat > vsSPIKES.$cell.dat    
#    awk '! /#/ {print $1}' ./50/vs.$cell.dat > vs.$cell.dat
    > vsSPIKES.$cell.dat;> vs.$cell.dat
    {
    for i in `ls  */vsSPIKES.$cell.dat|sort -n`
    do 
	   grep '^50' $i | awk '{print $2}'
	   awk '! /#/ {print $2}' $i > /tmp/vsSPIKES.dat
           cp vsSPIKES.$cell.dat /tmp/vsSPIKES.$cell.dat
	   paste /tmp/vsSPIKES.$cell.dat /tmp/vsSPIKES.dat > vsSPIKES.$cell.dat
    done 
    for i in `ls  */vs.$cell.dat|sort -n`
    do 
	   awk '! /#/ {print $2}' $i > /tmp/vs.dat
           cp vs.$cell.dat /tmp/vs.$cell.dat
	   paste /tmp/vs.$cell.dat /tmp/vs.dat > vs.$cell.dat
    done
    } > /tmp/vsSPIKES5810.dat
    paste freq.dat /tmp/vsSPIKES5810.dat > vsSPIKES_5810.$cell.dat
    rm -f /tmp/vsSPIKES* /tmp/vs.*
    sed -e 's/CELL/'$cell'/' -e 's/vs_5810/vsSPIKES_5810/' tMTF.gnu | gnuplot
    sed -e 's/CELL/'$cell'/' tMTFmap.gnu | gnuplot
    display tMTF_spikes.$cell.eps &
    display tMTF5810_spikes.$cell.eps &
done

for cell in `seq 0 1 3`; do
    {
    for i in `ls  */spctVS.$cell.dat|sort -n`
    do
        fm=${i%%/*t}
	   awk '! /#/ {x=$3/$2; print $1,'$fm',x}' $i 
        echo ""
    done
    } > spctVS.$cell.dat
    grep '^50' spctVS.$cell.dat | awk '{print $2,$3}' > spctVS_5810.$cell.dat
    sed -e 's/CELL/'$cell'/' \
        -e 's/DATAFILE/spctVS_5810/' \
        -e 's/EPSFILE/siMTF5810/' tMTF.gnu | gnuplot
    sed -e 's/CELL/'$cell'/' \
        -e 's/DATAFILE/spctVS/' \
        -e 's/EPSFILE/siMTF/' tMTFmap.gnu | gnuplot
    display siMTF5180.$cell.eps &
    display siMTF.$cell.eps &
done

# 
# 
# 
# 
# touch TS_vs.dat
# {
# for i in `ls  */vs.0.dat|sort -n`
# do 
# 	grep '^50' $i | awk '{print $2}'
# 	awk '{print $2}' $i > /tmp/vs.dat
# 	paste TS_vs.dat /tmp/vs.dat > TS_vs.dat
# done
# } > tstellate.dat
# paste freq.dat tstellate.dat > TS_vs_5810.dat
# rm -f tstellate.dat
# gnuplot TS_temporal.gnu
# 
# 
# touch DS_vs.dat
# {
# for i in `ls */vs.2.dat|sort -n`
# do 
# 	grep '^50' $i | awk '{print $2}'
# 	awk '{print $2}' $i > /tmp/vs.dat
# 	paste DS_vs.dat /tmp/vs.dat > DS_vs.dat
# done
# } > dstellate.dat
# paste freq.dat dstellate.dat > DS_vs_5810.dat
# rm -f dstellate.dat
# sed 's/TS/DS/' TS_temporal.gnu | gnuplot
# 
# touch TV_vs.dat
# {
# for i in `ls */vs.1.dat|sort -n`; 
# do 
#     grep '^50' $i | awk '{print $2}'; done
# 	awk '{print $2}' $i > /tmp/vs.dat
# 	paste TV_vs.dat /tmp/vs.dat > TV_vs.dat
# }> tuberculo.dat
# paste freq.dat tuberculo.dat > TV_vs_5810.dat
# rm -f tuberculo.dat
# sed 's/TS/TV/' TS_temporal.gnu | gnuplot
# 
# touch G_vs.dat
# {
# for i in `ls */vs.3.dat|sort -n`
# do 
# 	grep '^50' $i | awk '{print $2}'
# 	awk '{print $2}' $i > /tmp/vs.dat
# 	paste G_vs.dat /tmp/vs.dat > G_vs.dat
# done
# } > golgi.dat
# paste freq.dat golgi.dat > G_vs_5810.dat
# rm -f golgi.dat
# sed 's/TS/G/' TS_temporal.gnu | gnuplot
# 
# 
# 
# display TS_tMTF.eps &
# display DS_tMTF.eps &
# display TV_tMTF.eps &
# display G_tMTF.eps &


#rm -f *.dat
