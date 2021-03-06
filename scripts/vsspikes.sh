#!/bin/bash
# vsspikes.sh: Reorganise VS information 
#
#  - (C) Michael Eager (mick.eager@gmail.com)

set -eu

ls [0-9]*| grep :|sed 's/://'| sort -n> freq.dat
top=$(head -1 freq.dat)
awk '! /#/  {print $1,$2}' ${top}/rateplace.0.dat > cf.dat
[ ! -f response_area.0.dat ] && $(dirname $0)/response_area.sh
awk '{print $1}' response_area.0.dat > fm.dat

ls [0-9]*| grep :|sed 's/://'| sort -n> freq.dat


# X Y Z method with vs and rayleigh
for cell in $(seq 0 1 5); do
    > vsSPIKES.$cell.dat;   
    for i in $(ls  */vsSPIKES.$cell.dat|sort -n)
    do 
	## col 5 now contains pval calculated in NEURON
	   awk '! /#/ {R=$2*$4;z=0;if($4>0){z=(R^2)/$4};pval=exp(sqrt(1+(4*$4)+4*(($4)^2-z*$4))-(1+2*$4));print $2,$3,z,pval,$4,$5}' $i > vsSPIKES.dat
	   paste -d ' ' cf.dat vsSPIKES.dat >> vsSPIKES.$cell.dat
           echo "" >> vsSPIKES.$cell.dat
    done 
    paste -d ' ' fm.dat vsSPIKES.$cell.dat > vsSPIKES.dat
    mv vsSPIKES.dat vsSPIKES.$cell.dat

done
rm -f vsSPIKES.dat cf.dat fm.dat freq.dat



## MATRIX method

# for cell in $(seq 0 1 3); do
#     > vsSPIKES.$cell.dat;
#     for i in $(ls  */vsSPIKES.$cell.dat|sort -n)
#     do 
# 	   awk '! /#/ {print $2,$3}' $i > vsSPIKES.dat
#            cp vsSPIKES.$cell.dat vsSPIKES.$cell.dat
# 	   paste vsSPIKES.$cell.dat vsSPIKES.dat > vsSPIKES.$cell.dat
#     done 
#     rm -f vsSPIKES*    
# done
