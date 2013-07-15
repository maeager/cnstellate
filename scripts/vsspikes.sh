#!/bin/bash
# vsspikes.sh: Reorganise VS information 
#
#  - (C) Michael Eager (mick.eager@gmail.com)

set -eu

ls [0-9]*| grep :|sed 's/://'| sort -n> /tmp/freq.dat
top=$(head -1 /tmp/freq.dat)
awk '! /#/  {print $1,$2}' ${top}/rateplace.0.dat > /tmp/cf.dat
[ ! -f response_area.0.dat ] && $(dirname $0)/response_area.sh
awk '{print $1}' response_area.0.dat > /tmp/fm.dat

ls [0-9]*| grep :|sed 's/://'| sort -n> freq.dat


# X Y Z method with vs and rayleigh
for cell in $(seq 0 1 5); do
    > vsSPIKES.$cell.dat;   
    for i in $(ls  */vsSPIKES.$cell.dat|sort -n)
    do 
	   awk '! /#/ {R=$2*$4;z=0;if ($4>0){z=(R^2)/$4};pval=exp(sqrt(1+(4*$4)+4*(($4)^2-z*$4))-(1+2*$4));print $2,$3,z,pval}' $i > /tmp/vsSPIKES.dat
	   paste -d ' ' /tmp/cf.dat /tmp/vsSPIKES.dat >> vsSPIKES.$cell.dat
           echo "" >> vsSPIKES.$cell.dat
    done 
    paste -d ' ' /tmp/fm.dat vsSPIKES.$cell.dat > /tmp/vsSPIKES.dat
    mv /tmp/vsSPIKES.dat vsSPIKES.$cell.dat

done
rm -f /tmp/vsSPIKES.dat /tmp/cf.dat /tmp/fm.dat /tmp/freq.dat



## MATRIX method

# for cell in $(seq 0 1 3); do
#     > vsSPIKES.$cell.dat;
#     for i in $(ls  */vsSPIKES.$cell.dat|sort -n)
#     do 
# 	   awk '! /#/ {print $2,$3}' $i > /tmp/vsSPIKES.dat
#            cp vsSPIKES.$cell.dat /tmp/vsSPIKES.$cell.dat
# 	   paste /tmp/vsSPIKES.$cell.dat /tmp/vsSPIKES.dat > vsSPIKES.$cell.dat
#     done 
#     rm -f /tmp/vsSPIKES*    
# done
