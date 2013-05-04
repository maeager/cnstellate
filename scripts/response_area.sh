#!/bin/bash


rm -f response_area.[0-3].dat

for cell in $(seq 0 1 3); do
{
for i in `ls */rateplace.$cell.dat|sort -n`
do 
        spl=`echo $i| sed 's/\(.*\)\/.*/\1/'`
#	echo $spl
	grep '^[0-9]' $i | awk '{print '$spl',$1,$2,$3}'
	echo ""
done
}>> response_area.$cell.dat
done

# gnuplot response_area.gnu

# rm -f response_area.1.dat
# {
# for i in `ls */rateplace.1.dat|sort -n`
# do 
#         spl=`echo $i| sed 's/\(.*\)\/rateplace.1.dat/\1/'`
# #	echo $spl
# 	grep '^[0-9]' $i | awk '{print '$spl',$1,$2,$3}' 
# 	echo ""
# done
# }>> response_area.1.dat
# #sed 's/response_area.0/response_area.1/' response_area.gnu | gnuplot
# #display response_area-1.eps &

# rm -f response_area.2.dat
# {
# for i in `ls */rateplace.2.dat|sort -n`
# do 
#         spl=`echo $i| sed 's/\(.*\)\/rateplace.2.dat/\1/'`
# #	echo $spl
# 	grep '^[0-9]' $i | awk '{print '$spl',$2,$3}' 
# 	echo ""
# done
# }>> response_area.2.dat
# #sed 's/response_area.0/response_area.2/' response_area.gnu | gnuplot

# rm -f response_area.3.dat
# {
# for i in `ls */rateplace.3.dat|sort -n`
# do 
#         spl=`echo $i| sed 's/\(.*\)\/rateplace.3.dat/\1/'`
# #	echo $spl
# 	grep '^[0-9]' $i | awk '{print '$spl',$2,$3}' 
# 	echo ""
# done
# }>> response_area.3.dat
#sed 's/response_area.0/response_area.3/' response_area.gnu | gnuplot

#display response_area-0.eps &
#display response_area-1.eps &
#display response_area-2.eps &
#display response_area-3.eps &

#gnuplot response_area_redefined.gnu
#display response_area_redefined.eps &
