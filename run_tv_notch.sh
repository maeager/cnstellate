#!/bin/sh
set -e

./i686/special TV_Notch.hoc -c "NotchRun110()" -c "quit()"
./i686/special TV_Notch.hoc -c "NotchRun100()"  -c "quit()"
./i686/special TV_Notch.hoc -c "NotchRun90()"  -c "quit()"
./i686/special TV_Notch.hoc -c "NotchRun80()" -c "quit()"
./i686/special TV_Notch.hoc -c "NotchRun70()" -c "quit()"
./i686/special TV_Notch.hoc -c "NotchRun50()" -c "quit()"
./i686/special TV_Notch.hoc -c "NotchRunErev70()" -c "quit()"
./i686/special TV_Notch.hoc -c "NotchRunErev75()"  -c "quit()"
./i686/special TV_Notch.hoc -c "NotchRun_negativehsrtv()" -c "quit()"


# in rate_level folders collate cell responses to single file
ls |sort -n | xargs -n1 awk '{split(FILENAME,f,".");printf("%d\t%s\t%g\n",$2,f[1],$3)}' |sed '/^48000/G'> tv.dat
ls |sort -n | xargs -n1 awk '{split(FILENAME,f,".");printf("%d\t%s\t%g\n",$2,f[1],$4)}' |sed '/^48000/G'> ds.dat
ls |sort -n | xargs -n1 awk '{split(FILENAME,f,".");printf("%d\t%s\t%g\n",$2,f[1],$5)}' |sed '/^48000/G'> glg.dat
    ls |sort -n | xargs -n1 awk '{split(FILENAME,f,".");printf("%d\t%s\t%g\n",$2,f[1],$6)}' |sed '/^48000/G'> hsr.dat
    ls |sort -n | xargs -n1 awk '{split(FILENAME,f,".");printf("%d\t%s\t%g\n",$2,f[1],$7)}' |sed '/^48000/G'> lsr.dat
 
