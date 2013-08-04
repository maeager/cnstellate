#!/bin/bash

ResponsesScriptsPath=$(dirname $0)
ResponsesMPath=${ResponsesScriptsPath}/../mfiles
# cd  /media/c4bb64a6-7c5f-4dc1-9965-b0f4c1117b36/Work-archive/cnstellate-2.58/TStellate2_CS/ModulationTransferFunction/

for i in $(find -maxdepth 1 -type d | grep -v -e '^.$' | tr -d './' |sort -n)
do
   (cd $i; \
   $ResponsesScriptsPath/first_level_si_rate.sh)
done
