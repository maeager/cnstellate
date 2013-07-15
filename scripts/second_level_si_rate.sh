#!/bin/bash

ResponsesScriptsPath=$(dirname $0)

for i in $(find -maxdepth 1 -type d | grep -v -e '^.$' | tr -d './' |sort -n)
do
   (\
   cd $i; \
   ${ResponsesScriptsPath}/first_level_si_rate.sh \
   )
done
