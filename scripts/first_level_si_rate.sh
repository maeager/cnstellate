#!/bin/bash
# first_level_si_rate.sh: Organise SI and Rate data at first level 
#
#  - (C) Michael Eager (mick.eager@gmail.com)

set -eu

ResponsesScriptsPath=$(dirname $0)
ResponsesMPath=${ResponsesScriptsPath}/../mfiles

${ResponsesScriptsPath}/response_area.sh; 
${ResponsesScriptsPath}/vsspikes.sh; 
nice octave -f -q  --eval "tic;addpath('"$ResponsesMPath"');make_an_response_area(200,40000);toc;" 
   

