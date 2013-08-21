#!/bin/bash

ResponsesScriptsPath=$(dirname $0)
ResponsesMPath=${ResponsesScriptsPath}/../mfiles

${ResponsesScriptsPath}/response_area.sh;
nice octave -f -q  --eval "tic;addpath('"$ResponsesMPath"');make_an_response_area(200,40000);toc;" 
