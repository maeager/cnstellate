#!/bin/bash

ResponsesScriptsPath=${HOME}/Work/cnstellate/scripts
ResponsesMPath=${HOME}/Work/cnstellate/mfiles
# cd  /media/c4bb64a6-7c5f-4dc1-9965-b0f4c1117b36/Work-archive/cnstellate-2.58/TStellate2_CS/ModulationTransferFunction/

   ${ResponsesScriptsPath}/response_area.sh; 
   ${ResponsesScriptsPath}/vsspikes.sh; 
   nice octave -f -q  --eval "tic;addpath('"$ResponsesMPath"');make_an_response_area(200,40000);toc;" 
   

