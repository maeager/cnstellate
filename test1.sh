#!/bin/ bash
ARCH=$(arch)
NRNHOME="$HOME/src/neuron/nrnmpi"
NRNBIN="${NRNHOME}/$ARCH/bin/"
NRNIV="${NRNBIN}nrniv"
LNRNMECH="$(pwd)/${ARCH}/.libs/libnrnmech.so"
MPIRUN=$(which mpirun)

if [[ $(ps -u $USER | grep ssh-agent) == "" ]]
then
    nohup ssh-agent -s > ~/.ssh-agent
fi
source ~/.ssh-agent
ssh-add -t 900000

dorun() {
        echo "starting dorun $np"
        "${MPIRUN}" -np $1 "${NRNIV}" -dll "${LNRNMECH}" par_batch1.hoc >& temp
        rm -r -f result.$1
        mkdir result.$1
        cp vowel* result.$1
        sort -k 1n,1n -k 2n,2n out.dat > out.$1
        mv temp stdout.$1
        echo "completed dorun $np"
}



#for np in 25 20 15;
#  dorun $np
#done
dorun 3
