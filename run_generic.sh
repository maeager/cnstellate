#!/bin/sh

mpispecial -np $1 par_$2.hoc

mkdir $2
mv $2.* $2/
