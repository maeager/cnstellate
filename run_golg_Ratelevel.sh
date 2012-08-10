#!/bin/sh
set -o errexit


mpispecial -np 1 par_golgi_Ratelevel.hoc

mkdir golgi_Ratelevel
mv golgi_Ratelevel.* golgi_Ratelevel
