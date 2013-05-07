#!/bin/bash

sed 's/CELL/TS/' RL_NRL.gnu | gnuplot
sed 's/CELL/DS/' RL_NRL.gnu | gnuplot
sed 's/CELL/TV/' RL_NRL.gnu | gnuplot
sed 's/CELL/G/' RL_NRL.gnu | gnuplot
