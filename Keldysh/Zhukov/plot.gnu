#!/bin/bash

reset

set terminal postscript eps enhanced color font 'Helvetica, 24'
set output 'ValidationOfBicoloreAsymptotics.eps'

bicolor="800x1030nm/DLG800raEE2(1030)=0.dat"
monochrome="Monochrome/DLG800mono.dat"

plot \
bicolor u 1:2 w lp t 'Bicolor asympt. case', \
monochrome u 6:5 w l t 'Monochrome case'
