#!/bin/bash

# Produces the plot submitted to MetaNano conference proceeding 2018 with Prof. Bulgakova and A. Dostovalov. 

set -e
infile="PeriodWithThickness.csv"
if [ -e ${infile} ]; then
  echo "Input file ${infile} already exists. Plotting directly. "
else
  ./berini2000.py &> PeriodWithThickness.csv
fi
gnuplot plotThicknessPeriod.gnu --persist
