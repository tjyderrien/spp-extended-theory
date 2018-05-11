#!/bin/bash

# Produces the plot submitted to MetaNano conference proceeding 2018 with Prof. Bulgakova and A. Dostovalov. 

set -e
./dostovalov.py &> Dostovalov-PeriodWithThickness.csv
gnuplot plotThicknessPeriod.gnu --persist
