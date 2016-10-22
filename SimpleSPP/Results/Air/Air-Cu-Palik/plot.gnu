#!gnuplot

reset
set output 'SummaryOfAnalysis.eps' 
set terminal postscript eps enhanced color font 'Helvetica, 24'

file="ResultsZhenyaCu.csv"

set xlabel 'Laser wavelength (nm)'
set ylabel 'SPP period (nm)'

set key out above center

set xtics 200
set ytics 200

plot file u 1:3 w p lc 1 lw 5 t 'Simple SPP theory', \
file u 1:4 w p lc 3 lw 5 t  'Roughness coupled SPP model'
