#!gnuplot

reset

xscale=1e9
yscale=1e9

set log x 

set terminal postscript eps enhanced color font "Helvetica, 26"
set output '20140721-SPPperiodOfGoldFilmInAirOnSilicaSubstrate.eps'

set xlabel 'Au film thickness (nm)'
set ylabel 'SPP period (nm)'

plot "Result.tmp" u ($2*xscale):($3*xscale) w l t 'Period' 
