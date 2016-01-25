#!gnuplot

reset

set xlabel 'Mg thickness (m)'
set ylabel 'SPP decay depth (nm)'

set output '20150221-Air-Mg-W.eps'
set terminal postscript eps enhanced color font 'Helvetica, 26'

plot "500nm/Damping.awked.tmp" u 2:3 w l t '500 nm', \
"1000nm/Damping.awked.tmp" u 2:3 w l t '1000 nm'
