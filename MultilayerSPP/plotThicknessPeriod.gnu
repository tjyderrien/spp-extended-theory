#!gnuplot

set output 'plotThicknessPeriod.eps'
set terminal postscript eps enhanced color font "Helvetica, 24"
set log x
set format "%g"
set xlabel 'Thickness (m)'
set ylabel 'Period (m)'
plot "Dostovalov-PeriodWithThickness.csv" u 1:2 w p t '{/Symbol L}_{SPP}'

set output
set terminal x11
replot


