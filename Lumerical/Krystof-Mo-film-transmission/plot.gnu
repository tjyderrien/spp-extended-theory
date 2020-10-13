#!gnuplot

reset

set terminal postscript eps enhanced color font 'Helvetica, 24'
set output 'Wavelength-MoVsMoOnSi.eps'

set xlabel 'Wavelength (um)'
set ylabel 'Absorption'

p "< awk '{ print }' Mo-*.txt3.csv"     u 1:2 w l lc rgbcolor "red"  t 'Mo', \
  "< awk '{ print }' MoonSi-*.txt3.csv" u 1:2 w l lc rgbcolor "blue" t 'Mo|Si'



