#!gnuplot

reset

set terminal postscript eps enhanced color font 'Helvetica, 24'
set output 'Wavelength-AuNP.eps'

set xlabel 'Wavelength (um)'
set ylabel 'R, T'

set key top left

p "< awk '{ print }' 01-Substrates/SiOx-40nm-SiSubstrate-WavelengthRTA.txt1.csv" u 1:2 w l lt 1 lc rgbcolor "red"  t 'SiOx(40 nm)|Si, R', \
  "< awk '{ print }' 01-Substrates/SiOx-44nm-SiSubstrate-WavelengthRTA.txt1.csv" u 1:2 w l lt 2 lc rgbcolor "red" t 'SiOx(44 nm)|Si, R', \
  "< awk '{ print }' 01-Substrates/SiOx-40nm-SiSubstrate-WavelengthRTA.txt2.csv" u 1:2 w l lt 1 lc rgbcolor "blue"  t 'SiOx(40 nm)|Si, T', \
  "< awk '{ print }' 01-Substrates/SiOx-44nm-SiSubstrate-WavelengthRTA.txt2.csv" u 1:2 w l lt 2 lc rgbcolor "blue" t 'SiOx(44 nm)|Si, T'



