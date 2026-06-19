#!gnuplot

reset

set terminal postscript eps enhanced color font 'Helvetica, 24'
set output 'Wavelength-SiOxOnSiO2-EffectOfFilmThickness.eps'

set xlabel 'Wavelength (um)'
set ylabel 'R, T'

set key center right

p  "< awk '{ print }' 01-Substrates/SiOx-40nm-SiO2Substrate-WavelengthRTA.txt1.csv" u 1:2 w l lt 1 lc rgbcolor "red"  t 'SiOx(40 nm)|SiO_2, R', \
   "< awk '{ print }' 01-Substrates/SiOx-44nm-SiO2Substrate-WavelengthRTA.txt1.csv" u 1:2 w l lt 2 lc rgbcolor "red" t 'SiOx(44 nm)|SiO_2, R', \
   "< awk '{ print }' 01-Substrates/SiOx-40nm-SiO2Substrate-WavelengthRTA.txt2.csv" u 1:2 w l lt 1 lc rgbcolor "blue"  t 'SiOx(40 nm)|SiO_2, T', \
   "< awk '{ print }' 01-Substrates/SiOx-44nm-SiO2Substrate-WavelengthRTA.txt2.csv" u 1:2 w l lt 2 lc rgbcolor "blue" t 'SiOx(44 nm)|SiO_2, T'

# "< awk '{ print }' 01-Substrates/SiOx-40nm-SiSubstrate-WavelengthRTA.txt1.csv" u 1:2 w l lt 1 lc rgbcolor "red"  t 'SiOx(40 nm)|Si, R', \
#  "< awk '{ print }' 01-Substrates/SiOx-44nm-SiSubstrate-WavelengthRTA.txt1.csv" u 1:2 w l lt 2 lc rgbcolor "red" t 'SiOx(44 nm)|Si, R', \
#  "< awk '{ print }' 01-Substrates/SiOx-40nm-SiSubstrate-WavelengthRTA.txt2.csv" u 1:2 w l lt 1 lc rgbcolor "blue"  t 'SiOx(40 nm)|Si, T', \
#  "< awk '{ print }' 01-Substrates/SiOx-44nm-SiSubstrate-WavelengthRTA.txt2.csv" u 1:2 w l lt 2 lc rgbcolor "blue" t 'SiOx(44 nm)|Si, T', \



