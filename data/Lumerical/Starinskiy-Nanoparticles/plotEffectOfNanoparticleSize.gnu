#!gnuplot

reset

set terminal postscript eps enhanced color font 'Helvetica, 24'
set output 'Wavelength-EffectOfNPsize.eps'

folder="03-NanoparticlesInFilmDepositedOnSubstrates/"

set xlabel 'Wavelength (um)'
set ylabel 'R, T'

set key top left

p "< awk '{ print }' ".folder."AuNP-10nm-InsideSiOx-44nm-SiSubstrate-WavelengthRTA.txt1.csv" u 1:2 w l lt 1 lc rgbcolor "red"  t 'AuNP(10nm) in SiOx(44nm) | Si, R', \
  "< awk '{ print }' ".folder."AuNP-14nm-InsideSiOx-44nm-SiSubstrate-WavelengthRTA.txt1.csv" u 1:2 w l lt 2 lc rgbcolor "red" t 'AuNP(14nm) in SiOx(44nm) | Si, R', \
  "< awk '{ print }' ".folder."AuNP-10nm-InsideSiOx-44nm-SiSubstrate-WavelengthRTA.txt2.csv" u 1:2 w l lt 1 lc rgbcolor "blue"  t 'AuNP(10nm) in SiOx(44nm) | Si, T', \
  "< awk '{ print }' ".folder."AuNP-14nm-InsideSiOx-44nm-SiSubstrate-WavelengthRTA.txt2.csv" u 1:2 w l lt 2 lc rgbcolor "blue" t 'AuNP(14nm) in SiOx(44nm) | Si, T'


