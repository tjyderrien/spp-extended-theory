#!gnuplot

reset

set terminal postscript eps enhanced color font 'Helvetica, 24'
set output 'Wavelength-EffectOfBoundaryConditions.eps'

folder="04-RandomAuNPdistribEmbeddedInSiOx44nm-SiSubstrate/"

set xlabel 'Wavelength (um)'
set ylabel 'R, T'

set key top Left

set title 'Au NP random (x16), 14 nm radius in SiOx (44 nm) on Si'

p "< awk '{ print }' ".folder."UsualPeriodicConditions.txt1.csv" u 1:2 w l lt 2 lc rgbcolor "red" t 'Periodic, R', \
  "< awk '{ print }' ".folder."UsualPeriodicConditions.txt2.csv" u 1:2 w l lt 2 lc rgbcolor "blue" t 'Periodic, T', \
  "< awk '{ print }' ".folder."BlochPeriodicConditions.txt1.csv" u 1:2 w l lt 2 lc rgbcolor "red" t 'Bloch, R', \
  "< awk '{ print }' ".folder."BlochPeriodicConditions.txt2.csv" u 1:2 w l lt 2 lc rgbcolor "blue" t 'Bloch, T'

