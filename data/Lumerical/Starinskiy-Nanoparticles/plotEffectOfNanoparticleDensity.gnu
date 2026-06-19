#!gnuplot

reset

set terminal postscript eps enhanced color font 'Helvetica, 24'
set output 'Wavelength-EffectOf2NPdistanceInSiOxDepositedOnSi.eps'

folder="03-NanoparticlesInFilmDepositedOnSubstrates/"

set xlabel 'Wavelength (um)'
set ylabel 'R, T'

set key top Left

set title 'Au NP (x2), 14 nm radius gapped in SiOx (44 nm) on Si'

p "< awk '{ print }' ".folder."AuNPx2-14nmRadius-CentersSeparatedby40nm-InsideSiOx-44nm-SiSubstrate-WavelengthRTA.txt1.csv" u 1:2 w l lt 2 lc rgbcolor "red" t 'gap=40nm, R', \
  "< awk '{ print }' ".folder."AuNPx2-14nmRadius-CentersSeparatedby40nm-InsideSiOx-44nm-SiSubstrate-WavelengthRTA.txt2.csv" u 1:2 w l lt 2 lc rgbcolor "blue" t 'gap=40nm, T', \
  "< awk '{ print }' ".folder."AuNPx2-14nmRadius-CentersSeparatedby60nm-InsideSiOx-44nm-SiSubstrate-WavelengthRTA.txt1.csv" u 1:2 w l lt 2 lc rgbcolor "red" t 'gap=60nm, R', \
  "< awk '{ print }' ".folder."AuNPx2-14nmRadius-CentersSeparatedby60nm-InsideSiOx-44nm-SiSubstrate-WavelengthRTA.txt2.csv" u 1:2 w l lt 2 lc rgbcolor "blue" t 'gap=60nm, T'

