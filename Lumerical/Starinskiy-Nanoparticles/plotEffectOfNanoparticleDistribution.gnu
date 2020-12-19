#!gnuplot

reset

set terminal postscript eps enhanced color font 'Helvetica, 24'

filmthickness="40"
set output 'Wavelength-SiOx('.filmthickness.'nm)OnSiO2-EffectOfRandomDistribution.eps'

folder="02-NanoparticlesOnSubstrates/"
folder2="03-NanoparticlesInFilmDepositedOnSubstrates/"
T_exp_AuNPSiO2="MeasuredTransmissionAuNP-SiO2.csv"
T_exp_AuNPSiOxSiO2="MeasuredTransmission-SiOx-AuNP-SiO2.csv"

xscale=1E-3

set xlabel 'Wavelength (um)'
set ylabel 'T'

set key below

p "< awk '{ print }' ".folder."AuNPx413-RandomUniformDistributionOnSquareArea250nmSide-11.8nmPM3nmNormalDistributionOfDiameter-SiO2substrate-RTA-Periodic-Seed1234-EvenBetterMesh.txt2.csv" u 1:2 w l lt 1 lc rgbcolor "blue"  t 'AuNP-Random-1234|SiO_2, T', \
   "< awk '{ print }' ".folder."AuNPx423-RandomUniformDistributionOnSquareArea250nmSide-11.8nmPM3nmNormalDistributionOfDiameter-SiO2substrate-RTA-Periodic-Seed4321-EvenBetterMesh.txt2.csv" u 1:2 w l lt 2 lc rgbcolor "blue" t 'AuNP-Random-4321|SiO_2, T', \
   "< awk '{ print }' ".folder."AuNPx426-RandomUniformDistributionOnSquareArea250nmSide-11.8nmPM3nmNormalDistributionOfDiameter-SiO2substrate-RTA-Periodic-Seed1423-EvenBetterMesh-2.txt2.csv" u 1:2 w l lt 3 lc rgbcolor "blue" t 'AuNP-Random-1423|SiO_2, T', \
   "< awk '{ print }' ".folder2."AuNPx413-RandomUniformDistributionOnSquareArea250nmSide-11.8nmPM3nmNormalDistributionOfDiameter-InsideSiOxFilm".filmthickness."nm-SiO2substrate-RTA-Periodic-Seed1234-BetterMesh.txt2.csv" u 1:2 w l lc rgbcolor "red" lt 1 t "(AuNP-Random-1234)+SiO_x(".filmthickness."nm)|SiO_2, T", \
   "< awk '{ print }' ".folder2."AuNPx423-RandomUniformDistributionOnSquareArea250nmSide-11.8nmPM3nmNormalDistributionOfDiameter-InsideSiOxFilm".filmthickness."nm-SiO2substrate-RTA-Periodic-Seed4321-BetterMesh.txt2.csv" u 1:2 w l lc rgbcolor "red" lt 2 t "(AuNP-Random-4321)+SiO_x(".filmthickness."nm)|SiO_2, T", \
   "< awk '{ print }' ".folder2."AuNPx426-RandomUniformDistributionOnSquareArea250nmSide-11.8nmPM3nmNormalDistributionOfDiameter-InsideSiOxFilm".filmthickness."nm-SiO2substrate-RTA-Periodic-Seed1423-BetterMesh.txt2.csv" u 1:2 w l lc rgbcolor "red" lt 3 t "(AuNP-Random-1423)+SiO_x(".filmthickness."nm)|SiO_2, T", \
   T_exp_AuNPSiO2 u     (xscale*$1):2 w l lc rgbcolor "blue" lt 4 t 'Exp. (AuNP 11.8nm)|SiO_2', \
   T_exp_AuNPSiOxSiO2 u (xscale*$1):2 w l lc rgbcolor "red"  lt 4 t 'Exp. SiO_x(AuNP 11.8nm)|SiO_2'


# p "< awk '{ print }' ".folder."AuNPx413-RandomUniformDistributionOnSquareArea250nmSide-11.8nmPM3nmNormalDistributionOfDiameter-SiO2substrate-RTA-Periodic-Seed1234-EvenBetterMesh.txt1.csv" u 1:2 w l lt 1 lc rgbcolor "red" t 'AuNP-Random-1234|SiO_2, R', \
#   "< awk '{ print }' ".folder."AuNPx423-RandomUniformDistributionOnSquareArea250nmSide-11.8nmPM3nmNormalDistributionOfDiameter-SiO2substrate-RTA-Periodic-Seed4321-EvenBetterMesh.txt1.csv" u 1:2 w l lt 2 lc rgbcolor "red" t 'AuNP-Random-4321|SiO_2, R', \
# "< awk '{ print }' ".folder."AuNPx426-RandomUniformDistributionOnSquareArea250nmSide-11.8nmPM3nmNormalDistributionOfDiameter-SiO2substrate-RTA-Periodic-Seed1423-EvenBetterMesh-2.txt1.csv" u 1:2 w l lt 3 lc rgbcolor "red" t 'AuNP-Random-1423|SiO_2, R', \


# "< awk '{ print }' 01-Substrates/SiOx-40nm-SiSubstrate-WavelengthRTA.txt1.csv" u 1:2 w l lt 1 lc rgbcolor "red"  t 'SiOx(40 nm)|Si, R', \
#  "< awk '{ print }' 01-Substrates/SiOx-44nm-SiSubstrate-WavelengthRTA.txt1.csv" u 1:2 w l lt 2 lc rgbcolor "red" t 'SiOx(44 nm)|Si, R', \
#  "< awk '{ print }' 01-Substrates/SiOx-40nm-SiSubstrate-WavelengthRTA.txt2.csv" u 1:2 w l lt 1 lc rgbcolor "blue"  t 'SiOx(40 nm)|Si, T', \
#  "< awk '{ print }' 01-Substrates/SiOx-44nm-SiSubstrate-WavelengthRTA.txt2.csv" u 1:2 w l lt 2 lc rgbcolor "blue" t 'SiOx(44 nm)|Si, T', \



