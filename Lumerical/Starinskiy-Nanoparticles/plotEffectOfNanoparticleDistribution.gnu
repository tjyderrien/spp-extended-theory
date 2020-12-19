#!gnuplot

reset

set terminal postscript eps enhanced color font 'Helvetica, 18'

filmthickness="44"
NPdiameter="11.8"
set output 'Wavelength-SiOx('.filmthickness.'nm)OnSiO2-EffectOfRandomDistribution-VariousNPSizes.eps'

folder="02-NanoparticlesOnSubstrates/"
folder2="03-NanoparticlesInFilmDepositedOnSubstrates/"
T_exp_AgNPSiO2="MeasuredTransmissionAuNP-SiO2.csv"
T_exp_AgNPSiOxSiO2="MeasuredTransmission-SiOx-AuNP-SiO2.csv"

xscale=1E-3

set xlabel 'Wavelength (um)'
set ylabel 'T'

set key below

p  "< awk '{ print }'  ".folder."AuNPx426-RandomUniformDistributionOnSquareArea250nmSide-11.8nmPM3nmNormalDistributionOfDiameter-SiO2substrate-RTA-Periodic-Seed1423-EvenBetterMesh-2.txt2.csv" u 1:2 w l lt 1 lc rgbcolor "blue" t 'AuNP-Random-1423 (11.8 nm)|SiO_2, T', \
   "< awk '{ print }'  ".folder."AuNP707-RandomUniformDistributionOnSquareArea250nmSide-9.23nmPM3nmNormalDistributionOfDiameter-SiO2substrate-RTA-Periodic-Seed1423-BetterMesh0.15-AuJohnson-WideSpectrum.txt2.csv" u 1:2 w l lt 1 lc rgbcolor "red" t 'AuNP-Random-1423 (9.23 nm)|SiO_2, T', \
   "< awk '{ print }'  ".folder."AuNP521-RandomUniformDistributionOnSquareArea250nmSide-6.4nmPM3nmNormalDistributionOfDiameter-SiO2substrate-RTA-Periodic-Seed1423-BetterMesh0.15-AuJohnson-WideSpectrum.txt2.csv" u 1:2 w l lt 1 lc rgbcolor "black" t 'AuNP-Random-1423 (6.4 mn)|SiO_2, T', \
   "< awk '{ print }' ".folder2."AuNPx426-RandomUniformDistributionOnSquareArea250nmSide-11.8nmPM3nmNormalDistributionOfDiameter-InsideSiOxFilm".filmthickness."nm-SiO2substrate-RTA-Periodic-Seed1423-BetterMesh.txt2.csv" u 1:2 w l lc rgbcolor "blue" lt 1 t "AuNP-Random-1423 (9.23 nm)+SiO_x(".filmthickness."nm)|SiO_2, T", \
   "< awk '{ print }' ".folder2."AuNP707-RandomUniformDistributionOnSquareArea250nmSide-9.23nmPM3nmNormalDistributionOfDiameter-InsideSiOxFilm44nm-SiO2substrate-RTA-Periodic-Seed1423-BetterMesh-SiOxOfSergey-AuJohnson-WideSpectrum-Mesh1.5nm.txt2.csv" u 1:2 w l lc rgbcolor "red" lt 1 t "(AuNP-Random-1423)+SiO_x(".filmthickness."nm)|SiO_2, T", \
   "< awk '{ print }' ".folder2."AuNP521-RandomUniformDistributionOnSquareArea250nmSide-6.4nmPM3nmNormalDistributionOfDiameter-InsideSiOxFilm".filmthickness."nm-SiO2substrate-RTA-Periodic-Seed1423-BetterMesh-SiOxOfSergey-AuJohnson-WideSpectrum-Mesh1.5nm.txt2.csv" u 1:2 w l lc rgbcolor "black" lt 1 t "(AuNP-Random-1423)+SiO_x(".filmthickness."nm)|SiO_2, T", \
   T_exp_AgNPSiO2 u (xscale*$1):2 every 10 w p lc rgbcolor "blue"      ps 0.5 lt 1 t 'Exp. AuNP|SiO_2', \
   T_exp_AgNPSiO2 u (xscale*$1):3 every 10 w p lc rgbcolor "red"       ps 0.5 lt 1 notitle, \
   T_exp_AgNPSiO2 u (xscale*$1):4 every 10 w p lc rgbcolor "black"     ps 0.5 lt 1 notitle, \
   T_exp_AgNPSiOxSiO2 u (xscale*$1):2 every 100 w p lc rgbcolor "blue"  ps 0.5 lt 2 t 'Exp. SiO_x(AuNP)|SiO_2', \
   T_exp_AgNPSiOxSiO2 u (xscale*$1):3 every 100 w p lc rgbcolor "red"   ps 0.5 lt 2 notitle, \
   T_exp_AgNPSiOxSiO2 u (xscale*$1):4 every 100 w p lc rgbcolor "black" ps 0.5 lt 2 notitle

# Effect of random distribution

# p "< awk '{ print }' ".folder."AuNPx413-RandomUniformDistributionOnSquareArea250nmSide-11.8nmPM3nmNormalDistributionOfDiameter-SiO2substrate-RTA-Periodic-Seed1234-EvenBetterMesh.txt2.csv" u 1:2 w l lt 1 lc rgbcolor "blue"  t 'AuNP-Random-1234 (11.8 nm)|SiO_2, T', \
#   "< awk '{ print }' ".folder."AuNPx423-RandomUniformDistributionOnSquareArea250nmSide-11.8nmPM3nmNormalDistributionOfDiameter-SiO2substrate-RTA-Periodic-Seed4321-EvenBetterMesh.txt2.csv" u 1:2 w l lt 2 lc rgbcolor "blue" t 'AuNP-Random-4321 (11.8 nm)|SiO_2, T', \
#   "< awk '{ print }' ".folder2."AuNPx413-RandomUniformDistributionOnSquareArea250nmSide-11.8nmPM3nmNormalDistributionOfDiameter-InsideSiOxFilm".filmthickness."nm-SiO2substrate-RTA-Periodic-Seed1234-BetterMesh.txt2.csv" u 1:2 w l lc rgbcolor "red" lt 1 t "(AuNP-Random-1234)+SiO_x(".filmthickness."nm)|SiO_2, T", \
#   "< awk '{ print }' ".folder2."AuNPx423-RandomUniformDistributionOnSquareArea250nmSide-11.8nmPM3nmNormalDistributionOfDiameter-InsideSiOxFilm".filmthickness."nm-SiO2substrate-RTA-Periodic-Seed4321-BetterMesh.txt2.csv" u 1:2 w l lc rgbcolor "red" lt 2 t "(AuNP-Random-4321)+SiO_x(".filmthickness."nm)|SiO_2, T", \



# p "< awk '{ print }' ".folder."AuNPx413-RandomUniformDistributionOnSquareArea250nmSide-11.8nmPM3nmNormalDistributionOfDiameter-SiO2substrate-RTA-Periodic-Seed1234-EvenBetterMesh.txt1.csv" u 1:2 w l lt 1 lc rgbcolor "red" t 'AuNP-Random-1234|SiO_2, R', \
#   "< awk '{ print }' ".folder."AuNPx423-RandomUniformDistributionOnSquareArea250nmSide-11.8nmPM3nmNormalDistributionOfDiameter-SiO2substrate-RTA-Periodic-Seed4321-EvenBetterMesh.txt1.csv" u 1:2 w l lt 2 lc rgbcolor "red" t 'AuNP-Random-4321|SiO_2, R', \
# "< awk '{ print }' ".folder."AuNPx426-RandomUniformDistributionOnSquareArea250nmSide-11.8nmPM3nmNormalDistributionOfDiameter-SiO2substrate-RTA-Periodic-Seed1423-EvenBetterMesh-2.txt1.csv" u 1:2 w l lt 3 lc rgbcolor "red" t 'AuNP-Random-1423|SiO_2, R', \


# "< awk '{ print }' 01-Substrates/SiOx-40nm-SiSubstrate-WavelengthRTA.txt1.csv" u 1:2 w l lt 1 lc rgbcolor "red"  t 'SiOx(40 nm)|Si, R', \
#  "< awk '{ print }' 01-Substrates/SiOx-44nm-SiSubstrate-WavelengthRTA.txt1.csv" u 1:2 w l lt 2 lc rgbcolor "red" t 'SiOx(44 nm)|Si, R', \
#  "< awk '{ print }' 01-Substrates/SiOx-40nm-SiSubstrate-WavelengthRTA.txt2.csv" u 1:2 w l lt 1 lc rgbcolor "blue"  t 'SiOx(40 nm)|Si, T', \
#  "< awk '{ print }' 01-Substrates/SiOx-44nm-SiSubstrate-WavelengthRTA.txt2.csv" u 1:2 w l lt 2 lc rgbcolor "blue" t 'SiOx(44 nm)|Si, T', \



