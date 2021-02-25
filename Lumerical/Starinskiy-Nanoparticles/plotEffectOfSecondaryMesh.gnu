#!gnuplot

reset

set terminal postscript eps enhanced color font 'Helvetica, 24'

filmthickness="44"
NPsize="6.4"
seed="1423"

set output 'Wavelength-SiOx('.filmthickness.'nm)OnSiO2-EffectOfSecondMeshPrecision.eps'

folder="02-NanoparticlesOnSubstrates/"
folder2="03-NanoparticlesInFilmDepositedOnSubstrates/"
T_exp_AgNPSiO2="MeasuredTransmissionAuNP-SiO2.csv"
T_exp_AgNPSiOxSiO2="MeasuredTransmission-SiOx-AuNP-SiO2.csv"

xscale=1E-3

set xlabel 'Wavelength (um)'
set ylabel 'T'

set key below

p folder2."./AuNP521-RandomUniformDistributionOnSquareArea250nmSide-".NPsize."nmPM3nmNormalDistributionOfDiameter-InsideSiOxFilm".filmthickness."nm-SiO2substrate-RTA-Periodic-Seed".seed."-BetterMesh-SiOxOfSergey-AuJohnson-WideSpectrum.txt2.csv" u 1:2 w l t 'Mesh 2.5', \
  folder2."./AuNP521-RandomUniformDistributionOnSquareArea250nmSide-".NPsize."nmPM3nmNormalDistributionOfDiameter-InsideSiOxFilm".filmthickness."nm-SiO2substrate-RTA-Periodic-Seed".seed."-BetterMesh-SiOxOfSergey-AuJohnson-WideSpectrum-Mesh2.0nm.txt2.csv" u 1:2 w l t 'Mesh 2.0', \
  folder2."AuNP521-RandomUniformDistributionOnSquareArea250nmSide-".NPsize."nmPM3nmNormalDistributionOfDiameter-InsideSiOxFilm".filmthickness."nm-SiO2substrate-RTA-Periodic-Seed".seed."-BetterMesh-SiOxOfSergey-AuJohnson-WideSpectrum-Mesh1.5nm.txt2.csv" u 1:2 w l t 'Mesh 1.5', \
  folder2."AuNP521-RandomUniformDistributionOnSquareArea250nmSide-".NPsize."nmPM3nmNormalDistributionOfDiameter-InsideSiOxFilm".filmthickness."nm-SiO2substrate-RTA-Periodic-Seed".seed."-BetterMesh-SiOxOfSergey-AuJohnson-WideSpectrum-Mesh1.0nm.txt2.csv" u 1:2 w l t 'Mesh 1.0', \
  folder2."AuNP521-RandomUniformDistributionOnSquareArea250nmSide-".NPsize."nmPM3nmNormalDistributionOfDiameter-InsideSiOxFilm".filmthickness."nm-SiO2substrate-RTA-Periodic-Seed".seed."-BetterMesh-SiOxOfSergey-AuJohnson-WideSpectrum-NoMesh.txt2.csv" u 1:2  w l t 'No sec. mesh"


# p "< awk '{ print }' ".folder."AuNPx413-RandomUniformDistributionOnSquareArea250nmSide-11.8nmPM3nmNormalDistributionOfDiameter-SiO2substrate-RTA-Periodic-Seed1234-EvenBetterMesh.txt1.csv" u 1:2 w l lt 1 lc rgbcolor "red" t 'AuNP-Random-1234|SiO_2, R', \
#   "< awk '{ print }' ".folder."AuNPx423-RandomUniformDistributionOnSquareArea250nmSide-11.8nmPM3nmNormalDistributionOfDiameter-SiO2substrate-RTA-Periodic-Seed4321-EvenBetterMesh.txt1.csv" u 1:2 w l lt 2 lc rgbcolor "red" t 'AuNP-Random-4321|SiO_2, R', \
# "< awk '{ print }' ".folder."AuNPx426-RandomUniformDistributionOnSquareArea250nmSide-11.8nmPM3nmNormalDistributionOfDiameter-SiO2substrate-RTA-Periodic-Seed1423-EvenBetterMesh-2.txt1.csv" u 1:2 w l lt 3 lc rgbcolor "red" t 'AuNP-Random-1423|SiO_2, R', \


# "< awk '{ print }' 01-Substrates/SiOx-40nm-SiSubstrate-WavelengthRTA.txt1.csv" u 1:2 w l lt 1 lc rgbcolor "red"  t 'SiOx(40 nm)|Si, R', \
#  "< awk '{ print }' 01-Substrates/SiOx-44nm-SiSubstrate-WavelengthRTA.txt1.csv" u 1:2 w l lt 2 lc rgbcolor "red" t 'SiOx(44 nm)|Si, R', \
#  "< awk '{ print }' 01-Substrates/SiOx-40nm-SiSubstrate-WavelengthRTA.txt2.csv" u 1:2 w l lt 1 lc rgbcolor "blue"  t 'SiOx(40 nm)|Si, T', \
#  "< awk '{ print }' 01-Substrates/SiOx-44nm-SiSubstrate-WavelengthRTA.txt2.csv" u 1:2 w l lt 2 lc rgbcolor "blue" t 'SiOx(44 nm)|Si, T', \



