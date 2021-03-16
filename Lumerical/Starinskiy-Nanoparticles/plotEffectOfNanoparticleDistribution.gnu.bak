#!gnuplot

reset

set terminal postscript eps enhanced color font 'Helvetica, 18'

filmthickness="44"
NPdiameter="11.8"

dataname="RandomUniformDistributionOnSquareArea250nmSide"
substrate="SiO2substrate-RTA-Periodic"


folder="02-NanoparticlesOnSubstrates/"
folder2="03-NanoparticlesInFilmDepositedOnSubstrates/"
T_exp_AuNPSiO2="MeasuredTransmissionAuNP-SiO2.csv"
T_exp_AuNPSiOxSiO2="MeasuredTransmission-SiOx-AuNP-SiO2.csv"

xscale=1E-3

set xlabel 'Wavelength (um)'
set ylabel 'T'

# set title 'NP Au 11.8 nm'

set key below

set output 'EffectOfMesh-AuNP-SiO2.eps'

### Plotting the effect of mesh
p "./02-NanoparticlesOnSubstrates/AuNP521-RandomUniformDistributionOnSquareArea250nmSide-6.4nmPM3nmNormalDistributionOfDiameter-SiO2substrate-RTA-Periodic-Seed4321-BetterMesh1.0-AuJohnson-WideSpectrum.txt2.csv" w l lt 1 lc rgbcolor "black" t "NP|SiO_2, diam 6.4 nm, Mesh 1.0 nm", \
  "03-NanoparticlesInFilmDepositedOnSubstrates/AuNP521-RandomUniformDistributionOnSquareArea250nmSide-6.4nmPM3nmNormalDistributionOfDiameter-InsideSiOxFilm44nm-SiO2substrate-RTA-Periodic-Seed1423-BetterMesh-SiOxOfSergey-AuJohnson-WideSpectrum-Mesh2.0nm.txt2.csv" w l lt 1 lc rgbcolor "gray" t "SiO_x|NP|SiO_2, diam 6.4 nm, Mesh 2 nm", \
  "03-NanoparticlesInFilmDepositedOnSubstrates/AuNP521-RandomUniformDistributionOnSquareArea250nmSide-6.4nmPM3nmNormalDistributionOfDiameter-InsideSiOxFilm44nm-SiO2substrate-RTA-Periodic-Seed1423-BetterMesh-SiOxOfSergey-AuJohnson-WideSpectrum-Mesh1.5nm.txt2.csv" w l lt 2 lc rgbcolor "gray" t "SiO_x|NP|SiO_2, diam 6.4 nm, Mesh 1.5 nm", \
  "03-NanoparticlesInFilmDepositedOnSubstrates/AuNP521-RandomUniformDistributionOnSquareArea250nmSide-6.4nmPM3nmNormalDistributionOfDiameter-InsideSiOxFilm44nm-SiO2substrate-RTA-Periodic-Seed1423-BetterMesh-SiOxOfSergey-AuJohnson-WideSpectrum-Mesh1.0nm.txt2.csv" w l lt 3 lc rgbcolor "gray" t "SiO_x|NP|SiO_2, diam 6.4 nm, Mesh 1 nm", \
  "02-NanoparticlesOnSubstrates/AuNP571-RandomUniformDistributionOnSquareArea250nmSide-11.8nmPM3nmNormalDistributionOfDiameter-SiO2substrate-RTA-Periodic-Seed1423-BetterMesh1.5-AuJohnson-WideSpectrum.txt2.csv"  w l lt 1 lc rgbcolor "blue" t "NP|SiO_2, diam 11.8 nm, Mesh 1.5 nm", \
  "02-NanoparticlesOnSubstrates/AuNP571-RandomUniformDistributionOnSquareArea250nmSide-11.8nmPM4nmNormalDistributionOfDiameter-SiO2substrate-RTA-Periodic-Seed1423-BetterMesh0.15-AuJohnson-WideSpectrum.txt2.csv" w l lt 2 lc rgbcolor "blue" t "NP|SiO_2, diam 11.8 nm, Mesh 1.5 nm", \
  "02-NanoparticlesOnSubstrates/AuNP571-RandomUniformDistributionOnSquareArea250nmSide-11.8nmPM4nmNormalDistributionOfDiameter-SiO2substrate-RTA-Periodic-Seed1423-BetterMesh0.13-AuJohnson-WideSpectrum.txt2.csv" w l lt 3 lc rgbcolor "blue" t "NP|SiO_2, diam 11.8 nm, Mesh 1.3 nm", \
  "02-NanoparticlesOnSubstrates/AuNPxxx-RandomUniformDistributionOnSquareArea250nmSide-11.8nmPM3nmNormalDistributionOfDiameter-SiO2substrate-RTA-Periodic-Seed4321-BetterMesh1.0-AuJohnson-WideSpectrum.txt2.csv" w l lt 4 lc rgbcolor "blue" t "NP|SiO_2, diam 11.8 nm, Mesh 1.0 nm", \
  "03-NanoparticlesInFilmDepositedOnSubstrates/AuNP571-RandomUniformDistributionOnSquareArea250nmSide-11.8nmPM4nmNormalDistributionOfDiameter-InsideSiOxFilm44nm-SiO2substrate-RTA-Periodic-Seed1423-BetterMesh-SiOxOfSergey-AuJohnson-WideSpectrum-Mesh1.5nm.txt2.csv" w l lt 1 lc rgbcolor "red" t  "SiO_x|NP|SiO_2, diam 11.8 nm, Mesh 1.5 nm", \
  "03-NanoparticlesInFilmDepositedOnSubstrates/AuNP571-RandomUniformDistributionOnSquareArea250nmSide-11.8nmPM4nmNormalDistributionOfDiameter-InsideSiOxFilm44nm-SiO2substrate-RTA-Periodic-Seed1423-BetterMesh-SiOxOfSergey-AuJohnson-WideSpectrum-Mesh1.3nm.txt2.csv" w l lt 2 lc rgbcolor "red" t  "SiO_x|NP|SiO_2, diam 11.8 nm, Mesh 1.3 nm", \
  "03-NanoparticlesInFilmDepositedOnSubstrates/AuNP384-RandomUniformDistributionOnSquareArea250nmSide-11.8nmPM3nmNormalDistributionOfDiameter-InsideSiOxFilm44nm-SiO2substrate-RTA-Periodic-Seed1423-BetterMesh-SiOxOfSergey-AuJohnson-WideSpectrum-Mesh1.5nm-FDTD200nm.txt2.csv" w l lt 3 lc rgbcolor "red" t  "SiO_x|384 NP|SiO_2, diam 11.8 nm, Mesh 1.5 nm, FDTD 200 nm", \
  "03-NanoparticlesInFilmDepositedOnSubstrates/AuNP384-RandomUniformDistributionOnSquareArea250nmSide-11.8nmPM3nmNormalDistributionOfDiameter-InsideSiOxFilm44nm-SiO2substrate-RTA-Periodic-Seed1423-BetterMesh-SiOxOfSergey-AuJohnson-WideSpectrum-Mesh1.3nm.txt2.csv" w l lt 4 lc rgbcolor "red" t  "SiO_x|384 NP|SiO_2, diam 11.8 nm, Mesh 1.3 nm", \
  "03-NanoparticlesInFilmDepositedOnSubstrates/AuNPxxx-RandomUniformDistributionOnSquareArea250nmSide-11.8nmPM3nmNormalDistributionOfDiameter-InsideSiOxFilm44nm-SiO2substrate-RTA-Periodic-Seed4321-BetterMesh-SiOxOfSergey-AuJohnson-WideSpectrum-Mesh1.0nm-FDTD250nm.txt2.csv" w l lt 5 lc rgbcolor "red" t  "SiO_x|384 NP|SiO_2, diam 11.8 nm, Mesh 1.0 nm"
  
exit

# Running again for other seeds, using same mesh. 
set output 'EffectOfSeed-AuNP-SiO2.eps'

p "02-NanoparticlesOnSubstrates/AuNP571-RandomUniformDistributionOnSquareArea250nmSide-11.8nmPM4nmNormalDistributionOfDiameter-SiO2substrate-RTA-Periodic-Seed1423-BetterMesh0.15-AuJohnson-WideSpectrum.txt2.csv" w l lt 1 lc rgbcolor "blue" t "NP1423|SiO_2, diam 11.8 nm, Mesh 1.5 nm", \
  "02-NanoparticlesOnSubstrates/AuNPxxx-RandomUniformDistributionOnSquareArea250nmSide-11.8nmPM3nmNormalDistributionOfDiameter-SiO2substrate-RTA-Periodic-Seed1234-BetterMesh1.5-AuJohnson-WideSpectrum.txt2.csv" w l lt 2 lc rgbcolor "blue" t "NP1234|SiO_2, diam 11.8 nm, Mesh 1.5 nm", \
  "02-NanoparticlesOnSubstrates/AuNPxxx-RandomUniformDistributionOnSquareArea250nmSide-11.8nmPM3nmNormalDistributionOfDiameter-SiO2substrate-RTA-Periodic-Seed4321-BetterMesh1.5-AuJohnson-WideSpectrum.txt2.csv" w l lt 3 lc rgbcolor "blue" t "NP4321|SiO_2, diam 11.8 nm, Mesh 1.5 nm", \
  "03-NanoparticlesInFilmDepositedOnSubstrates/AuNP384-RandomUniformDistributionOnSquareArea250nmSide-11.8nmPM3nmNormalDistributionOfDiameter-InsideSiOxFilm44nm-SiO2substrate-RTA-Periodic-Seed1423-BetterMesh-SiOxOfSergey-AuJohnson-WideSpectrum-Mesh1.5nm-FDTD200nm.txt2.csv" w l lt 1 lc rgbcolor "red" t "SiO_x|NP1423|SiO_2, diam 11.8 nm, Mesh 1.5 nm", \
  "03-NanoparticlesInFilmDepositedOnSubstrates/AuNPxxx-RandomUniformDistributionOnSquareArea250nmSide-11.8nmPM3nmNormalDistributionOfDiameter-InsideSiOxFilm44nm-SiO2substrate-RTA-Periodic-Seed1234-BetterMesh-SiOxOfSergey-AuJohnson-WideSpectrum-Mesh1.5nm-FDTD250nm.txt2.csv" w l lt 2 lc rgbcolor "red" t "SiO_x|NP1234|SiO_2, diam 11.8 nm, Mesh 1.5 nm", \
  "03-NanoparticlesInFilmDepositedOnSubstrates/AuNPxxx-RandomUniformDistributionOnSquareArea250nmSide-11.8nmPM3nmNormalDistributionOfDiameter-InsideSiOxFilm44nm-SiO2substrate-RTA-Periodic-Seed4321-BetterMesh-SiOxOfSergey-AuJohnson-WideSpectrum-Mesh1.5nm-FDTD250nm.txt2.csv" w l lt 3 lc rgbcolor "red" t "SiO_x|NP4321|SiO_2, diam 11.8 nm, Mesh 1.5 nm"

exit

# Effect of the NP size can be studied only now, once we have converged mesh and checked effect of the seed. 
set output 'EffectOfNPSize-AuNP-SiO2.eps'
# Effect of NP diameter
p  folder."AuNP521-".dataname."-6.4nmPM3nmNormalDistributionOfDiameter-".substrate."-Seed1423-BetterMesh0.15-AuJohnson-WideSpectrum.txt2.csv" u 1:2 w l lt 1 lc rgbcolor "black" t 'AuNP-Random-1423 (6.4 mn)|SiO_2, T', \
   folder."AuNP707-".dataname."-9.23nmPM3nmNormalDistributionOfDiameter-".substrate."-Seed1423-BetterMesh0.15-AuJohnson-WideSpectrum.txt2.csv" u 1:2 w l lt 1 lc rgbcolor "red" t 'AuNP-Random-1423 (9.23 nm)|SiO_2, T', \
   folder2."AuNPx426-".dataname."-11.8nmPM3nmNormalDistributionOfDiameter-InsideSiOxFilm".filmthickness."nm-".substrate."-Seed1423-BetterMesh.txt2.csv" u 1:2 w l lc rgbcolor "blue" lt 1 t "AuNP-Random-1423 (9.23 nm)+SiO_x(".filmthickness."nm)|SiO_2, T", \
   folder2."AuNP707-".dataname."-9.23nmPM3nmNormalDistributionOfDiameter-InsideSiOxFilm44nm-".substrate."-Seed1423-BetterMesh-SiOxOfSergey-AuJohnson-WideSpectrum-Mesh1.5nm.txt2.csv" u 1:2 w l lc rgbcolor "red" lt 1 t "(AuNP-Random-1423)+SiO_x(".filmthickness."nm)|SiO_2, T", \
   folder2."AuNP521-".dataname."-6.4nmPM3nmNormalDistributionOfDiameter-InsideSiOxFilm".filmthickness."nm-".substrate."-Seed1423-BetterMesh-SiOxOfSergey-AuJohnson-WideSpectrum-Mesh1.5nm.txt2.csv" u 1:2 w l lc rgbcolor "black" lt 1 t "(AuNP-Random-1423)+SiO_x(".filmthickness."nm)|SiO_2, T", \
   T_exp_AgNPSiO2 u (xscale*$1):2 every 10 w p lc rgbcolor "blue"      ps 0.5 lt 1 t 'Exp. AuNP|SiO_2', \
   T_exp_AgNPSiO2 u (xscale*$1):3 every 10 w p lc rgbcolor "red"       ps 0.5 lt 1 notitle, \
   T_exp_AgNPSiO2 u (xscale*$1):4 every 10 w p lc rgbcolor "black"     ps 0.5 lt 1 notitle, \
   T_exp_AgNPSiOxSiO2 u (xscale*$1):2 every 100 w p lc rgbcolor "blue"  ps 0.5 lt 2 t 'Exp. SiO_x(AuNP)|SiO_2', \
   T_exp_AgNPSiOxSiO2 u (xscale*$1):3 every 100 w p lc rgbcolor "red"   ps 0.5 lt 2 notitle, \
   T_exp_AgNPSiOxSiO2 u (xscale*$1):4 every 100 w p lc rgbcolor "black" ps 0.5 lt 2 notitle

# Effect of random distribution
p "< awk '{ print }' ".folder."AuNPx413-".dataname."-11.8nmPM3nmNormalDistributionOfDiameter-".substrate."-Seed1234-EvenBetterMesh.txt2.csv" u 1:2 w l lt 1 lc rgbcolor "blue"  t 'AuNP-Random-1234 (11.8 nm)|SiO_2, T', \
  "< awk '{ print }' ".folder."AuNPx423-".dataname."-11.8nmPM3nmNormalDistributionOfDiameter-".substrate."-Seed4321-EvenBetterMesh.txt2.csv" u 1:2 w l lt 2 lc rgbcolor "blue" t 'AuNP-Random-4321 (11.8 nm)|SiO_2, T', \
  "< awk '{ print }' ".folder2."AuNPx413-".dataname."-11.8nmPM3nmNormalDistributionOfDiameter-InsideSiOxFilm".filmthickness."nm-".substrate."-Seed1234-BetterMesh.txt2.csv" u 1:2 w l lc rgbcolor "red" lt 1 t "(AuNP-Random-1234)+SiO_x(".filmthickness."nm)|SiO_2, T", \
  "< awk '{ print }' ".folder2."AuNPx423-".dataname."-11.8nmPM3nmNormalDistributionOfDiameter-InsideSiOxFilm".filmthickness."nm-".substrate."-Seed4321-BetterMesh.txt2.csv" u 1:2 w l lc rgbcolor "red" lt 2 t "(AuNP-Random-4321)+SiO_x(".filmthickness."nm)|SiO_2, T", \

### Plotting the effect of seed for randomization
p "< awk '{ print }' ".folder."AuNPx413-".dataname."-11.8nmPM3nmNormalDistributionOfDiameter-".substrate."-Seed1234-EvenBetterMesh.txt2.csv" u 1:2 w l lt 1 lc rgbcolor "blue"  t 'AuNP-Random-1234|SiO_2, T', \
   "< awk '{ print }' ".folder."AuNPx423-".dataname."-11.8nmPM3nmNormalDistributionOfDiameter-".substrate."-Seed4321-EvenBetterMesh.txt2.csv" u 1:2 w l lt 2 lc rgbcolor "blue" t 'AuNP-Random-4321|SiO_2, T', \
   "< awk '{ print }' ".folder."AuNPx426-".dataname."-11.8nmPM3nmNormalDistributionOfDiameter-".substrate."-Seed1423-EvenBetterMesh-2.txt2.csv" u 1:2 w l lt 3 lc rgbcolor "blue" t 'AuNP-Random-1423|SiO_2, T', \
   "< awk '{ print }' ".folder."AuNPx423-".dataname."-11.8nmPM3nmNormalDistributionOfDiameter-".substrate."-Seed4321-EvenBetterMesh.txt2.csv" u 1:2 w l lt 4 lc rgbcolor "blue" t 'AuNP-Random-4321|SiO_2, T'

### Plotting the effect of NP dispersion
p  "< awk '{ print }' ".folder."AuNPx426-".dataname."-11.8nmPM3nmNormalDistributionOfDiameter-".substrate."-Seed1423-EvenBetterMesh-2.txt2.csv" u 1:2 w l lt 1 lc rgbcolor "blue" t 'AuNP-Random-1423|SiO_2, T', \
   "< awk '{ print }' ".folder."AuNP571-".dataname."-11.8nmPM3nmNormalDistributionOfDiameter-".substrate."-Seed1423-BetterMesh0.13-AuJohnson-WideSpectrum.txt2.csv" u 1:2 w l lt 2 lc rgbcolor "blue" t 'AuNP-Random-1423, Mesh 0.13, diam. pm 3 nm', \
   "< awk '{ print }' ".folder."AuNP571-".dataname."-11.8nmPM4nmNormalDistributionOfDiameter-".substrate."-Seed1423-BetterMesh0.13-AuJohnson-WideSpectrum.txt2.csv" u 1:2 w l lt 3 lc rgbcolor "blue" t 'AuNP-Random-1423, Mesh 0.13, diam. pm 4 nm'

### Plotting effect of the embedding
p  "< awk '{ print }' ".folder."AuNP571-".dataname."-11.8nmPM3nmNormalDistributionOfDiameter-".substrate."-Seed1423-BetterMesh0.13-AuJohnson-WideSpectrum.txt2.csv" u 1:2 w l lt 4 lc rgbcolor "blue" t 'AuNP-Random-1423, Mesh 0.13, diam. pm 3 nm', \
   "< awk '{ print }' ".folder."AuNP571-".dataname."-11.8nmPM4nmNormalDistributionOfDiameter-".substrate."-Seed1423-BetterMesh0.13-AuJohnson-WideSpectrum.txt2.csv" u 1:2 w l lt 5 lc rgbcolor "blue" t 'AuNP-Random-1423, Mesh 0.13, diam. pm 4 nm', \
   "< awk '{ print }' ".folder."AuNP571-".dataname."-11.8nmPM4nmNormalDistributionOfDiameter-".substrate."-Seed1423-BetterMesh0.15-AuJohnson-WideSpectrum.txt2.csv" u 1:2 w l lt 6 lc rgbcolor "blue" t 'AuNP-Random-1423, Mesh 0.15, diam. pm 4 nm', \
   "< awk '{ print }' ".folder2."AuNPx413-".dataname."-11.8nmPM3nmNormalDistributionOfDiameter-InsideSiOxFilm".filmthickness."nm-".substrate."-Seed1234-BetterMesh.txt2.csv" u 1:2 w l lc rgbcolor "red" lt 1 t "(AuNP-Random-1234)+SiO_x(".filmthickness."nm)|SiO_2, T", \
   "< awk '{ print }' ".folder2."AuNPx423-".dataname."-11.8nmPM3nmNormalDistributionOfDiameter-InsideSiOxFilm".filmthickness."nm-".substrate."-Seed4321-BetterMesh.txt2.csv" u 1:2 w l lc rgbcolor "red" lt 2 t "(AuNP-Random-4321)+SiO_x(".filmthickness."nm)|SiO_2, T", \
   "< awk '{ print }' ".folder2."AuNPx426-".dataname."-11.8nmPM3nmNormalDistributionOfDiameter-InsideSiOxFilm".filmthickness."nm-".substrate."-Seed1423-BetterMesh.txt2.csv" u 1:2 w l lc rgbcolor "red" lt 3 t "(AuNP-Random-1423)+SiO_x(".filmthickness."nm)|SiO_2, T", \
   "< awk '{ print }' ".folder2."AuNP571-".dataname."-11.8nmPM4nmNormalDistributionOfDiameter-InsideSiOxFilm".filmthickness."nm-".substrate."-Seed1423-BetterMesh-SiOxOfSergey-AuJohnson-WideSpectrum-Mesh1.5nm.txt2.csv" u 1:2 w l lc rgbcolor "red" lt 4 t "(AuNP-Random-1423)+SiO_x(".filmthickness."nm)|SiO_2, T, Mesh 1.5, Diam. pm 4 nm", \
   "< awk '{ print }' ".folder2."AuNP571-".dataname."-11.8nmPM4nmNormalDistributionOfDiameter-InsideSiOxFilm".filmthickness."nm-".substrate."-Seed1423-BetterMesh-SiOxOfSergey-AuJohnson-WideSpectrum-Mesh1.3nm.txt2.csv" u 1:2 w l lc rgbcolor "red" lt 5 t "(AuNP-Random-1423)+SiO_x(".filmthickness."nm)|SiO_2, T, Mesh 1.3, Diam. pm 4 nm", \
   "< awk '{ print }' ".folder2."AuNP571-".dataname."-11.8nmPM3nmNormalDistributionOfDiameter-InsideSiOxFilm".filmthickness."nm-".substrate."-Seed1423-BetterMesh-SiOxOfSergey-AuJohnson-WideSpectrum-Mesh1.3nm.txt2.csv" u 1:2 w l lc rgbcolor "red" lt 6 t "(AuNP-Random-1423)+SiO_x(".filmthickness."nm)|SiO_2, T, Mesh 1.3, Diam. pm 3 nm", \
   T_exp_AuNPSiO2 u     (xscale*$1):2 w l lc rgbcolor "blue" lt 4 t 'Exp. (AuNP 11.8nm)|SiO_2', \
   T_exp_AuNPSiOxSiO2 u (xscale*$1):2 w l lc rgbcolor "red"  lt 4 t 'Exp. SiO_x(AuNP 11.8nm)|SiO_2'

# Effect of randomization seed (reflectivity)
# p "< awk '{ print }' ".folder."AuNPx413-".dataname."-11.8nmPM3nmNormalDistributionOfDiameter-".substrate."-Seed1234-EvenBetterMesh.txt1.csv" u 1:2 w l lt 1 lc rgbcolor "red" t 'AuNP-Random-1234|SiO_2, R', \
#   "< awk '{ print }' ".folder."AuNPx423-".dataname."-11.8nmPM3nmNormalDistributionOfDiameter-".substrate."-Seed4321-EvenBetterMesh.txt1.csv" u 1:2 w l lt 2 lc rgbcolor "red" t 'AuNP-Random-4321|SiO_2, R', \
# "< awk '{ print }' ".folder."AuNPx426-".dataname."-11.8nmPM3nmNormalDistributionOfDiameter-".substrate."-Seed1423-EvenBetterMesh-2.txt1.csv" u 1:2 w l lt 3 lc rgbcolor "red" t 'AuNP-Random-1423|SiO_2, R', \

# Effect of embedding film thickness (reflectivity, transmission)
# "< awk '{ print }' 01-Substrates/SiOx-40nm-SiSubstrate-WavelengthRTA.txt1.csv" u 1:2 w l lt 1 lc rgbcolor "red"  t 'SiOx(40 nm)|Si, R', \
#  "< awk '{ print }' 01-Substrates/SiOx-44nm-SiSubstrate-WavelengthRTA.txt1.csv" u 1:2 w l lt 2 lc rgbcolor "red" t 'SiOx(44 nm)|Si, R', \
#  "< awk '{ print }' 01-Substrates/SiOx-40nm-SiSubstrate-WavelengthRTA.txt2.csv" u 1:2 w l lt 1 lc rgbcolor "blue"  t 'SiOx(40 nm)|Si, T', \
#  "< awk '{ print }' 01-Substrates/SiOx-44nm-SiSubstrate-WavelengthRTA.txt2.csv" u 1:2 w l lt 2 lc rgbcolor "blue" t 'SiOx(44 nm)|Si, T', \



