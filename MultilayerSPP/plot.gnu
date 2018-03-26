#!gnuplot

reset

# FIG 1
set terminal postscript eps enhanced color font "Helvetica, 24"
set output 'EffectOfThickness.eps'

fileThickness="/home/thibault/Documents/spp-extended-theory.lipss_database/LIPSS-database/Experiments/Dostovalov/Cr/SummaryOfResults.csv"
fileMaterials="/home/thibault/Documents/spp-extended-theory.lipss_database/LIPSS-database/Experiments/Dostovalov/Cr/SummaryOfResults-Substrate.csv"

set title 'Laser wavelength 1026 nm'

unset key
# set key out bottom center horizontal Left
# set key out right center Left

set log x
# set xtics 0.1
# set ytics 0.1

set xrange [10:400]
set yrange [10:1200]

set xlabel 'Film thickness (nm)'
set ylabel 'Period (nm)'

wavelength=1026
nCr = {-0.67,24.87}**0.5 #Cr
nCrO2 = {4.9173, 0.1784}**0.5 #CrO2
print real(nCrO2)
Reduced(n)=wavelength/real(n)

plot fileThickness u 4:9 w l lw 6 lc 1 lt 2 t 'TF-SPP {/Symbol l}', \
fileThickness u 4:8 w l lw 6 lc 3 lt 1 t 'TF-SPP coupled', \
fileThickness u 4:10 w l lw 6 lc 7 lt 3 t 'TF-SPP {/Symbol l}/n[Cr]', \
wavelength w l lw 3 lc 7 lt 2 t '{/Symbol l}', \
Reduced(nCr) w l lw 3 lc 7 lt 4 t '{/Symbol l}/n[Cr]', \
fileThickness u 4:6 w p lc 1 lw 10 pt 4 t 'Measured LSFL ||', \
fileThickness u 4:5 w p lc 3 lw 15 pt 7 t 'Measured HSFL perp.'

# Reduced(nCrO2) w l lw 3 lc 7 lt 3 t '{/Symbol l}/n[CrO_2]', \

# FIG 2

# fileMaterials u 4:6 w p lw 4 lc 5 t 'Other subs. LSFL', \
# fileMaterials u 4:7 w p lw 4 lc 8 t 'Other subs. HSFL', \
# fileThickness u 4:11 w l lw 4 lc 2 lt 4 t 'TF-SPP Maple', \

set output 'EffectOfSubstrate.eps'
plot wavelength w l lw 3 lc 7 lt 2 t '{/Symbol l}', \
Reduced(nCrO2) w l lw 3 lc 7 lt 3 t '{/Symbol l}/n[CrO_2]', \
Reduced(nCr) w l lw 3 lc 7 lt 4 t '{/Symbol l}/n[Cr]', \
fileMaterials u 4:9 w l lw 4 lc 1 lt 2 t 'TF-SPP {/Symbol l}', \
fileMaterials u 4:10 w l lw 4 lc 7 lt 3 t 'TF-SPP {/Symbol l}/n[Cr]', \
fileMaterials u 4:8 w l lw 4 lc 1 lt 1 t 'TF-SPP coupled', \
fileMaterials u 4:11 w l lw 4 lc 2 lt 4 t 'TF-SPP Maple', \
fileMaterials u 4:6 w p lc 1 lw 5 t 'Exp. LSFL ||', \
fileMaterials u 4:5 w p lc 3 lw 8 t 'Exp. HSFL perp.'
