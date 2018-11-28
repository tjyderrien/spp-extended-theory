#!gnuplot

set output 'CrCr2O3CrO2-SPP-result.eps'
set terminal postscript eps enhanced color font 'Helvetica, 24'

set size 1.3, 1.0

wavelength=1026e-9

set colors podo

filename =  "Dostovalov-SPPmodes-CrOxideCompound.csv" # summary of everything
filename0 = "Dostovalov-SPPmodes-CrOxideCompound-Cr-branch0.csv" #branch (k2,k3) --
filename1 = "Dostovalov-SPPmodes-CrOxideCompound-Cr-branch1.csv" # branch -+
filename2 = "Dostovalov-SPPmodes-CrOxideCompound-Cr-branch2.csv" # branch +-
filename3 = "Dostovalov-SPPmodes-CrOxideCompound-Cr-branch3.csv" # branch ++
Bonus     = "../../Cr-Cr2O3/Dostovalov-SPPmodes-Cr-Cr2O3-branch0.csv" #branch -- from another calculation
Waveguide  = "Dostovalov-SPPmodes-CrOxideCompound-Cr-WaveGuideAnalysis.csv" #lambda/n

set xlabel 'Fraction of oxide (%)'
set ylabel 'SPP period (nm)'

#unset key 
set key outside right Left reverse

yscale = 1E9
size_of_points = 1.3

set log y
plot filename0 u ($1*100):(yscale*$6) w p lc rgb "red"    ps 1.3 pt 12 t 'Branch -, -', \
     filename1 u ($1*100):(yscale*$6) w p lc rgb "blue"   ps 1.0 pt 9 t 'Branch -, +', \
     filename2 u ($1*100):(yscale*$6) w p lc rgb "web-green"  ps 0.9 pt 5 t 'Branch +, -', \
     filename3 u ($1*100):(yscale*$6) w p lc rgb "orange" ps size_of_points-0.5 pt 7 t 'Branch +, +', \
     Bonus     u ($1*100):(yscale*$6) w p lc rgb "black"  ps size_of_points-0.5 pt 11 t 'Branch -, - (Cr+Cr_2O_3)', \
     Waveguide u ($1*100):(yscale*$6) w l lc rgb "black" lt 3 t 'Waveguiding mode {/Symbol l}/n (nm)',\
     wavelength*yscale w l lc rgb "black" lw 2 dt 3 notitle #t 'Laser wavelength {/Symbol l} (nm)'
unset log y     
set output 'CrCr2O3CrO2-SPP-subwavelength-LessDense.eps'
set yrange [0:1200]
replot

### PRINTING EPSILON(oxide)

set output 'CrCr2O3CrO2-epsilon.eps'
set ylabel 'Dielectric permittivity'
unset yrange
plot filename u ($1*100):($2) w l lc rgb "blue" lw 2 t 'Re({/Symbol e})', \
     filename u ($1*100):($3) w l lc rgb "red"  lw 2 t 'Im({/Symbol e})'
