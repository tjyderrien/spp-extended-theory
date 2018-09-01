#!gnuplot

set output 'CrCr2O3CrO2-SPP-result.eps'
set terminal postscript eps enhanced color font 'Helvetica, 24'

wavelength=1026e-9

set colors podo

filename =  "Dostovalov-SPPmodes-CrOxideCompound.csv" # summary of everything
filename0 = "Dostovalov-SPPmodes-CrOxideCompound-Cr-branch0.csv" #branch (k2,k3) --
filename1 = "Dostovalov-SPPmodes-CrOxideCompound-Cr-branch1.csv" # branch -+
filename2 = "Dostovalov-SPPmodes-CrOxideCompound-Cr-branch2.csv" # branch +-
filename3 = "Dostovalov-SPPmodes-CrOxideCompound-Cr-branch3.csv" # branch ++

set xlabel 'Fraction of oxide (%)'
set ylabel 'SPP period (nm)'

set key right outside Left

yscale = 1E9

plot filename0 u ($1*100):(yscale*$6) w p lc 1 ps 0.4 pt 1 t 'Branch -, -', \
     filename1 u ($1*100):(yscale*$6) w p lc 2 ps 0.6 pt 9 t 'Branch -, +', \
     filename2 u ($1*100):(yscale*$6) w p lc 3 ps 0.6 pt 5 t 'Branch +, -', \
     filename3 u ($1*100):(yscale*$6) w p lc 4 ps 0.6 pt 7 t 'Branch +, +', \
     wavelength*yscale w l lc 1 lw 2 dt 3 notitle #t 'Laser wavelength {/Symbol l} (nm)'
     
set output 'CrCr2O3CrO2-SPP-subwavelength.eps'
set yrange [0:1200]
replot

### PRINTING EPSILON(oxide)

set output 'CrCr2O3CrO2-epsilon.eps'
set ylabel 'Dielectric permittivity'
unset yrange
plot filename u ($1*100):($2) w l lc 1 lw 2 t 'Re({/Symbol e})', \
     filename u ($1*100):($3) w l lc 2 lw 2 t 'Im({/Symbol e})'
