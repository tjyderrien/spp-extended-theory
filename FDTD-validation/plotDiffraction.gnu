#!gnuplot

# PRE-TREATMENT
# system "bash ../../calcDiffractionPatternSphere.sh"	

reset


# plot pm3d with Diffraction.dat
set terminal x11 2
set xlabel 'theta index'
set ylabel 'phi index'
splot "Diffraction.dat" u 1:2:3 w p t 'Ex', \
"Diffraction.dat" u 1:2:4 w p t 'Ey', \
"Diffraction.dat" u 1:2:5 w p t 'Ez'


reset
set terminal x11 1

# calculating total scattered energy (summed on all calculated points 
system "echo 'Total intensity out of the simulation box'"
system "echo '/!\\ must contain enough points in resolution to be accurate'"
system "awk -f '../../sumLines.awk' xnff_set001_*"


## plot r(t) for various angles
# plotting brut data

set yrange [1e-10:1]
unset log y

plot "xnff_set001_000_000" u 1:2 w d t '{/Symbol t}=0, {/Symbol f}=0', \
"xnff_set001_000_001" u 1:2 w d t '{/Symbol t}=0, {/Symbol f}=90', \
"xnff_set001_000_002" u 1:2 w d t '{/Symbol t}=0, {/Symbol f}=180', \
"xnff_set001_000_003" u 1:2 w d t '{/Symbol t}=0, {/Symbol f}=270'

