#!gnuplot

reset

# set output 'MeasuredArea.eps'
# set terminal postscript eps enhanced color font 'Helvetica, 24'

set xlabel 'Angle {/Symbol t} (deg)'
set ylabel 'Intensity (u.a.)'

set xrange [-50:50]
set grid
# set key out center top
unset key

Cyan="CyanZone.csv"
Yellow="YellowZone.csv"
Green="GreenZone.csv"
Norm1=9000
Norm2=9000
Norm3=38000

plot Cyan u 1:($2/Norm1) w l lt 1 lc 7 lw 3 notitle, \
Yellow u 1:($2/Norm2) w l lt 2 lc 1 lw 3 notitle, \
Green u 1:($2/Norm3) w l lt 1 lc rgbcolor "green" lw 3 notitle
