#!gnuplot

DataFile="ExperimentalError.csv"
Data10k="< awk '{ if($2==10000) print }' ".DataFile
Data2k="< awk '{ if($2==2000) print }' ".DataFile

DataFileJMo="ExperimentalError-OrientationJ-Mo.csv"
DataFileJTi="ExperimentalError-OrientationJ-Ti.csv"

DataMo10kJ="< awk '{ if($2==10000) print }' ".DataFileJMo
DataTi10kJ="< awk '{ if($2==10000) print }' ".DataFileJTi
DataTi2kJ="< awk '{ if($2==2000) print }' ".DataFileJTi

reset
set xrange [1:55]
unset log x
set output 'ErrorBar.eps'
set terminal postscript eps enhanced color font 'Helvetica, 24'

set xlabel 'Number of periodic structures in the probed area'
set ylabel 'Orientation dispersion {/Symbol q} (deg)'

plot Data10k u 3:4:5  w yerrorbars lw 3 lc 1 t '10,000x', \
Data2k u 3:4:5 w yerrorbars lw 3 lc 3 t '2,000x'
#, \
#Data10k u 3:4:5 w l lw 3 lc 1 notitle smooth csplines, \
#Data2k u 3:4:5 w l lw 3 lc 3 notitle smooth csplines 

########################################################

reset

print "Plotting experimental measurement using ImageJ>OrientationJ"

set key bottom right

set terminal postscript eps enhanced color font 'Helvetica, 24'
set output "ErrorBar-OrientationJ.eps"

set xlabel 'Side of the measured square ({/Symbol m}m)'
set ylabel 'Orientation dispersion {/Symbol q} (deg)'

set xrange [1:30]

xscale=1E6

plot DataMo10kJ u ($3*xscale):5 w p lw 3 lc 1 lt 1 notitle, \
DataMo10kJ u ($3*xscale):5 w l smooth sbezier lw 3 lc 1 lt 1 t 'Mo: 10,000x', \
DataTi10kJ u ($3*xscale):5 w p lw 3 lc 1 lt 2 notitle, \
DataTi10kJ u ($3*xscale):5 w l smooth sbezier lw 3 lc 1 lt 2 t 'Ti: 10,000x', \
DataTi2kJ u ($3*xscale):5 w p lw 3 lc 3 lt 2 notitle, \
DataTi2kJ u ($3*xscale):5 w l smooth sbezier lw 3 lc 3 lt 2 t 'Ti: 2,000x'
