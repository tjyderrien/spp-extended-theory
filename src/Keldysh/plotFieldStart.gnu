#!gnuplot
reset

file="StarkEffect-Si-4bands-2.56eV.dat"

set log x

plot file u 1:2 w p notitle, \
     file u 1:3 w p notitle, \
     file u 1:4 w p notitle, \
     file u 1:5 w p notitle, \
     file u 1:6 w p notitle, \
     file u 1:7 w p notitle, \
     file u 1:8 w p notitle
