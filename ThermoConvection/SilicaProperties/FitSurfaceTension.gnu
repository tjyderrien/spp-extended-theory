#!gnuplot

reset

file="SurfaceTension.csv"

sigma(T) = a*T+b

fit sigma(x) file u 1:2 via a, b

set xlabel 'T [K]'
set ylabel '{/Symbol s} [N.m]'
plot file u 1:2 w p t 'Data', \
sigma(x) w l t 'Fit'

print a, b
