#!gnuplot

reset

file="Wray1959-Fig3.csv"

kappa(T) = a3*T**3+a2*T**2+a1*T+a0

fit kappa(x) file u 1:2 via a3, a2, a1, a0
set xrange [1:2000]
set xlabel 'T [K]'
set ylabel '{/Symbol k} [W/m/K]'
plot file u 1:2 w p t 'Data', \
kappa(x) w l t 'Fit'

print a3, a2, a1, a0
