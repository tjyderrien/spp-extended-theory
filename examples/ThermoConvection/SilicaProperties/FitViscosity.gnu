#!gnuplot
reset
file="Urbain1982-SiO2-viscosity.csv"

lneta(T) = a*T+b

fit lneta(x) file u 1:2 via a, b

a = 6.23888379570223e0
b = -14.6668118241314e0


set terminal x11 1
plot file u 1:2 w p t 'Data', \
lneta(x) w l t 'Fit'

print a, b


###### Now we plot the real viscosity and convert to SI system
set terminal x11 2
eta(T) = exp(1E4*a/T+b)
unitY = 1E1 #10 Poise = 1 Pa.s. 


set xrange [1300:2000] 
set log y
set xlabel 'T [K]'
plot eta(x)/unitY w l t '\Symbol{e}(T) [Pa.s]'
