#!gnuplot

reset

a0 = 0.031
a1 = -48.08e-9
a2 = 29713e-18

e=1.6e-19

Ed = 14.5
Eg = 2.57
E0 = 5.4

n(E) = (1.+E0*Ed/(E0**2-E**2))**0.5
k(lambda)=a0+a1/lambda+a2/lambda**2

c = 3E8
h = 6.63e-34
energy(lambda) = h*c/lambda/e #energy: eV, lambda: m
lambda(energy) = h*c/energy/e #energy: eV, lambda: m 

## PLOT ( energy )
# set xlabel 'Energy (eV)' 
# set xrange [0:]

# plot n(x) w l t 'n', \
# k(lambda(x)) w l t 'k'

## PLOT ( wavelength )

wavelengthmin = 300e-9
wavelengthmax = 2300e-9

set log x

set xlabel 'wavelength (m)'
set xrange [wavelengthmin:wavelengthmax]

# unset log xy

plot n(energy(x)) w l t 'n', \
k(x) w l t 'k'

print "Extrapolation to 1026 nm:"
wavelength0=1026e-9
ncomplexRe = n(energy(wavelength0))
ncomplexIm = k(wavelength0)

ncomplex = {1.,0.}*ncomplexRe + {0.,1.}*ncomplexIm

print "n=", ncomplexRe, ", k=", ncomplexIm
print "epsilon complex:", ncomplex**2
