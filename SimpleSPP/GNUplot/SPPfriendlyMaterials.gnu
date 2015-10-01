#!gnuplot

reset
unset multiplot 

load '/media/thibault/PHD/Travail/Berlin/Calculs/SingleSPPfile.gnu'

epsilonMachine=1e-2

epsilonAir={1e0,0e0}
epsilonAg800={-31.0213e0,0.4095e0}
epsilonAg400={-4.4222e0,0.2103e0}
epsilonTi800={-6.2067e0,25.2004e0}
epsilonTi400={-4.3620e0,12.3621e0}
epsilonSi800={13.64e0,0.048e0}
epsilonSi400={30.8469e0,4.2994e0}
epsilonSi343={18.8384e0,31.6251e0}
epsilonGe800={22.0765e0,3.0347e0}
epsilonGe400={12.2406e0,18.3387e0}
epsilonW800={4.4655e0,19.905e0}
epsilonW400={4.0653e0,16.453e0}
epsilonMg400={0.17730e0,3.5274}**2
epsilonMg800={0.73837,7.6743}**2
############### SPP friendly materials ###################
set output '20141119-SPPfriendlymaterials.eps'
set terminal postscript eps enhanced color font 'Times, 8' size 8cm, 18cm

set ylabel 'Im({/Symbol e})'
set xlabel ''

# set ytics 10
# set xtics 10
# set xrange [-30:30]
set yrange [epsilonMachine:]

unset log x
set log y

set samples 1000000


ImEps2(epsilon1,ReEps2)=(-real(epsilon1)/(imag(epsilon1))) * ReEps2
ImEps22(epsilon1,ReEps2)=(imag(epsilon1)/(real(epsilon1))) * ReEps2

set key left


set multiplot layout 5,2

set title 'Ag, 800 nm'
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -31.0213e0 0.4095e0 800" u 4:($5+epsilonMachine) w p lw 6 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -31.0213e0 0.4095e0 800" u 4:($5+epsilonMachine):1 w labels center offset 0.0,1.0 notitle, \
ImEps2(epsilonAir, x) w l notitle lc 7 lw 2, \
ImEps22(epsilonAir, x) w l notitle lc 3 lw 2, \
ImEps2(epsilonAg800, x) w l notitle lc 1 lw 3, \
ImEps22(epsilonAg800, x) w l notitle lc 4 lw 3
# 

set title 'Ag, 400 nm'
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -4.4222e0 0.2103e0 400" u 4:($5+epsilonMachine) w p lw 6 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -4.4222e0 0.2103e0 400" u 4:($5+epsilonMachine):1 w labels center offset 0.0,1.0 notitle, \
ImEps2(epsilonAir, x) w l notitle lc 7 lw 2, \
ImEps22(epsilonAir, x) w l notitle lc 3 lw 2, \
ImEps2(epsilonAg400, x) w l notitle lc 1 lw 3, \
ImEps22(epsilonAg400, x) w l notitle lc 4 lw 3

set title 'Ti, 800 nm'
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -6.2067e0 25.2004e0 800" u 4:($5+epsilonMachine) w p lw 6 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -6.2067e0 25.2004e0 800" u 4:($5+epsilonMachine):1 w labels center offset 0.0,1.0 notitle, \
ImEps2(epsilonAir, x) w l notitle lc 7 lw 2, \
ImEps22(epsilonAir, x) w l notitle lc 3 lw 2, \
ImEps2(epsilonTi800, x) w l notitle lc 1 lw 3, \
ImEps22(epsilonTi800, x) w l notitle lc 4 lw 3

set title 'Ti, 400 nm'
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -4.3620e0 12.3621e0 400" u 4:($5+epsilonMachine) w p lw 6 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -4.3620e0 12.3621e0 400" u 4:($5+epsilonMachine):1 w labels center offset 0.0,1.0 notitle, \
ImEps2(epsilonAir, x) w l notitle lc 7 lw 2, \
ImEps22(epsilonAir, x) w l notitle lc 3 lw 2, \
ImEps2(epsilonTi400, x) w l notitle lc 1 lw 3, \
ImEps22(epsilonTi400, x) w l notitle lc 4 lw 3

set title 'W, 800 nm'
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 4.4655e0 19.905e0 800" u 4:($5+epsilonMachine) w p lw 6 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 4.4655e0 19.905e0 800" u 4:($5+epsilonMachine):1 w labels center offset 0.0,1.0 notitle, \
ImEps2(epsilonAir, x) w l notitle lc 7 lw 2, \
ImEps22(epsilonAir, x) w l notitle lc 3 lw 2, \
ImEps2(epsilonW800, x) w l notitle lc 1 lw 3, \
ImEps22(epsilonW800, x) w l notitle lc 4 lw 3

set title 'W, 400 nm'
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 4.0653e0 16.453e0 400" u 4:($5+epsilonMachine) w p lw 6 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 4.0653e0 16.453e0 400" u 4:($5+epsilonMachine):1 w labels center offset 0.0,1.0 notitle, \
ImEps2(epsilonAir, x) w l notitle lc 7 lw 2, \
ImEps22(epsilonAir, x) w l notitle lc 3 lw 2, \
ImEps2(epsilonW400, x) w l notitle lc 1 lw 3, \
ImEps22(epsilonW400, x) w l notitle lc 4 lw 3

set xlabel 'Re({/Symbol e})' #Only for last plot

set yrange [epsilonMachine:]

set title 'Ge, 800 nm'
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 22.0765e0 3.0347e0 800" u 4:($5+epsilonMachine) w p lw 6 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 22.0765e0 3.0347e0 800" u 4:($5+epsilonMachine):1 w labels center offset 0.0,1.0 notitle, \
ImEps2(epsilonAir, x) w l notitle lc 7 lw 2, \
ImEps22(epsilonAir, x) w l notitle lc 3 lw 2, \
ImEps2(epsilonGe800, x) w l notitle lc 1 lw 3, \
ImEps22(epsilonGe800, x) w l notitle lc 4 lw 3

set title 'Ge, 400 nm'
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 12.2406e0 18.3387e0 400" u 4:($5+epsilonMachine) w p lw 6 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 12.2406e0 18.3387e0 400" u 4:($5+epsilonMachine):1 w labels center offset 0.0,1.0 notitle, \
ImEps2(epsilonAir, x) w l notitle lc 7 lw 2, \
ImEps22(epsilonAir, x) w l notitle lc 3 lw 2, \
ImEps2(epsilonGe400, x) w l notitle lc 1 lw 3, \
ImEps22(epsilonGe400, x) w l notitle lc 4 lw 3

set title 'Si, 800 nm'
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 13.64e0 0.048e0 800" u 4:($5+epsilonMachine) w p lw 6 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 13.64e0 0.048e0 800" u 4:($5+epsilonMachine):1 w labels center offset 0.0,1.0 notitle, \
ImEps2(epsilonAir, x) w l notitle lc 7 lw 2, \
ImEps22(epsilonAir, x) w l notitle lc 3 lw 2, \
ImEps2(epsilonSi800, x) w l notitle lc 1 lw 3, \
ImEps2(epsilonSi800, x) w l notitle lc 4 lw 3

set title 'Si, 400 nm'
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 30.8469e0 4.2994e0 400" u 4:($5+epsilonMachine) w p lw 6 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 30.8469e0 4.2994e0 400" u 4:($5+epsilonMachine):1 w labels center offset 0.0,1.0 notitle, \
ImEps2(epsilonAir, x) w l notitle lc 7 lw 2, \
ImEps22(epsilonAir, x) w l notitle lc 3 lw 2, \
ImEps2(epsilonSi400, x) w l notitle lc 1 lw 3, \
ImEps22(epsilonSi400, x) w l notitle lc 4 lw 3

unset multiplot



################### SPP period ##############
# unset multiplot
# 
# set output '20141208-SPPfriendlyMaterials-Period.eps'
# set terminal postscript eps enhanced color font 'Times, 8' size 8cm,18cm
# 
# set multiplot layout 4,2
# 
# unset log y
# set xlabel 'Re({/Symbol e})'
# set ylabel 'Period (nm)'
# 
# set title 'Ag, 800 nm'
# plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -31.0213e0 0.4095e0 800" u 4:(1e9*period(beta(800e-9,{-31.0213e0,0.4095e0},$4*Unit+$5*Imaginary))) w p lw 4 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -31.0213e0 0.4095e0 800" u 4:(1e9*period(beta(800e-9,{-31.0213e0,0.4095e0},$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle
# 
# set title 'Ag, 400 nm'
# plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -4.4222e0 0.2103e0 400" u 4:(1e9*period(beta(400e-9,{-4.4222e0,0.2103e0},$4*Unit+$5*Imaginary))) w p lw 4 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -4.4222e0 0.2103e0 400" u 4:(1e9*period(beta(400e-9,{-4.4222e0,0.2103e0},$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle
# 
# set title 'Ti, 800 nm'
# plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -6.2067e0 25.2004e0 800" u 4:(1e9*period(beta(800e-9,{-6.2067e0,25.2004e0},$4*Unit+$5*Imaginary))) w p lw 4 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -6.2067e0 25.2004e0 800" u 4:(1e9*period(beta(800e-9,{-6.2067e0,25.2004e0},$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle
# 
# set title 'Ti, 400 nm'
# plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -4.3620e0 12.3621e0 400" u 4:(1e9*period(beta(400e-9,{-4.3620e0,12.3621e0},$4*Unit+$5*Imaginary))) w p lw 4 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -4.3620e0 12.3621e0 400" u 4:(1e9*period(beta(400e-9,{-4.3620e0,12.3621e0},$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle
# 
# set title 'W, 800 nm'
# plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 4.4655e0 19.905e0 800" u 4:(1e9*period(beta(800e-9,{4.4655e0,19.905e0},$4*Unit+$5*Imaginary))) w p lw 4 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 4.4655e0 19.905e0 800" u 4:(1e9*period(beta(800e-9,{4.4655e0,19.905e0},$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle
# 
# set title 'W, 400 nm'
# plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 4.0653e0 16.453e0 400" u 4:(1e9*period(beta(400e-9,{4.0653e0,16.453e0},$4*Unit+$5*Imaginary))) w p lw 4 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 4.0653e0 16.453e0 400" u 4:(1e9*period(beta(400e-9,{4.0653e0,16.453e0},$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle
# 
# set title 'Ge, 800 nm'
# plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 22.0765e0 3.0347e0 800" u 4:(1e9*period(beta(800e-9,{22.0765e0,3.0347e0},$4*Unit+$5*Imaginary))) w p lw 4 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 22.0765e0 3.0347e0 800" u 4:(1e9*period(beta(800e-9,{22.0765e0,3.0347e0},$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle
# 
# set title 'Ge, 400 nm'
# plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 12.2406e0 18.3387e0 400" u 4:(1e9*period(beta(400e-9,{12.2406e0,18.3387e0},$4*Unit+$5*Imaginary))) w p lw 4 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 12.2406e0 18.3387e0 400" u 4:(1e9*period(beta(400e-9,{12.2406e0,18.3387e0},$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle
# 
# unset multiplot

################ SPP decay depth ###########
# 
# unset multiplot
# 
# set output '20141208-SPPfriendlyMaterials-DecayDepth.eps'
# set terminal postscript eps enhanced color font 'Times, 10' size 8cm,18cm
# 
# set multiplot layout 4,2
# 
# unset log y
# set xlabel 'Re({/Symbol e})'
# set ylabel 'Decay depth (nm)'
# 
# set title 'Ag, 800 nm'
# plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -31.0213e0 0.4095e0 800" u 4:(1e9*DecayDepth(kzSPP(800e-9,epsilonAg800,$4*Unit+$5*Imaginary))) w p lw 4 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -31.0213e0 0.4095e0 800" u 4:(1e9*DecayDepth(kzSPP(800e-9,epsilonAg800,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle
# 
# set title 'Ag, 400 nm'
# plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -4.4222e0 0.2103e0 400" u 4:(1e9*DecayDepth(kzSPP(400e-9,epsilonAg400,$4*Unit+$5*Imaginary))) w p lw 4 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -4.4222e0 0.2103e0 400" u 4:(1e9*DecayDepth(kzSPP(400e-9,epsilonAg400,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle
# 
# set title 'Ti, 800 nm'
# plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -6.2067e0 25.2004e0 800" u 4:(1e9*DecayDepth(kzSPP(800e-9,epsilonTi800,$4*Unit+$5*Imaginary))) w p lw 4 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -6.2067e0 25.2004e0 800" u 4:(1e9*DecayDepth(kzSPP(800e-9,epsilonTi800,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle
# 
# set title 'Ti, 400 nm'
# plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -4.3620e0 12.3621e0 400" u 4:(1e9*DecayDepth(kzSPP(400e-9,epsilonTi400,$4*Unit+$5*Imaginary))) w p lw 4 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -4.3620e0 12.3621e0 400" u 4:(1e9*DecayDepth(kzSPP(400e-9,epsilonTi400,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle
# 
# set title 'W, 800 nm'
# plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 4.4655e0 19.905e0 800" u 4:(1e9*DecayDepth(kzSPP(800e-9,epsilonW800,$4*Unit+$5*Imaginary))) w p lw 4 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 4.4655e0 19.905e0 800" u 4:(1e9*DecayDepth(kzSPP(800e-9,epsilonW800,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle
# 
# set title 'W, 400 nm'
# plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 4.0653e0 16.453e0 400" u 4:(1e9*DecayDepth(kzSPP(400e-9,epsilonW400,$4*Unit+$5*Imaginary))) w p lw 4 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 4.0653e0 16.453e0 400" u 4:(1e9*DecayDepth(kzSPP(400e-9,epsilonW400,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle
# 
# set title 'Ge, 800 nm'
# plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 22.0765e0 3.0347e0 800" u 4:(1e9*DecayDepth(kzSPP(800e-9,epsilonGe800,$4*Unit+$5*Imaginary))) w p lw 4 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 22.0765e0 3.0347e0 800" u 4:(1e9*DecayDepth(kzSPP(800e-9,epsilonGe800,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle
# 
# set title 'Ge, 400 nm'
# plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 12.2406e0 18.3387e0 400" u 4:(1e9*DecayDepth(kzSPP(800e-9,epsilonGe400,$4*Unit+$5*Imaginary))) w p lw 4 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 12.2406e0 18.3387e0 400" u 4:(1e9*DecayDepth(kzSPP(800e-9,epsilonGe400,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle
# 
# unset multiplot
