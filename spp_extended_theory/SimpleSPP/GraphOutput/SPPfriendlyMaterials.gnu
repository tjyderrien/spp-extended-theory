#!gnuplot

reset
unset multiplot 

load '/media/thibault/PHD/Travail/Berlin/Calculs/SingleSPPfile.gnu'

epsilonMachine=1e-2

epsilonAir={1e0,0e0}
epsilonAu800={-24.0620e0,1.5069e0}
epsilonAu400={-1.6580e0, 5.7356e0}
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
epsilonSiO2800={2.1121e0, 0e0}
epsilonSiO2400={2.1612e0, 0e0}
epsilonTiO2800={2.5197e,0e0}**2
epsilonTiO2400={}**2

# ############### SPP friendly materials ###################
# ### HERE WE PLOT IM(EPSILON) AS FUNCTION OF REAL(EPSILON). 
# ### BUT IT IS HARDLY VISIBLE, THOUGH MEANINGFUL. 
# ##########################################################
# 
# 
# ImEps2(epsilon1,ReEps2)=(-real(epsilon1)/(imag(epsilon1))) * ReEps2
# # ImEps22(epsilon1,ReEps2)=(imag(epsilon1)/(real(epsilon1))) * ReEps2
# 
# 
# set yrange [epsilonMachine:]
# 
# unset log x
# set log y
# set log y2
# 
# set samples 1000000
# 
# set key left
# 
# set xlabel 'Re({/Symbol e}_2)'
# set ylabel 'Im({/Symbol e}_2), {/Symbol l}=800 nm' textcolor rgbcolor "red"
# set y2label 'Im({/Symbol e}_2), {/Symbol l}=400 nm' textcolor rgbcolor "blue"
# 
# set ytics nomirror format "%g" textcolor rgbcolor "red"
# set y2tics nomirror format "%g" textcolor rgbcolor "blue"
# 
# # set multiplot layout 2,3
# 
# set output '20141119-SPPfriendlymaterials-Au.eps'
# set terminal postscript eps enhanced color font 'Times, 20' #size 18cm, 8cm
# 
# set title 'Au'
# plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -24.0620e0 1.5069e0 800" u 4:($5+epsilonMachine) w p lw 6 lc 1 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -24.0620e0 1.5069e0 800" u 4:($5+epsilonMachine):1 w labels center offset 0.0,1.0 textcolor rgbcolor 1 notitle, \
# ImEps2(epsilonAu800, x) w l notitle lc 1 lw 3, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -1.6580e0 5.7356e0 400" u 4:($5+epsilonMachine) w p lw 5 lc 3 notitle axis x1y2, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -1.6580e0 5.7356e0 400" u 4:($5+epsilonMachine):1 w labels center offset 0.0,1.0 notitle axis x1y2, \
# ImEps2(epsilonAu400, x) w l notitle lc 3 lw 3 axis x1y2
# 
# 
# set output '20141119-SPPfriendlymaterials-Au-ImagZero.eps'
# set terminal postscript eps enhanced color font 'Times, 20' #size 18cm, 8cm
# 
# set title 'Au'
# plot "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv -24.0620e0 1.5069e0 800" u 4:($5+epsilonMachine) w p lw 6 lc 1 notitle, \
# "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv -24.0620e0 1.5069e0 800" u 4:($5+epsilonMachine):1 w labels center offset 0.0,1.0 textcolor rgbcolor 1 notitle, \
# ImEps2(epsilonAu800, x) w l notitle lc 1 lw 3
# # "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv -1.6580e0 5.7356e0 400" u 4:($5+epsilonMachine) w p lw 5 lc 3 notitle axis x1y2, \
# # "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv -1.6580e0 5.7356e0 400" u 4:($5+epsilonMachine):1 w labels center offset 0.0,1.0 notitle axis x1y2, \
# # ImEps2(epsilonAu400, x) w l notitle lc 3 lw 3 axis x1y2
# 
# 
# set output '20141119-SPPfriendlymaterials-Ag.eps'
# set terminal postscript eps enhanced color font 'Times, 20' #size 18cm, 8cm
# 
# set title 'Ag'
# plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -31.0213e0 0.4095e0 800" u 4:($5+epsilonMachine) w p lw 6 lc 1 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -31.0213e0 0.4095e0 800" u 4:($5+epsilonMachine):1 w labels center offset 0.0,1.0 textcolor rgbcolor 1 notitle, \
# ImEps2(epsilonAg800, x) w l notitle lc 1 lw 3, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -4.4222e0 0.2103e0 400" u 4:($5+epsilonMachine) w p lw 5 lc 3 notitle axis x1y2, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -4.4222e0 0.2103e0 400" u 4:($5+epsilonMachine):1 w labels center offset 0.0,1.0 notitle axis x1y2, \
# ImEps2(epsilonAg400, x) w l notitle lc 3 lw 3 axis x1y2
# 
# 
# set output '20141119-SPPfriendlymaterials-Ag-ImagZero.eps'
# set terminal postscript eps enhanced color font 'Times, 20' #size 18cm, 8cm
# 
# set title 'Ag'
# plot "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv -31.0213e0 0.4095e0 800" u 4:($5+epsilonMachine) w p lw 6 lc 1 notitle, \
# "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv -31.0213e0 0.4095e0 800" u 4:($5+epsilonMachine):1 w labels center offset 0.0,1.0 textcolor rgbcolor 1 notitle, \
# ImEps2(epsilonAg800, x) w l notitle lc 1 lw 3, \
# "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv -4.4222e0 0.2103e0 400" u 4:($5+epsilonMachine) w p lw 5 lc 3 notitle axis x1y2, \
# "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv -4.4222e0 0.2103e0 400" u 4:($5+epsilonMachine):1 w labels center offset 0.0,1.0 notitle axis x1y2, \
# ImEps2(epsilonAg400, x) w l notitle lc 3 lw 3 axis x1y2
# 
# 
# set output '20141119-SPPfriendlymaterials-Ti.eps'
# set terminal postscript eps enhanced color font 'Times, 20' #size 18cm, 8cm
# 
# set title 'Ti'
# plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -6.2067e0 25.2004e0 800" u 4:($5+epsilonMachine) w p lw 6 lc 1 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -6.2067e0 25.2004e0 800" u 4:($5+epsilonMachine):1 w labels center offset 0.0,1.0 notitle, \
# ImEps2(epsilonTi800, x) w l notitle lc 1 lw 3, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -4.3620e0 12.3621e0 400" u 4:($5+epsilonMachine) w p lw 5 lc 3 notitle axis x1y2, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -4.3620e0 12.3621e0 400" u 4:($5+epsilonMachine):1 w labels center offset 0.0,1.0 notitle axis x1y2, \
# ImEps2(epsilonTi400, x) w l notitle lc 3 lw 3 axis x1y2
# 
# 
# set output '20141119-SPPfriendlymaterials-Ti-ImagZero.eps'
# set terminal postscript eps enhanced color font 'Times, 20' #size 18cm, 8cm
# 
# set title 'Ti'
# plot "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv -6.2067e0 25.2004e0 800" u 4:($5+epsilonMachine) w p lw 6 lc 1 notitle, \
# "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv -6.2067e0 25.2004e0 800" u 4:($5+epsilonMachine):1 w labels center offset 0.0,1.0 notitle, \
# ImEps2(epsilonTi800, x) w l notitle lc 1 lw 3, \
# "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv -4.3620e0 12.3621e0 400" u 4:($5+epsilonMachine) w p lw 5 lc 3 notitle axis x1y2, \
# "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv -4.3620e0 12.3621e0 400" u 4:($5+epsilonMachine):1 w labels center offset 0.0,1.0 notitle axis x1y2, \
# ImEps2(epsilonTi400, x) w l notitle lc 3 lw 3 axis x1y2
# 
# 
# set output '20141119-SPPfriendlymaterials-W.eps'
# set terminal postscript eps enhanced color font 'Times, 20' #size 18cm, 8cm
# 
# set title 'W'
# plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 4.4655e0 19.905e0 800" u 4:($5+epsilonMachine) w p lw 6 lc 1 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 4.4655e0 19.905e0 800" u 4:($5+epsilonMachine):1 w labels center offset 0.0,1.0 notitle, \
# ImEps2(epsilonW800, x) w l notitle lc 1 lw 3, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 4.0653e0 16.453e0 400" u 4:($5+epsilonMachine) w p lw 5 lc 3 notitle axis x1y2, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 4.0653e0 16.453e0 400" u 4:($5+epsilonMachine):1 w labels center offset 0.0,1.0 notitle axis x1y2, \
# ImEps2(epsilonW400, x) w l notitle lc 3 lw 3 axis x1y2
#  
# 
# set output '20141119-SPPfriendlymaterials-W-ImagZero.eps'
# set terminal postscript eps enhanced color font 'Times, 20' #size 18cm, 8cm
# 
# set title 'W'
# plot "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv 4.4655e0 19.905e0 800" u 4:($5+epsilonMachine) w p lw 6 lc 1 notitle, \
# "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv 4.4655e0 19.905e0 800" u 4:($5+epsilonMachine):1 w labels center offset 0.0,1.0 notitle, \
# ImEps2(epsilonW800, x) w l notitle lc 1 lw 3, \
# "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv 4.0653e0 16.453e0 400" u 4:($5+epsilonMachine) w p lw 5 lc 3 notitle axis x1y2, \
# "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv 4.0653e0 16.453e0 400" u 4:($5+epsilonMachine):1 w labels center offset 0.0,1.0 notitle axis x1y2, \
# ImEps2(epsilonW400, x) w l notitle lc 3 lw 3 axis x1y2
#  
# 
# set output '20141119-SPPfriendlymaterials-Ge.eps'
# set terminal postscript eps enhanced color font 'Times, 20' #size 18cm, 8cm
# 
# set title 'Ge'
# plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 22.0765e0 3.0347e0 800" u 4:($5+epsilonMachine) w p lw 6 lc 1 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 22.0765e0 3.0347e0 800" u 4:($5+epsilonMachine):1 w labels center offset 0.0,1.0 notitle, \
# ImEps2(epsilonGe800, x) w l lc 1 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 12.2406e0 18.3387e0 400" u 4:($5+epsilonMachine) w p lw 5 lc 3 notitle axis x1y2, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 12.2406e0 18.3387e0 400" u 4:($5+epsilonMachine):1 w labels center offset 0.0,1.0 notitle axis x1y2, \
# ImEps2(epsilonGe400, x) w l notitle lw 3 lc 3 axis x1y2
# 
# 
# set output '20141119-SPPfriendlymaterials-Ge-ImagZero.eps'
# set terminal postscript eps enhanced color font 'Times, 20' #size 18cm, 8cm
# 
# set title 'Ge'
# plot "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv 22.0765e0 3.0347e0 800" u 4:($5+epsilonMachine) w p lw 6 lc 1 notitle, \
# "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv 22.0765e0 3.0347e0 800" u 4:($5+epsilonMachine):1 w labels center offset 0.0,1.0 notitle, \
# ImEps2(epsilonGe800, x) w l lc 1 notitle, \
# "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv 12.2406e0 18.3387e0 400" u 4:($5+epsilonMachine) w p lw 5 lc 3 notitle axis x1y2, \
# "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv 12.2406e0 18.3387e0 400" u 4:($5+epsilonMachine):1 w labels center offset 0.0,1.0 notitle axis x1y2, \
# ImEps2(epsilonGe400, x) w l notitle lw 3 lc 3 axis x1y2
# 
# 
# # 
# # unset multiplot
# # 
# # set output '20141119-SPPfriendlymaterials2.eps'
# # set terminal postscript eps enhanced color font 'Times, 8' size 8cm, 8cm
# # 
# # set multiplot layout 2,3
# # 
# # set yrange [epsilonMachine:]
# # set ylabel '' #Im({/Symbol e})'
# 
# set output '20141119-SPPfriendlymaterials-Si.eps'
# set terminal postscript eps enhanced color font 'Times, 20' #size 18cm, 8cm
# 
# set title 'Si'
# plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 13.64e0 0.048e0 800" u 4:($5+epsilonMachine) w p lw 6 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 13.64e0 0.048e0 800" u 4:($5+epsilonMachine):1 w labels center offset 0.0,1.0 notitle, \
# ImEps2(epsilonSi800, x) w l notitle lc 1 lw 3, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 30.8469e0 4.2994e0 400" u 4:($5+epsilonMachine) w p lw 6 lc 3 notitle axis x1y2, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 30.8469e0 4.2994e0 400" u 4:($5+epsilonMachine):1 w labels center offset 0.0,1.0 notitle axis x1y2, \
# ImEps2(epsilonSi400, x) w l notitle lc 3 lw 3 axis x1y2
# 
# 
# set output '20141119-SPPfriendlymaterials-Si-ImagZero.eps'
# set terminal postscript eps enhanced color font 'Times, 20' #size 18cm, 8cm
# 
# set title 'Si'
# plot "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv 13.64e0 0.048e0 800" u 4:($5+epsilonMachine) w p lw 6 notitle, \
# "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv 13.64e0 0.048e0 800" u 4:($5+epsilonMachine):1 w labels center offset 0.0,1.0 notitle, \
# ImEps2(epsilonSi800, x) w l notitle lc 1 lw 3
# # "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv 30.8469e0 4.2994e0 400" u 4:($5+epsilonMachine) w p lw 6 lc 3 notitle axis x1y2, \
# # "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv 30.8469e0 4.2994e0 400" u 4:($5+epsilonMachine):1 w labels center offset 0.0,1.0 notitle axis x1y2, \
# # ImEps2(epsilonSi400, x) w l notitle lc 3 lw 3 axis x1y2
# 
# 
# set output '20141119-SPPfriendlymaterials-SiO2.eps'
# set terminal postscript eps enhanced color font 'Times, 20' #size 18cm, 8cm
# 
# set xrange [-8:0]
# 
# set title 'SiO_2'
# plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 2.1121e0 0.000e0 800" u 4:($5+epsilonMachine) w p lw 6 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 2.1121e0 0.000e0 800" u 4:($5+epsilonMachine):1 w labels center offset 0.0,1.0 notitle, \
# ImEps2(epsilonSiO2800, x) w l notitle lc 1 lw 3, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 2.1612e0 0e0 400" u 4:($5+epsilonMachine) w p lw 6 lc 3 notitle axis x1y2, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 2.1612e0 0e0 400" u 4:($5+epsilonMachine):1 w labels center offset 0.0,1.0 notitle axis x1y2, \
# ImEps2(epsilonSiO2400, x) w l notitle lc 3 lw 3 axis x1y2
# 
# set output '20141119-SPPfriendlymaterials-SiO2-ImagZero.eps'
# set terminal postscript eps enhanced color font 'Times, 20' #size 18cm, 8cm
# 
# set title 'SiO_2'
# plot "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv 2.1121e0 0.000e0 800" u 4:($5+epsilonMachine) w p lw 6 notitle, \
# "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv 2.1121e0 0.000e0 800" u 4:($5+epsilonMachine):1 w labels center offset 0.0,1.0 notitle, \
# ImEps2(epsilonSiO2800, x) w l notitle lc 1 lw 3, \
# "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv 2.1612e0 0e0 400" u 4:($5+epsilonMachine) w p lw 6 lc 3 notitle axis x1y2, \
# "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv 2.1612e0 0e0 400" u 4:($5+epsilonMachine):1 w labels center offset 0.0,1.0 notitle axis x1y2, \
# ImEps2(epsilonSiO2400, x) w l notitle lc 3 lw 3 axis x1y2
# 
# 
# set output '20141119-SPPfriendlymaterials-Air.eps'
# set terminal postscript eps enhanced color font 'Times, 20' # size 8cm, 8cm
# 
# set title 'Air'
# plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 1e0 0.000e0 800" u 4:($5+epsilonMachine) w p lw 6 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 1e0 0.000e0 800" u 4:($5+epsilonMachine):1 w labels center offset 0.0,1.0 notitle, \
# ImEps2(epsilonAir, x) w l notitle lc 1 lw 3, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 1e0 0e0 400" u 4:($5+epsilonMachine) w p lw 6 lc 3 notitle axis x1y2, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 1e0 0e0 400" u 4:($5+epsilonMachine):1 w labels center offset 0.0,1.0 notitle axis x1y2, \
# ImEps2(epsilonAir, x) w l notitle lc 3 lw 3 axis x1y2
# 
# 
# set output '20141119-SPPfriendlymaterials-Air-ImagZero.eps'
# set terminal postscript eps enhanced color font 'Times, 20' # size 8cm, 8cm
# 
# set title 'Air'
# plot "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv 1e0 0.000e0 800" u 4:($5+epsilonMachine) w p lw 6 notitle, \
# "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv 1e0 0.000e0 800" u 4:($5+epsilonMachine):1 w labels center offset 0.0,1.0 notitle, \
# "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv 1e0 0e0 400" u 4:($5+epsilonMachine) w p lw 6 lc 3 notitle axis x1y2, \
# "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv 1e0 0e0 400" u 4:($5+epsilonMachine):1 w labels center offset 0.0,1.0 notitle axis x1y2
# set xrange [*:*]
# # unset multiplot
# 
# ########################################################
# ### SAME BUT WITH ATOMIC NUMBERS AND SPP CONDITION VALUE
# ### MORE VISIBLE, THOUGH PHYSICAL MEANING IS LESS CLEAR
# ########################################################
# ############### SPP friendly materials ###################
# set grid
# 
# set xlabel 'Atomic mass (uam)' 
# set ylabel 'f({/Symbol e}_1, {/Symbol e}_2), {/Symbol l}=800 nm' textcolor rgbcolor "red"
# set y2label 'f({/Symbol e}_1, {/Symbol e}_2), {/Symbol l}=400 nm' textcolor rgbcolor "blue"
# 
# unset log xy
# unset log y2
# 
# set xrange [*:*]
# set yrange [*:*]
# set y2range [*:*]
# 
# set key left
# # set xlabel 'Re({/Symbol e})' 
# 
# # set multiplot layout 2,3
# 
# set output '20141119-SPPfriendlymaterials-ConditionWithAtomicNumber-Au.eps'
# set terminal postscript eps enhanced color font 'Times, 20' # size 8cm, 8cm
# 
# set title 'Au'
# plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -24.0620e0 1.5069e0 800" u 10:($11) w p lw 6 lc 1 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -24.0620e0 1.5069e0 800" u 10:($11):1 w labels center offset 0.0,1.0 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -1.6580e0 5.7356e0 400" u 10:($11) w p lw 5 lc 3 axis x1y2 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -1.6580e0 5.7356e0 400" u 10:($11):1 w labels center offset 0.0,1.0 axis x1y2 notitle
# 
# set output '20141119-SPPfriendlymaterials-ConditionWithAtomicNumber-Ag.eps'
# set terminal postscript eps enhanced color font 'Times, 20' # size 8cm, 8cm
# 
# set title 'Ag'
# plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -31.0213e0 0.4095e0 800" u 10:($11) w p lw 6 lc 1 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -31.0213e0 0.4095e0 800" u 10:($11):1 w labels center offset 0.0,1.0 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -4.4222e0 0.2103e0 400" u 10:($11) w p lw 5 lc 3 axis x1y2 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -4.4222e0 0.2103e0 400" u 10:($11):1 w labels center offset 0.0,1.0 axis x1y2 notitle
# 
# set output '20141119-SPPfriendlymaterials-ConditionWithAtomicNumber-Ti.eps'
# set terminal postscript eps enhanced color font 'Times, 20' # size 8cm, 8cm
# 
# 
# set title 'Ti'
# plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -6.2067e0 25.2004e0 800" u 10:($11) w p lw 6 lc 1 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -6.2067e0 25.2004e0 800" u 10:($11):1 w labels center offset 0.0,1.0 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -4.3620e0 12.3621e0 400" u 10:($11) w p lw 5 lc 3 notitle axis x1y2, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -4.3620e0 12.3621e0 400" u 10:($11):1 w labels center offset 0.0,1.0 axis x1y2 notitle
# # 
# # unset multiplot
# # 
# # set output '20141119-SPPfriendlymaterials-ConditionWithAtomicNumber2.eps'
# # set terminal postscript eps enhanced color font 'Times, 8' size 8cm, 8cm
# # 
# # set multiplot layout 2,3
# 
# set output '20141119-SPPfriendlymaterials-ConditionWithAtomicNumber-W.eps'
# set terminal postscript eps enhanced color font 'Times, 20' # size 8cm, 8cm
# 
# set title 'W'
# plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 4.4655e0 19.905e0 800" u 10:($11) w p lw 6 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 4.4655e0 19.905e0 800" u 10:($11):1 w labels center offset 0.0,1.0 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 4.0653e0 16.453e0 400" u 10:($11) w p lw 5 lc 3 notitle axis x1y2, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 4.0653e0 16.453e0 400" u 10:($11):1 w labels center offset 0.0,1.0 notitle axis x1y2
# 
# set output '20141119-SPPfriendlymaterials-ConditionWithAtomicNumber-Ge.eps'
# set terminal postscript eps enhanced color font 'Times, 20' # size 8cm, 8cm
# 
# set title 'Ge'
# plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 22.0765e0 3.0347e0 800" u 10:($11) w p lw 6 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 22.0765e0 3.0347e0 800" u 10:($11):1 w labels center offset 0.0,1.0 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 12.2406e0 18.3387e0 400" u 10:($11) w p lw 5 lc 3 notitle axis x1y2, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 12.2406e0 18.3387e0 400" u 10:($11):1 w labels center offset 0.0,1.0 axis x1y2 notitle
# 
# set output '20141119-SPPfriendlymaterials-ConditionWithAtomicNumber-Si.eps'
# set terminal postscript eps enhanced color font 'Times, 20' # size 8cm, 8cm
# 
# set title 'Si'
# plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 13.64e0 0.048e0 800" u 10:($11) w p lw 6 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 13.64e0 0.048e0 800" u 10:($11):1 w labels center offset 0.0,1.0 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 30.8469e0 4.2994e0 400" u 10:($11) w p lw 6 lc 3 notitle axis x1y2, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 30.8469e0 4.2994e0 400" u 10:($11):1 w labels center offset 0.0,1.0 notitle axis x1y2
# 
# set output '20141119-SPPfriendlymaterials-ConditionWithAtomicNumber-SiO2.eps'
# set terminal postscript eps enhanced color font 'Times, 20' # size 8cm, 8cm
# 
# set title 'SiO_2'
# plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 2.1121e0 0.000e0 800" u 10:($11) w p lw 6 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 2.1121e0 0.000e0 800" u 10:($11):1 w labels center offset 0.0,1.0 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 2.1612e0 0e0 400" u 10:($11) w p lw 6 lc 3 notitle axis x1y2, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 2.1612e0 0e0 400" u 10:($11):1 w labels center offset 0.0,1.0 notitle axis x1y2
# 
# set output '20141119-SPPfriendlymaterials-ConditionWithAtomicNumber-Air.eps'
# set terminal postscript eps enhanced color font 'Times, 20' # size 8cm, 8cm
# 
# set title 'Air'
# plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 1e0 0.000e0 800" u 10:($11) w p lw 6 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 1e0 0.000e0 800" u 10:($11):1 w labels center offset 0.0,1.0 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 1e0 0e0 400" u 10:($11) w p lw 6 lc 3 notitle axis x1y2, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 1e0 0e0 400" u 10:($11):1 w labels center offset 0.0,1.0 notitle axis x1y2
# 
# # unset multiplot
# 
# 
# ####################################################################
# ### SAME BUT WITH BAND GAP OR FERMI ENERGY 
# ### MORE VISIBLE, PHYSICALLY MEANING, but some materials are missing
# ####################################################################
# set grid
# 
# set xlabel 'Fermi Energy (eV < 0), Band gap energy (> 0)' 
# set ylabel 'f({/Symbol e}_1, {/Symbol e}_2), {/Symbol l}=800 nm' textcolor rgbcolor "red"
# set y2label 'f({/Symbol e}_1, {/Symbol e}_2), {/Symbol l}=400 nm' textcolor rgbcolor "blue"
# 
# unset log xy
# unset log y2
# 
# set xrange [*:*]
# set yrange [*:*]
# set y2range [*:*]
# 
# set key left
# # set xlabel 'Re({/Symbol e})' 
# 
# # set multiplot layout 2,3
# 
# set output '20141119-SPPfriendlymaterials-ConditionWithBandGap-Au.eps'
# set terminal postscript eps enhanced color font 'Times, 20' # size 8cm, 8cm
# 
# set title 'Au'
# plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -24.0620e0 1.5069e0 800" u 2:($11) w p lw 6 lc 1 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -24.0620e0 1.5069e0 800" u 2:($11):1 w labels center offset 0.0,1.0 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -1.6580e0 5.7356e0 400" u 2:($11) w p lw 5 lc 3 axis x1y2 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -1.6580e0 5.7356e0 400" u 2:($11):1 w labels center offset 0.0,1.0 axis x1y2 notitle
# 
# set output '20141119-SPPfriendlymaterials-ConditionWithBandGap-Ag.eps'
# set terminal postscript eps enhanced color font 'Times, 20' # size 8cm, 8cm
# 
# set title 'Ag'
# plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -31.0213e0 0.4095e0 800" u 2:($11) w p lw 6 lc 1 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -31.0213e0 0.4095e0 800" u 2:($11):1 w labels center offset 0.0,1.0 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -4.4222e0 0.2103e0 400" u 2:($11) w p lw 5 lc 3 axis x1y2 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -4.4222e0 0.2103e0 400" u 2:($11):1 w labels center offset 0.0,1.0 axis x1y2 notitle
# 
# set output '20141119-SPPfriendlymaterials-ConditionWithBandGap-Ti.eps'
# set terminal postscript eps enhanced color font 'Times, 20' # size 8cm, 8cm
# 
# 
# set title 'Ti'
# plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -6.2067e0 25.2004e0 800" u 2:($11) w p lw 6 lc 1 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -6.2067e0 25.2004e0 800" u 2:($11):1 w labels center offset 0.0,1.0 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -4.3620e0 12.3621e0 400" u 2:($11) w p lw 5 lc 3 notitle axis x1y2, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -4.3620e0 12.3621e0 400" u 2:($11):1 w labels center offset 0.0,1.0 axis x1y2 notitle
# # 
# # unset multiplot
# # 
# # set output '20141119-SPPfriendlymaterials-ConditionWithBandGap2.eps'
# # set terminal postscript eps enhanced color font 'Times, 8' size 8cm, 8cm
# # 
# # set multiplot layout 2,3
# 
# set output '20141119-SPPfriendlymaterials-ConditionWithBandGap-W.eps'
# set terminal postscript eps enhanced color font 'Times, 20' # size 8cm, 8cm
# 
# set title 'W'
# plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 4.4655e0 19.905e0 800" u 2:($11) w p lw 6 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 4.4655e0 19.905e0 800" u 2:($11):1 w labels center offset 0.0,1.0 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 4.0653e0 16.453e0 400" u 2:($11) w p lw 5 lc 3 notitle axis x1y2, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 4.0653e0 16.453e0 400" u 2:($11):1 w labels center offset 0.0,1.0 notitle axis x1y2
# 
# set output '20141119-SPPfriendlymaterials-ConditionWithBandGap-Ge.eps'
# set terminal postscript eps enhanced color font 'Times, 20' # size 8cm, 8cm
# 
# set title 'Ge'
# plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 22.0765e0 3.0347e0 800" u 2:($11) w p lw 6 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 22.0765e0 3.0347e0 800" u 2:($11):1 w labels center offset 0.0,1.0 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 12.2406e0 18.3387e0 400" u 2:($11) w p lw 5 lc 3 notitle axis x1y2, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 12.2406e0 18.3387e0 400" u 2:($11):1 w labels center offset 0.0,1.0 axis x1y2 notitle
# 
# set output '20141119-SPPfriendlymaterials-ConditionWithBandGap-Si.eps'
# set terminal postscript eps enhanced color font 'Times, 20' # size 8cm, 8cm
# 
# set title 'Si'
# plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 13.64e0 0.048e0 800" u 2:($11) w p lw 6 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 13.64e0 0.048e0 800" u 2:($11):1 w labels center offset 0.0,1.0 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 30.8469e0 4.2994e0 400" u 2:($11) w p lw 6 lc 3 notitle axis x1y2, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 30.8469e0 4.2994e0 400" u 2:($11):1 w labels center offset 0.0,1.0 notitle axis x1y2
# 
# set output '20141119-SPPfriendlymaterials-ConditionWithBandGap-SiO2.eps'
# set terminal postscript eps enhanced color font 'Times, 20' # size 8cm, 8cm
# 
# set title 'SiO_2'
# plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 2.1121e0 0.000e0 800" u 2:($11) w p lw 6 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 2.1121e0 0.000e0 800" u 2:($11):1 w labels center offset 0.0,1.0 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 2.1612e0 0e0 400" u 2:($11) w p lw 6 lc 3 notitle axis x1y2, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 2.1612e0 0e0 400" u 2:($11):1 w labels center offset 0.0,1.0 notitle axis x1y2
# 
# 
# set output '20141119-SPPfriendlymaterials-ConditionWithBandGap-Air.eps'
# set terminal postscript eps enhanced color font 'Times, 20' # size 8cm, 8cm
# 
# set title 'Air'
# plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 1e0 0.000e0 800" u 2:($11) w p lw 6 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 1e0 0.000e0 800" u 2:($11):1 w labels center offset 0.0,1.0 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 1e0 0e0 400" u 2:($11) w p lw 6 lc 3 notitle axis x1y2, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 1e0 0e0 400" u 2:($11):1 w labels center offset 0.0,1.0 notitle axis x1y2
# 
# 
# ####################################################################
# ### SAME BUT USING Neglecting Imag parts CONDITIONs, as a function of BAND GAP OR FERMI ENERGY 
# ####################################################################
# set grid
# 
# set xlabel 'Fermi Energy (eV < 0), Band gap energy (> 0)' 
# set ylabel 'f, Im({/Symbol e}_i)=0, {/Symbol l}=800 nm' textcolor rgbcolor "red"
# set y2label 'f, Im({/Symbol e}_i)=0, {/Symbol e}_2), {/Symbol l}=400 nm' textcolor rgbcolor "blue"
# 
# unset log xy
# unset log y2
# 
# set xrange [*:*]
# set yrange [*:*]
# set y2range [*:*]
# 
# set key left
# # set xlabel 'Re({/Symbol e})' 
# 
# # set multiplot layout 2,3
# 
# set output '20141119-SPPfriendlymaterials-ConditionWithBandGap-Au-ImagZero.eps'
# set terminal postscript eps enhanced color font 'Times, 20' # size 8cm, 8cm
# 
# set title 'Au'
# plot "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv -24.0620e0 1.5069e0 800" u 2:($11) w p lw 6 lc 1 notitle, \
# "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv -24.0620e0 1.5069e0 800" u 2:($11):1 w labels center offset 0.0,1.0 notitle
# # "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv -1.6580e0 5.7356e0 400" u 2:($11) w p lw 5 lc 3 axis x1y2 notitle, \
# # "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv -1.6580e0 5.7356e0 400" u 2:($11):1 w labels center offset 0.0,1.0 axis x1y2 notitle
# 
# set output '20141119-SPPfriendlymaterials-ConditionWithBandGap-Ag-ImagZero.eps'
# set terminal postscript eps enhanced color font 'Times, 20' # size 8cm, 8cm
# 
# set title 'Ag'
# plot "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv -31.0213e0 0.4095e0 800" u 2:($11) w p lw 6 lc 1 notitle, \
# "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv -31.0213e0 0.4095e0 800" u 2:($11):1 w labels center offset 0.0,1.0 notitle, \
# "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv -4.4222e0 0.2103e0 400" u 2:($11) w p lw 5 lc 3 axis x1y2 notitle, \
# "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv -4.4222e0 0.2103e0 400" u 2:($11):1 w labels center offset 0.0,1.0 axis x1y2 notitle
# 
# set output '20141119-SPPfriendlymaterials-ConditionWithBandGap-Ti-ImagZero.eps'
# set terminal postscript eps enhanced color font 'Times, 20' # size 8cm, 8cm
# 
# set title 'Ti'
# plot "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv -6.2067e0 25.2004e0 800" u 2:($11) w p lw 6 lc 1 notitle, \
# "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv -6.2067e0 25.2004e0 800" u 2:($11):1 w labels center offset 0.0,1.0 notitle, \
# "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv -4.3620e0 12.3621e0 400" u 2:($11) w p lw 5 lc 3 notitle axis x1y2, \
# "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv -4.3620e0 12.3621e0 400" u 2:($11):1 w labels center offset 0.0,1.0 axis x1y2 notitle
# # 
# # unset multiplot
# # 
# # set output '20141119-SPPfriendlymaterials-ConditionWithBandGap2.eps'
# # set terminal postscript eps enhanced color font 'Times, 8' size 8cm, 8cm
# # 
# # set multiplot layout 2,3
# 
# set output '20141119-SPPfriendlymaterials-ConditionWithBandGap-W-ImagZero.eps'
# set terminal postscript eps enhanced color font 'Times, 20' # size 8cm, 8cm
# 
# set title 'W'
# plot "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv 4.4655e0 19.905e0 800" u 2:($11) w p lw 6 notitle, \
# "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv 4.4655e0 19.905e0 800" u 2:($11):1 w labels center offset 0.0,1.0 notitle, \
# "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv 4.0653e0 16.453e0 400" u 2:($11) w p lw 5 lc 3 notitle axis x1y2, \
# "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv 4.0653e0 16.453e0 400" u 2:($11):1 w labels center offset 0.0,1.0 notitle axis x1y2
# 
# set output '20141119-SPPfriendlymaterials-ConditionWithBandGap-Ge-ImagZero.eps'
# set terminal postscript eps enhanced color font 'Times, 20' # size 8cm, 8cm
# 
# set title 'Ge'
# plot "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv 22.0765e0 3.0347e0 800" u 2:($11) w p lw 6 notitle, \
# "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv 22.0765e0 3.0347e0 800" u 2:($11):1 w labels center offset 0.0,1.0 notitle, \
# "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv 12.2406e0 18.3387e0 400" u 2:($11) w p lw 5 lc 3 notitle axis x1y2, \
# "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv 12.2406e0 18.3387e0 400" u 2:($11):1 w labels center offset 0.0,1.0 axis x1y2 notitle
# 
# set output '20141119-SPPfriendlymaterials-ConditionWithBandGap-Si-ImagZero.eps'
# set terminal postscript eps enhanced color font 'Times, 20' # size 8cm, 8cm
# 
# set title 'Si'
# plot "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv 13.64e0 0.048e0 800" u 2:($11) w p lw 6 notitle, \
# "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv 13.64e0 0.048e0 800" u 2:($11):1 w labels center offset 0.0,1.0 notitle
# # "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv 30.8469e0 0.048e0 400" u 2:($11) w p lw 6 lc 3 notitle axis x1y2, \
# # "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv 30.8469e0 0.048e0 400" u 2:($11):1 w labels center offset 0.0,1.0 notitle axis x1y2
# 
# set output '20141119-SPPfriendlymaterials-ConditionWithBandGap-SiO2-ImagZero.eps'
# set terminal postscript eps enhanced color font 'Times, 20' # size 8cm, 8cm
# 
# set title 'SiO_2'
# plot "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv 2.1121e0 0.000e0 800" u 2:($11) w p lw 6 notitle, \
# "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv 2.1121e0 0.000e0 800" u 2:($11):1 w labels center offset 0.0,1.0 notitle, \
# "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv 2.1612e0 0e0 400" u 2:($11) w p lw 6 lc 3 notitle axis x1y2, \
# "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv 2.1612e0 0e0 400" u 2:($11):1 w labels center offset 0.0,1.0 notitle axis x1y2
# 
# 
# set output '20141119-SPPfriendlymaterials-ConditionWithBandGap-Air-ImagZero.eps'
# set terminal postscript eps enhanced color font 'Times, 20' # size 8cm, 8cm
# 
# set title 'Air'
# plot "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv 1e0 0.000e0 800" u 2:($11) w p lw 6 notitle, \
# "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv 1e0 0.000e0 800" u 2:($11):1 w labels center offset 0.0,1.0 notitle, \
# "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv 1e0 0e0 400" u 2:($11) w p lw 6 lc 3 notitle axis x1y2, \
# "< sh SPPfriendlySelect-ImagZero.sh MaterialOpticalDatabaseForPlasmonics.csv 1e0 0e0 400" u 2:($11):1 w labels center offset 0.0,1.0 notitle axis x1y2
# 

#############################################
################### SPP period ##############
####### Plotting as function of Re(eps) #####
unset y2tics 

set ytics format "%g" textcolor rgbcolor "black"
set ylabel textcolor rgbcolor "black"
set y2label ''

set yrange [0:900]

unset multiplot
unset grid
# set key top right
# set multiplot layout 2,2

unset log y
set xlabel 'Re({/Symbol e})'
set ylabel 'Field period (nm)'
# 
# set output '20141208-SPPfriendlyMaterials-Period-Au.eps'
# set terminal postscript eps enhanced color font 'Helvetica, 22' size 10cm, 8cm
# 
# set xrange [0:35]
# 
# set title 'Au'
# 
# set key top right
# plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -24.0620e0 1.5069e0 800" u 4:(1e9*period(beta(800e-9,epsilonAu800,$4*Unit+$5*Imaginary))) w p lc 1 lw 4 t '{/Symbol l}=800 nm', \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -24.0620e0 1.5069e0 800" u 4:(1e9*period(beta(800e-9,epsilonAu800,$4*Unit+$5*Imaginary))):1 w labels center offset 0.0,1.0 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -1.6580e0 5.7356e0 400" u 4:(1e9*period(beta(400e-9,epsilonAu400,$4*Unit+$5*Imaginary))) w p lw 4 lc 3 t '{/Symbol l}=400 nm', \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -1.6580e0 5.7356e0 400" u 4:(1e9*period(beta(400e-9,epsilonAu400,$4*Unit+$5*Imaginary))):1 w labels center offset 0.0,1.0 notitle
# 
# set output '20141208-SPPfriendlyMaterials-Period-Au-ImagZero.eps'
# set terminal postscript eps enhanced color font 'Helvetica, 22' size 10cm, 8cm
# 
# set xrange [0:35]
# 
# set title 'Au'
# plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -24.0620e0 1.5069e0 800" u 4:(1e9*period(beta(800e-9,real(epsilonAu800),$4))) w p lc 1 lw 4 t '{/Symbol l}=800 nm', \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -24.0620e0 1.5069e0 800" u 4:(1e9*period(beta(800e-9,real(epsilonAu800),$4))):1 w labels center offset 0.0,1.0 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -1.6580e0 5.7356e0 400" u 4:(1e9*period(beta(400e-9,real(epsilonAu400),$4))) w p lw 4 lc 3 t '{/Symbol l}=400 nm', \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -1.6580e0 5.7356e0 400" u 4:(1e9*period(beta(400e-9,real(epsilonAu400),$4))):1 w labels center offset 0.0,1.0 notitle
# 
# 
# set output '20141208-SPPfriendlyMaterials-Period-Ag.eps'
# set terminal postscript eps enhanced color font 'Helvetica, 22' size 10cm, 8cm
# 
# set key top left
# 
# set title 'Ag'
# plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -31.0213e0 0.4095e0 800" u 4:(1e9*period(beta(800e-9,epsilonAg800,$4*Unit+$5*Imaginary))) w p lc 1 lw 4 t '{/Symbol l}=800 nm', \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -31.0213e0 0.4095e0 800" u 4:(1e9*period(beta(800e-9,epsilonAg800,$4*Unit+$5*Imaginary))):1 w labels center offset 0.0,1.0 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -4.4222e0 0.2103e0 400" u 4:(1e9*period(beta(400e-9,epsilonAg400,$4*Unit+$5*Imaginary))) w p lw 4 lc 3 t '{/Symbol l}=400 nm', \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -4.4222e0 0.2103e0 400" u 4:(1e9*period(beta(400e-9,epsilonAg400,$4*Unit+$5*Imaginary))):1 w labels center offset 0.0,1.0 notitle
# 
# 
# set output '20141208-SPPfriendlyMaterials-Period-Ag-ImagZero.eps'
# set terminal postscript eps enhanced color font 'Helvetica, 22' size 10cm, 8cm
# 
# set key top left
# 
# set title 'Ag'
# plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -31.0213e0 0.4095e0 800" u 4:(1e9*period(beta(800e-9,real(epsilonAu400),$4))) w p lc 1 lw 4 t '{/Symbol l}=800 nm', \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -31.0213e0 0.4095e0 800" u 4:(1e9*period(beta(800e-9,real(epsilonAu400),$4))):1 w labels center offset 0.0,1.0 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -4.4222e0 0.2103e0 400" u 4:(1e9*period(beta(400e-9,real(epsilonAg400),$4))) w p lw 4 lc 3 t '{/Symbol l}=400 nm', \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -4.4222e0 0.2103e0 400" u 4:(1e9*period(beta(400e-9,real(epsilonAg400),$4))):1 w labels center offset 0.0,1.0 notitle
# 
# set key top right
# 
# set yrange [0:900]
# 
# set output '20141208-SPPfriendlyMaterials-Period-Ti.eps'
# set terminal postscript eps enhanced color font 'Helvetica, 22' size 10cm, 8cm
# 
# set title 'Ti'
# plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -6.2067e0 25.2004e0 800" u 4:(1e9*period(beta(800e-9,epsilonTi800,$4*Unit+$5*Imaginary))) w p lc 1 lw 4 t '{/Symbol l}=800 nm', \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -6.2067e0 25.2004e0 800" u 4:(1e9*period(beta(800e-9,epsilonTi800,$4*Unit+$5*Imaginary))):1 w labels center offset 0.0,1.0 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -4.3620e0 12.3621e0 400" u 4:(1e9*period(beta(400e-9,epsilonTi400,$4*Unit+$5*Imaginary))) w p lw 4 lc 3 t '{/Symbol l}=400 nm', \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -4.3620e0 12.3621e0 400" u 4:(1e9*period(beta(400e-9,epsilonTi400,$4*Unit+$5*Imaginary))):1 w labels center offset 0.0,1.0 notitle
# 
# set output '20141208-SPPfriendlyMaterials-Period-Ti-ImagZero.eps'
# set terminal postscript eps enhanced color font 'Helvetica, 22' size 10cm, 8cm
# 
# set title 'Ti'
# plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -6.2067e0 25.2004e0 800" u 4:(1e9*period(beta(800e-9,real(epsilonTi800),$4))) w p lc 1 lw 4 t '{/Symbol l}=800 nm', \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -6.2067e0 25.2004e0 800" u 4:(1e9*period(beta(800e-9,real(epsilonTi800),$4))):1 w labels center offset 0.0,1.0 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -4.3620e0 12.3621e0 400" u 4:(1e9*period(beta(400e-9,real(epsilonTi400),$4))) w p lw 4 lc 3 t '{/Symbol l}=400 nm', \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -4.3620e0 12.3621e0 400" u 4:(1e9*period(beta(400e-9,real(epsilonTi400),$4))):1 w labels center offset 0.0,1.0 notitle
# 
# 
# # unset multiplot
# 
# set output '20141208-SPPfriendlyMaterials-Period-W.eps'
# set terminal postscript eps enhanced color font 'Helvetica, 22' size 10cm, 8cm
# 
# set yrange [10:900]
# 
# set title 'W'
# set log y
# plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 4.4655e0 19.905e0 800" u 4:(1e9*period(beta(800e-9,epsilonW800,$4*Unit+$5*Imaginary))) w p lw 4 t '{/Symbol l}=800 nm', \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 4.4655e0 19.905e0 800" u 4:(1e9*period(beta(800e-9,epsilonW800,$4*Unit+$5*Imaginary))):1 w labels center offset 0.0,1.0 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 4.0653e0 16.453e0 400" u 4:(1e9*period(beta(400e-9,epsilonW400,$4*Unit+$5*Imaginary))) w p lw 4 lc 3 t '{/Symbol l}=400 nm', \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 4.0653e0 16.453e0 400" u 4:(1e9*period(beta(400e-9,epsilonW400,$4*Unit+$5*Imaginary))):1 w labels center offset 0.0,1.0 notitle
# unset log y
# 
# set output '20141208-SPPfriendlyMaterials-Period-W-ImagZero.eps'
# set terminal postscript eps enhanced color font 'Helvetica, 22' size 10cm, 8cm
# 
# set yrange [10:900]
# 
# set title 'W'
# set log y
# plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 4.4655e0 19.905e0 800" u 4:(1e9*period(beta(800e-9,real(epsilonW800),$4))) w p lw 4 t '{/Symbol l}=800 nm', \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 4.4655e0 19.905e0 800" u 4:(1e9*period(beta(800e-9,real(epsilonW800),$4))):1 w labels center offset 0.0,1.0 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 4.0653e0 16.453e0 400" u 4:(1e9*period(beta(400e-9,real(epsilonW400),$4))) w p lw 4 lc 3 t '{/Symbol l}=400 nm', \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 4.0653e0 16.453e0 400" u 4:(1e9*period(beta(400e-9,real(epsilonW400),$4))):1 w labels center offset 0.0,1.0 notitle
# unset log y
# 
# set output '20141208-SPPfriendlyMaterials-Period-Ge.eps'
# set terminal postscript eps enhanced color font 'Helvetica, 22' size 10cm, 8cm
# 
# set title 'Ge'
# plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 22.0765e0 3.0347e0 800" u 4:(1e9*period(beta(800e-9,epsilonGe800,$4*Unit+$5*Imaginary))) w p lw 4 t '{/Symbol l}=800 nm', \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 22.0765e0 3.0347e0 800" u 4:(1e9*period(beta(800e-9,epsilonGe800,$4*Unit+$5*Imaginary))):1 w labels center offset 0.0,1.0 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 12.2406e0 18.3387e0 400" u 4:(1e9*period(beta(400e-9,epsilonGe400,$4*Unit+$5*Imaginary))) w p lw 4 lc 3 t '{/Symbol l}=400 nm', \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 12.2406e0 18.3387e0 400" u 4:(1e9*period(beta(400e-9,epsilonGe400,$4*Unit+$5*Imaginary))):1 w labels center offset 0.0,1.0 notitle
# 
# set output '20141208-SPPfriendlyMaterials-Period-Ge-ImagZero.eps'
# set terminal postscript eps enhanced color font 'Helvetica, 22' size 10cm, 8cm
# 
# set title 'Ge'
# plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 22.0765e0 3.0347e0 800" u 4:(1e9*period(beta(800e-9,real(epsilonGe800),$4))) w p lw 4 t '{/Symbol l}=800 nm', \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 22.0765e0 3.0347e0 800" u 4:(1e9*period(beta(800e-9,real(epsilonGe800),$4))):1 w labels center offset 0.0,1.0 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 12.2406e0 18.3387e0 400" u 4:(1e9*period(beta(400e-9,real(epsilonGe400),$4))) w p lw 4 lc 3 t '{/Symbol l}=400 nm', \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 12.2406e0 18.3387e0 400" u 4:(1e9*period(beta(400e-9,real(epsilonGe400),$4))):1 w labels center offset 0.0,1.0 notitle
# 
# set output '20141208-SPPfriendlyMaterials-Period-Si.eps'
# set terminal postscript eps enhanced color font 'Helvetica, 22' size 10cm, 8cm
# 
# set title 'Si'
# set log y
# plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 13.64e0 0.048e0 800" u 4:(1e9*period(beta(800e-9,epsilonSi800,$4*Unit+$5*Imaginary))) w p lw 4 t '{/Symbol l}=800 nm', \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 13.64e0 0.048e0 800" u 4:(1e9*period(beta(800e-9,epsilonSi800,$4*Unit+$5*Imaginary))):1 w labels center offset 0.0,1.0 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 30.8469e0 4.2994e0 400" u 4:(1e9*period(beta(400e-9,epsilonSi400,$4*Unit+$5*Imaginary))) w p lw 4 lc 3 t '{/Symbol l}=400 nm', \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 30.8469e0 4.2994e0 400" u 4:(1e9*period(beta(400e-9,epsilonSi400,$4*Unit+$5*Imaginary))):1 w labels center offset 0.0,1.0 notitle
# unset log y
# 
# set output '20141208-SPPfriendlyMaterials-Period-Si-ImagZero.eps'
# set terminal postscript eps enhanced color font 'Helvetica, 22' size 10cm, 8cm
# 
# set title 'Si'
# set log y
# plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 13.64e0 0.048e0 800" u 4:(1e9*period(beta(800e-9,real(epsilonSi800),$4))) w p lw 4 t '{/Symbol l}=800 nm', \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 13.64e0 0.048e0 800" u 4:(1e9*period(beta(800e-9,real(epsilonSi800),$4))):1 w labels center offset 0.0,1.0 notitle, \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 30.8469e0 4.2994e0 400" u 4:(1e9*period(beta(400e-9,real(epsilonSi400),$4))) w p lw 4 lc 3 t '{/Symbol l}=400 nm', \
# "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 30.8469e0 4.2994e0 400" u 4:(1e9*period(beta(400e-9,real(epsilonSi400),$4))):1 w labels center offset 0.0,1.0 notitle
# unset log y

set output '20141208-SPPfriendlyMaterials-Period-SiO2.eps'
set terminal postscript eps enhanced color font 'Helvetica, 22' size 10cm, 8cm

set yrange [0:600]
set xrange [-70:0]
set title 'SiO_2'
set key bottom left
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 2.1121e0 0e0 800" u 4:(1e9*period(beta(800e-9,epsilonSiO2800,$4*Unit+$5*Imaginary))) w p lw 4 t '{/Symbol l}=800 nm', \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 2.1121e0 0e0 800" u 4:(1e9*period(beta(800e-9,epsilonSiO2800,$4*Unit+$5*Imaginary))):1 w labels center offset 0.0,1.0 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 2.1612e0 0e0 400" u 4:(1e9*period(beta(400e-9,epsilonSiO2400,$4*Unit+$5*Imaginary))) w p lw 4 lc 3 t '{/Symbol l}=400 nm', \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 2.1612e0 0e0 400" u 4:(1e9*period(beta(400e-9,epsilonSiO2400,$4*Unit+$5*Imaginary))):1 w labels center offset 0.0,1.0 notitle

set output '20141208-SPPfriendlyMaterials-Period-SiO2-ImagZero.eps'
set terminal postscript eps enhanced color font 'Helvetica, 22' size 10cm, 8cm

set yrange [0:600]
set xrange [-70:0]

set title 'SiO_2'
set key bottom left
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 2.1121e0 0e0 800" u 4:(1e9*period(beta(800e-9,real(epsilonSiO2800),$4))) w p lw 4 t '{/Symbol l}=800 nm', \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 2.1121e0 0e0 800" u 4:(1e9*period(beta(800e-9,real(epsilonSiO2800),$4))):1 w labels center offset 0.0,1.0 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 2.1612e0 0e0 400" u 4:(1e9*period(beta(400e-9,real(epsilonSiO2400),$4))) w p lw 4 lc 3 t '{/Symbol l}=400 nm', \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 2.1612e0 0e0 400" u 4:(1e9*period(beta(400e-9,real(epsilonSiO2400),$4))):1 w labels center offset 0.0,1.0 notitle

set output '20141208-SPPfriendlyMaterials-Period-Air.eps'
set terminal postscript eps enhanced color font 'Helvetica, 22' size 10cm, 8cm

set key bottom left
set yrange [200:900]

set title 'Air'
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 1e0 0e0 800" u 4:(1e9*period(beta(800e-9,epsilonAir,$4*Unit+$5*Imaginary))) w p lw 4 t '{/Symbol l}=800 nm', \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 1e0 0e0 800" u 4:(1e9*period(beta(800e-9,epsilonAir,$4*Unit+$5*Imaginary))):1 w labels center offset 0.0,1.0 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 1e0 0e0 400" u 4:(1e9*period(beta(400e-9,epsilonAir,$4*Unit+$5*Imaginary))) w p lw 4 lc 3 t '{/Symbol l}=400 nm', \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 1e0 0e0 400" u 4:(1e9*period(beta(400e-9,epsilonAir,$4*Unit+$5*Imaginary))):1 w labels center offset 0.0,1.0 notitle

set output '20141208-SPPfriendlyMaterials-Period-Air-ImagZero.eps'
set terminal postscript eps enhanced color font 'Helvetica, 22' size 10cm, 8cm

set key bottom left

set yrange [200:900]

set title 'Air'
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 1e0 0e0 800" u 4:(1e9*period(beta(800e-9,epsilonAir,$4))) w p lw 4 t '{/Symbol l}=800 nm', \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 1e0 0e0 800" u 4:(1e9*period(beta(800e-9,epsilonAir,$4))):1 w labels center offset 0.0,1.0 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 1e0 0e0 400" u 4:(1e9*period(beta(400e-9,epsilonAir,$4))) w p lw 4 lc 3 t '{/Symbol l}=400 nm', \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 1e0 0e0 400" u 4:(1e9*period(beta(400e-9,epsilonAir,$4))):1 w labels center offset 0.0,1.0 notitle



###########################################################
################### SPP period ############################
####### Plotting as function of band gap / fermi energy #####
set grid

set xlabel 'Fermi Energy (eV < 0), Band gap energy (> 0)' 

unset y2tics 

set ytics format "%g" textcolor rgbcolor "black"
set ylabel textcolor rgbcolor "black"
set y2label ''

unset log xy
unset log y2

set xrange [*:*]
set yrange [*:*]
set y2range [*:*]

unset multiplot

set key top right
# set multiplot layout 2,2

unset log y
set ylabel 'Period (nm)'

set output '20141208-SPPfriendlyMaterials-Period-ConditionWithBandGap-Au.eps'
set terminal postscript eps enhanced color font 'Times, 20' #size 8cm,8cm

set title 'Au'
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -24.0620e0 1.5069e0 800" u 2:(1e9*period(beta(800e-9,epsilonAu800,$4*Unit+$5*Imaginary))) w p lc 1 lw 4 t '{/Symbol l}=800 nm', \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -24.0620e0 1.5069e0 800" u 2:(1e9*period(beta(800e-9,epsilonAu800,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -1.6580e0 5.7356e0 400" u 2:(1e9*period(beta(400e-9,epsilonAu400,$4*Unit+$5*Imaginary))) w p lw 4 lc 3 t '{/Symbol l}=400 nm', \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -1.6580e0 5.7356e0 400" u 2:(1e9*period(beta(400e-9,epsilonAu400,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle


set output '20141208-SPPfriendlyMaterials-Period-ConditionWithBandGap-Ag.eps'
set terminal postscript eps enhanced color font 'Times, 20' #size 8cm,8cm

set key top left
set log y
set title 'Ag'
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -31.0213e0 0.4095e0 800" u 2:(1e9*period(beta(800e-9,epsilonAg800,$4*Unit+$5*Imaginary))) w p lc 1 lw 4 t '{/Symbol l}=800 nm', \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -31.0213e0 0.4095e0 800" u 2:(1e9*period(beta(800e-9,epsilonAg800,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -4.4222e0 0.2103e0 400" u 2:(1e9*period(beta(400e-9,epsilonAg400,$4*Unit+$5*Imaginary))) w p lw 4 lc 3 t '{/Symbol l}=400 nm', \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -4.4222e0 0.2103e0 400" u 2:(1e9*period(beta(400e-9,epsilonAg400,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle
unset log y
set output '20141208-SPPfriendlyMaterials-Period-ConditionWithBandGap-Ti.eps'
set terminal postscript eps enhanced color font 'Times, 20' #size 8cm,8cm

set title 'Ti'
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -6.2067e0 25.2004e0 800" u 2:(1e9*period(beta(800e-9,epsilonTi800,$4*Unit+$5*Imaginary))) w p lc 1 lw 4 t '{/Symbol l}=800 nm', \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -6.2067e0 25.2004e0 800" u 2:(1e9*period(beta(800e-9,epsilonTi800,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -4.3620e0 12.3621e0 400" u 2:(1e9*period(beta(400e-9,epsilonTi400,$4*Unit+$5*Imaginary))) w p lw 4 lc 3 t '{/Symbol l}=400 nm', \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -4.3620e0 12.3621e0 400" u 2:(1e9*period(beta(400e-9,epsilonTi400,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle

# unset multiplot

set output '20141208-SPPfriendlyMaterials-Period-ConditionWithBandGap-W.eps'
set terminal postscript eps enhanced color font 'Times, 20' #size 8cm,8cm

set title 'W'
set log y
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 4.4655e0 19.905e0 800" u 2:(1e9*period(beta(800e-9,epsilonW800,$4*Unit+$5*Imaginary))) w p lw 4 t '{/Symbol l}=800 nm', \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 4.4655e0 19.905e0 800" u 2:(1e9*period(beta(800e-9,epsilonW800,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 4.0653e0 16.453e0 400" u 2:(1e9*period(beta(400e-9,epsilonW400,$4*Unit+$5*Imaginary))) w p lw 4 lc 3 t '{/Symbol l}=400 nm', \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 4.0653e0 16.453e0 400" u 2:(1e9*period(beta(400e-9,epsilonW400,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle
unset log y

set output '20141208-SPPfriendlyMaterials-Period-ConditionWithBandGap-Ge.eps'
set terminal postscript eps enhanced color font 'Times, 20' #size 8cm,8cm

set title 'Ge'
set log y
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 22.0765e0 3.0347e0 800" u 2:(1e9*period(beta(800e-9,epsilonGe800,$4*Unit+$5*Imaginary))) w p lw 4 t '{/Symbol l}=800 nm', \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 22.0765e0 3.0347e0 800" u 2:(1e9*period(beta(800e-9,epsilonGe800,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 12.2406e0 18.3387e0 400" u 2:(1e9*period(beta(400e-9,epsilonGe400,$4*Unit+$5*Imaginary))) w p lw 4 lc 3 t '{/Symbol l}=400 nm', \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 12.2406e0 18.3387e0 400" u 2:(1e9*period(beta(400e-9,epsilonGe400,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle
unset log y

set output '20141208-SPPfriendlyMaterials-Period-ConditionWithBandGap-Si.eps'
set terminal postscript eps enhanced color font 'Times, 20' #size 8cm,8cm

set title 'Si'
set log y
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 13.64e0 0.048e0 800" u 2:(1e9*period(beta(800e-9,epsilonSi800,$4*Unit+$5*Imaginary))) w p lw 4 t '{/Symbol l}=800 nm', \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 13.64e0 0.048e0 800" u 2:(1e9*period(beta(800e-9,epsilonSi800,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 30.8469e0 4.2994e0 400" u 2:(1e9*period(beta(400e-9,epsilonSi400,$4*Unit+$5*Imaginary))) w p lw 4 lc 3 t '{/Symbol l}=400 nm', \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 30.8469e0 4.2994e0 400" u 2:(1e9*period(beta(400e-9,epsilonSi400,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle
unset log y

set output '20141208-SPPfriendlyMaterials-Period-ConditionWithBandGap-SiO2.eps'
set terminal postscript eps enhanced color font 'Times, 20' #size 8cm,8cm

set title 'SiO2'
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 2.1121e0 0e0 800" u 2:(1e9*period(beta(800e-9,epsilonSiO2800,$4*Unit+$5*Imaginary))) w p lw 4 t '{/Symbol l}=800 nm', \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 2.1121e0 0e0 800" u 2:(1e9*period(beta(800e-9,epsilonSiO2800,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 2.1612e0 0e0 400" u 2:(1e9*period(beta(400e-9,epsilonSiO2400,$4*Unit+$5*Imaginary))) w p lw 4 lc 3 t '{/Symbol l}=400 nm', \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 2.1612e0 0e0 400" u 2:(1e9*period(beta(400e-9,epsilonSiO2400,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle


set output '20141208-SPPfriendlyMaterials-Period-ConditionWithBandGap-Air.eps'
set terminal postscript eps enhanced color font 'Times, 20' #size 8cm,8cm

set title 'Air'
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 1e0 0e0 800" u 2:(1e9*period(beta(800e-9,epsilonAir,$4*Unit+$5*Imaginary))) w p lw 4 t '{/Symbol l}=800 nm', \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 1e0 0e0 800" u 2:(1e9*period(beta(800e-9,epsilonAir,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 1e0 0e0 400" u 2:(1e9*period(beta(400e-9,epsilonAir,$4*Unit+$5*Imaginary))) w p lw 4 lc 3 t '{/Symbol l}=400 nm', \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 1e0 0e0 400" u 2:(1e9*period(beta(400e-9,epsilonAir,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle

# unset multiplot


################ SPP decay depth ###########
##### expressed with Re(epsilon) ###########
unset multiplot

# set multiplot layout 2,2

set log y
# set log y2

set xlabel 'Re({/Symbol e})'
set ylabel 'Decay depth (nm)'
# set y2label 'Decay depth (nm)'

set ytics #nomirror textcolor rgbcolor "red"
# set y2tics nomirror textcolor rgbcolor "blue"


set output '20141208-SPPfriendlyMaterials-DecayDepth-Au800.eps'
set terminal postscript eps enhanced color font 'Times, 20' #size 8cm,8cm

set title 'Au, 800 nm'
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -24.0620e0 0.4095e0 800" u 4:(1e9*DecayDepth(kzSPP(800e-9,epsilonAg800,$4*Unit+$5*Imaginary))) w p lw 4 lc 1 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -24.0620e0 0.4095e0 800" u 4:(1e9*DecayDepth(kzSPP(800e-9,epsilonAu800,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -24.0620e0 0.4095e0 800" u 4:(1e9*DecayDepth(kzSPP(800e-9,$4*Unit+$5*Imaginary, epsilonAu800))) w p lw 4 lc 3 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -24.0620e0 0.4095e0 800" u 4:(1e9*DecayDepth(kzSPP(800e-9,$4*Unit+$5*Imaginary, epsilonAu800))):1 w labels center offset 3.0,0.0 notitle

set output '20141208-SPPfriendlyMaterials-DecayDepth-Au400.eps'
set terminal postscript eps enhanced color font 'Times, 20' #size 8cm,8cm

set title 'Au, 400 nm'
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -1.6580e0 5.7356e0 400" u 4:(1e9*DecayDepth(kzSPP(400e-9,epsilonAu400,$4*Unit+$5*Imaginary))) w p lw 4 lc 1 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -1.6580e0 5.7356e0 400" u 4:(1e9*DecayDepth(kzSPP(400e-9,epsilonAu400,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -1.6580e0 5.7356e0 400" u 4:(1e9*DecayDepth(kzSPP(400e-9,$4*Unit+$5*Imaginary, epsilonAu400))) w p lw 4 lc 3 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -1.6580e0 5.7356e0 400" u 4:(1e9*DecayDepth(kzSPP(400e-9,$4*Unit+$5*Imaginary, epsilonAu400))):1 w labels center offset 3.0,0.0 notitle

set output '20141208-SPPfriendlyMaterials-DecayDepth-Ag800.eps'
set terminal postscript eps enhanced color font 'Times, 20' #size 8cm,8cm

set title 'Ag, 800 nm'
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -31.0213e0 0.4095e0 800" u 4:(1e9*DecayDepth(kzSPP(800e-9,epsilonAg800,$4*Unit+$5*Imaginary))) w p lw 4 lc 1 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -31.0213e0 0.4095e0 800" u 4:(1e9*DecayDepth(kzSPP(800e-9,epsilonAg800,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -31.0213e0 0.4095e0 800" u 4:(1e9*DecayDepth(kzSPP(800e-9,$4*Unit+$5*Imaginary, epsilonAg800))) w p lw 4 lc 3 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -31.0213e0 0.4095e0 800" u 4:(1e9*DecayDepth(kzSPP(800e-9,$4*Unit+$5*Imaginary, epsilonAg800))):1 w labels center offset 3.0,0.0 notitle

set output '20141208-SPPfriendlyMaterials-DecayDepth-Ag400.eps'
set terminal postscript eps enhanced color font 'Times, 20' #size 8cm,8cm

set title 'Ag, 400 nm'
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -4.4222e0 0.2103e0 400" u 4:(1e9*DecayDepth(kzSPP(400e-9,epsilonAg400,$4*Unit+$5*Imaginary))) w p lc 1 lw 4 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -4.4222e0 0.2103e0 400" u 4:(1e9*DecayDepth(kzSPP(400e-9,epsilonAg400,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -4.4222e0 0.2103e0 400" u 4:(1e9*DecayDepth(kzSPP(400e-9,$4*Unit+$5*Imaginary,epsilonAg400))) w p lw 4 lc 3 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -4.4222e0 0.2103e0 400" u 4:(1e9*DecayDepth(kzSPP(400e-9,$4*Unit+$5*Imaginary,epsilonAg400))):1 w labels center offset 3.0,0.0 notitle

set output '20141208-SPPfriendlyMaterials-DecayDepth-Ti800.eps'
set terminal postscript eps enhanced color font 'Times, 20' #size 8cm,8cm

set title 'Ti, 800 nm'
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -6.2067e0 25.2004e0 800" u 4:(1e9*DecayDepth(kzSPP(800e-9,epsilonTi800,$4*Unit+$5*Imaginary))) w p lw 4 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -6.2067e0 25.2004e0 800" u 4:(1e9*DecayDepth(kzSPP(800e-9,epsilonTi800,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -6.2067e0 25.2004e0 800" u 4:(1e9*DecayDepth(kzSPP(800e-9,$4*Unit+$5*Imaginary,epsilonTi800))) w p lw 4 lc 3 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -6.2067e0 25.2004e0 800" u 4:(1e9*DecayDepth(kzSPP(800e-9,$4*Unit+$5*Imaginary,epsilonTi800))):1 w labels center offset 3.0,0.0 notitle

set output '20141208-SPPfriendlyMaterials-DecayDepth-Ti400.eps'
set terminal postscript eps enhanced color font 'Times, 20' #size 8cm,8cm

set title 'Ti, 400 nm'
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -4.3620e0 12.3621e0 400" u 4:(1e9*DecayDepth(kzSPP(400e-9,epsilonTi400,$4*Unit+$5*Imaginary))) w p lw 4 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -4.3620e0 12.3621e0 400" u 4:(1e9*DecayDepth(kzSPP(400e-9,epsilonTi400,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -4.3620e0 12.3621e0 400" u 4:(1e9*DecayDepth(kzSPP(400e-9,$4*Unit+$5*Imaginary,epsilonTi400))) w p lw 4 lc 3 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -4.3620e0 12.3621e0 400" u 4:(1e9*DecayDepth(kzSPP(400e-9,$4*Unit+$5*Imaginary,epsilonTi400))):1 w labels center offset 3.0,0.0 notitle


set output '20141208-SPPfriendlyMaterials-DecayDepth-W800.eps'
set terminal postscript eps enhanced color font 'Times, 20' #size 8cm,8cm

set title 'W, 800 nm'
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 4.4655e0 19.905e0 800" u 4:(1e9*DecayDepth(kzSPP(800e-9,epsilonW800,$4*Unit+$5*Imaginary))) w p lw 4 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 4.4655e0 19.905e0 800" u 4:(1e9*DecayDepth(kzSPP(800e-9,epsilonW800,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 4.4655e0 19.905e0 800" u 4:(1e9*DecayDepth(kzSPP(800e-9,$4*Unit+$5*Imaginary,epsilonW800))) w p lw 4 lc 3 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 4.4655e0 19.905e0 800" u 4:(1e9*DecayDepth(kzSPP(800e-9,$4*Unit+$5*Imaginary,epsilonW800))):1 w labels center offset 3.0,0.0 notitle

set output '20141208-SPPfriendlyMaterials-DecayDepth-W400.eps'
set terminal postscript eps enhanced color font 'Times, 20' #size 8cm,8cm

set title 'W, 400 nm'
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 4.0653e0 16.453e0 400" u 4:(1e9*DecayDepth(kzSPP(400e-9,epsilonW400,$4*Unit+$5*Imaginary))) w p lw 4 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 4.0653e0 16.453e0 400" u 4:(1e9*DecayDepth(kzSPP(400e-9,epsilonW400,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 4.0653e0 16.453e0 400" u 4:(1e9*DecayDepth(kzSPP(400e-9,$4*Unit+$5*Imaginary,epsilonW400))) w p lw 4 lc 3 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 4.0653e0 16.453e0 400" u 4:(1e9*DecayDepth(kzSPP(400e-9,$4*Unit+$5*Imaginary,epsilonW400))):1 w labels center offset 3.0,0.0 notitle

set output '20141208-SPPfriendlyMaterials-DecayDepth-Ge800.eps'
set terminal postscript eps enhanced color font 'Times, 20' #size 8cm,8cm

set title 'Ge, 800 nm'
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 22.0765e0 3.0347e0 800" u 4:(1e9*DecayDepth(kzSPP(800e-9,epsilonGe800,$4*Unit+$5*Imaginary))) w p lw 4 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 22.0765e0 3.0347e0 800" u 4:(1e9*DecayDepth(kzSPP(800e-9,epsilonGe800,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 22.0765e0 3.0347e0 800" u 4:(1e9*DecayDepth(kzSPP(800e-9,$4*Unit+$5*Imaginary,epsilonGe800))) w p lw 4 lc 3 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 22.0765e0 3.0347e0 800" u 4:(1e9*DecayDepth(kzSPP(800e-9,$4*Unit+$5*Imaginary,epsilonGe800))):1 w labels center offset 3.0,0.0 notitle

set output '20141208-SPPfriendlyMaterials-DecayDepth-Ge400.eps'
set terminal postscript eps enhanced color font 'Times, 20' #size 8cm,8cm

set title 'Ge, 400 nm'
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 12.2406e0 18.3387e0 400" u 4:(1e9*DecayDepth(kzSPP(800e-9,epsilonGe400,$4*Unit+$5*Imaginary))) w p lw 4 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 12.2406e0 18.3387e0 400" u 4:(1e9*DecayDepth(kzSPP(800e-9,epsilonGe400,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 12.2406e0 18.3387e0 400" u 4:(1e9*DecayDepth(kzSPP(800e-9,$4*Unit+$5*Imaginary,epsilonGe400))) w p lw 4 lc 3 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 12.2406e0 18.3387e0 400" u 4:(1e9*DecayDepth(kzSPP(800e-9,$4*Unit+$5*Imaginary,epsilonGe400))):1 w labels center offset 3.0,0.0 notitle


set output '20141208-SPPfriendlyMaterials-DecayDepth-Si800.eps'
set terminal postscript eps enhanced color font 'Times, 20' #size 8cm,8cm

set title 'Si, 800 nm'
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 13.64e0 0.048e0 800" u 4:(1e9*DecayDepth(kzSPP(800e-9,epsilonSi800,$4*Unit+$5*Imaginary))) w p lw 4 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 13.64e0 0.048e0 800" u 4:(1e9*DecayDepth(kzSPP(800e-9,epsilonSi800,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 13.64e0 0.048e0 800" u 4:(1e9*DecayDepth(kzSPP(800e-9,$4*Unit+$5*Imaginary,epsilonSi800))) w p lw 4 lc 3 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 13.64e0 0.048e0 800" u 4:(1e9*DecayDepth(kzSPP(800e-9,$4*Unit+$5*Imaginary,epsilonSi800))):1 w labels center offset 3.0,0.0 notitle


set output '20141208-SPPfriendlyMaterials-DecayDepth-Si400.eps'
set terminal postscript eps enhanced color font 'Times, 20' #size 8cm,8cm

set title 'Si, 400 nm'
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 30.8469e0 4.2994e0 400" u 4:(1e9*DecayDepth(kzSPP(400e-9,epsilonSi400,$4*Unit+$5*Imaginary))) w p lw 4 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 30.8469e0 4.2994e0 400" u 4:(1e9*DecayDepth(kzSPP(400e-9,epsilonSi400,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 30.8469e0 4.2994e0 400" u 4:(1e9*DecayDepth(kzSPP(400e-9,$4*Unit+$5*Imaginary,epsilonSi400))) w p lw 4 lc 3 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 30.8469e0 4.2994e0 400" u 4:(1e9*DecayDepth(kzSPP(400e-9,$4*Unit+$5*Imaginary,epsilonSi400))):1 w labels center offset 3.0,0.0 notitle


set output '20141208-SPPfriendlyMaterials-DecayDepth-SiO2800.eps'
set terminal postscript eps enhanced color font 'Times, 20' #size 8cm,8cm

set title 'SiO_2, 800 nm'
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 2.1121e0 0e0 800" u 4:(1e9*DecayDepth(kzSPP(800e-9,epsilonSiO2800,$4*Unit+$5*Imaginary))) w p lw 4 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 2.1121e0 0e0 800" u 4:(1e9*DecayDepth(kzSPP(800e-9,epsilonSiO2800,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 2.1121e0 0e0 800" u 4:(1e9*DecayDepth(kzSPP(800e-9,$4*Unit+$5*Imaginary,epsilonSiO2800))) w p lw 4 lc 3 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 2.1121e0 0e0 800" u 4:(1e9*DecayDepth(kzSPP(800e-9,$4*Unit+$5*Imaginary,epsilonSiO2800))):1 w labels center offset 3.0,0.0 notitle

set output '20141208-SPPfriendlyMaterials-DecayDepth-SiO2400.eps'
set terminal postscript eps enhanced color font 'Times, 20' #size 8cm,8cm

set title 'SiO_2, 400 nm'
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 2.1612e0 0e0 400" u 4:(1e9*DecayDepth(kzSPP(400e-9,epsilonSiO2400,$4*Unit+$5*Imaginary))) w p lw 4 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 2.1612e0 0e0 400" u 4:(1e9*DecayDepth(kzSPP(400e-9,epsilonSiO2400,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 2.1612e0 0e0 400" u 4:(1e9*DecayDepth(kzSPP(400e-9,$4*Unit+$5*Imaginary,epsilonSiO2400))) w p lw 4 lc 3 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 2.1612e0 0e0 400" u 4:(1e9*DecayDepth(kzSPP(400e-9,$4*Unit+$5*Imaginary,epsilonSiO2400))):1 w labels center offset 3.0,0.0 notitle

unset multiplot

set output '20141208-SPPfriendlyMaterials-DecayDepth-Air800.eps'
set terminal postscript eps enhanced color font 'Times, 20' #size 8cm,8cm

set title 'Air, 800 nm'
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 1e0 0e0 800" u 4:(1e9*DecayDepth(kzSPP(800e-9,epsilonAir,$4*Unit+$5*Imaginary))) w p lw 4 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 1e0 0e0 800" u 4:(1e9*DecayDepth(kzSPP(800e-9,epsilonAir,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 1e0 0e0 800" u 4:(1e9*DecayDepth(kzSPP(800e-9,$4*Unit+$5*Imaginary,epsilonAir))) w p lw 4 lc 3 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 1e0 0e0 800" u 4:(1e9*DecayDepth(kzSPP(800e-9,$4*Unit+$5*Imaginary,epsilonAir))):1 w labels center offset 3.0,0.0 notitle

set output '20141208-SPPfriendlyMaterials-DecayDepth-Air400.eps'
set terminal postscript eps enhanced color font 'Times, 20' #size 8cm,8cm

set title 'Air, 400 nm'
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 1e0 0e0 400" u 4:(1e9*DecayDepth(kzSPP(400e-9,epsilonAir,$4*Unit+$5*Imaginary))) w p lw 4 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 1e0 0e0 400" u 4:(1e9*DecayDepth(kzSPP(400e-9,epsilonAir,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 1e0 0e0 400" u 4:(1e9*DecayDepth(kzSPP(400e-9,$4*Unit+$5*Imaginary,epsilonAir))) w p lw 4 lc 3 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 1e0 0e0 400" u 4:(1e9*DecayDepth(kzSPP(400e-9,$4*Unit+$5*Imaginary,epsilonAir))):1 w labels center offset 3.0,0.0 notitle

unset multiplot

################ SPP decay depth ###########
##### expressed with Band gap ###########
unset multiplot

# set multiplot layout 2,2

set log y
# set log y2

set xlabel 'Fermi Energy (eV < 0), Band gap energy (> 0)' 
set ylabel 'Decay depth (nm)'
# set y2label 'Decay depth in substrate (nm)'

set ytics # nomirror textcolor rgbcolor "red"
# set y2tics nomirror textcolor rgbcolor "blue"


set output '20141208-SPPfriendlyMaterials-DecayDepth-BandGap-Au800.eps'
set terminal postscript eps enhanced color font 'Times, 20' #size 8cm,8cm

set title 'Au, 800 nm'
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -24.0620e0 0.4095e0 800" u 2:(1e9*DecayDepth(kzSPP(800e-9,epsilonAu800,$4*Unit+$5*Imaginary))) w p lw 4 lc 1 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -24.0620e0 0.4095e0 800" u 2:(1e9*DecayDepth(kzSPP(800e-9,epsilonAu800,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -24.0620e0 0.4095e0 800" u 2:(1e9*DecayDepth(kzSPP(800e-9,$4*Unit+$5*Imaginary, epsilonAu800))) w p lw 4 lc 3 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -24.0620e0 0.4095e0 800" u 2:(1e9*DecayDepth(kzSPP(800e-9,$4*Unit+$5*Imaginary, epsilonAu800))):1 w labels center offset 3.0,0.0 notitle

set output '20141208-SPPfriendlyMaterials-DecayDepth-BandGap-Au400.eps'
set terminal postscript eps enhanced color font 'Times, 20' #size 8cm,8cm

set title 'Au, 400 nm'
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -1.6580e0 5.7356e0 400" u 2:(1e9*DecayDepth(kzSPP(400e-9,epsilonAu400,$4*Unit+$5*Imaginary))) w p lw 4 lc 1 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -1.6580e0 5.7356e0 400" u 2:(1e9*DecayDepth(kzSPP(400e-9,epsilonAu400,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -1.6580e0 5.7356e0 400" u 2:(1e9*DecayDepth(kzSPP(400e-9,$4*Unit+$5*Imaginary, epsilonAu400))) w p lw 4 lc 3 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -1.6580e0 5.7356e0 400" u 2:(1e9*DecayDepth(kzSPP(400e-9,$4*Unit+$5*Imaginary, epsilonAu400))):1 w labels center offset 3.0,0.0 notitle 

set output '20141208-SPPfriendlyMaterials-DecayDepth-BandGap-Ag800.eps'
set terminal postscript eps enhanced color font 'Times, 20' #size 8cm,8cm

set title 'Ag, 800 nm'
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -31.0213e0 0.4095e0 800" u 2:(1e9*DecayDepth(kzSPP(800e-9,epsilonAg800,$4*Unit+$5*Imaginary))) w p lw 4 lc 1 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -31.0213e0 0.4095e0 800" u 2:(1e9*DecayDepth(kzSPP(800e-9,epsilonAg800,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -31.0213e0 0.4095e0 800" u 2:(1e9*DecayDepth(kzSPP(800e-9,$4*Unit+$5*Imaginary, epsilonAg800))) w p lw 4 lc 3 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -31.0213e0 0.4095e0 800" u 2:(1e9*DecayDepth(kzSPP(800e-9,$4*Unit+$5*Imaginary, epsilonAg800))):1 w labels center offset 3.0,0.0 notitle

set output '20141208-SPPfriendlyMaterials-DecayDepth-BandGap-Ag400.eps'
set terminal postscript eps enhanced color font 'Times, 20' #size 8cm,8cm

set title 'Ag, 400 nm'
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -4.4222e0 0.2103e0 400" u 2:(1e9*DecayDepth(kzSPP(400e-9,epsilonAg400,$4*Unit+$5*Imaginary))) w p lc 1 lw 4 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -4.4222e0 0.2103e0 400" u 2:(1e9*DecayDepth(kzSPP(400e-9,epsilonAg400,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -4.4222e0 0.2103e0 400" u 2:(1e9*DecayDepth(kzSPP(400e-9,$4*Unit+$5*Imaginary,epsilonAg400))) w p lw 4 lc 3 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -4.4222e0 0.2103e0 400" u 2:(1e9*DecayDepth(kzSPP(400e-9,$4*Unit+$5*Imaginary,epsilonAg400))):1 w labels center offset 3.0,0.0 notitle

set output '20141208-SPPfriendlyMaterials-DecayDepth-BandGap-Ti800.eps'
set terminal postscript eps enhanced color font 'Times, 20' #size 8cm,8cm

set title 'Ti, 800 nm'
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -6.2067e0 25.2004e0 800" u 2:(1e9*DecayDepth(kzSPP(800e-9,epsilonTi800,$4*Unit+$5*Imaginary))) w p lw 4 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -6.2067e0 25.2004e0 800" u 2:(1e9*DecayDepth(kzSPP(800e-9,epsilonTi800,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -6.2067e0 25.2004e0 800" u 2:(1e9*DecayDepth(kzSPP(800e-9,$4*Unit+$5*Imaginary,epsilonTi800))) w p lw 4 lc 3 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -6.2067e0 25.2004e0 800" u 2:(1e9*DecayDepth(kzSPP(800e-9,$4*Unit+$5*Imaginary,epsilonTi800))):1 w labels center offset 3.0,0.0 notitle

set output '20141208-SPPfriendlyMaterials-DecayDepth-BandGap-Ti400.eps'
set terminal postscript eps enhanced color font 'Times, 20' #size 8cm,8cm

set title 'Ti, 400 nm'
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -4.3620e0 12.3621e0 400" u 2:(1e9*DecayDepth(kzSPP(400e-9,epsilonTi400,$4*Unit+$5*Imaginary))) w p lw 4 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -4.3620e0 12.3621e0 400" u 2:(1e9*DecayDepth(kzSPP(400e-9,epsilonTi400,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -4.3620e0 12.3621e0 400" u 2:(1e9*DecayDepth(kzSPP(400e-9,$4*Unit+$5*Imaginary,epsilonTi400))) w p lw 4 lc 3 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv -4.3620e0 12.3621e0 400" u 2:(1e9*DecayDepth(kzSPP(400e-9,$4*Unit+$5*Imaginary,epsilonTi400))):1 w labels center offset 3.0,0.0 notitle


set output '20141208-SPPfriendlyMaterials-DecayDepth-BandGap-W800.eps'
set terminal postscript eps enhanced color font 'Times, 20' #size 8cm,8cm

set title 'W, 800 nm'
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 4.4655e0 19.905e0 800" u 2:(1e9*DecayDepth(kzSPP(800e-9,epsilonW800,$4*Unit+$5*Imaginary))) w p lw 4 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 4.4655e0 19.905e0 800" u 2:(1e9*DecayDepth(kzSPP(800e-9,epsilonW800,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 4.4655e0 19.905e0 800" u 2:(1e9*DecayDepth(kzSPP(800e-9,$4*Unit+$5*Imaginary,epsilonW800))) w p lw 4 lc 3 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 4.4655e0 19.905e0 800" u 2:(1e9*DecayDepth(kzSPP(800e-9,$4*Unit+$5*Imaginary,epsilonW800))):1 w labels center offset 3.0,0.0 notitle

set output '20141208-SPPfriendlyMaterials-DecayDepth-BandGap-W400.eps'
set terminal postscript eps enhanced color font 'Times, 20' #size 8cm,8cm

set title 'W, 400 nm'
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 4.0653e0 16.453e0 400" u 2:(1e9*DecayDepth(kzSPP(400e-9,epsilonW400,$4*Unit+$5*Imaginary))) w p lw 4 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 4.0653e0 16.453e0 400" u 2:(1e9*DecayDepth(kzSPP(400e-9,epsilonW400,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 4.0653e0 16.453e0 400" u 2:(1e9*DecayDepth(kzSPP(400e-9,$4*Unit+$5*Imaginary,epsilonW400))) w p lw 4 lc 3 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 4.0653e0 16.453e0 400" u 2:(1e9*DecayDepth(kzSPP(400e-9,$4*Unit+$5*Imaginary,epsilonW400))):1 w labels center offset 3.0,0.0 notitle

set output '20141208-SPPfriendlyMaterials-DecayDepth-BandGap-Ge800.eps'
set terminal postscript eps enhanced color font 'Times, 20' #size 8cm,8cm

set title 'Ge, 800 nm'
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 22.0765e0 3.0347e0 800" u 2:(1e9*DecayDepth(kzSPP(800e-9,epsilonGe800,$4*Unit+$5*Imaginary))) w p lw 4 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 22.0765e0 3.0347e0 800" u 2:(1e9*DecayDepth(kzSPP(800e-9,epsilonGe800,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 22.0765e0 3.0347e0 800" u 2:(1e9*DecayDepth(kzSPP(800e-9,$4*Unit+$5*Imaginary,epsilonGe800))) w p lw 4 lc 3 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 22.0765e0 3.0347e0 800" u 2:(1e9*DecayDepth(kzSPP(800e-9,$4*Unit+$5*Imaginary,epsilonGe800))):1 w labels center offset 3.0,0.0 notitle

set output '20141208-SPPfriendlyMaterials-DecayDepth-BandGap-Ge400.eps'
set terminal postscript eps enhanced color font 'Times, 20' #size 8cm,8cm

set title 'Ge, 400 nm'
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 12.2406e0 18.3387e0 400" u 2:(1e9*DecayDepth(kzSPP(800e-9,epsilonGe400,$4*Unit+$5*Imaginary))) w p lw 4 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 12.2406e0 18.3387e0 400" u 2:(1e9*DecayDepth(kzSPP(800e-9,epsilonGe400,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 12.2406e0 18.3387e0 400" u 2:(1e9*DecayDepth(kzSPP(800e-9,$4*Unit+$5*Imaginary,epsilonGe400))) w p lw 4 lc 3 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 12.2406e0 18.3387e0 400" u 2:(1e9*DecayDepth(kzSPP(800e-9,$4*Unit+$5*Imaginary,epsilonGe400))):1 w labels center offset 3.0,0.0 notitle


set output '20141208-SPPfriendlyMaterials-DecayDepth-BandGap-Si800.eps'
set terminal postscript eps enhanced color font 'Times, 20' #size 8cm,8cm

set title 'Si, 800 nm'
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 13.64e0 0.048e0 800" u 2:(1e9*DecayDepth(kzSPP(800e-9,epsilonSi800,$4*Unit+$5*Imaginary))) w p lw 4 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 13.64e0 0.048e0 800" u 2:(1e9*DecayDepth(kzSPP(800e-9,epsilonSi800,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 13.64e0 0.048e0 800" u 2:(1e9*DecayDepth(kzSPP(800e-9,$4*Unit+$5*Imaginary,epsilonSi800))) w p lw 4 lc 3 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 13.64e0 0.048e0 800" u 2:(1e9*DecayDepth(kzSPP(800e-9,$4*Unit+$5*Imaginary,epsilonSi800))):1 w labels center offset 3.0,0.0 notitle


set output '20141208-SPPfriendlyMaterials-DecayDepth-BandGap-Si400.eps'
set terminal postscript eps enhanced color font 'Times, 20' #size 8cm,8cm

set title 'Si, 400 nm'
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 30.8469e0 4.2994e0 400" u 2:(1e9*DecayDepth(kzSPP(400e-9,epsilonSi400,$4*Unit+$5*Imaginary))) w p lw 4 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 30.8469e0 4.2994e0 400" u 2:(1e9*DecayDepth(kzSPP(400e-9,epsilonSi400,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 30.8469e0 4.2994e0 400" u 2:(1e9*DecayDepth(kzSPP(400e-9,$4*Unit+$5*Imaginary,epsilonSi400))) w p lw 4 lc 3 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 30.8469e0 4.2994e0 400" u 2:(1e9*DecayDepth(kzSPP(400e-9,$4*Unit+$5*Imaginary,epsilonSi400))):1 w labels center offset 3.0,0.0 notitle


set output '20141208-SPPfriendlyMaterials-DecayDepth-BandGap-SiO2800.eps'
set terminal postscript eps enhanced color font 'Times, 20' #size 8cm,8cm

set title 'SiO_2, 800 nm'
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 2.1121e0 0e0 800" u 2:(1e9*DecayDepth(kzSPP(800e-9,epsilonSiO2800,$4*Unit+$5*Imaginary))) w p lw 4 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 2.1121e0 0e0 800" u 2:(1e9*DecayDepth(kzSPP(800e-9,epsilonSiO2800,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 2.1121e0 0e0 800" u 2:(1e9*DecayDepth(kzSPP(800e-9,$4*Unit+$5*Imaginary,epsilonSiO2800))) w p lw 4 lc 3 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 2.1121e0 0e0 800" u 2:(1e9*DecayDepth(kzSPP(800e-9,$4*Unit+$5*Imaginary,epsilonSiO2800))):1 w labels center offset 3.0,0.0 notitle

set output '20141208-SPPfriendlyMaterials-DecayDepth-BandGap-SiO2400.eps'
set terminal postscript eps enhanced color font 'Times, 20' #size 8cm,8cm

set title 'SiO_2, 400 nm'
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 2.1612e0 0e0 400" u 2:(1e9*DecayDepth(kzSPP(400e-9,epsilonSiO2400,$4*Unit+$5*Imaginary))) w p lw 4 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 2.1612e0 0e0 400" u 2:(1e9*DecayDepth(kzSPP(400e-9,epsilonSiO2400,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 2.1612e0 0e0 400" u 2:(1e9*DecayDepth(kzSPP(400e-9,$4*Unit+$5*Imaginary,epsilonSiO2400))) w p lw 4 lc 3 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 2.1612e0 0e0 400" u 2:(1e9*DecayDepth(kzSPP(400e-9,$4*Unit+$5*Imaginary,epsilonSiO2400))):1 w labels center offset 3.0,0.0 notitle


set output '20141208-SPPfriendlyMaterials-DecayDepth-BandGap-Air800.eps'
set terminal postscript eps enhanced color font 'Times, 20' #size 8cm,8cm

set title 'Air, 800 nm'
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 1e0 0e0 800" u 2:(1e9*DecayDepth(kzSPP(800e-9,epsilonAir,$4*Unit+$5*Imaginary))) w p lw 4 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 1e0 0e0 800" u 2:(1e9*DecayDepth(kzSPP(800e-9,epsilonAir,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 1e0 0e0 800" u 2:(1e9*DecayDepth(kzSPP(800e-9,$4*Unit+$5*Imaginary,epsilonAir))) w p lw 4 lc 3 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 1e0 0e0 800" u 2:(1e9*DecayDepth(kzSPP(800e-9,$4*Unit+$5*Imaginary,epsilonAir))):1 w labels center offset 3.0,0.0 notitle

set output '20141208-SPPfriendlyMaterials-DecayDepth-BandGap-Air400.eps'
set terminal postscript eps enhanced color font 'Times, 20' #size 8cm,8cm

set title 'Air, 400 nm'
plot "< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 1e0 0e0 400" u 2:(1e9*DecayDepth(kzSPP(400e-9,epsilonAir,$4*Unit+$5*Imaginary))) w p lw 4 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 1e0 0e0 400" u 2:(1e9*DecayDepth(kzSPP(400e-9,epsilonAir,$4*Unit+$5*Imaginary))):1 w labels center offset 3.0,0.0 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 1e0 0e0 400" u 2:(1e9*DecayDepth(kzSPP(400e-9,$4*Unit+$5*Imaginary,epsilonAir))) w p lw 4 lc 3 notitle, \
"< sh SPPfriendlySelect.sh MaterialOpticalDatabaseForPlasmonics.csv 1e0 0e0 400" u 2:(1e9*DecayDepth(kzSPP(400e-9,$4*Unit+$5*Imaginary,epsilonAir))):1 w labels center offset 3.0,0.0 notitle


unset multiplot
