#! gnuplot

unset multiplot

reset
set terminal postscript eps enhanced color font "Helvetica, 24" size 3.5, 2.62


# ---- Params ----



#		# --------------------------------------#
#     	  	# --------------------------------------#
#     	  	#					#
#     		#                FILES         		#
#     		#					#
#     		# --------------------------------------#
#     		# --------------------------------------#



DLOA_Ti_OrigSEM  = "/home/yo/WORK/HiLASE/STUDIES/SPP-Extended-Theory/spp-extended-theory/LIPSS-database/Experiments/Iaroslav/OriginalSEM/OrientationComparison/Al_DLOA_SizeSquaredArea.dat"

DLOA_Al_OrigSEM  = "/home/yo/WORK/HiLASE/STUDIES/SPP-Extended-Theory/spp-extended-theory/LIPSS-database/Experiments/Iaroslav/OriginalSEM/OrientationComparison/Ti_DLOA_SizeSquaredArea.dat"



# ---- Distribution of angles from picts




AngleDistrib_Al_5um      = "/home/yo/WORK/HiLASE/STUDIES/SPP-Extended-Theory/spp-extended-theory/LIPSS-database/Experiments/Iaroslav/OriginalSEM/Al/20161110_Al_Square_05.03um.txt"
AngleDistrib_Al_7p5um    = "/home/yo/WORK/HiLASE/STUDIES/SPP-Extended-Theory/spp-extended-theory/LIPSS-database/Experiments/Iaroslav/OriginalSEM/Al/20161110_Al_Square_07.51um.txt"
AngleDistrib_Al_10um     = "/home/yo/WORK/HiLASE/STUDIES/SPP-Extended-Theory/spp-extended-theory/LIPSS-database/Experiments/Iaroslav/OriginalSEM/Al/20161110_Al_Square_10.00um.txt"
AngleDistrib_Al_12p5um   = "/home/yo/WORK/HiLASE/STUDIES/SPP-Extended-Theory/spp-extended-theory/LIPSS-database/Experiments/Iaroslav/OriginalSEM/Al/20161110_Al_Square_12.54um.txt"
AngleDistrib_Al_15um     = "/home/yo/WORK/HiLASE/STUDIES/SPP-Extended-Theory/spp-extended-theory/LIPSS-database/Experiments/Iaroslav/OriginalSEM/Al/20161110_Al_Square_15.03um.txt"
AngleDistrib_Al_17p5um   = "/home/yo/WORK/HiLASE/STUDIES/SPP-Extended-Theory/spp-extended-theory/LIPSS-database/Experiments/Iaroslav/OriginalSEM/Al/20161110_Al_Square_17.46um.txt"
AngleDistrib_Al_20um     = "/home/yo/WORK/HiLASE/STUDIES/SPP-Extended-Theory/spp-extended-theory/LIPSS-database/Experiments/Iaroslav/OriginalSEM/Al/20161110_Al_Square_20.00um.txt"
AngleDistrib_Al_25um     = "/home/yo/WORK/HiLASE/STUDIES/SPP-Extended-Theory/spp-extended-theory/LIPSS-database/Experiments/Iaroslav/OriginalSEM/Al/20161110_Al_Square_25.03um.txt"
AngleDistrib_Al_30um     = "/home/yo/WORK/HiLASE/STUDIES/SPP-Extended-Theory/spp-extended-theory/LIPSS-database/Experiments/Iaroslav/OriginalSEM/Al/20161110_Al_Square_30.00um.txt"
AngleDistrib_Al_35um     = "/home/yo/WORK/HiLASE/STUDIES/SPP-Extended-Theory/spp-extended-theory/LIPSS-database/Experiments/Iaroslav/OriginalSEM/Al/20161110_Al_Square_34.97um.txt"
AngleDistrib_Al_40um     = "/home/yo/WORK/HiLASE/STUDIES/SPP-Extended-Theory/spp-extended-theory/LIPSS-database/Experiments/Iaroslav/OriginalSEM/Al/20161110_Al_Square_40.00um.txt"
AngleDistrib_Al_45um     = "/home/yo/WORK/HiLASE/STUDIES/SPP-Extended-Theory/spp-extended-theory/LIPSS-database/Experiments/Iaroslav/OriginalSEM/Al/20161110_Al_Square_44.97um.txt"
AngleDistrib_Al_50um     = "/home/yo/WORK/HiLASE/STUDIES/SPP-Extended-Theory/spp-extended-theory/LIPSS-database/Experiments/Iaroslav/OriginalSEM/Al/20161110_Al_Square_50.05um.txt"





#		# --------------------------------------#
#     	  	# --------------------------------------#
#     	  	#					#
#     		#                PLOTS         		#
#     		#					#
#     		# --------------------------------------#
#     		# --------------------------------------#


#   ----   DLOA fctn of the size of the squared area   ----   #
#   ----   DLOA fctn of the size of the squared area   ----   #
#   ----   DLOA fctn of the size of the squared area   ----   #



set output "FIG_DLOA_fctn_SizeSquaredArea_OrigSEM_TiAl_20161110_2.eps"

set xlabel "Side length a ({/Symbol m}m)"

set ylabel "DLOA {/Symbol dq} (degree)"

plot \
DLOA_Ti_OrigSEM  u 1:2 w p  pointsize 1.5 lw 2 lt 5 lc 3 title "Ti", \
DLOA_Ti_OrigSEM  u 1:2 w l lw 1 lt 1 lc 3 notitle "Ti", \
DLOA_Al_OrigSEM  u 1:2 w p  pointsize 1.5 lw 2 lt 7 lc 2 title "Al", \
DLOA_Al_OrigSEM  u 1:2 w l lw 1 lt 1 lc 2 notitle "Al", \




#   ----   Angular Orientation distribution   ----   #
#   ----   Angular Orientation distribution   ----   #
#   ----   Angular Orientation distribution   ----   #

set output "FIG_AngularDistrib_OrigSEM_Al_20161110_2.eps"

set key left top

set xlabel "{/Symbol q} (degree)"
set xrange  [-90:90]

set ylabel "Counts (a.u.)"

plot \
AngleDistrib_Al_10um  u 1:4 w l lw 1 lt 1 lc 1 title "10 {/Symbol m}m", \
AngleDistrib_Al_20um  u 1:4 w l lw 3 lt 2 lc 4 title "20 {/Symbol m}m", \
AngleDistrib_Al_50um  u 1:4 w l lw 4 lt 1 lc 0 title "50 {/Symbol m}m", \
0.5 w l lw 2 lt 2 lc 9 notitle, \
#  AngleDistrib_Al_15um  u 1:4 w l lw 2 lt 2 lc 3 title "15 {/Symbol m}m", \
#  AngleDistrib_Al_30um  u 1:4 w l lw 6 lt 0 lc 2 title "30 {/Symbol m}m", \
#  AngleDistrib_Al_40um  u 1:4 w l lw 6 lt 0 lc 8 title "40 {/Symbol m}m", \

