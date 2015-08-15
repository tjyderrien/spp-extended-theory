#!/bin/bash

SPPexttheoryPath="/home/thibault/Documents/Codes/SPP-extended-theory/FDTD-validation/GsVit"

rm Diffraction.dat
echo "# theta, phi, Ix, Iy, Iz" >> Diffraction.dat
for i in `ls xnff* -1` 
do
	# collect the theta, phi coordinates from filename
	theta=`echo $i | tr "_" " " | awk '{ print $3}' `
	phi=`echo $i | tr "_" " " | awk '{ print $4}' `
#	echo theta=$theta, phi=$phi
	
	# calculate the total field for each coordinate
	value=`awk -f ${SPPexttheoryPath}/sumLines.awk $i`
#	echo $value
	
	#export to a file
	echo "$theta $phi $value i" >> Diffraction.dat
done

# calculate the total field intensity (use it to normalize !!)
echo "Reflectivity per axis"
indexTheta=`tail -n1 Diffraction.dat | awk '{ print $1 }' `
Rx=`awk -v var1=$indexTheta '{ if($1==var1 && $2=="000") print $3 }' Diffraction.dat | head -n1 `

Ry=`awk -v var1=$indexTheta '{ if($1==var1 && $2=="000") print $4 }' Diffraction.dat | head -n1 `

Rz=`awk -v var1=$indexTheta '{ if($1==var1 && $2=="000") print $5 }' Diffraction.dat | head -n1 `

echo $Rx, $Ry, $Rz

# let's sum every diffraction angles
echo "Total intensity per axis"
Total=`awk 'BEGIN {sumX=0e0; sumY=0e0; sumZ=0e0; } { if($2=="000") { sumX=sumX+$3; sumY=sumY+$4; sumZ=sumZ+$5;} } END {print sumX, sumY, sumZ}' Diffraction.dat`
TotalX=`echo $Total | awk '{ print $1 }' `
TotalY=`echo $Total | awk '{ print $2 }' `
TotalZ=`echo $Total | awk '{ print $3 }' `

echo $TotalX, $TotalY, $TotalZ
