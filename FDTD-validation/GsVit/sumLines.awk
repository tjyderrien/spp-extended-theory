#!gawk

BEGIN {sumX=0e0; sumY=0e0; sumZ=0e0; }
{

# Field intensity is calcualted. Integral on space and time is calculated.
# How to calculate reflectivity actually ? 

sumX=sumX+$2**2; 
sumY=sumY+$3**2; 
sumZ=sumZ+$4**2;

# We can normalize to the sum of energy we have into whole diffracted energy

}
END {print sumX, sumY, sumZ}
