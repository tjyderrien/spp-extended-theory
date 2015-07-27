#!gawk

BEGIN {sumX=0e0; sumY=0e0; sumZ=0e0; }
{sumX=sumX+$2**2; sumY=sumY+$3**2; sumZ=sumZ+$4**2;}
END {print sumX, sumY, sumZ}
