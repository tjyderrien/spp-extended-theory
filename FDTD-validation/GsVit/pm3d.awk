#! gawk

BEGIN {
lambda=515e-9
pi=3.141592654
Xvalue=0
Xrange=40e0
Trange=10e-9
Yvalue=50e0
# Xrange=2*lambda
# Yrange=Xrange #BonseSipe fourier space
Yrange=20e-6 #real space
}

{
####### CONSTRUCTION DE CARTES X Z
# On parcourt le fichier : quel que soit X, on parcourt Z
# on voudrait sauter une ligne lorsque X change pour faire du PM3D
# if($2!=Xvalue && $2<Xrange && $2>-Xrange && $3<Yrange) {
if($1!=Yvalue) { print ""; }
# }
#definir la valeur de X
Yvalue=$1; 

# afficher la ligne dans tous les cas que l'on veut afficher
# if($3<Xrange && $3>-Xrange && $4<Yrange && $4>-Yrange) print $3,$4,$8/lambda,"i"; 
# if($2<Xrange && $2>-Xrange && $3<Yrange) {
	print $0;
# }

}
