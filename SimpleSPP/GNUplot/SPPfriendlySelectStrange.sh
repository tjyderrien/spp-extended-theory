#! /bin/gawk
gawk -v eps1re=$2 -v eps1im=$3 -v lambda=$4 '
{
# $4 will be eps2re, $5 will be eps2im

if( eps1im=0e0 && $5==0e0 ) { #case of perfect materials
  if( (eps1re*$4 < 0) && ($3==lambda) && (eps1re < - $4)) {print $0, eps1re*$4+eps1im*$5}
}
else { #case of at least one absorbing materials
  if( (eps1re*$4+eps1im*$5 < 0) && ($3==lambda) ) {print $0, eps1re*$4+eps1im*$5}
  
}
# condition 2 also
# if( (eps1re*$4+eps1im*$5 < 0) && ($3==lambda) && ($4*eps1im-eps1re*$5 < 0) ) {print $0}
}
' $1