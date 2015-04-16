#! /bin/gawk
gawk -v eps1re=$2 -v eps1im=$3 -v lambda=$4 '
{
if( (eps1re*$4+eps1im*$5 < 0) && ($3==lambda) && ($4*eps1im-eps1re*$5 < 0) ) {print $0}
}
' $1