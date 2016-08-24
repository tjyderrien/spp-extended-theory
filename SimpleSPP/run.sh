#!/bin/bash

echo "Multi-material multi-wavelength plotting script."
echo "Usage: modify the script, and ./run.sh."

Materials="Au Ag Ti"
Sources="Palik Johnson"

for iMaterials in $Materials
do
 for iSource in $Sources
 do
  ./plotMultiwavelength.py $iMaterials $iSource --no-show
 done
done
