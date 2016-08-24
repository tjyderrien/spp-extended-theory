#!/bin/bash

echo "Multi-material multi-wavelength plotting script."
echo "Usage: modify the script, and ./run.sh."
echo "Workload is parallelized using bash."

Materials="Au Ag Ti"
Sources="Palik Johnson"

echo "Materials: $Materials"
echo "Sources: $Sources"

for iMaterials in $Materials
do
 for iSource in $Sources
 do
  ./plotMultiwavelength.py $iMaterials $iSource 1E-9 --no-show &
 done
done
