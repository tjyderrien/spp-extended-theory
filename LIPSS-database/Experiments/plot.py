#!/usr/bin/env python
#-*- coding: utf-8 -*-
from SimpleSPProutines import *
#reset

print "==== Plot the exp. LIPSS regularity as function of theoretical Lspp... ===="

print "| Reading the experimental file..."
database = "Iaroslav/ExperimentalData_1030nm_GaussianFit.csv"
ExpData = loadtxt(database, dtype='str', delimiter='\t', skiprows=4)

#print ExpData
print "" 

print "| Selecting the listed materials..."
Materials = ExpData[:,0];
Sources = ExpData[:,10];
print "|| Materials = "+str(Materials)
print "|| Sources = "+str(Sources)
print "" 

#TODO: make 2 object: (i) Literature (ii) Iaroslav points. 

print "| Extracting the dispersion of available exp. data..."
DispersionAngle = np.asfarray(ExpData[:,1]);
DispersionAngleError = np.asfarray(ExpData[:,2]);
print DispersionAngle
print "" 

print "| Selecting the laser wavelengths..."
Wavelengths = ExpData[:,3]
print "|| Wavelengths = "+str(Wavelengths)
print "" 

print "| Grabbing the available theoretical database..."

unit = 1E9

# Select database
database="../../SimpleSPP/MaterialOpticalDatabaseForPlasmonics.csv"
dbarray = loadtxt(database, dtype='str', delimiter='\t')

plt.figure()
ax = plt.subplot(111)
plt.title('Origin of LIPSS regularity')
plt.xlabel(r'$L_{SPP}$ ($\mu$m)')
plt.ylabel(r'$\delta \theta$ (deg)')

#print "Size="+str(Materials.size)
Lspp = np.zeros(0)
DeltaLspp = np.zeros(0)
errorX = np.zeros(len(Materials))

for i in np.arange(0, Materials.size, 1):
  wavelength = float(Wavelengths[i])
  material = Materials[i]
  source = Sources[i]
  print "|| Treating (wavelength, material) = "+str(wavelength)+", "+str(material)+"."

  #=====
  #print "|| Selecting the material of interfaces..."
  DataMaterial1 = FilterDatabaseContains(dbarray, "Air", 0)
  print "** Material1 filtering #"+str(i)+" returned "+str(len(DataMaterial1))+" entries."
  
  print "|| Selecting the wavelength..."
  DataMaterial1 = FilterDatabase(DataMaterial1, str(int(wavelength)), 2)
  print "** Material1 filtering #"+str(i)+" returned "+str(len(DataMaterial1))+" entries."
  print DataMaterial1
  print ""
  
  DataMaterial2 = FilterDatabaseContains(dbarray, str(material), 0)
  print "** Material2 ("+str(material)+") filtering #"+str(i)+" returned "+str(len(DataMaterial2))+" entries."
  #print DataMaterial2
  print ""
  
  DataMaterial2 = FilterDatabase(DataMaterial2, str(int(wavelength)), 2)
  print "** Wavelength ("+str(int(wavelength))+") filtering returned "+str(len(DataMaterial2))+" entries."
  #print DataMaterial2
  # Check that number of solutions is one for each research. 
  
  #DataMaterial2 = FilterDatabaseContains(DataMaterial2, str(source), 0)
  #print "** Source ("+str(source)+") filtering #"+str(i)+" returned "+str(len(DataMaterial2))+" entries."
  #print DataMaterial2
  #print ""
  
  #print "|| Material 2 data = "+str(DataMaterial2)

  print "|| Extracting epsilons..."
  #print DataMaterial1
  Material1=DataMaterial1[0,0]
  BandGap1=DataMaterial1[0,1] 
  wavelength1=DataMaterial1[0,2]
  RealEps1=DataMaterial1[0,3]
  ImagEps1= DataMaterial1[0,4]
  
  Material2=DataMaterial2[0,0];
  BandGap2=DataMaterial2[0,1];
  wavelength2=DataMaterial2[0,2];
  RealEps2=DataMaterial2[0,3];
  ImagEps2=DataMaterial2[0,4];

  #print RealEps1, ImagEps1
  #====
  #print "Building theoretical Lspp..."
  eps1 = float(RealEps1) + 1.0j*float(ImagEps1)
  eps2 = float(RealEps2) + 1.0j*float(ImagEps2)
  
  #print "Debug: wavelength2 = ", wavelength2
  
  NewLspp = DecayLengthSPP(betaSPP(float(wavelength)/unit, eps1, eps2))
  NewDeltaLspp = deltaLspp(float(wavelength)/unit, eps1, eps2, 0e0, 0e0, eta, eta)
  Lspp = np.append(Lspp, NewLspp)
  DeltaLspp = np.append(DeltaLspp, NewDeltaLspp)
  print "[New entry] Material: "+str(Material2)+", eps2="+str(eps2)+", Lspp="+str(NewLspp)+", dLspp="+str(NewDeltaLspp)

print "" 
print "====== FINAL RESULTS ======"
print "Dispersion angle: "+str(DispersionAngle)
print "Lspp: "+str(Lspp)
print "dLspp: "+str(DeltaLspp)
print ""
print "| Plotting the results..."
plt.errorbar(1E6*Lspp, DispersionAngle, yerr=DispersionAngleError, fmt='.')
#plt.errorbar(1E6*Lspp, DispersionAngle, xerr=1e6*DeltaLspp, yerr=DispersionAngleError, fmt='.')
ax.set_xscale('log')
#ax.view([])
#ax.set_yscale('log')
#plt.show()
plt.savefig('OriginOfRegularity.eps')
print "Figure OriginOfRegularity.eps was saved successfully."