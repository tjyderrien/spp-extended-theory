#!/usr/bin/env python3
#-*- coding: utf-8 -*-

# Copyright (C) 2013-2024 T. J.-Y. Derrien
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>

## @Packagee plotMultiwavelength 
# Prepares plots presenting the properties of SPP as function of numerous materials. 

# IMPORT LIBRARIES
from spp_extended_theory.Libs.libMaterials import *
import sys
from spp_extended_theory.Libs.libSPP import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline, interp1d

precision = 1E-6
interpolation_method_list={0: "spline linear",
                           1: "interp1d nearest", #works well on a small interval
                           2: "interp1d linear", #works well on a small interval
                           3: "interp1d cubic", #works well on a small interval
                           4: "spline cubic" #can generates singularities on wide spectral range"spline cubic", #can generates singularities on wide spectral range
                           }
interpolation_method = interpolation_method_list[0]
# Now, we shall construct database for SPP lifetimes. Actually, SPP lifetime require the knowledge of all spectrum of response to be known.

# =======================================================
if(len(sys.argv)<=1):
  print("Usage: ./plotMultiwavelength.py           \ ")
  print("    <Name of the substrate (Air, Be, Au, ...)> \ ")
  print("    <Source for data: Palik or name of the 1st author> \ ")
  print("    [<precision>: 1E-9 by default>] \ ")
  print("    [--no-show]")
  print("Example: ./plotMultiwavelength.py Au Johnson")
  sys.exit()

#============ Manage the command line input =================
query = sys.argv[1]
#wavelength = 1E-9*float(sys.argv[2])
try:
  source = "-"+sys.argv[2]
except:
  source = ""

try:
  precision = float(sys.argv[3])
except:
  print("** Warning: precision was choosed by default: 1E-9 m")

# Command line option for avoiding the visual plotting...
try:
  if(sys.argv[4] == "--no-show"):
    ShowPictures = False
  else:
    ShowPictures = True
except: 
  print("** Warning: No-show command was not defined.")
  ShowPictures = True

MaterialFile2 = query+source

# ============= BYPASSING COMMAND LINE =============

SpecialWavelength=800e-9 #this wavelength will be indicated on the graphs
xmin=0   #nm
xmax=2E3 #nm
MaterialFolder="Database"

MaterialFile1="Air"
#MaterialFile1="Al2O3-Palik"
#MaterialFile1="SiO2-Palik"
#MaterialFile1="TiO2-Palik"
#MaterialFile1="Si-Palik"
#MaterialFile1="ZnO-Bond"

#MaterialFile1="Ag-Johnson"
#MaterialFile1="Ag-Palik"

#MaterialFile2="Al-Palik"
#MaterialFile2="Ti-Palik"
#MaterialFile2="Ag-Palik"
#MaterialFile2="Au-Palik"
#MaterialFile2="Au-Johnson"

#MaterialFile2="SiO2-Palik"
#MaterialFile2="Ti-Johnson"
#MaterialFile2="Mo-Palik"
#MaterialFile2="Si-Palik"

# Managing source files units (quite artificial...)
UnitMat1=1e10
if(source == "-Palik"):
  UnitMat2=1e10
  print("** Info: Palik optical data selected...")
elif(source == "-Johnson"):
  UnitMat2=1e6
  print("** Info: Johnson optical data selected...")
elif(source == "-GoriAndBond"):
  UnitMat2=1e9
  print("** Info: GoriAndBond optical data was selected...")  
else:
  print("** Warning: Rare source of optical data was selected... ")
  UnitMat2=1e6

#================ Loading Material dielectric complex permittivity into arrays ====================
# try:
MaterialArray2 = np.loadtxt(MaterialFolder+'/'+MaterialFile2, delimiter='\t', skiprows=4)
"""except: 
	print(("** Warning: Could not read "+MaterialFile2+" database."))
	print("** Warning: Attempting second method...")
	try:
		MaterialArray2 = np.loadtxt(MaterialFolder+'/'+MaterialFile2, delimiter=' ', skiprows=4)
		print("Success.")
	except:
		print("** Error: Also failed reading of database... Exiting.")
		sys.exit()"""
		
#wavelengths2 = MaterialArray2[:,0]
nlines, ncols = MaterialArray2.shape

#print nlines, ncols

try:
  MaterialArray1 = np.loadtxt(MaterialFolder+'/'+MaterialFile1, delimiter='\t', skiprows=4)
except:
  try:
    MaterialArray1 = np.loadtxt(MaterialFolder+'/'+MaterialFile1, delimiter=' ', skiprows=4)
  except:
    print(("** Warning: Material 1 ("+MaterialFile1+") was not found in "+MaterialFolder+"."))
    print("** Warning: Material 1 was replaced by Air.")
    
    MaterialArray1 = np.zeros((nlines, 3))
    MaterialArray1[:,1] = np.ones(nlines) #Air index is 1. 
    MaterialArray1[:,0] = MaterialArray2[:,0] #Copy the table of wavelengths

#print "Shape of MaterialArray1 is "+str(MaterialArray1.shape)
#print "Shape of MaterialArray2 is "+str(MaterialArray2.shape)

# Let's build eps1 and eps2
wavelengths1=MaterialArray1[:,0]
wavelengths1=np.multiply(wavelengths1,UnitMat1**(-1)) #converting wavelength to emters
n=MaterialArray1[:,1]; k=MaterialArray1[:,2]
eps1=np.power(np.add(n,np.multiply(1j, k)),2) #conversion to epsilon
del n, k
#print wavelengths1

wavelengths2=MaterialArray2[:,0]
wavelengths2=np.multiply(wavelengths2,UnitMat2**(-1)) #converting wavelength to meters
n=MaterialArray2[:,1]; k=MaterialArray2[:,2]

if(source=="-GoriAndBond"):
  eps2=np.add(n,np.multiply(1j, k)) #the value read in the file is ALREADY an epsilon!
else: 
  eps2=np.power(np.add(n,np.multiply(1j, k)),2)
del n, k
#print eps2

# Before doing calculations, we shall interpolate the most dense mesh on the second, and take their intersection. 
print(("Wavelength mesh size 1="+str(eps1.size)))
print(("Wavelength mesh size 2="+str(eps2.size)))

# dense mesh generation
# TODO: Error on results when taking Palik data on wide spectrum! 
# wavelengths = np.arange(,,precision)

# Si - laser
# wavelengths_min = 0.21E-6 #np.amin(wavelengths2)
# wavelengths_max = 0.8E-6  #np.amax(wavelengths2)

# Water - laser
# wavelengths_min = 70e-9 #np.amin(wavelengths2)
# wavelengths_max = 2E-6  #np.amax(wavelengths2)

# Water - micro-waves
wavelengths_min = 1E-3 #np.amin(wavelengths2)
wavelengths_max = 0.1E0  #1E0  #0.8E-6  #np.amax(wavelengths2)

try:
  wavelengths = np.arange( wavelengths_min, wavelengths_max, 1e-10)
except MemoryError:
  wavelengths = np.linspace( wavelengths_min, wavelengths_max, 1_000_000)
print(wavelengths)
print("Checking if wavelength range is reasonable...")
if(wavelengths.size > 1E5 ):
  print("** Warning: interpolation may be very long to perform...")
  print("**        Attempting logscale.")
  wavelengths_log = np.arange(
    np.log10(wavelengths_min),
    np.log10(wavelengths_max),
    10)
  wavelengths     = np.power(10,wavelengths_log)
# sys.exit()"""


# if(interpolation_method in interpolation_method_list):
if(interpolation_method=="spline linear"):
  order=1
  print("** INFO: selected spline linear interpolation")
  feps1r=InterpolatedUnivariateSpline(wavelengths1, eps1.real, k=order)
  feps1i=InterpolatedUnivariateSpline(wavelengths1, eps1.imag, k=order)
  feps2r=InterpolatedUnivariateSpline(wavelengths2, eps2.real, k=order)
  feps2i=InterpolatedUnivariateSpline(wavelengths2, eps2.imag, k=order)
elif(interpolation_method=="spline cubic"):
  order=2
  print("** INFO: selected spline cubic interpolation")
  feps1r=InterpolatedUnivariateSpline(wavelengths1, eps1.real, k=order)
  feps1i=InterpolatedUnivariateSpline(wavelengths1, eps1.imag, k=order)
  feps2r=InterpolatedUnivariateSpline(wavelengths2, eps2.real, k=order)
  feps2i=InterpolatedUnivariateSpline(wavelengths2, eps2.imag, k=order)
elif(interpolation_method=="interp1d nearest"):
  print("** INFO: selected interp1d nearest interpolation")
  feps1r = interp1d(wavelengths1, eps1.real, kind="nearest") #, k=order)
  feps1i = interp1d(wavelengths1, eps1.imag, kind="nearest") #, k=order)
  feps2r = interp1d(wavelengths2, eps2.real, kind="nearest") #, k=order)
  feps2i = interp1d(wavelengths2, eps2.imag, kind="nearest") #, k=order)
elif(interpolation_method=="interp1d linear"):
  print("** INFO: selected interp1d linear interpolation")
  feps1r = interp1d(wavelengths1, eps1.real, kind="linear") #, k=order)
  feps1i = interp1d(wavelengths1, eps1.imag, kind="linear") #, k=order)
  feps2r = interp1d(wavelengths2, eps2.real, kind="linear") #, k=order)
  feps2i = interp1d(wavelengths2, eps2.imag, kind="linear") #, k=order)
elif(interpolation_method=="interp1d cubic"):
  print("** INFO: selected interp1d cubic interpolation")
  feps1r = interp1d(wavelengths1, eps1.real, kind="cubic") #, k=order)
  feps1i = interp1d(wavelengths1, eps1.imag, kind="cubic") #, k=order)
  feps2r = interp1d(wavelengths2, eps2.real, kind="cubic") #, k=order)
  feps2i = interp1d(wavelengths2, eps2.imag, kind="cubic") #, k=order)
"""else:
  print("** Error: interpolation method "+str(interpolation_method)+" is not taken into account. ")
  sys.exit()"""

print(("Interpolating on Wavelength mesh size = "+str(wavelengths.size)))

eps1new=np.add(feps1r(wavelengths),np.multiply(1.0j, feps1i(wavelengths)))
eps2new=np.add(feps2r(wavelengths),np.multiply(1.0j, feps2i(wavelengths)))

#========= Multiwavelength data: DATA ARE NOW READY ======

print("Checking quality of interpolation for the dielectric function...")
plt.figure()
plt.xlabel('Wavelength (m)')
plt.ylabel('epsilon')
plt.semilogx(wavelengths1, eps1.real, 'b+', label=r'Re('+MaterialFile1+')')
plt.semilogx(wavelengths, eps1new.real, 'b-', label='interp Re('+MaterialFile1+')')
plt.semilogx(wavelengths1, eps1.imag, 'k+', label='Im('+MaterialFile1+')')
plt.semilogx(wavelengths, eps1new.imag, 'k-', label='interp Im('+MaterialFile1+')')
plt.semilogx(wavelengths2, eps2.real, 'r+', label='Re('+MaterialFile2+')')
plt.semilogx(wavelengths, eps2new.real, 'r-', label='interp Re('+MaterialFile2+')')
plt.semilogx(wavelengths2, eps2.imag, 'g+', label='Im('+MaterialFile2+')')
plt.semilogx(wavelengths, eps2new.imag, 'g-', label='interp Im('+MaterialFile2+')')
plt.legend(loc=2)
plt.title(r'$\varepsilon(\lambda)$')
plt.savefig(MaterialFile1+MaterialFile2+'epsilon.png')
plt.show()

print("Plot the SPP dispersion relation...")

betaSPP = np.vectorize(betaSPP)
omega = np.vectorize(omega)

kspp = betaSPP(wavelengths, eps1new, eps2new)
omegaspp = omega(wavelengths) #we here set resonance omega_spp = omega_laser
#omegaspp = np.sort(omegaspp)

plt.figure()
plt.xlabel(r'k $(m^{-1})$')
plt.ylabel(r'$\omega$ ($s^{-1}$)')
plt.plot(kspp.real, omegaspp, label='SPP')
plt.plot(omega(wavelengths)/c, omega(wavelengths), label='Light line')
plt.plot(omega(wavelengths)/c, omega(np.add(np.multiply(wavelengths,0e0), SpecialWavelength)), label='Laser '+str(round(SpecialWavelength*1E9,0))+' nm')
plt.title(r'Dispersion at $'+MaterialFile1+r'$/$'+MaterialFile2+r'$ interface')
plt.legend(loc=2)
plt.savefig(MaterialFile1+MaterialFile2+'Dispersion.eps')
#plt.show()

print("Plot the SPP period with wavelength...")

plt.figure()
plt.xlabel(r'Wavelength $\lambda$ $(nm)$')
plt.ylabel(r'$\Lambda_{\mbox{SPP}}$ ($nm$)')
plt.loglog(1e9*wavelengths, 1e9 * (2e0*pi/kspp.real)) #label='Near-field period'
#plt.title('Period of field at $'+MaterialFile1+'$/$'+MaterialFile2+'$ interface')
plt.legend(loc=2)
plt.grid(True)
plt.savefig(MaterialFile1+MaterialFile2+'Period.eps')
#plt.show()

print("Plot the SPP mean-free path with wavelength...")

plt.figure()
plt.xlabel(r'Wavelength $\lambda$ $(nm)$')
plt.ylabel(r'SPP mean-free-path $L_{\mbox{SPP}}$ (\mbox{$\mu$m})')
plt.plot(1e9*wavelengths, 1e6 * (0.5E0/kspp.imag)) #label=MaterialFile1+'/'+MaterialFile2
#plt.title('Period of field at $'+MaterialFile1+'$/$'+MaterialFile2+'$ interface')
plt.legend(loc=2)
plt.grid(True)
plt.savefig(MaterialFile1+MaterialFile2+'MeanFreePath.svg')
#plt.show()
plt.plot(1e9*wavelengths, 1e6 * (0.5E0/kspp.imag), label=MaterialFile1+'/'+MaterialFile2)
plt.savefig(MaterialFile1+MaterialFile2+'MeanFreePath-LogLog.eps')

print("Plot the lifetime with wavelength...")
# RealDerivativeByComplex = np.vectorize(RealDerivativeByComplex)
SPPgroupVelocity = RealDerivativeByComplex(omegaspp, kspp)
SPPphaseVelocity = np.divide(omegaspp,kspp)

SPPgroupVelocityHohenau = HohenauGroupVelocity(omegaspp, kspp) #uses Hohenau / Jackson definition of group velocity for non-absorbing materials

print((np.shape(omegaspp)))
print((np.shape(kspp)))
print((np.shape(SPPgroupVelocity)))
print((np.shape(SPPgroupVelocityHohenau)))

SPPgroupVelocityPlot = np.clip(SPPgroupVelocity.real, 0, 1000E8)
SPPgroupVelocityHohenau = np.clip(SPPgroupVelocityHohenau, 0, 1000E8)

#print SPPgroupVelocity.real
#TODO: Group velocity can be negative, it may designates gain medium?
#See [Hohenau and JR Krenn, PRB 78, 155405 (2008)]
plt.figure()
plt.xlabel(r'Wavelength $\lambda$ (m)')
plt.ylabel(r'Velocity $v$ (m/s)')
plt.semilogx(2*pi*c/omegaspp[1:], SPPgroupVelocityPlot, 'b-', label=r'$Re(v_g)$')
plt.semilogx(2*pi*c/omegaspp[1:], SPPgroupVelocityHohenau, 'k-', label=r'$v_g[Re(\beta)]$')
plt.semilogx(2*pi*c/omegaspp, SPPphaseVelocity, 'r--', label=r'$v_{\phi}$')
#plt.plot(1e9*2*pi*c/omegaspp, c, label=r'$c$')
#plt.plot(omega(wavelengths)/c, omega(wavelengths), label='Light line')
#plt.plot(omega(wavelengths)/c, omega(np.add(np.multiply(wavelengths,0e0), SpecialWavelength)), label=r'$c$')
plt.title(MaterialFile1+'/'+MaterialFile2+' interface')
plt.legend(loc=4)
# plt.xticks(np.arange(0, 3500, 500))
# plt.axis([xmin,xmax,-300,300])
plt.savefig(MaterialFile1+MaterialFile2+'Velocities.eps')
plt.savefig(MaterialFile1+MaterialFile2+'Velocities.png')
plt.show()

print("Plot SPP lifetime with wavelength...")
LifeTimeOld = LifeTimeRaether(kspp[1:], eps1new[1:], eps2new[1:])
LifeTimeNew = LifeTimeVg(kspp[1:], (SPPgroupVelocity.real))
LifeTimeRe  = LifeTimeVg(kspp[1:], (SPPgroupVelocityHohenau))
LifeTimePhase = LifeTimeVg(kspp[1:], (SPPphaseVelocity[1:].real))
LifeTimeApprox = 2e0*(eps2new.real)**2 /( omegaspp * eps1new.real**2 * eps2new.imag) #TODO: origin of this formula ?

print("Lifetime approx.")
print(LifeTimeApprox)

# Let's clip Lifetime where they are negative (negative group velocity...). 
LifeTimeOld = np.clip(LifeTimeOld, 0, 1)
LifeTimeNew = np.clip(LifeTimeNew, 0, 1)
LifeTimeRe  = np.clip(LifeTimeRe, 0, 1)
LifeTimePhase  = np.clip(LifeTimePhase, 0, 1)
LifeTimeApprox = np.clip(LifeTimeApprox, 0, 1)

# Plotting the SPP lifetime
HohenauLabel="Hohenau formula"#17
RaetherLabel="Raether formula" #18
plt.figure()
plt.xlabel('Wavelength $\lambda$ (m)')
plt.ylabel(r'SPP lifetime $\tau_{\textrm{SPP}}$ (s)')
line2,=plt.loglog(2*pi*c/omegaspp[1:], LifeTimeRe, 'r--', label=r'$\omega$ real, $k$ complex')
line1,=plt.loglog(2*pi*c/omegaspp[1:], LifeTimeOld, 'k-', label=r'$\omega$ complex, $k$ real')
#line3,=plt.plot(1e9*2*pi*c/omegaspp[1:], 1E12*LifeTimeNew, 'r--', label=r'$\tau_{SPP}=L_{SPP}/v_g$, Eq. (Wirtinger), $\omega \in \mathbb{R}$, $\beta \in \mathbb{C}$')
#plt.plot(1e9*2*pi*c/omegaspp[1:], 1E12*LifeTimePhase, 'r--', label=r'Phase, real $\omega$, complex $\beta$')
#plt.plot(1e9*2*pi*c/omegaspp, 1E12*LifeTimeApprox, 'b--', label=r'Phase, approx.')
#plt.title(MaterialFile1+'/'+MaterialFile2+' interface')
# change style of lines
plt.setp(line1, linewidth=3); plt.setp(line2, linewidth=3); #plt.setp(line3, linewidth=3); 
plt.legend(loc="best")
#plt.xticks(np.linspace(0, 3500e-9, 500))
if(query == "Ti"): #used only for changing the visual scale!
  maxLifeTime=0.05
elif (query == "Ag"):
  maxLifeTime=7
elif (query == "Au"):
  maxLifeTime=10
else:
  maxLifeTime=10
#plt.axis([xmin,xmax,0,maxLifeTime])
plt.tight_layout()
plt.savefig(MaterialFile1+MaterialFile2+'Lifetime.eps')
plt.savefig(MaterialFile1+MaterialFile2+'Lifetime.png')
plt.show()


# 
print("Plotting the SPP decay depth and optical penetration depth...")

SPPdecayDepth=DecayDepth(kzSPP(wavelengths, eps1new, eps2new))
SPPdecayDepth2=DecayDepth(kzSPP(wavelengths, eps2new, eps1new))
OpticalPenetrationDepth = 2e0*omega(wavelengths)/c*(eps2new**0.5e0)
OpticalPenetrationDepth = (OpticalPenetrationDepth.imag)**(-1e0)

plt.figure()
plt.xlabel(r'Wavelength $\lambda$ (nm)')
plt.ylabel('Decay depth (m)')
plt.loglog(1E9*2.*pi*c/omegaspp, SPPdecayDepth, 'r-', label='Plasmonic, medium 1')
plt.loglog(1E9*2.*pi*c/omegaspp, SPPdecayDepth2, 'r--', label='Plasmonic, medium 2')
plt.loglog(1E9*2.*pi*c/omegaspp, OpticalPenetrationDepth, 'bx', label='Optical, medium 2')
plt.legend(loc=1)
#plt.axis([xmin,xmax,0,10e0])
plt.grid()
plt.savefig(MaterialFile1+MaterialFile2+'SppDecayDepth.eps')
plt.savefig(MaterialFile1+MaterialFile2+'SppDecayDepth.png')

#fraction = 70./100.
#EffectivePermittivity = MaxwellGarnett2(eps1new, eps2new, 1.-fraction)
#print "Plotting Maxwell-Garnett 2-material mixing."
#plt.figure()
#plt.title('Mixture: '+MaterialFile1+'('+str(int(100.*fraction))+' perc.)'+'/'+MaterialFile2)
#plt.xlabel(r'Wavelength $\lambda$ (nm)')
#plt.plot(1e9*wavelengths, EffectivePermittivity.real, 'b-', label=r'$Re(\varepsilon_{eff})$')
#plt.plot(1e9*wavelengths, EffectivePermittivity.imag, 'r-', label=r'$Im(\varepsilon_{eff})$')
#plt.legend(loc=1)
#plt.xlim((200.,2000.))
#filename=MaterialFile1+'fraction'+str(fraction)+MaterialFile2+'fraction'+str(1.-fraction)+'Garnett'
#plt.savefig(filename+'.eps')
#plt.savefig(filename+'.png')

if(ShowPictures):
  plt.show()

