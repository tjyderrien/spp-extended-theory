#!/usr/bin/env python
#-*- coding: utf-8 -*-

#TODO: To verify coupling efficiency factor from Sipe theory: plot the efficiency factor maximum as function of the laser wavelength, and correlate with papers such as Endriz and Spicer, PRB 4, 4144 (1971); Benneth and Porteus, JOSA 51, 123 (1961)

# IMPORT LIBRARIES
from libSPP import *
#from plotGraph import *

def G(s):
	return 0.5*(cmath.sqrt(s**2+4e0)+s)-cmath.sqrt(s**2+1e0)

def F(s):
	return cmath.sqrt(s**2+1)-s

def R(eps):
	return (eps-1)/(eps+1)



def gammaz(epsilon, f, s):
	return 0.25*(epsilon-1e0)/pi/(epsilon-(1-f)*(epsilon-1)*(F(s)+R(epsilon)*G(s)))

def gammat(epsilon, f, s):
	return 0.25*(epsilon-1e0)/pi/(1e0+0.5e0*(1e0-f)*(epsilon-1)*(F(s)+R(epsilon)*G(s)))

def tz(epsilon, theta):
	return 2e0*cmath.sin(theta)/(epsilon*abs(cmath.cos(theta))+(epsilon-cmath.sin(theta)**2)**(0.5))

def tx(epsilon, theta):
	return 2e0*(epsilon-cmath.sin(theta)**2)**(0.5e0)/(epsilon*abs(cmath.cos(theta))+(epsilon-cmath.sin(theta)**2)**(0.5))

def ts(epsilon, theta):
	return 2*abs(cmath.cos(theta))/(abs(cmath.cos(theta))+(epsilon-cmath.sin(theta)**2)**(0.5))

def hzz(epsilon, kappa):
	return (2*1j)*kappa**2/(epsilon*cmath.sqrt(1-kappa**2)+cmath.sqrt(epsilon-kappa**2))

def hzk(epsilon, kappa):
	return (2*1j)*kappa*cmath.sqrt(1-kappa**2)/(epsilon*cmath.sqrt(1-kappa**2)+cmath.sqrt(epsilon-kappa**2))

def hkz(epsilon, kappa):
	return (2*1j)*kappa*cmath.sqrt(epsilon-kappa**2)/(epsilon*cmath.sqrt(1-kappa**2)+cmath.sqrt(epsilon-kappa**2))

def hkk(epsilon, kappa):
	return (2*1j)*cmath.sqrt((epsilon-kappa**2)*(1-kappa**2))/(epsilon*cmath.sqrt(1-kappa**2)+cmath.sqrt(epsilon-kappa**2))

def hss(epsilon, kappa):
	return (2*1j)/(cmath.sqrt(1-kappa**2)+cmath.sqrt(epsilon-kappa**2))

def kappapn(kappa):
	return cmath.sqrt(kappa[0]**2+(cmath.sin(theta)+kappa[1])**2)

def kappamn(kappa):
	return cmath.sqrt(kappa[0]**2+(cmath.sin(theta)-kappa[1])**2)

def kpDotY(theta, kappa):
	return (cmath.sin(theta)+kappa[1])/(kappapn(kappa))

def kmDotY(theta, kappa):
	return (cmath.sin(theta)-kappa[1])/(kappamn(kappa))

def kpDotX(kappa):
	return kappa[0]/kappapn(kappa)

def kmDotX(kappa):
	return kappa[0]/kappamn(kappa)

def kpn(kappap): 
	return cmath.sqrt(kappap[0]**2+kappap[1]**2)

def kmn(kappam):
	return cmath.sqrt(kappam[0]**2+kappam[1]**2)

def vsp(theta, f, s, epsilon, kappa, kappap):
	return (hss(epsilon, kpn(kappap)) * kpDotY(theta, kappa)**2 + hkk(epsilon, kpn(kappap))*kpDotX(kappa)**2) *gammat(epsilon, f, s)*abs(ts(epsilon, theta))**2

def vsm(theta, f, s, epsilon, kappa, kappam):
	return (hss(epsilon, kmn(kappam)) * kmDotY(theta, kappa)**2 + hkk(epsilon, kmn(kappam))*kmDotX(kappa)**2) *gammat(epsilon, f, s)*abs(ts(epsilon, theta))**2

def vpp(theta, f, s, epsilon, kappa, kappap):
	return (hss(epsilon, kpn(kappap))*kpDotX(kappa)**2+hkk(epsilon, kpn(kappap))*kpDotY(theta, kappa)**2)*gammat(epsilon, f, s)*abs(tx(epsilon, theta))**2+hkz(epsilon, kpn(kappa))*kpDotY(theta, kappa)*gammaz(epsilon, f, s)*epsilon*(tx(epsilon, theta).conjugate())*tz(epsilon, theta)+hzk(epsilon, kpn(kappap))*kpDotY(theta, kappa)*gammat(epsilon, f, s)*tx(epsilon, theta)*(tz(epsilon, theta).conjugate())+hzz(epsilon, kpn(kappap))*gammaz(epsilon, f, s)*epsilon*abs(tz(epsilon, theta))**2

def vpm(theta, f, s, epsilon, kappa, kappam):
	return (hss(epsilon, kmn(kappam))*kmDotX(kappa)**2+hkk(epsilon, kmn(kappam))*kmDotY(theta, kappa)**2)*gammat(epsilon, f, s)*abs(tx(epsilon, theta))**2+hkz(epsilon, kmn(kappam))*kmDotY(theta, kappa)*gammaz(epsilon, f, s)*epsilon*(tx(epsilon, theta).conjugate())*tz(epsilon, theta)+hzk(epsilon, kmn(kappam))*kmDotY(theta, kappa)*gammat(epsilon, f, s)*tx(epsilon, theta)*(tz(epsilon, theta).conjugate())+hzz(epsilon, kmn(kappam))*gammaz(epsilon, f, s)*epsilon*abs(tz(epsilon, theta))**2

def etas(theta, f, s, epsilon, kappa, kappap, kappam):
	return 2*pi*abs(vsp(theta, f, s, epsilon, kappa, kappap)+(vsm(theta, f, s, epsilon, kappa, kappam).conjugate()))

def etap(theta, f, s, epsilon, kappa, kappap, kappam):
	return 2*pi*abs(vpp(theta, f, s, epsilon, kappa, kappap)+(vpm(theta, f, s, epsilon, kappa, kappam).conjugate()))


#==== Attempting a 1D plot
# Known quantities
theta = 0e0 #Single value here, but we can vectorize functions easily later
f = 0.1e0 #Filling factor: taken from Bonse et al, JAP (2009)
s = 0.4e0 #Shape factor: taken from Bonse et al, JAP (2009)
#wavelength = 800e-9
#epsilon = 12.80259+0.00109j
#epsilon = -97.593456+25.2698472743j

# Meshes for solution
#kappax = np.arange(0, 4, 0.1)
#kappay = np.arange(0, 4, 0.1)
#kappa = np.array([wavelength * 1, wavelength * 0]); #test values
kappax = 4e0 ; kappay = 0e0*kappax; #test values

#for kappax in meshkappa:

#ftab = np.arange(0, 1, 0.1)
print "Info: Generating the mesh..."
kapparange = np.arange(0.1,4,0.1)
#for wavelength in wavelengths
#for f in ftab:
#idtab = 0
#etaresult = np.zeros(kapparange.shape)
#kappax = 0e0
#for kappay in kapparange:
  
  ### Defining simple quantities for Sipe model
  #kappa = np.array([kappax, kappay])
  #kappai = np.array([-cmath.sin(theta), 0])
  #kappap = kappai + kappa; kappam = kappai - kappa

  ##print "kappax = "+str(kappax)
  ##print idtab
  #etaresult[idtab] = etas(theta, f, s, epsilon, kappa, kappap, kappam)
  ##print etaresult[idtab]
  ##print "eta = "+str(etaresult)
  #idtab = idtab+1

##=========== Make a 1D plot

#plt.figure()
#plt.xlabel(r'$\kappa_x$')
#plt.ylabel(r'$\eta$')
#print kapparange.shape, etaresult.shape
#plt.plot(kapparange, etaresult, '-', label='Sipe')
#print etaresult
#plt.savefig('SipeEtaKappaX.eps')
##exit()

#=========== 2D plot

query = 'Air'
#query2= 'InP (Bonse 2005)'
#query2= 'Mo (Ordal 1988)'
query2= 'Cu (Palik)'
print "Caution: the expression must be exactly the one of MaterialDatabase.csv."

wavelength = 355
select = str(wavelength)
unit = 1E-9
wavelength = wavelength * unit
print "Wavelength = "+str(wavelength/unit)+" nm."

## Generate the database
SPPdb = GenerateDatabase() #Generate from MaterialDatabase.csv
print "SPP database has "+str(len(SPPdb))+" entries."

print "Full Database:"
print SPPdb

# Select the material of interface 1
SPPdb = FilterDatabase(SPPdb, query, 0)
print "Filter on materials: SPP database has now "+str(len(SPPdb))+" entries."

#print SPPdb #works well

# Filter database on wavelength
try: 
	title = select+' nm'
	SPPdb = FilterDatabase(SPPdb, select+".0", 2)
	print "Filter on wavelength: SPP database "+title+" has "+str(len(SPPdb))+" entries."
except:
	print "Exception: no optical data is available for "+query+" at "+title+"."
	exit()
  
#print SPPdb

# Filter database on materials

try: 
	title = query2
	SPPdb = FilterDatabase(SPPdb, query2, 1)
	print "Filter on material: SPP database "+title+" has "+str(len(SPPdb))+" entries."
except:
	print "Exception: no optical data is available for "+query+" at "+title+"."
	exit()

if(len(SPPdb)==0):
  print "SPP database returned 0 matching result."
  exit()
  
# Extract materials from database
Material1, Material2, Wavelength, OldSPPactiveBool, NewSPPactiveBool, SPPperiod, SPPperiodError, SPPdecayDepth1, SPPdecayDepth2, Reflectivity, OpticalPenetration1, OpticalPenetration2, SPPdecayLength, eps1rM, eps1cM, eps2rM, eps2cM, k1imag, k2imag, DeltaLsppValue = ExtractDataDb(SPPdb)

# Calculation of refractive index
eps1rM=np.asfarray(eps1rM)
eps1cM=np.asfarray(eps1cM)
eps2rM=np.asfarray(eps2rM)
eps2cM=np.asfarray(eps2cM)

epsilon1 = np.add(eps1rM,np.multiply(1e0j, eps1cM))
epsilon2 = np.add(eps2rM,np.multiply(1e0j, eps2cM))
  
print "Mesh generation..."
precision = 2.5e-2
kx = np.arange(-4e0,4e0,precision)
ky = np.arange(-4e0,4e0,precision)
kxx, kyy = np.meshgrid(ky, kx)

# calculating Sipe efficiency for many materials
print "Calculating efficiency for all (kx, ky) values at wavelength "+title+"."

print "kxx shape = "+str(kxx.shape)+"."
etaSipe = np.zeros(kxx.shape)

materialIndex = 0
print "Preparing 2D figure for material "+str(Material2[materialIndex])
print epsilon2[materialIndex]

for m in np.arange(0,(kx.size),1):
	#idy=0
	for n in np.arange(0,ky.size,1):
		#print "[Debug]"+str(m)+", "+str(n)
		kappa = np.array([kx[m], ky[n]])
		kappai = np.array([-cmath.sin(theta), 0])
		kappap = kappai + kappa; kappam = kappai - kappa
		etaSipe[m,n] = etap(theta, f, s, epsilon2[materialIndex], kappa, kappap, kappam)
		#idy=idy+1
	#idx=idx+1	

print etaSipe
 
numberlevels = 8

maximum = np.amax(etaSipe)
print maximum
# Plot the graph
plt.figure()
levels = np.arange(0,maximum,maximum/numberlevels)
CS = plt.contourf(kxx, kyy, etaSipe, levels=levels, cmap=plt.cm.Blues)
plt.xlabel(r'$\kappa_x$')
plt.ylabel(r'$\kappa_y$')
plt.colorbar(CS)
plt.savefig('Sipe2d'+str(wavelength/unit)+'nm-'+query2+'.eps')
plt.show()


# TODO 
# - Automatize the inverse Fourier transform to check regularity and pattern shape: see formula in my thesis. 
# - Automatic calculation of orientation angle precision
# - Can be great to plot directly precision angle as a function of materials. 