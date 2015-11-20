#!/usr/bin/env python
#-*- coding: utf-8 -*-

# IMPORT LIBRARIES
from SimpleSPProutines import *
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

# Known quantities
theta = 0 #Single value here, but we can vectorize functions easily later
f = 0.5 #Filling factor
s = 0.4 #Shape factor
wavelength = 1030e-9
epsilon = 12.80259+0.00109j

# Meshes for solution
#kappax = np.arange(0, 4, 0.1)
#kappay = np.arange(0, 4, 0.1)
#kappa = np.array([wavelength * 1, wavelength * 0]); #test values
kappax = 4 ; kappay = 0*kappax; #test values


#G = np.vectorize(G)
#F = np.vectorize(F)
#R = np.vectorize(R)
#gammaz = np.vectorize(gammaz)
#gammat = np.vectorize(gammat)
#tz = np.vectorize(tz)
#tx = np.vectorize(tx)
#ts = np.vectorize(ts)
#hzz = np.vectorize(hzz)
#hzk = np.vectorize(hzk)
#hkz = np.vectorize(hkz)
#hkk = np.vectorize(hkk)
#hss = np.vectorize(hss)
#kappapn = np.vectorize(kappapn)
#kappamn = np.vectorize(kappamn)
#kpDotY = np.vectorize(kpDotY)
#kpDotX = np.vectorize(kpDotX)
#kmDotY = np.vectorize(kmDotY)
#kmDotX = np.vectorize(kmDotX)
#kpn = np.vectorize(kpn)
#kmn = np.vectorize(kmn)
#vsp = np.vectorize(vsp)
#vsm = np.vectorize(vsm)
#vpp = np.vectorize(vpp) 
#vpm = np.vectorize(vpm)
#etas = np.vectorize(etas)
#etap = np.vectorize(etap)

#for kappax in meshkappa:

ftab = np.arange(0, 1, 0.1)
kapparange = np.arange(0.1,4,0.1)
#for wavelength in wavelengths
#for f in ftab:
for kappax in kapparange:
  
  ## Defining simple quantities for Sipe model
  kappa = np.array([kappax, kappay])
  kappai = np.array([-cmath.sin(theta), 0])
  kappap = kappai + kappa; kappam = kappai - kappa

  print "kappax = "+str(kappax)
  etaresult = etas(theta, f, s, epsilon, kappa, kappap, kappam)
  print "eta = "+str(etaresult)



