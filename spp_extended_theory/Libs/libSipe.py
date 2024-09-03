## @package libSipe
# Module libSipe provides the basic functions to describe coupling efficiency factor as function of wavenumber kappa. 
# The system is a semi-infinite medium described by an homogeneous complex-valued dielectric permittivity. 
# The model was strictly taken from [Bonse, J. et al, J. Appl. Phys. 97, 013538 (2005)], which summarizes
# the Sipe model given in [Sipe, J. E. et al. Phys. Rev. B 27, 1141-1154 (1983)].
# NOTE: module was only validated for normal incidence. 

import cmath

from scipy.constants import pi

def G(s): #validated on Bonse et al 2005.
	return 0.5*(cmath.sqrt(s**2+4e0)+s)-cmath.sqrt(s**2+1e0)

def F(s): #validated on Bonse et al 2005.
	return cmath.sqrt(s**2+1)-s

## Surface reflectivity computed with oversimplified Fresnel formula at normal incidence. 
def R(eps): #validated on Bonse et al 2005
	return (eps-1.)/(eps+1.)

def gammaz(epsilon, f, s): #validated on Bonse et al 2005
	return 0.25*(epsilon-1e0)/pi/(epsilon-(1.-f)*(epsilon-1.)*(F(s)+R(epsilon)*G(s)))

def gammat(epsilon, f, s): #validated on Bonse et al 2005
	return 0.25*(epsilon-1e0)/pi/(1e0+0.5e0*(1e0-f)*(epsilon-1.)*(F(s)-R(epsilon)*G(s)))

def tz(epsilon, theta): #validated on Bonse et al 2005 - ERRATUM 2024.
	return 2e0*cmath.sin(theta)*abs(cmath.cos(theta))/(epsilon*abs(cmath.cos(theta))+(epsilon-cmath.sin(theta)**2)**(0.5))

def tx(epsilon, theta): #validated on Bonse et al 2005 - ERRATUM 2024.
	return 2e0*(epsilon-cmath.sin(theta)**2)**(0.5e0) * abs(cmath.cos(theta)) / (epsilon*abs(cmath.cos(theta))+(epsilon-cmath.sin(theta)**2)**(0.5))

def ts(epsilon, theta): #validated on Bonse et al 2005
	return 2e0*abs(cmath.cos(theta))/(abs(cmath.cos(theta))+(epsilon-cmath.sin(theta)**2)**(0.5))

def hzz(epsilon, kappa): #validated on Bonse et al 2005
	return (2.*1j) * kappa**2 / (epsilon*cmath.sqrt(1.-kappa**2)+cmath.sqrt(epsilon-kappa**2))

def hzk(epsilon, kappa): #validated on Bonse et al 2005
	return (2.*1j)*kappa*cmath.sqrt(1.-kappa**2)/(epsilon*cmath.sqrt(1.-kappa**2)+cmath.sqrt(epsilon-kappa**2))

def hkz(epsilon, kappa): #validated on Bonse et al 2005
	return (2.*1j)*kappa*cmath.sqrt(epsilon-kappa**2)/(epsilon*cmath.sqrt(1.-kappa**2)+cmath.sqrt(epsilon-kappa**2))

def hkk(epsilon, kappa): #validated on Bonse et al 2005 - ERRATUM 2024.
	return (2.*1j)*cmath.sqrt(epsilon-kappa**2)*cmath.sqrt(1.-kappa**2)/(epsilon*cmath.sqrt(1.-kappa**2)+cmath.sqrt(epsilon-kappa**2))

def hss(epsilon, kappa): #validated on Bonse et al 2005
	return (2.*1j)/(cmath.sqrt(1.-kappa**2)+cmath.sqrt(epsilon-kappa**2))

def kappapn(kappa,theta=0.): #validated on Bonse et al 2005
	return cmath.sqrt(kappa[0]**2+(cmath.sin(theta)+kappa[1])**2)

def kappamn(kappa,theta=0.): #validated on Bonse et al 2005
	return cmath.sqrt(kappa[0]**2+(cmath.sin(theta)-kappa[1])**2)

def kpDotY(theta, kappa): #validated on Bonse et al 2005
	return (cmath.sin(theta)+kappa[1])/(kappapn(kappa))

def kmDotY(theta, kappa): #validated on Bonse et al 2005
	return (cmath.sin(theta)-kappa[1])/(kappamn(kappa))

def kpDotX(kappa): #validated on Bonse et al 2005
	return kappa[0]/kappapn(kappa)

def kmDotX(kappa): #validated on Bonse et al 2005
	return -kappa[0]/kappamn(kappa)

def kpn(kappap):
	return cmath.sqrt(kappap[0]**2+kappap[1]**2)

def kmn(kappam):
	return cmath.sqrt(kappam[0]**2+kappam[1]**2)

def vsp(theta, f, s, epsilon, kappa, kappap): #validated on Bonse et al 2005
	return (hss(epsilon, kpn(kappap)) * kpDotY(theta, kappa)**2 + hkk(epsilon, kpn(kappap))*kpDotX(kappa)**2) *gammat(epsilon, f, s)*abs(ts(epsilon, theta))**2

def vsm(theta, f, s, epsilon, kappa, kappam): #validated on Bonse et al 2005
	return (hss(epsilon, kmn(kappam)) * kmDotY(theta, kappa)**2 + hkk(epsilon, kmn(kappam))*kmDotX(kappa)**2) *gammat(epsilon, f, s)*abs(ts(epsilon, theta))**2

def vpp(theta, f, s, epsilon, kappa, kappap): #validated on Bonse et al 2005
	return (hss(epsilon, kpn(kappap))*kpDotX(kappa)**2+hkk(epsilon, kpn(kappap))*kpDotY(theta, kappa)**2)*gammat(epsilon, f, s)\
		   *abs(tx(epsilon, theta))**2\
		   +hkz(epsilon, kpn(kappap))*kpDotY(theta, kappa)*gammaz(epsilon, f, s)*epsilon*(tx(epsilon, theta).conjugate())*tz(epsilon, theta)\
		   +hzk(epsilon, kpn(kappap))*kpDotY(theta, kappa)*gammat(epsilon, f, s)*tx(epsilon, theta)*(tz(epsilon, theta).conjugate())\
		   +hzz(epsilon, kpn(kappap))*gammaz(epsilon, f, s)*epsilon*abs(tz(epsilon, theta))**2

def vpm(theta, f, s, epsilon, kappa, kappam): #validated on Bonse et al 2005
	return (hss(epsilon, kmn(kappam))*kmDotX(kappa)**2+hkk(epsilon, kmn(kappam))*kmDotY(theta, kappa)**2)*gammat(epsilon, f, s)*abs(tx(epsilon, theta))**2\
		   +hkz(epsilon, kmn(kappam))*kmDotY(theta, kappa)*gammaz(epsilon, f, s)*epsilon*(tx(epsilon, theta).conjugate())*tz(epsilon, theta)\
		   +hzk(epsilon, kmn(kappam))*kmDotY(theta, kappa)*gammat(epsilon, f, s)*tx(epsilon, theta)*(tz(epsilon, theta).conjugate())\
		   +hzz(epsilon, kmn(kappam))*gammaz(epsilon, f, s)*epsilon*abs(tz(epsilon, theta))**2

def etas(theta, f, s, epsilon, kappa, kappap, kappam): #validated on Bonse et al 2005
	return 2.*pi*abs(vsp(theta, f, s, epsilon, kappa, kappap)+(vsm(theta, f, s, epsilon, kappa, kappam).conjugate()))

def etap(theta, f, s, epsilon, kappa, kappap, kappam): #validated on Bonse et al 2005
	return 2.*pi*abs(vpp(theta, f, s, epsilon, kappa, kappap)+(vpm(theta, f, s, epsilon, kappa, kappam).conjugate()))

## Unfortunately, the way that functions were coded is not possible to vectorize automatically. 
# We have to use loops then :-( 
#G      = np.vectorize(G)
#F      = np.vectorize(F)
#R      = np.vectorize(R)
#gammaz = np.vectorize(gammaz)
#gammat = np.vectorize(gammat)
#tz     = np.vectorize(tz)
#tx     = np.vectorize(tx)
#ts     = np.vectorize(ts)
#hzz    = np.vectorize(hzz)
#hzk    = np.vectorize(hzk)
#hkz    = np.vectorize(hkz)
#hkk    = np.vectorize(hkk)
#hss    = np.vectorize(hss)
#kappapn= np.vectorize(kappapn)
#kappamn= np.vectorize(kappamn)
#kpDotY = np.vectorize(kpDotY)
#kmDotY = np.vectorize(kmDotY)
#kpDotX = np.vectorize(kpDotX)
#kmDotX = np.vectorize(kmDotX)
#kpn    = np.vectorize(kpn)
#kmn    = np.vectorize(kmn)
#vsp    = np.vectorize(vsp)
#vsm    = np.vectorize(vsm)
#vpp    = np.vectorize(vpp)
#vpm    = np.vectorize(vpm)
#etas   = np.vectorize(etas)
#etap   = np.vectorize(etap)
