#!gnuplot

# single SPP functions files


### Dielectric functions of different materials 
# lambda=800e-9; c=2.99792458e8; ec=1.602176487e-19; me=9.1e-31; epsilon0=8.854187e-12;
c=2.99792458e8; ec=1.60217646e-19; me=9.10938188e-31; epsilon0=8.85418781762e-12;
h=6.63e-34; hbar=h/(2e0*pi)

omega(l)=2e0*pi*c/l

meffSi=0.18e0; 
meffaSi=0.2e0;
meffWater=0.5e0;
meffTiO2=1e0;
meffSiO2=0.49e0; #Temnov PhD thesis
meffSiC=0.35e0; #Garry Harris book average value

muSiC(ne)=1e-4*1.02e8*(ne*1e-6)**(-0.326e0) #GL Harris, HS Henry, A Jackson, Carrier mobilities and concentrations in SiC, input: Ne in cm^-3, output: V/m^2/s^1

nuSi=(1.1e-15)**(-1e0); 
nuaSi=(1.0e-15)**(-1e0);
nuWater=(1.7e-15)**(-1e0);
nuTiO2=(1.0e-15)**(-1e0);
nuSiO2=(0.4e-15)**(-1e0); #Temnov PhD thesis
# nuSiC2(ne)=(58e-15)**(-1e0); #calculated from e- mobility
nuSiC(ne)=ec/(meffSiC*me*muSiC(ne)) #Input: Ne in m^-3

Unit={1e0,0e0}
Imaginary={0e0,1e0}

epsilonAir=1e0
# epsilonSi0={13.64e0,0.048e0} #Si palik 800 nm
epsilonSi0={3.692e0,0.0065}**2 #Si palik 800 nm Bonse 2009
epsilonSi0_1030nm={12.8e0,0.001414418e0} #Si palik 1030 nm my code Flaps2D
# epsilonSi0={5.57e0,0.387e0}**2 #Si palik 400 nm
epsilonaSi0={14.9e0,0.627e0} #a-Si palik
epsilonWater={1.326e0,0e0}**2
epsilonTi={-2.85e0,19.1e0}
epsilonTiO2={2.51975e0,0.003e0}**2 #ordinary mode
epsilonTiO22={2.7912e0,0e0}**2 #extraordinary mode
epsilonSiO2={1.453e0,0e0}**2 #Temnov PhD thesis
epsilonSiC0={2.598e0,0e0}**2 #Palik

sigma1Si=1.1e5 #Bonse 2009 value
sigma2Si=6.8e-11 #Bonse 2009 value

Ncr(lambda,meff,nu,epsilonInf)=me*meff*epsilon0*real(epsilonInf)/ec**2 * (omega(lambda)**2 + nu**2)
# Ncr(lambda,meff,nu,epsilon)=real(epsilon)*epsilon0*me*meff*(1e0+omega(lambda)**2e0*nu**-2e0)/(ec*nu**-1e0)**2e0 #Jorn formula
Nspp(lambda,meff,nu,epsilonInf)=me*meff*epsilon0*real(epsilonInf+1e0)/ec**2 * (omega(lambda)**2 + nu**2)

epsilon(lambda, ne,epsilonInf,meff,nu)=epsilonInf-((ne*ec**2/(me*meff*epsilon0))/omega(lambda)**2) * 1e0/(1e0+{0e0,1e0}*nu/omega(lambda))
absorption(lambda, epsilon)=(2e0*omega(lambda)/c * imag(sqrt(epsilon)))
### SPP function at single interface
beta(lambda,epsilon1,epsilon2)=omega(lambda)/c*sqrt(epsilon1*epsilon2/(epsilon1+epsilon2))
period(beta)=2e0*pi/real(beta)

sakabe(lambda,ne,meff)=lambda*(1e0+(4e0*pi*c**2*1e-7*ne*ec**2/(me*meff) / (omega(lambda))**2e0)*(-1e0))**(-0.5e0)

reflectivity(epsilon1, epsilon2)=abs(((epsilon1**0.5e0-epsilon2**0.5e0)/(epsilon1**0.5e0+epsilon2**0.5e0))**2)

### SPP decay depth
DecayDepth(kzSPP)=2e0*pi/real(kzSPP)
kzSPP(lambda,epsilon1,epsilon2)=sqrt(beta(lambda,epsilon1,epsilon2)**2-epsilon2*(omega(lambda)**2/c**2))

# DecayDepthRatio(lambda,epsilon1,epsilon2)=DecayDepth(kzSPP(lambda, epsilon1, epsilon2))/absorption(lambda, epsilon2)

### SPP EXCITATION CONDITION FUNCTIONS

### ConditionPerfect: complex, complex --> bool
ConditionPerfect(epsilon1, epsilon2)=((real(epsilon1)*real(epsilon2) < 0e0) && (real(epsilon1) + real(epsilon2) < 0e0))

## ConditionPerfect: complex, complex --> bool
ConditionAbs(epsilon1,epsilon2) = (real(epsilon1)*real(epsilon2) + imag(epsilon1)*imag(epsilon2) < 0e0)

#### Analytical model for Ne calculation
R=reflectivity(epsilonAir,epsilonSi0)
sigma1Si(lambda)=absorption(lambda,epsilonSi0)

# Ne(lambda,epsilonSi0,intensity)=sigma1Si/(hbar*omega(lambda))*(1e0-R)*intensity*tau+sigma2Si/(2e0*hbar*omega(lambda))*(1e0-R)**2*intensity**2*tau


Ne(lambda,epsilonSi0,fluence)=fluence*(1e0-R)*(sigma1Si(lambda)+sigma2Si*fluence*(1e0-R)/(2e0*sqrt(2e0*pi)*tau))/(hbar*omega(lambda))

