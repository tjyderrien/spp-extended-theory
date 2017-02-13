#!/usr/bin/env python
#-*- coding: utf-8 -*-

# IMPORT LIBRARIES
from libSPP import *
#from plotGraph import *

# PHYSICAL INPUT
#wavelength = 800e-9
#epsAir=1e0
#epsSi0=13.64+0.048j
#meffe=0.18
#nuSi=(1.1e-15)**-1

#example=period(betaSPP(wavelength, epsAir, Drude(wavelength, 1e28, epsSi0, nuSi)))

#print example

print "Generating SPP database..."
SPPdb = GenerateDatabase()

SppOutput = 'SPPactiveInterfaces.dat'
print "Exporting to "+SppOutput+"..."
print 
ExportToTxt(SPPdb, SppOutput)
# Writing table caption
f=open(SppOutput, "a")
f.write("#Material1\tMaterial2\tWavelength\tOldSPPactiveBool\tNewSPPactiveBool\tSPPperiod\tSPPperiodError\tSPPdecayDepth1\tSPPdecayDepth2\tReflectivity\tOpticalPenetration1\tOpticalPenetration2\tSPPdecayLength\teps1.real\teps1.imag\teps2.real\teps2.imag\tSPPdepthImagk1\tSPPdepthImagk2")
f.close()
print "Exported. Please open file "+SppOutput+"."

deltaBetaSPP = np.vectorize(deltaBetaSPP)
deltaPeriodSPP = np.vectorize(deltaPeriodSPP)
deltaLspp = np.vectorize(deltaLspp)

##  plot precision of Lambda over precision of epsilon

#eta = np.arange(1e-3, 1e0, 1e-3)

#wavelength = 1030e-9
#epsAir = 1+0j
#epsAu = -26.154188586+1.8503881331j
#epsTi = -4.2656+27.277j
#epsMo = -11.728+20.297j
#epsAl = -90.19+27.70j

#eps1 = epsAir
#eps2 = epsAl

#DeltaPeriod = deltaPeriodSPP(wavelength, eps1, eps2, eta, eta, eta, eta)
#SPPperiod = period(betaSPP(wavelength, eps1, eps2))
#DeltaPeriodRel = DeltaPeriod / SPPperiod

#print "Wavelength (nm)"
#print "Precision on period: "+str(1e9*DeltaPeriod)+" nm"
#print "Relative precision on period "+str(100*DeltaPeriodRel)+"%"

#fig = plt.figure()
#plt.title(r'Period uncertainty at $\lambda=$'+str(wavelength*1E9)+' nm')
#ax1 = fig.add_subplot(111)
#ax1.semilogx(eta, 1e2*DeltaPeriodRel, '-r', label='abs')
#ax1.set_xlabel(r'$\eta$')
#ax1.set_ylabel('Relative uncertainty (\%)', color='r')
#for tl in ax1.get_yticklabels():
    #tl.set_color('r')
#ax2 = ax1.twinx()
#ax2.semilogx(eta, 1e9*DeltaPeriod, '-b', label='rel')
#ax2.set_ylabel(r'$\delta \Lambda_{SPP}$ (nm)', color='b')
#for tl in ax2.get_yticklabels():
    #tl.set_color('b')
#plt.grid()
#plt.show()
#plt.savefig('PeriodError.eps')

#print 1e9*deltaPeriodSPP(800e-9, 1+0j, -26.154188586+1.8503881331j, 1e-4, 1e-4, 1e-4, 1e-4)
#print deltaBetaSPP(800e-9, 1e0+0e0j, -26.154188586e0+1j*1.8503881331e0, 1e-4, 1e-4, 1e-4, 1e-4)


""" TODO: interface with HTML for publication on the web. 
1. Put results into a NP.array.
2. Use a converter to HTML, CSV and PDF maybe. 
"""

# SPPactiveInterfacesArray contains all data we need, just remove lines starting with #. 

