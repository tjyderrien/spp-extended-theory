#!/usr/bin/python2
from libMaterials import reflectivity, Drude
import matplotlib.pyplot as plt
import numpy as np

#Reflectivity for Sokolowski-Tinten 2000

wavelength=625E-9

LogNexc=np.linspace(19+6, np.log10(1E28), 100)
Nexc=np.power(10, LogNexc)
plt.xlabel(r'Excited electron density (m$^{-3}$)')
plt.ylabel('Surface reflectivity (%)')

R1=reflectivity(1., Drude(wavelength, Nexc, 15.1737532077e0+0.171398142447e0j, 0.5E-15**-1, 0.18), 49E0*np.pi/180E0, "P")
R2=reflectivity(1., Drude(wavelength, Nexc, 15.1737532077e0+0.171398142447e0j, 1.0E-15**-1, 0.18), 49E0*np.pi/180E0, "P")
R3=reflectivity(1., Drude(wavelength, Nexc, 15.1737532077e0+0.171398142447e0j, 1.1E-15**-1, 0.18), 49E0*np.pi/180E0, "P")

plt.semilogx(Nexc, R1, "-", label=r"$\nu^{-1}=0.5$ fs")
plt.semilogx(Nexc, R2, "--", label=r"$\nu^{-1}=1.0$ fs")
plt.semilogx(Nexc, R3, "-.", label=r"$\nu^{-1}=1.1$ fs")

plt.legend(loc="best")
plt.savefig("Reflectivity625nm-Si"+".eps")
plt.show()
