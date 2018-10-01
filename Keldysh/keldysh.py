#!/usr/bin/env python2.7
#-*- coding: utf-8 -*-
## @package Keldysh
## Computes the Keldysh excitation rate of quasi-free electrons
# This module aims to calculate the density of excited electrons as function of laser parameters. 
# Two types of usage are planned :
# * Generating tables to use directly into simulation codes
# * Outputing density in certain conditions. 

from libKeldysh       import * 
from libKeldyshPulses import *
from libKeldyshZhukov import *
#from libKeldyshUlrich import *

Header="[keldysh] "

print ""
print "** Welcome to SPP-extended-theory suite."
print "** Author(s): T.J.-Y. Derrien"
print ""
print "** Loaded Keldysh module [Keldysh, Sov. J. Exp. Th. Phys. 47, 5 (1964)]..."
print "** Loaded Gruzdev formula [Gruzdev, Optical Engineering 53, 122515 (2014)]"

print "** Info: this file contains examples how to use the Keldysh library. "
print "         It also contains validation cases of the present theory on Si and known references. "

## Test the Keldysh model using silicon band gap given by the LDA functionals. 
def SiliconLDAbandGap(): #{{{
  print "Defining Si material parameters..."

  Egap = 2.56e0*e; #LDA band gap of Si: 2.58 eV. #1.12e0*e for indirect band gap; 
  meff=0.2226e0; #Effective mass of Si
  Ntotal=1.*5E28

  wavelength = 800e-9;
  tau=10e-15; dt = 1E-17; CEP=0e0
  PeakFluence = 0.01*1E4 #J/cm2 * 1E4 = J/m2
  PeakField   = np.sqrt(2e0 * PeakFluence / (tau * c * epsilon_0))

  t0=0. #defines the instant 0.
  Delay = 20e-15 #delay between maxima of the pulses
  tmin=-1.*tau + t0; tmax=1.*tau + Delay + t0

  instants = np.arange(tmin, tmax, dt)
  #print "Time range: "+str(instants.min())+", "+str(instants.max())+"."

  PeakField2  = PeakField
  CEP2        = pi/3.
  wavelength2 = wavelength / 2.

  print Header+"** Test: building single pulse centered on 0..."
  FieldEnvelope1, RealField1 = PulseSquaredSinTemporalShape(instants, tau, PeakField, wavelength, CEP, t0, 0.)

  print Header+"** Test: We build a second pulse with a delay..."
  FieldEnvelope2, RealField2 = PulseSquaredSinTemporalShape(instants, tau, PeakField2, wavelength2, CEP2, t0, Delay)

  print Header+"** Test: building a bicolor double pulse"

  FieldEnvelopeTot, RealFieldTot = PulseSquaredSinTemporalShapeDoublePulse(instants, tau, tau, PeakField, PeakField2, wavelength, wavelength2, CEP, CEP2, t0, Delay)

  plt.plot(instants, RealField1.real, 'r-')
  plt.plot(instants, FieldEnvelope1.real, 'r--')
  plt.plot(instants, RealField2.real, 'b-')
  plt.plot(instants, FieldEnvelope2.real, 'b--')
  plt.plot(instants, RealFieldTot.real, 'k-')
  plt.plot(instants, FieldEnvelopeTot.real, 'k--')
  plt.xlabel('')
  filename = 'PulseEnvelopes-'+str(wavelength*1E9)+'nm-'+str(wavelength2*1E9)+'nm-CEP2'+str(CEP/pi * 180.)+'+delay'+str(Delay*1E15)+'fs'
  plt.savefig(filename+'.eps')
  plt.savefig(filename+'.png')
  #plt.show()

  print Header+"** Info: PulseEnvelope.EPS and PNG were written in the current folder. "

  order = 50
  ShowPlot = True

  print Header+"** Test 0: Convergence test using the Keldysh-Gruzdev formulas..."
  plotPulseToDensity(Egap, meff, wavelength, tau, PeakFluence, 1E-17, order, ShowPlot)
  plotPulseToDensity(Egap, meff, wavelength, tau, PeakFluence, 5E-17, order, ShowPlot)
  plotPulseToDensity(Egap, meff, wavelength, tau, PeakFluence, 1E-16, order, ShowPlot)
  print "Checking dt convergence..."
  print ""
  print "Checking order convergence..."
  plotPulseToDensity(Egap, meff, wavelength, tau, PeakFluence, dt, 10, ShowPlot)
  plotPulseToDensity(Egap, meff, wavelength, tau, PeakFluence, dt, 20, ShowPlot)
  plotPulseToDensity(Egap, meff, wavelength, tau, PeakFluence, dt, 30, ShowPlot)
  plotPulseToDensity(Egap, meff, wavelength, tau, PeakFluence, dt, 40, ShowPlot)
  plotPulseToDensity(Egap, meff, wavelength, tau, PeakFluence, dt, 50, ShowPlot)

  print Header+"** Test 2: computing the W_PI values from self-coded and Gruzdev theory..."
  timeKeldysh, N_Keldysh_SI, N_Gruzdev_SI = plotPulseToDensity(Egap, meff, wavelength, tau, FieldEnvelope.real, dt, order, ShowPlot, 0e0, Ntotal)

  print Header+"** Test 3: computing the W_PI values from Vladimir Zhukov tables..."
  wPI_Zhukov = VZ_generateWpiTables(FieldEnvelope1, FieldEnvelope2, wavelength, wavelength2, CEP, CEP2, Egap, meff, tau, tau, Delay, dt, Ntotal, t0)
#}}}

## Repeats the results obtained in Gulley, Opt. Eng. 51, 121805 (2012). 
def SilicaGulley2012(): #{{{
  print "Defining SiO2 material parameters from [Gulley 2012]..."

  Egap = 9e0*e; #band gap of SiO2
  meff = 1e0; #Effective mass of SiO2
  Ntotal=10.*5E28; #valence band electron density #to avoid limitation

  wavelength = 800e-9; wavelength2 = 800e-9
  tau=10e-15; dt = 1E-17; CEP=0e0
  PeakFluence = 1E19*10e-15 #setting by the peak intensity #0.01*1E4 #J/cm2 * 1E4 = J/m2
  PeakField   = np.sqrt(2e0 * PeakFluence / (tau * c * epsilon_0))

  t0=0. #defines the instant 0.
  Delay = 0. #delay between maxima of the pulses
  tmin=-1.*tau + t0; tmax=1.*tau + Delay + t0

  instants = np.arange(tmin, tmax, dt)
  #print "Time range: "+str(instants.min())+", "+str(instants.max())+"."

  PeakField2  = 0. #/ 2.
  CEP2        = 0. #pi/3.
  #wavelength2 = wavelength
  
  print Header+"** Test: building single pulse centered on 0..."
  FieldEnvelope1, RealField1 = PulseSquaredSinTemporalShape(instants, tau, PeakField, wavelength, CEP, t0, 0.)

  print Header+"** Test: We build a second pulse with a delay..."
  FieldEnvelope2, RealField2 = PulseSquaredSinTemporalShape(instants, tau, PeakField2, wavelength2, CEP, t0, Delay)

  print Header+"** Test: building a bicolor double pulse"

  FieldEnvelopeTot, RealFieldTot = PulseSquaredSinTemporalShapeDoublePulse(instants, tau, tau, PeakField, PeakField, wavelength, wavelength2, CEP, CEP2, t0, Delay)

  plt.plot(instants, RealField1.real, '-')
  plt.plot(instants, FieldEnvelope1.real, '--')
  plt.plot(instants, RealField2.real, '-')
  plt.plot(instants, FieldEnvelope2.real, '--')
  plt.plot(instants, RealFieldTot.real, '-')
  plt.plot(instants, FieldEnvelopeTot.real, '--')
  plt.xlabel('')
  plt.savefig('PulseEnvelopes.eps')
  plt.savefig('PulseEnvelopes.png')
  #plt.show()

  print Header+"** Info: PulseEnvelope.EPS and PNG were written in the current folder. "

  order = 50
  ShowPlot = True

  print Header+"** Test 1: computing the W_PI values from self-coded and validated Gruzdev theory..."
  timeKeldysh, N_Keldysh_SI, N_Gruzdev_SI = plotPulseToDensity(Egap, meff, wavelength, tau, FieldEnvelope1.real, dt, order, ShowPlot, 0e0, Ntotal)
  
  #print Header+"** Test 0: Convergence test using the Keldysh-Gruzdev formulas..."
  #plotPulseToDensity(Egap, meff, wavelength, tau, PeakFluence, 1E-17, order, ShowPlot)
  #plotPulseToDensity(Egap, meff, wavelength, tau, PeakFluence, 5E-17, order, ShowPlot)
  #plotPulseToDensity(Egap, meff, wavelength, tau, PeakFluence, 1E-16, order, ShowPlot)
  
  #print "Checking dt convergence..."
  #print ""
  #print "Checking order convergence..."
  
  #plotPulseToDensity(Egap, meff, wavelength, tau, PeakFluence, dt, 10, ShowPlot)
  #plotPulseToDensity(Egap, meff, wavelength, tau, PeakFluence, dt, 20, ShowPlot)
  #plotPulseToDensity(Egap, meff, wavelength, tau, PeakFluence, dt, 30, ShowPlot)
  #plotPulseToDensity(Egap, meff, wavelength, tau, PeakFluence, dt, 40, ShowPlot)
  #plotPulseToDensity(Egap, meff, wavelength, tau, PeakFluence, dt, 50, ShowPlot)

  #print Header+"** Test 1: computing the W_PI values from self-coded and validated Gruzdev theory..."
  #timeKeldysh, N_Keldysh_SI, N_Gruzdev_SI = plotPulseToDensity(Egap, meff, wavelength, tau, FieldEnvelope1.real, dt, order, ShowPlot, 0e0, Ntotal)

  #print Header+"** Test 2: computing the W_PI values from Vladimir Zhukov tables..."
  #print Header+"           WE DONT HAVE THEM FOR THIS BAND GAP. Contact zukov@ict.nsc.ru."
  #wPI_Zhukov = VZ_generateWpiTables(FieldEnvelope1, FieldEnvelope2, wavelength, wavelength2, CEP, CEP2, Egap, meff, tau, tau, Delay, dt, Ntotal, t0)
  
  return 0
#}}}

## Plots the silicon tunneling case
def SiliconTunneling(): #{{{
    EfieldLog = np.linspace(9, 11, 100)
    Efield = np.power(10.,EfieldLog)
    wavelength=800e-9
    meff=0.226
    Egap = 2.56*e                            
    wTunnel = KeldyshTunnelingLimit(Egap, meff, wavelength, Efield)
    gamma = gammaKeldysh(Egap, meff, Efield, wavelength)
    plt.figure()
    ax1 = plt.subplot(111)
    ax1.loglog(Efield, wTunnel, 'r-+', label=r"$w_{tunnel}$")
    ax1.set_xlabel(r"$E_{peak}$ (V/m)")
    ax1.set_ylabel(r"$w_{tunnel}$ (m$^{-1}$)")
    ax12 = ax1.twinx()
    ax12.loglog(Efield, gamma, 'k--', label=r"$\gamma$")
    ax12.set_ylabel(r"$\gamma$")
    plt.tight_layout()
    filename="Si-800nm-SiliconLDAbandGap"
    plt.savefig(filename+".eps")
    plt.savefig(filename+".png")
    plt.show()
    
#}}}                                
    

SiliconTunneling()


#SilicaGulley2012()
#SiliconLDAbandGap()
#SilicaGraef2017()



### Build the famous mapping of N_exc(intensity) from Keldysh theory. 

#StephaneGraf = "SiO2"
#HamedMerdji = "ZnO"

#choice = "SiO2"

#if(choice == StephaneGraf):
  #fluencies = 1E4*np.arange(0.1, 10, 0.5) #array([1E10, 1E11, 1E12, 1E13])*1E4 #W/m2
  #tau = 300e-15
  #intensities = fluencies #warning, it's a trick! 

#elif(choice == HamedMerdji): 
  ## Hamed Merdji group case
  #intensities = np.power(10., 4.+np.arange(10., 13., 0.1)) #array([1E10, 1E11, 1E12, 1E13])*1E4 #W/m2
#else: 
  #print Header+"Please define a new set of laser parameters in libKeldysh.py."

#print intensities 

#count = 0
#Nexc = np.zeros(intensities.size)
#for intensity in intensities:
  #if(choice == HamedMerdji): 
    #Nexc[count] = ZnOMerdji2017(intensity)
  #elif(choice == StephaneGraf):
    #Nexc[count] = SilicaGraef2017(intensity)
  #else: 
    #print Header+"** Error in Keldysh.py when computing N_exc. "
  #count += 1
  
#print Nexc

#plt.figure()
#plt.loglog(1E-4*intensities, 1E-6*Nexc)
#plt.ylabel(r'$N_{exc}^{max}$, $cm^{-3}$')
##plt.title(r"Wavelength $\lambda = $"+str(wavelength*1E6)+r" $\mu$m.")

#if(choice == HamedMerdji):
  #plt.xlabel(r"$I_{max}$, $W/cm^2$")
  #plt.savefig('Keldysh-NexcOfIntensity-Merdji-ZnO-3200nm-100fs.eps')

#elif(choice == StephaneGraf): 
  #plt.xlabel(r"$\phi_0$, $J/cm^2$")
  #plt.savefig('Keldysh-NexcOfIntensity-Graf-SiO2-1025nm-300fs.eps')
  
#else: 
  #print Header+"Error when plotting the final figure."
