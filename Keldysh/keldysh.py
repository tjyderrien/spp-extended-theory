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
from libStark         import *
#from libKeldyshUlrich import *

Header="[keldysh] "

print ""
print "** Welcome to SPP-extended-theory suite."
print "** Author(s): T.J.-Y. Derrien"
print ""
print "** Loaded Keldysh module [Keldysh, Sov. J. Exp. Th. Phys. 47, 5 (1964)]..."
print "** Loaded Gruzdev formula [Gruzdev, Optical Engineering 53, 122515 (2014)]"
print "** Loading Stark module [De Giovannini, U.; Hubener, H. & Rubio, A., Nano Letters, 16, 7993-7998 (2016)]"

print "** Info: this file contains examples how to use the Keldysh library. "
print "         It also contains validation cases of the present theory on Si and known references. "

## Test the Keldysh model using silicon band gap given by the LDA functionals. 
def SiliconLDAbandGap(): #{{{
  print "Defining Si material parameters..."

  Egap = 2.56e0*e; #3.4e0*e #LDA band gap of Si: 2.58 eV. #1.12e0*e for indirect band gap; 
  meff=0.2226e0; #Effective mass of Si
  Ntotal=1.*5E28
  order=200
  
  wavelength = 800e-9;
  tau=7e-15; dt = 1E-17; CEP=0e0
  #PeakFluence = 1.0*1E4 #J/cm2 * 1E4 = J/m2
  #PeakField   = np.sqrt(2e0 * PeakFluence / (tau * c * epsilon_0))
  PeakField = 5E9
  print Header+"Peak field ="+str(PeakField/1E9)+" V/nm"
  
  print "=== SIMPLE QUANTITIES ==="
  
  gamma = gammaKeldysh(Egap, meff, PeakField, wavelength)
  print "Adiabadicity parameter:"+str(gamma)
  
  k1 = Keldysh1phi(gamma); k2 = Keldysh2theta(gamma)
  print "Keldysh1 phi(gamma)   = "+str(k1)
  print "Keldysh2 theta(gamma) = "+str(k2)
  
  EgapEff = EffectiveGap(Egap, k1, k2)
  print "Effective gap: "+str(EgapEff/e)+" eV."

  KeldyshFunctionResult = KeldyshFunction( k1, k2, EgapEff, order, wavelength )
  print "KeldyshFunction_Keldysh: "+str(KeldyshFunctionResult)
  
  KeldyshFunctionResultG = KeldyshFunction_Gruzdev( k1, k2, EgapEff, order, wavelength )
  print "KeldyshFunction_Gruzdev: "+str(KeldyshFunctionResult)
  
  KeldyshFunctionResultGulley = KeldyshFunction_Gulley(k1, k2, EgapEff, order, wavelength)
  print "KeldyshFunction_Gulley: "+str(KeldyshFunctionResultGulley)
  
  wPI = IonizationRate(k1, k2, KeldyshFunctionResult, EgapEff, wavelength, meff)
  
  wPIgulley = IonizationRate(k1, k2, KeldyshFunctionResultGulley, EgapEff, wavelength, meff)
  
  wPIg = IonizationRate_Gruzdev(k1, k2, KeldyshFunctionResultG, EgapEff, wavelength, meff)
  print "wPI(Keldysh)="+str(wPI)
  print "wPI(Gruzdev)="+str(wPIg)
  print "wPI(Gulley) ="+str(wPIgulley)

  print "===== COMPARING THE EFFECTIVE GAPS using Stark effect ===="
  
  print "== Preparing Giovannini et al model... =="
  DME = 1. #from the paper #-1+2j #arbitrary!
  Efield_SI_log = np.linspace(8,11, 50)
  Efield_SI     = np.power(10.,Efield_SI_log)

  E_gap_SI      = Egap
  omega_SI      = 2.*np.pi * c / wavelength

  print "Intuitive Floquet energy shift (eV): "+str(omega_SI*hbar/e)

  Efield_AU = Field_SI_to_AU(Efield_SI)
  E_gap_AU  = Energy_eV_to_Hartree(E_gap_SI/e)
  omega_AU  = Energy_eV_to_Hartree(omega_SI*hbar/e)

  
  print "== 4x4 original numerical attempt"
  Enumerical_min = Stark4bandsEnergyShift_notcorrected_numerical(0., omega_AU, E_gap_AU, DME, 50)
  
  #Enumerical_max = Stark6bandsEnergyShift_notcorrected_numerical(np.max(Efield_AU), omega_AU, E_gap_AU, DME, 50)
  
  # Removing numerical degeneracies
  Enumerical_min = np.round(Enumerical_min, 8)
  Enumerical_min_set = set(Enumerical_min.flatten())
  #Enumerical_max = np.round(Enumerical_max, 8)
  #Enumerical_max_set = set(Enumerical_max.flatten())
  
  Enumerical_min_eV = np.round(Energy_Hartree_to_eV(Enumerical_min), 8)
  Enumerical_min_set_eV = set(Enumerical_min_eV.flatten())
  #Enumerical_max_eV = np.round(Energy_Hartree_to_eV(Enumerical_max), 8)
  #Enumerical_max_set_eV = set(Enumerical_max_eV.flatten())
  
  #Enumerical_min_set = [set(v) for v in Enumerical_min]
  #Enumerical_max_set = map(np.unique, Enumerical_max)
  
  # For comparison, we compute the 4x4 original. 
  ENC1, ENC2, ENC3, ENC4, EgapShift_AU = Stark4bandsEnergyShift_notcorrected_exact(0., omega_AU, E_gap_AU, DME)
  
  print "4x4 original - exact values"
  print "eV: "
  print Energy_Hartree_to_eV(ENC1)
  print Energy_Hartree_to_eV(ENC2)
  print Energy_Hartree_to_eV(ENC3)
  print Energy_Hartree_to_eV(ENC4)
  print ""
  print "4x4 original - numerical attempt"
  print "eV:"
  #print Enumerical_min_set
  print np.array(Enumerical_min_set_eV)

  ## NOTE: Numerical solver works well for 4x4 original case. 
  
  # Now we have corrected the 6x6 and found exact solution
  E1, E2, E3, E4, E5, E6 = Stark6bandsEnergyShift_modified_exact(0., omega_AU, E_gap_AU, DME)
  
## Generalizing to many fields

  print "== Preparing Keldysh model of Stark effect... =="
  gamma_t   = gammaKeldysh(E_gap_SI, 1.0, Efield_SI, wavelength)
  k1_t      = Keldysh1phi(gamma_t); k2_t = Keldysh2theta(gamma_t)
  EgapEff_t = EffectiveGap(E_gap_SI, k1_t, k2_t)
  
  plt.figure()
  plt.semilogx(Efield_SI, EgapEff_t/e, 'r-', label=r'$E_g^{eff}$, Keldysh-Stark (1964)')
  #plt.semilogx(Efield_SI, Energy_Hartree_to_eV(EgapShift_AU), 'b-', label=r'Floquet $E_g$ (2016)')
  plt.semilogx(Efield_SI, Energy_Hartree_to_eV(E1-E2), 'b-', label=r"$E_g + \Delta E_{Stark}=E_1-E_2$")
  plt.semilogx(Efield_SI, Energy_Hartree_to_eV(E1), 'k--', label="Floquet bands")
  plt.semilogx(Efield_SI, Energy_Hartree_to_eV(E2), 'k--')
  plt.semilogx(Efield_SI, Energy_Hartree_to_eV(E3), 'k--')
  plt.semilogx(Efield_SI, Energy_Hartree_to_eV(E4), 'k--')
  plt.semilogx(Efield_SI, Energy_Hartree_to_eV(E5), 'k--')
  plt.semilogx(Efield_SI, Energy_Hartree_to_eV(E6), 'k--')
  #plt.semilogx(Efield_SI, Energy_Hartree_to_eV(ENC1), 'b--', label="Floquet bands")
  #plt.semilogx(Efield_SI, Energy_Hartree_to_eV(ENC2), 'b--')
  #plt.semilogx(Efield_SI, Energy_Hartree_to_eV(ENC3), 'b--')
  #plt.semilogx(Efield_SI, Energy_Hartree_to_eV(ENC4), 'b--')
  plt.xlabel("Field amplitude (V/m)")
  plt.ylabel("Band energy level (eV)")
  plt.legend(loc='best')
  plt.tight_layout()
  plt.show()

  print "==== EXTRAPOLATION TO TEMPORAL ASPECTS ====="
  t0=0. #defines the instant 0.
  Delay = 0e-15 #delay between maxima of the pulses
  tmin=-2.*tau + t0; tmax=2.*tau + Delay + t0

  instants = np.arange(tmin, tmax, dt)
  #print "Time range: "+str(instants.min())+", "+str(instants.max())+"."

  PeakField2  = 0. #PeakField
  CEP2        = 0.*pi
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

  order = 100
  ShowPlot = True

  print Header+"** Test 0: Convergence test using the Keldysh-Gruzdev formulas..."
  print Header+"===== Checking order convergence... ======"
  #plotPulseToDensity(Egap, meff, wavelength, tau, FieldEnvelopeTot.real, dt, 10, ShowPlot)
  #plotPulseToDensity(Egap, meff, wavelength, tau, FieldEnvelopeTot.real, dt, 20, ShowPlot)
  #plotPulseToDensity(Egap, meff, wavelength, tau, FieldEnvelopeTot.real, dt, 30, ShowPlot)
  #plotPulseToDensity(Egap, meff, wavelength, tau, FieldEnvelopeTot.real, dt, 40, ShowPlot)
  #plotPulseToDensity(Egap, meff, wavelength, tau, FieldEnvelopeTot.real, dt, 50, ShowPlot)
  #plotPulseToDensity(Egap, meff, wavelength, tau, FieldEnvelopeTot.real, dt, 100, ShowPlot)
  #plotPulseToDensity(Egap, meff, wavelength, tau, FieldEnvelopeTot.real, dt, 150, ShowPlot)
  #plotPulseToDensity(Egap, meff, wavelength, tau, FieldEnvelopeTot.real, dt, 200, ShowPlot)
  print "======= Checking dt convergence... ======"
  plotPulseToDensity(Egap, meff, wavelength, tau, FieldEnvelopeTot.real, dt/2., order, ShowPlot)
  plotPulseToDensity(Egap, meff, wavelength, tau, FieldEnvelopeTot.real, dt/5., order, ShowPlot)
  plotPulseToDensity(Egap, meff, wavelength, tau, FieldEnvelopeTot.real, dt/10., order, ShowPlot)
  print ""

  print Header+"** Test 2: computing the W_PI values from self-coded and Gruzdev theory..."
  timeKeldysh, N_Keldysh_SI, N_Gruzdev_SI = plotPulseToDensity(Egap, meff, wavelength, tau, FieldEnvelopeTot.real, dt, order, ShowPlot, 0e0, Ntotal)

  print Header+"** Test 3: computing the W_PI values from Vladimir Zhukov tables..."
  wPI_Zhukov = VZ_generateWpiTables(FieldEnvelope1, FieldEnvelope2, wavelength, wavelength2, CEP, CEP2, Egap, meff, tau, tau, Delay, dt, Ntotal, t0)
#}}}

## Repeats the results obtained in Gulley, Opt. Eng. 51, 121805 (2012). 
def SilicaGulley2012(): #{{{
  print "Defining SiO2 material parameters from [Gulley 2012]..."
  numpoints=1000
  
  Egap = 9e0*e; #band gap of SiO2
  meff = 1e0; #Effective mass of SiO2
  N_total=10.*5E28; #valence band electron density #to avoid limitation

  wavelength = 800e-9; wavelength2 = 800e-9
  tau=140e-15; dt = 1E-17; CEP=0e0
  PeakIntensity_log = np.linspace(np.log10(1e11*1e4), np.log10(1e15*1e4), numpoints)
  PeakIntensity = np.power(10., PeakIntensity_log)
  PeakField = np.sqrt(2. * PeakIntensity / c / epsilon_0) #I = 0.5 c n0 eps0 E²
  print Header+"Peak field (min, max): "+str(np.min(PeakField)/1E9)+" V/nm, "+str(np.max(PeakField/1E9))+" V/nm."
  
  
  t0=0. #defines the instant 0.
  Delay = 0. #delay between maxima of the pulses
  tmin=-4.*tau + t0; tmax=4.*tau + Delay + t0

  instants = np.arange(tmin, tmax, dt)
  #print "Time range: "+str(instants.min())+", "+str(instants.max())+"."

  PeakField2  = 0. #/ 2.
  CEP2        = 0. #pi/3.
  #wavelength2 = wavelength
  
  #print Header+"** Test: building single pulse centered on 0..."
  #FieldEnvelope1, RealField1 = PulseSquaredSinTemporalShape(instants, tau, PeakField, wavelength, CEP, t0, 0.)

  #print Header+"** Test: We build a second pulse with a delay..."
  #FieldEnvelope2, RealField2 = PulseSquaredSinTemporalShape(instants, tau, PeakField2, wavelength2, CEP, t0, Delay)

  #print Header+"** Test: building a bicolor double pulse"

  #FieldEnvelopeTot, RealFieldTot = PulseSquaredSinTemporalShapeDoublePulse(instants, tau, tau, PeakField, PeakField, wavelength, wavelength2, CEP, CEP2, t0, Delay)

  #plt.plot(instants, RealField1.real, '-')
  #plt.plot(instants, FieldEnvelope1.real, '--')
  #plt.plot(instants, RealField2.real, '-')
  #plt.plot(instants, FieldEnvelope2.real, '--')
  #plt.plot(instants, RealFieldTot.real, '-')
  #plt.plot(instants, FieldEnvelopeTot.real, '--')
  #plt.xlabel('')
  #plt.savefig('PulseEnvelopes.eps')
  #plt.savefig('PulseEnvelopes.png')
  ##plt.show()

  #print Header+"** Info: PulseEnvelope.EPS and PNG were written in the current folder. "

  order = 100
  ShowPlot = True

  timeKeldysh, N_Keldysh_SI, N_Gruzdev_SI, gamma, wPI, wPIg, N_excited_Keldysh_trapz, N_excited_Gruzdev_trapz = generateWpiTables(Egap, meff, wavelength, tau, PeakField, dt, order, N_total, t0, False)
  
  ### plot w_PI(intensity)
  print Header+"Importing Gruzdev [2014] data..."
  #try:
  Gulley2012=np.loadtxt("Results/Gulley/Gulley-Fig2.csv", dtype='float', delimiter='\t')
      #Gruzdev2014=np.loadtxt("Gruzdev2014-Fig1.csv", dtype='float', delimiter=',')
      #print Gruzdev2014[:,0]
  #except: 
      #print Header+"** Warning: failed to import Gruzdev2014 data table..."
    
    
  print Header+"Plotting as function of laser field intensity ..."
  plt.figure()
  plt.xlabel("Intensity (W/cm$^{2}$)")
  plt.ylabel("$w_{PI}$ (cm$^{-3}$ fs$^{-1}$)")
  #plt.loglog(1e-4*FieldToIntensity(PeakField.real), 1e-6*1e-15*wPI,  linestyle="-", color="r", label=r"$w_{PI}$ "+ShortRefKeldysh)
  plt.loglog(1e-4*FieldToIntensity(PeakField.real), 1e-6*wPI, linestyle="-", color="b", label=r"$w_{PI}$ "+ShortRefGulley)
  plt.loglog(1e-4*Gulley2012[:,0], 1e-6*Gulley2012[:,1], linestyle="-", color="k", label="Data from "+ShortRefGulley) #JUST FOR VALIDATION. 
  plt.grid()
  plt.legend(loc='best')
  #plt.xlim((1E10, 1E14))
  #plt.ylim((1E20*1E6,1E40*1E6))
  plt.tight_layout()
  plt.savefig("Keldysh-Intensity-Wpi-Gulley2012.eps")
  plt.show()
  
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
    

def SilicaGruzdev2014(): #{{{
  print "Defining SiO2 material parameters..."

  Egap = 8.97*e; 
  meff=0.6e0; #Effective mass of Si
  N_total=1.*5E28
  numpoints = 1000
  optical_index = 1.5356
  
  wavelength = 800e-9;
  tau=35e-15; dt = 1E-17; CEP=0e0
  #PeakFluence = 1.0*1E4 #J/cm2 * 1E4 = J/m2
  #PeakField   = np.sqrt(2e0 * PeakFluence / (tau * c * epsilon_0))
  PeakIntensity_log = np.linspace(np.log10(1e14), np.log10(1e18), numpoints)
  PeakIntensity = np.power(10., PeakIntensity_log)
  PeakField = np.sqrt(2. * PeakIntensity / c / epsilon_0) #I = 0.5 c eps0 E²
  # I_inside_matter = 0.5 c eps0 n0 E**2: in matter, pulse is compressed in space. Therefore, intensity is stronger. As photon energy do not change, the time frequency does not change either.  
  LocalPeakField = np.sqrt(2. * PeakIntensity / optical_index/ c / epsilon_0) #I = 0.5 c n0 eps0 E²
  # NOTE: it looks the field is not affected by the refractive index, but the intensity is. This may originate from the change of cycle durations inside matter. Although the field may not change (without accounting for the induced fields). 
  
  print Header+"Peak field (min, max): "+str(np.min(PeakField)/1E9)+" V/nm, "+str(np.max(PeakField/1E9))+" V/nm."

  t0=0. #defines the instant 0.
  
  #print Header+"** Info: PulseEnvelope.EPS and PNG were written in the current folder. "

  order = 100
  ShowPlot = True

  print Header+"** Computing the W_PI values from Keldysh and Gruzdev theories..."
  
  timeKeldysh, N_Keldysh_SI, N_Gruzdev_SI, gamma, wPI, wPIg, N_excited_Keldysh_trapz, N_excited_Gruzdev_trapz = generateWpiTables(Egap, meff, wavelength, tau, LocalPeakField, dt, order, N_total, t0, False, optical_index)
  
  #timeKeldysh, N_Keldysh_SI, N_Gruzdev_SI, gamma, wPI, LocalwPIg, N_excited_Keldysh_trapz, N_excited_Gruzdev_trapz = generateWpiTables(Egap, meff, wavelength, tau, LocalPeakField, dt, order, N_total, t0, False)
  
  ### plot w_PI(intensity)
  print Header+"Importing Gruzdev [2014] data..."
  try:
      #Gruzdev2014=np.loadtxt("Gulley-Fig2.csv", dtype='float', delimiter='\t')
      Gruzdev2014=np.loadtxt("Gruzdev2014-Fig1.csv", dtype='float', delimiter=',')
      #print Gruzdev2014[:,0]
  except: 
      print Header+"** Warning: failed to import Gruzdev2014 data table..."
    
    
  print Header+"Plotting as function of laser field intensity ..."
  plt.figure()
  plt.xlabel("Intensity (W/cm$^{2}$)")
  plt.ylabel("$w_{PI}$ (cm$^{-3}$ fs$^{-1}$)")
  #plt.loglog(1e-4*FieldToIntensity(PeakField.real), 1e-6*1e-15*wPI,  linestyle="-", color="r", label=r"$w_{PI}$ "+ShortRefKeldysh)
  plt.loglog(1e-4*FieldToIntensity(PeakField.real), 1e-6*1e-15*wPIg, linestyle="-", color="b", label=r"$w_{PI}$ "+ShortRefGruzdev)
  #plt.loglog(1e-4*FieldToIntensity(PeakField.real), 1e-6*1e-15*wPIg, linestyle="-", color="g", label=r"$w_{PI}$ "+ShortRefGruzdev)
  plt.loglog(Gruzdev2014[:,0], Gruzdev2014[:,1], linestyle="--", color="k", label="Data from "+ShortRefGruzdev) #JUST FOR VALIDATION. 
  plt.grid()
  plt.legend(loc='best')
  plt.xlim((1E10, 1E14))
  #plt.ylim((1E20*1E6,1E40*1E6))
  plt.tight_layout()
  plt.savefig("Keldysh-Intensity-Wpi-Gruzdev2014.eps")
  plt.show()
  
  #timeKeldysh, N_Keldysh_SI, N_Gruzdev_SI = plotPulseToDensity(Egap, meff, wavelength, tau, PeakField, dt, order, ShowPlot, 0e0, N_total)
  
  print("")

#}}}

#SiliconTunneling()


SiliconLDAbandGap()
#SilicaGulley2012()
#SilicaGruzdev2014()
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
