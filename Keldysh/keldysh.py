#!/usr/bin/env python2.7
#-*- coding: utf-8 -*-
## @package Keldysh
## Computes the Keldysh excitation rate of quasi-free electrons
# This module aims to calculate the density of excited electrons as function of laser parameters. 
# Two types of usage are planned :
# * Generating tables to use directly into simulation codes
# * Outputing density in certain conditions. 

from libKeldysh       import *
from libKeldyshZhukov import *
from libStark         import *
#from libKeldyshUlrich import *

Header="[keldysh] "

if __name__ == "__main__":
    print("")
    print("** Welcome to spp_extended_theory suite.")
    print("** Author(s): T.J.-Y. Derrien")
    print("")
    print("** Loaded Keldysh module [Keldysh, Sov. J. Exp. Th. Phys. 47, 5 (1964)]...")
    print("** Loaded Gruzdev formula [Gruzdev, Optical Engineering 53, 122515 (2014)]")
    print("** Loading Stark module [De Giovannini, U.; Hubener, H. & Rubio, A., Nano Letters, 16, 7993-7998 (2016)]")

    print("** Info: this file contains examples how to use the Keldysh library. ")
    print("         It also contains validation cases of the present theory on Si and known references. ")

## Applies Keldysh model to the solid state parameters passed in argument. 
def Test_Keldysh(Egap, meff, PeakField, wavelength, order): #{{{
  print(Header+"Peak field ="+str(PeakField/1E9)+" V/nm")
  
  print("=== SIMPLE QUANTITIES ===")
  
  gamma = gammaKeldysh(Egap, meff, PeakField, wavelength)
  print("Adiabadicity parameter:"+str(gamma))
  
  k1 = Keldysh1phi(gamma); k2 = Keldysh2theta(gamma)
  print("Keldysh1 (Keldysh|Gruzdev) phi(gamma)   = "+str(k1))
  print("Keldysh2 (Keldysh|Gruzdev) theta(gamma) = "+str(k2))
  
  k1Gulley = Keldysh1phiGulley(gamma); k2Gulley = Keldysh2thetaGulley(k1Gulley)
  print("Keldysh1 (Gulley) phi(gamma)   = "+str(k1Gulley))
  print("Keldysh2 (Gulley) theta(gamma) = "+str(k2Gulley))
  
  EgapEff       = EffectiveGap(Egap, k1, k2)
  print("Effective gap (Keldysh|Gruzdev): "+str(EgapEff/e)+" eV.")
  EgapEffGulley = EffectiveGapGulley(Egap, PeakField, meff, wavelength)
  print("Effective gap (Gulley): "+str(EgapEffGulley/e)+" eV.")
  
  xGulley = GulleyX(Egap, gamma, k2Gulley, wavelength)
  print("Gulley X parameter: "+str(xGulley))
  
  print("Checking Gulley elliptics: ")
  print(Gulley_Compute_Elliptics(k1Gulley, k2Gulley))
  
  KeldyshFunctionResult = KeldyshFunction( k1, k2, EgapEff, order, wavelength )
  print("KeldyshFunction_Keldysh: "+str(KeldyshFunctionResult))
  
  KeldyshFunctionResultG = KeldyshFunction_Gruzdev( k1, k2, EgapEff, order, wavelength )
  print("KeldyshFunction_Gruzdev: "+str(KeldyshFunctionResultG))
  
  KeldyshFunctionResultGulley = KeldyshFunctionGulley(k1Gulley, k2Gulley, xGulley, gamma, order, wavelength)
  print("KeldyshFunction_Gulley: "+str(KeldyshFunctionResultGulley))
  
  wPI = IonizationRate(k1, k2, KeldyshFunctionResult, EgapEff, wavelength, meff)
  
  wPIgulley = IonizationRate_Gulley(k1Gulley, k2Gulley, KeldyshFunctionResultGulley, xGulley, wavelength, meff)
  
  wPIg = IonizationRate_Gruzdev(k1, k2, KeldyshFunctionResultG, EgapEff, wavelength, meff)
  print("wPI(Keldysh)="+str(wPI))
  print("wPI(Gruzdev)="+str(wPIg))
  print("wPI(Gulley) ="+str(wPIgulley))
#}}}

## Test the Keldysh model using silicon band gap given by the LDA functionals. Also plots Stark shift as function of the laser intensity. 
def SiliconLDAbandGap(): #{{{
  print("Defining Si material parameters...")
  unit = 1E-4 #W/m2 to W/cm2.
  Egap = 2.56*e; #SiO2 #2.56e0*e; #3.4e0*e #LDA band gap of Si: 2.58 eV. #1.12e0*e for indirect band gap; 
  meff=0.2226e0; #Effective mass of Si
  Ntotal=1.*5E28
  order=200
  Efield_min      = 0E0 
  Efield_max      = 4E10 #V/m
  num_fields=1200
  
  LogScale        =False 
  ShowKeldyshStark=False #Keldysh-Stark is not applicable in tunneling regime. Therefore, it is better to remove it, as it is misleading. NMB: the E_eff in Keldysh theory is applicable only from multiphotonic case.
  ShowBandGap     =True  #Add dots on figs to indicate band gap and replicates. 
  GSpointSize     = 15
  
  wavelength = 800e-9 #
  tau=7e-15; dt = 1E-17; CEP=0e0
  #PeakFluence = 1.0*1E4 #J/cm2 * 1E4 = J/m2
  #PeakField   = np.sqrt(2e0 * PeakFluence / (tau * c * epsilon_0))
  PeakIntensity = 3.2E14 #Zhukov example at 1.6 um
  RefractiveIndex = OpticalIndex[wavelength] #1: for far field, # 1.45: for near field in silica
  PeakField   = np.sqrt(2e0 * PeakIntensity * RefractiveIndex / (c * epsilon_0))
  #PeakField = 3E9
  
  Test_Keldysh(Egap, meff, PeakField, wavelength, order)

  #exit()
  print("===== COMPARING THE EFFECTIVE GAPS using Stark effect ====")
  exit()

  
  print("== Preparing Giovannini et al model... ==")
  DME = 1. #from the paper #-1+2j #arbitrary!
  if(LogScale): 
      Efield_SI_log = np.linspace(np.log10(Efield_min),np.log10(Efield_max), num_fields)
      Efield_SI     = np.power(10.,Efield_SI_log)
  else:
      Efield_SI     = np.linspace(Efield_min, Efield_max, num_fields)
  #Efield_SI     = 1e9  

  E_gap_SI      = Egap
  omega_SI      = 2.*np.pi * c / wavelength

  print("Replicas energy shifts (eV): "+str(omega_SI*hbar/e))

  Efield_AU = Field_SI_to_AU(Efield_SI)
  E_gap_AU  = Energy_eV_to_Hartree(E_gap_SI/e)
  omega_AU  = Energy_eV_to_Hartree(omega_SI*hbar/e)
  
  #def Test():
        
    
    #print("== USELESS: 2-band 1-photon (4x4) original numerical attempt")
    #Enumerical_min = Stark2bands1photon_EnergyShift_notcorrected_numerical(Efield_AU, omega_AU, E_gap_AU, DME, 50)
    
    ##Enumerical_max = Stark2bands1photon_EnergyShift_notcorrected_numerical(np.max(Efield_AU), omega_AU, E_gap_AU, DME, 50)
    
    ## Removing numerical degeneracies
    #Enumerical_min = np.round(Enumerical_min, 8)
    #Enumerical_min_set = set(Enumerical_min.flatten())
    ##Enumerical_max = np.round(Enumerical_max, 8)
    ##Enumerical_max_set = set(Enumerical_max.flatten())
    
    #Enumerical_min_eV = np.round(Energy_Hartree_to_eV(Enumerical_min), 8)
    #Enumerical_min_set_eV = set(Enumerical_min_eV.flatten())
    ##Enumerical_max_eV = np.round(Energy_Hartree_to_eV(Enumerical_max), 8)
    ##Enumerical_max_set_eV = set(Enumerical_max_eV.flatten())
    
    ##Enumerical_min_set = [set(v) for v in Enumerical_min]
    ##Enumerical_max_set = map(np.unique, Enumerical_max)
    
    ## For comparison, we compute the 4x4 original. 
    #ENC1, ENC2, ENC3, ENC4, EgapShift_AU = Stark2bands1photon_Cropped_EnergyShift_notcorrected_exact(Efield_AU, omega_AU, E_gap_AU, DME)
    
    #print("USELESS: 4x4 original - exact values")
    #print "eV: "
    #print Energy_Hartree_to_eV(ENC1)
    #print Energy_Hartree_to_eV(ENC2)
    #print Energy_Hartree_to_eV(ENC3)
    #print Energy_Hartree_to_eV(ENC4)
    #print ""
    #print "USELESS: 4x4 original - numerical attempt"
    #print "eV:"
    ##print Enumerical_min_set
    #print np.array(Enumerical_min_set_eV)

    ### NOTE: Numerical solver works well for 4x4 original case. 
    
    ## Now we have corrected the 6x6 and found exact solution
    #E1, E2, E3, E4, E5, E6 = Stark2bands1photon_EnergyShift_modified_exact(Efield_AU, omega_AU, E_gap_AU, DME)
    #print ""
    #print "VALID: 2 BANDS - 1 PHOTON - 6x6 modified exact"
    #print "eV: "
    #print Energy_Hartree_to_eV(E1)
    #print Energy_Hartree_to_eV(E2)
    #print Energy_Hartree_to_eV(E3)
    #print Energy_Hartree_to_eV(E4)
    #print Energy_Hartree_to_eV(E5)
    #print Energy_Hartree_to_eV(E6)
    
    ### Let's go for 12x12 matrix, only numerical. 
    #print ""
    #print "4 BANDS - 1 PHOTON - 12x12 - numerical approach"
    #FourBandsOnePhoton_E     = Stark4bands1photon_EnergyShift_notcorrected_numerical(Efield_AU, omega_AU, E_gap_AU, DME, 50)
    #FourBandsOnePhoton_Eigen = Stark4bands1photon_EnergyShift_eigen(Efield_AU, omega_AU, E_gap_AU, DME)
    
    #FourBandsOnePhoton_E_round = np.round(FourBandsOnePhoton_E, 8)
    #FourBandsOnePhoton_E_round_set = set(FourBandsOnePhoton_E_round.flatten())
    #FourBandsOnePhoton_Eigen_round = np.round(FourBandsOnePhoton_Eigen, 8)
    #FourBandsOnePhoton_Eigen_round_set = set(FourBandsOnePhoton_Eigen_round.flatten())
    ##Enumerical_max = np.round(Enumerical_max, 8)
    ##Enumerical_max_set = set(Enumerical_max.flatten())
    
    #FourBandsOnePhoton_E_round_eV = np.round(Energy_Hartree_to_eV(FourBandsOnePhoton_E_round), 5)
    #FourBandsOnePhoton_E_round_set_eV = set(FourBandsOnePhoton_E_round_eV.flatten())
    #FourBandsOnePhoton_Eigen_round_eV = np.round(Energy_Hartree_to_eV(FourBandsOnePhoton_Eigen_round), 5)
    #FourBandsOnePhoton_Eigen_round_set_eV = set(FourBandsOnePhoton_Eigen_round_eV.flatten())
    
    #print "eV:"
    #print FourBandsOnePhoton_E_round_set_eV
    #print "eV: eigensolver: "
    #print FourBandsOnePhoton_Eigen_round_eV
    
    #print ""
    #print "2 BANDS - 2 PHOTONS - 10x10 numerical (eigen solver)"
    #TwoBandsTwoPhotons_Eigen       = Stark2bands2photons_EnergyShift_eigen(Efield_AU, omega_AU, E_gap_AU, DME)
    #TwoBandsTwoPhotons_Eigen_eV    = Energy_Hartree_to_eV(TwoBandsTwoPhotons_Eigen)
    #TwoBandsTwoPhotons_Eigen_round     = np.round(TwoBandsTwoPhotons_Eigen, 8)
    #TwoBandsTwoPhotons_Eigen_round_set = set(TwoBandsTwoPhotons_Eigen_round.flatten())
    #TwoBandsTwoPhotons_Eigen_round_eV = np.round(Energy_Hartree_to_eV(TwoBandsTwoPhotons_Eigen_round), 5)
    #TwoBandsTwoPhotons_Eigen_round_set_eV = set(TwoBandsTwoPhotons_Eigen_round_eV.flatten())
    
    #print "eV: eigensolver: "
    #print Energy_Hartree_to_eV(TwoBandsTwoPhotons_Eigen)
    ##print TwoBandsTwoPhotons_Eigen_round_eV.flatten()
    #print np.shape(TwoBandsTwoPhotons_Eigen_round_eV)
  
  ##exit()

## Generalizing to many fields
  print("== Preparing Stark shift as function of field intensity... ==")
  
  ## Keldysh Stark shift, 2 levels, 1 photon. 
  gamma_t   = gammaKeldysh(E_gap_SI, 1.0, Efield_SI*np.sqrt(RefractiveIndex), wavelength)
  k1_t      = Keldysh1phi(gamma_t); k2_t = Keldysh2theta(gamma_t)
  EgapEff_t = EffectiveGap(E_gap_SI, k1_t, k2_t)
  
  print("Plotting scattered graph...")
  print(np.shape(Efield_AU))
  
  filename=str(round(E_gap_SI/e))+"eV-"+str(wavelength*1E9)+"nm"
  
  # PLOTTING THE 2-bands 1-photon Stark effect
  ax1=plt.figure()
  #plt.xlabel("Field amplitude (V/m)")
  plt.xlabel(r"Laser intensity (W/cm$^2$)")
  plt.ylabel("Band energy level (eV)")
  
  # Preparing the 2-levels 2-photon, and 4-levels 1-photon Stark shifts. 
  for element in Efield_SI:
      print(element)
      E1, E2, E3, E4, E5, E6         = Stark2bands1photon_EnergyShift_modified_exact(Field_SI_to_AU(element*np.sqrt(RefractiveIndex)), omega_AU, E_gap_AU, DME)
      TwoBandsOnePhoton_Eigen        = [E1, E2, E3, E4, E5, E6]
      TwoBandsOnePhoton_Eigen_eV    = Energy_Hartree_to_eV(TwoBandsOnePhoton_Eigen)
      print(TwoBandsOnePhoton_Eigen_eV)
      Intensity_el  = 0.5*c*epsilon_0*element**2*RefractiveIndex
      #plt.scatter(np.ones(np.size(TwoBandsOnePhoton_Eigen_eV))*element, TwoBandsOnePhoton_Eigen_eV, c="black", s=1)
      ## THIS IS THE MOST IMPORTANT PART.
      plt.scatter(unit * np.ones(np.size(TwoBandsOnePhoton_Eigen_eV))*Intensity_el, TwoBandsOnePhoton_Eigen_eV, c="black", s=1)
  Intensity_SI  = unit * 0.5*c*epsilon_0*Efield_SI**2*RefractiveIndex
  if(ShowKeldyshStark):
    plt.plot(Efield_SI,  0.5*EgapEff_t/e, 'r-',    label=r'$E_g^{eff}$, Keldysh-Stark (1964)')
    plt.plot(Efield_SI, -0.5*EgapEff_t/e, 'r-')
    plt.plot(Intensity_SI,  0.5*EgapEff_t/e, 'r-', label=r'$E_g^{eff}$, Keldysh-Stark (1964)')
    plt.plot(Intensity_SI, -0.5*EgapEff_t/e, 'r-')
  if(ShowBandGap):
      MPI_number=1 #here
      Intensity_min  = 0.5*c*epsilon_0*Efield_min**2*RefractiveIndex
      plt.scatter(np.ones(np.size(E_gap_SI/e))*Intensity_min, E_gap_SI/e/2, c="black", s=GSpointSize, label="Gap")
      plt.scatter(np.ones(np.size(E_gap_SI/e))*Intensity_min,-E_gap_SI/e/2, c="black", s=GSpointSize)
      for MPI in np.arange(1,MPI_number+1):
            plt.scatter(np.ones(np.size(E_gap_SI/e))*Intensity_min, np.real(E_gap_SI/e/2.)+MPI*h*c/wavelength/e, c="red",    s=GSpointSize, label="non-interacting replicates")
            plt.scatter(np.ones(np.size(E_gap_SI/e))*Intensity_min, np.real(E_gap_SI/e/2.)-MPI*h*c/wavelength/e, c="blue",   s=GSpointSize, label="non-interacting replicates")
            plt.scatter(np.ones(np.size(-E_gap_SI/e))*Intensity_min, np.real(-E_gap_SI/e/2.)+MPI*h*c/wavelength/e, c="red",  s=GSpointSize, label="non-interacting replicates")
            plt.scatter(np.ones(np.size(-E_gap_SI/e))*Intensity_min, np.real(-E_gap_SI/e/2.)-MPI*h*c/wavelength/e, c="blue", s=GSpointSize, label="non-interacting replicates")
  plt.title("2 bands, 1 photon transition")
  #plt.legend(loc='best')
  plt.tight_layout()
  plt.savefig(filename+"_TwoBandsOnePhoton.eps")
  #plt.show()
  
  
  # PLOTTING THE 4-bands 1-photon Stark effect
  plt.figure()
  #plt.xlabel("Field amplitude (V/m)")
  plt.xlabel(r"Laser intensity (W/cm$^2$)")
  plt.ylabel("Band energy level (eV)")
  
  for element in Efield_SI:
      print(element)
      FourBandsOnePhoton_Eigen       = Stark4bands1photon_EnergyShift_eigen(Field_SI_to_AU(element*np.sqrt(RefractiveIndex)), omega_AU, E_gap_AU, DME)
      FourBandsOnePhoton_Eigen_eV    = Energy_Hartree_to_eV(FourBandsOnePhoton_Eigen)
      print(FourBandsOnePhoton_Eigen_eV)
      Intensity_el  = 0.5*c*epsilon_0*element**2*RefractiveIndex
      #plt.scatter(np.ones(np.size(FourBandsOnePhoton_Eigen_eV))*element, FourBandsOnePhoton_Eigen_eV, c="black", s=1)
      plt.scatter(unit * np.ones(np.size(FourBandsOnePhoton_Eigen_eV))*Intensity_el, FourBandsOnePhoton_Eigen_eV, c="black", s=1)
  Intensity_SI  = unit * 0.5*c*epsilon_0*Efield_SI**2*RefractiveIndex
  #plt.plot(Efield_SI, 0.5*EgapEff_t/e, 'r-', label=r'$E_g^{eff}$, Keldysh-Stark (1964)')
  #plt.plot(Efield_SI, -0.5*EgapEff_t/e, 'r-')
  if(ShowKeldyshStark):
    plt.plot(Intensity_SI, 0.5*EgapEff_t/e, 'r-', label=r'$E_g^{eff}$, Keldysh-Stark (1964)')
    plt.plot(Intensity_SI, -0.5*EgapEff_t/e, 'r-')
  if(ShowBandGap):
      MPI_number=1 #here
      Intensity_min  = 0.5*c*epsilon_0*Efield_min**2*RefractiveIndex
      plt.scatter(np.ones(np.size(E_gap_SI/e))*Intensity_min, E_gap_SI/e/2, c="black", s=10, label="Gap")
      plt.scatter(np.ones(np.size(E_gap_SI/e))*Intensity_min,-E_gap_SI/e/2, c="black", s=10)
      
      for MPI in np.arange(1,MPI_number+1):
            plt.scatter(np.ones(np.size(E_gap_SI/e))*Intensity_min, np.real(E_gap_SI/e/2.)+MPI*h*c/wavelength/e, c="red",    s=GSpointSize, label="non-interacting replicates")
            plt.scatter(np.ones(np.size(E_gap_SI/e))*Intensity_min, np.real(E_gap_SI/e/2.)-MPI*h*c/wavelength/e, c="blue",   s=GSpointSize, label="non-interacting replicates")
            plt.scatter(np.ones(np.size(-E_gap_SI/e))*Intensity_min, np.real(-E_gap_SI/e/2.)+MPI*h*c/wavelength/e, c="red",  s=GSpointSize, label="non-interacting replicates")
            plt.scatter(np.ones(np.size(-E_gap_SI/e))*Intensity_min, np.real(-E_gap_SI/e/2.)-MPI*h*c/wavelength/e, c="blue", s=GSpointSize, label="non-interacting replicates")
  
  plt.title("4 bands (deg. 2), 1 photon transition")
  #plt.legend(loc='best')
  plt.tight_layout()
  plt.savefig(filename+"_FourBandsOnePhoton.eps")
  #plt.show()
  
  # PLOTTING THE 2-bands 2-photon Stark effect
  plt.figure()
  #plt.xlabel("Field amplitude (V/m)")
  plt.xlabel(r"Laser intensity (W/cm$^2$)")
  plt.ylabel("Band energy level (eV)")
  #print Efield_SI
  for element in Efield_SI: 
      print(element)
      TwoBandsTwoPhotons_Eigen       = Stark2bands2photons_EnergyShift_eigen(Field_SI_to_AU(element*np.sqrt(RefractiveIndex)), omega_AU, E_gap_AU, DME)
      TwoBandsTwoPhotons_Eigen_eV    = Energy_Hartree_to_eV(TwoBandsTwoPhotons_Eigen)
      print(TwoBandsTwoPhotons_Eigen_eV)
      Intensity_el  = 0.5*c*epsilon_0*element**2*RefractiveIndex
      #plt.scatter(np.ones(np.size(TwoBandsTwoPhotons_Eigen_eV))*element, TwoBandsTwoPhotons_Eigen_eV, c="black", s=1)
      plt.scatter(unit * np.ones(np.size(TwoBandsTwoPhotons_Eigen_eV))*Intensity_el, TwoBandsTwoPhotons_Eigen_eV, c="black", s=1) 
  #plt.scatter(Field_AU_to_SI(Efield_AU), TwoBandsTwoPhotons_Eigen_round_eV, s=1, c=(1,1,1))
  Intensity_SI  = unit * 0.5*c*epsilon_0*Efield_SI**2*RefractiveIndex
  #plt.plot(Efield_SI, 0.5*EgapEff_t/e, 'r-', label=r'$E_g^{eff}$, Keldysh-Stark (1964)')
  #plt.plot(Efield_SI, -0.5*EgapEff_t/e, 'r-')
  if(ShowKeldyshStark):
      plt.plot(Intensity_SI, 0.5*EgapEff_t/e, 'r-', label=r'$E_g^{eff}$, Keldysh-Stark (1964)')
      plt.plot(Intensity_SI, -0.5*EgapEff_t/e, 'r-')
  if(ShowBandGap):
      MPI_number=2 #here
      Intensity_min  = 0.5*c*epsilon_0*Efield_min**2*RefractiveIndex
      plt.scatter(np.ones(np.size(E_gap_SI/e))*Intensity_min, E_gap_SI/e/2, c="black", s=10, label="Gap")
      plt.scatter(np.ones(np.size(E_gap_SI/e))*Intensity_min,-E_gap_SI/e/2, c="black", s=10)
      for MPI in np.arange(1,MPI_number+1):
          plt.scatter(np.ones(np.size(E_gap_SI/e))*Intensity_min, np.real(E_gap_SI/e/2.)+MPI*h*c/wavelength/e, c="red",    s=GSpointSize, label="non-interacting replicates")
          plt.scatter(np.ones(np.size(E_gap_SI/e))*Intensity_min, np.real(E_gap_SI/e/2.)-MPI*h*c/wavelength/e, c="blue",   s=GSpointSize, label="non-interacting replicates")
          plt.scatter(np.ones(np.size(-E_gap_SI/e))*Intensity_min, np.real(-E_gap_SI/e/2.)+MPI*h*c/wavelength/e, c="red",  s=GSpointSize, label="non-interacting replicates")
          plt.scatter(np.ones(np.size(-E_gap_SI/e))*Intensity_min, np.real(-E_gap_SI/e/2.)-MPI*h*c/wavelength/e, c="blue", s=GSpointSize, label="non-interacting replicates")
  plt.title("2 bands, 2 photons transition")
  #plt.legend(loc='best')
  plt.tight_layout()
  #plt.xscale('log')
  plt.savefig(filename+"_TwoBandsTwoPhotons.eps")
  plt.show()
  
  #plt.semilogx(Efield_SI, Energy_Hartree_to_eV(EgapShift_AU), 'b-', label=r'Floquet $E_g$ (2016)')
  #plt.semilogx(Efield_SI, Energy_Hartree_to_eV(E1-E2), 'b-', label=r"$E_g + \Delta E_{Stark}=E_1-E_2$")
  
  
  #plt.semilogx(Efield_SI, Energy_Hartree_to_eV(E1), 'k--', label="Floquet bands")
  #plt.semilogx(Efield_SI, Energy_Hartree_to_eV(E2), 'k--')
  #plt.semilogx(Efield_SI, Energy_Hartree_to_eV(E3), 'k--')
  #plt.semilogx(Efield_SI, Energy_Hartree_to_eV(E4), 'k--')
  #plt.semilogx(Efield_SI, Energy_Hartree_to_eV(E5), 'k--')
  #plt.semilogx(Efield_SI, Energy_Hartree_to_eV(E6), 'k--')
  
  #plt.semilogx(Efield_SI, Energy_Hartree_to_eV(ENC1), 'b--', label="Floquet bands")
  #plt.semilogx(Efield_SI, Energy_Hartree_to_eV(ENC2), 'b--')
  #plt.semilogx(Efield_SI, Energy_Hartree_to_eV(ENC3), 'b--')
  #plt.semilogx(Efield_SI, Energy_Hartree_to_eV(ENC4), 'b--')
  
  print("##### PREPARING NUMERICAL INTEGRATION OF KELDYSH CONTOUR by VP Zhukov.")
  plt.figure()
  Intensity_SI  = 0.5*c*epsilon_0*Efield_SI**2*RefractiveIndex
  wPI_Zhukov = VZ_generateWpiTables(Efield_SI*np.sqrt(RefractiveIndex), 0., wavelength, 800e-9, 0., 0., Egap, meff, Ntotal)
  ExportToTxt(wPI_Zhukov, "Zhukov_Wpi"+str(wavelength*1E9)+"nm.csv")
  plt.xlabel(r"Laser intensity (W/m$^2$)")
  plt.ylabel(r"Excitation rate $w_{\mathrm{PI}}$ (m$^{-3}s^{-1}$)")
  plt.loglog(unit * Intensity_SI, wPI_Zhukov,'.', label="Num. int. Keldysh")
  plt.legend(loc='best')
  plt.tight_layout()
  plt.grid()
  filename="Zhukov_NumInt-"+str(round(Egap/e))+"eV-"+str(1E9*wavelength)+".eps"
  plt.savefig(filename)
  print(Header+"** Info:"+filename+"was created.")
  plt.show()
  exit()
#}}}

## Repeats results from Gulley2012, but could not be repeated so far. 
def KeldyshGulley():

  #print "==== EXTRAPOLATION TO TEMPORAL ASPECTS ====="
  #t0=0. #defines the instant 0.
  #Delay = 0e-15 #delay between maxima of the pulses
  #tmin=-2.*tau + t0; tmax=2.*tau + Delay + t0

  #instants = np.arange(tmin, tmax, dt)
  ##print "Time range: "+str(instants.min())+", "+str(instants.max())+"."

  #PeakField2  = 0. #PeakField
  #CEP2        = 0.*pi
  #wavelength2 = wavelength / 2.

  #print Header+"** Test: building single pulse centered on 0..."
  #FieldEnvelope1, RealField1 = PulseSquaredSinTemporalShape(instants, tau, PeakField, wavelength, CEP, t0, 0.)

  #print Header+"** Test: We build a second pulse with a delay..."
  #FieldEnvelope2, RealField2 = PulseSquaredSinTemporalShape(instants, tau, PeakField2, wavelength2, CEP2, t0, Delay)

  #print Header+"** Test: building a bicolor double pulse"

  #FieldEnvelopeTot, RealFieldTot = PulseSquaredSinTemporalShapeDoublePulse(instants, tau, tau, PeakField, PeakField2, wavelength, wavelength2, CEP, CEP2, t0, Delay)

  #plt.plot(instants, RealField1.real, 'r-')
  #plt.plot(instants, FieldEnvelope1.real, 'r--')
  #plt.plot(instants, RealField2.real, 'b-')
  #plt.plot(instants, FieldEnvelope2.real, 'b--')
  #plt.plot(instants, RealFieldTot.real, 'k-')
  #plt.plot(instants, FieldEnvelopeTot.real, 'k--')
  #plt.xlabel('')
  #filename = 'PulseEnvelopes-'+str(wavelength*1E9)+'nm-'+str(wavelength2*1E9)+'nm-CEP2'+str(CEP/pi * 180.)+'+delay'+str(Delay*1E15)+'fs'
  #plt.savefig(filename+'.eps')
  #plt.savefig(filename+'.png')
  ##plt.show()

  print(Header+"** Info: PulseEnvelope.EPS and PNG were written in the current folder. ")

  order = 100
  ShowPlot = True

  print(Header+"** Test 0: Convergence test using the Keldysh-Gruzdev formulas...")
  print(Header+"===== Checking order convergence... ======")
  #plotPulseToDensity(Egap, meff, wavelength, tau, FieldEnvelopeTot.real, dt, 10, ShowPlot)
  #plotPulseToDensity(Egap, meff, wavelength, tau, FieldEnvelopeTot.real, dt, 20, ShowPlot)
  #plotPulseToDensity(Egap, meff, wavelength, tau, FieldEnvelopeTot.real, dt, 30, ShowPlot)
  #plotPulseToDensity(Egap, meff, wavelength, tau, FieldEnvelopeTot.real, dt, 40, ShowPlot)
  #plotPulseToDensity(Egap, meff, wavelength, tau, FieldEnvelopeTot.real, dt, 50, ShowPlot)
  #plotPulseToDensity(Egap, meff, wavelength, tau, FieldEnvelopeTot.real, dt, 100, ShowPlot)
  #plotPulseToDensity(Egap, meff, wavelength, tau, FieldEnvelopeTot.real, dt, 150, ShowPlot)
  #plotPulseToDensity(Egap, meff, wavelength, tau, FieldEnvelopeTot.real, dt, 200, ShowPlot)
  print("======= Checking dt convergence... ======")
  plotPulseToDensity(Egap, meff, wavelength, tau, FieldEnvelopeTot.real, dt/2., order, ShowPlot)
  plotPulseToDensity(Egap, meff, wavelength, tau, FieldEnvelopeTot.real, dt/5., order, ShowPlot)
  plotPulseToDensity(Egap, meff, wavelength, tau, FieldEnvelopeTot.real, dt/10., order, ShowPlot)
  print("")

  print(Header+"** Test 2: computing the W_PI values from self-coded and Gruzdev theory...")
  timeKeldysh, N_Keldysh_SI, N_Gruzdev_SI = plotPulseToDensity(Egap, meff, wavelength, tau, FieldEnvelopeTot.real, dt, order, ShowPlot, 0e0, Ntotal)

  print(Header+"** Test 3: computing the W_PI values from Vladimir Zhukov tables...")
  wPI_Zhukov = VZ_generateWpiTables(FieldEnvelope1, FieldEnvelope2, wavelength, wavelength2, CEP, CEP2, Egap, meff)
#}}}

## Repeats the results obtained in Gulley, Opt. Eng. 51, 121805 (2012). 
def SilicaGulley2012(): #{{{
  print("Defining SiO2 material parameters from [Gulley 2012]...")
  numpoints=1000
  
  order = 50 
  
  Egap = 9e0*e; #band gap of SiO2
  meff = 1e0; #Effective mass of SiO2
  N_total=10.*5E28; #valence band electron density #to avoid limitation

  wavelength = 800e-9; wavelength2 = 800e-9
  tau=140e-15; dt = 1E-17; CEP=0e0
  PeakIntensity_log = np.linspace(np.log10(1e11*1e4), np.log10(1e15*1e4), numpoints)
  PeakIntensity = np.power(10., PeakIntensity_log)
  PeakField = np.sqrt(2. * PeakIntensity / c / epsilon_0) #I = 0.5 c n0 eps0 E²
  print(Header+"Peak field (min, max): "+str(np.min(PeakField)/1E9)+" V/nm, "+str(np.max(PeakField/1E9))+" V/nm.")
  
  #PeakField = 1E9
  
  Test_Keldysh(Egap, meff, PeakField, wavelength, order)
  #exit()
  
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

  #order = 100
  ShowPlot = True

  #timeKeldysh, N_Keldysh_SI, N_Gruzdev_SI, gamma, wPI, wPIg, N_excited_Keldysh_trapz, N_excited_Gruzdev_trapz = generateWpiTables(Egap, meff, wavelength, tau, PeakField, dt, order, N_total, t0, False)
  
  wPI, QGulley, xGulley = GenerateKeldyshGulleyDatabase(Egap, meff, wavelength, PeakField, order)
  
  ### plot w_PI(intensity)
  print(Header+"Importing Gruzdev [2014] data...")
  #try:
  Gulley2012=np.loadtxt("Results/Gulley/Gulley-Fig2.csv", dtype='float', delimiter='\t')
      #Gruzdev2014=np.loadtxt("Gruzdev2014-Fig1.csv", dtype='float', delimiter=',')
      #print Gruzdev2014[:,0]
  #except: 
      #print Header+"** Warning: failed to import Gruzdev2014 data table..."
    
    
  print(Header+"Plotting as function of laser field intensity ...")
  plt.figure()
  plt.xlabel("Intensity (W/cm$^{2}$)")
  plt.ylabel("$w_{PI}$ (cm$^{-3}$ s$^{-1}$)")
  #plt.loglog(1e-4*FieldToIntensity(PeakField.real), 1e-6*1e-15*wPI,  linestyle="-", color="r", label=r"$w_{PI}$ "+ShortRefKeldysh)
  plt.loglog(1e-4*FieldToIntensity(PeakField.real), 1e-6*wPI, linestyle="-", color="b", label=r"$w_{PI}$ "+ShortRefGulley)
  plt.loglog(1e-4*FieldToIntensity(PeakField.real), QGulley, linestyle="-", color="r", label=r"$Q(\gamma,x)$ "+ShortRefGulley)
  plt.loglog(1e-4*FieldToIntensity(PeakField.real), xGulley, linestyle="--", color="r", label=r"$x$ "+ShortRefGulley)
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
    
## Plots the wPI of Keldysh-Gruzdev model for SiO2 as function of intensity. 
def SilicaGruzdev2014(): #{{{
  print("Defining SiO2 material parameters...")

  Egap = 8.97*e; 
  meff=0.6e0; #Effective mass of Si
  N_total=1.*5E28
  numpoints = 1000
  wavelength = 1030e-9;
  if(wavelength==800e-9): 
    optical_index = 1.5356 #800nm
  elif (wavelength==1030E-9):
    optical_index = 1.45 #1030nm
  else: 
    print(Header+"** Warning: optical refractive index was taken equal to 1. ")
    optical_index=1. 
  
  tau=250e-15; dt = 1E-17; CEP=0e0
  
  #PeakFluence = 1.0*1E4 #J/cm2 * 1E4 = J/m2
  #PeakField   = np.sqrt(2e0 * PeakFluence / (tau * c * epsilon_0))
  PeakIntensity_log = np.linspace(np.log10(1e14), np.log10(1e18), numpoints)
  PeakIntensity = np.power(10., PeakIntensity_log)
  PeakField = np.sqrt(2. * PeakIntensity / c / epsilon_0) #I = 0.5 c eps0 E²
  # I_inside_matter = 0.5 c eps0 n0 E**2: in matter, pulse is compressed in space. Therefore, intensity is stronger. As photon energy do not change, the time frequency does not change either.  
  LocalPeakField = np.sqrt(2. * PeakIntensity / optical_index/ c / epsilon_0) #I = 0.5 c n0 eps0 E²
  # NOTE: it looks the field is not affected by the refractive index, but the intensity is. This may originate from the change of cycle durations inside matter. Although the field may not change (without accounting for the induced fields). 
  
  print(Header+"Peak field (min, max): "+str(np.min(PeakField)/1E9)+" V/nm, "+str(np.max(PeakField/1E9))+" V/nm.")

  t0=0. #defines the instant 0.
  
  #print Header+"** Info: PulseEnvelope.EPS and PNG were written in the current folder. "

  order = 100
  ShowPlot = True

  print(Header+"** Computing the W_PI values from Keldysh and Gruzdev theories...")
  
  timeKeldysh, N_Keldysh_SI, N_Gruzdev_SI, gamma, wPI, wPIg, N_excited_Keldysh_trapz, N_excited_Gruzdev_trapz = generateWpiTables(Egap, meff, wavelength, tau, LocalPeakField, dt, order, N_total, t0, False, optical_index)
  
  #timeKeldysh, N_Keldysh_SI, N_Gruzdev_SI, gamma, wPI, LocalwPIg, N_excited_Keldysh_trapz, N_excited_Gruzdev_trapz = generateWpiTables(Egap, meff, wavelength, tau, LocalPeakField, dt, order, N_total, t0, False)
  
  ### plot w_PI(intensity)
  print(Header+"Importing Gruzdev [2014] data...")
  try:
      #Gruzdev2014=np.loadtxt("Gulley-Fig2.csv", dtype='float', delimiter='\t')
      Gruzdev2014=np.loadtxt("Gruzdev2014-Fig1.csv", dtype='float', delimiter=',')
      #print Gruzdev2014[:,0]
  except: 
      print(Header+"** Warning: failed to import Gruzdev2014 data table...")
    
    
  print(Header+"Plotting as function of laser field intensity ...")
  plt.figure()
  plt.xlabel("Intensity (W/cm$^{2}$)")
  plt.ylabel("$w_{PI}$ (cm$^{-3}$ fs$^{-1}$)")
  #plt.loglog(1e-4*FieldToIntensity(PeakField.real), 1e-6*1e-15*wPI,  linestyle="-", color="r", label=r"$w_{PI}$ "+ShortRefKeldysh)
  plt.loglog(1e-4*FieldToIntensity(PeakField.real), 1e-6*1e-15*wPIg, linestyle="-", color="b", label=r"$w_{PI}$ "+ShortRefGruzdev+r", $\lambda=$"+str(wavelength*1E9)+" nm")
  #plt.loglog(1e-4*FieldToIntensity(PeakField.real), 1e-6*1e-15*wPIg, linestyle="-", color="g", label=r"$w_{PI}$ "+ShortRefGruzdev)
  #plt.loglog(Gruzdev2014[:,0], Gruzdev2014[:,1], linestyle="--", color="k", label="Data from "+ShortRefGruzdev) #JUST FOR VALIDATION. 
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

if __name__ == "__main__":
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
