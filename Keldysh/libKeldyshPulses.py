#!/usr/bin/env python2
#-*- coding: utf-8 -*-

# Copyright (C) 2013-2017 T. J.-Y. Derrien
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

## @package libKeldyshPulses
# This module defines various types of laser pulses and belongs to the libKeldysh package. 
# Several flavors of the Keldysh theory are available: 
# - Keldysh original paper in solid, for Kane band structure [compared with td-dft]
# - Keldysh paper with few terms corrected by Gruzdev [compared with td-dft]
# - Keldysh-Zhukov tables, where Keldysh theory was computed numerically without using the saddle point method
# - Keldysh-Shcheblanov model, improving rigor on the analytical integration [https://arxiv.org/abs/1706.07303]
# - Keldysh-Corkum model, allowing for analytical treatment of mulltiwavelength fields [Physical Review Letters, 2017, 118, 173601]

from scipy.constants import c, epsilon_0, mu_0, pi, e, m_e, h, hbar
import numpy as np

## Computes a step Heaviside function. 
def step(x):
    return 1.0 * (x > 0.0)


## Single pulse shape [table of peak_intensity(time)] evolution with time
#
# Defines the temporal shape of the laser pulse using a Gaussian law. 
# Polarzation is linear and constant
# @param t: instant to output (can be a table)
# @param tau: pulse duration (s)
# @param PeakIntensity: peak intensity (W/m^2)
# @param t0: instant for the peak intensity (t0=0 by default)
def PulseGaussianTemporalShape(t, tau, PeakIntensity, t0=0.):
  sigmaTau = sigmaFWHM(tau)
  #PeakIntensity = fluence/tau 
  #TODO: Missing coefficient on peak intensity ? 
  intensity = PeakIntensity * np.exp(-0.5 * ((t-t0)/(sigmaTau))**2 )
  return intensity

## Single pulse shape [table of peak_intensity(time)] evolution with time
# Defines the temporal shape of the laser pulse using a squared sinus law. 
# Polarization is linear and constant. 
# Outputs: <Array of real-valued field envelope, array of complex electric field>
# @param t: instants to output (can be a table)
# @param tau: pulse duration FWHM (s)
# @param PeakField: peak of the electric field envelope (V/m) (scalar only)
# @param t0: central instant for the laser pulse (t0=0 by default)
# @param PulseDelay: temporal delay between 2 pulses (in seconds)
def PulseSquaredSinTemporalShape(t, tau, PeakField, wavelength, CEP=0., t0=0., PulseDelay=0.):
  t1 = t0 + PulseDelay
  omega = 2e0*pi*c/wavelength
  H1 = step(t - t1 + tau) #! theer could be a mistake in pulse duration here!
  H2 = step(t - t1 - tau)
  Envelope = PeakField*np.sin(pi*(t-t1-tau)/(2e0*tau))**2 * H1 * (1.-H2)
  Phase = np.exp(1e0j*(omega*t+CEP))
  Field = Envelope * Phase
  return Envelope, Field

## Bi-color double pulse [table of TotalEnvelope(time), TotalField(time)] evolution with time (POLARIZATION IS FOR NOW NEGLECTED!)
# Output: Total envelope <array>, total field <array> at a given space point. 
# Construct the temporal shape of two-color laser pulses mixed together using a squared sinus law and a time delay. 
# Pulses CAN be of different wavelengths! 
# @param t: instants to output (can be a table)
# @param tau1: pulse 1 duration FWHM (s)
# @param tau2: pulse 2 duration FWHM (s)
# @param Efield1: pulse 1, (scalar) peak field of the envelope (V/m)
# @param Efield2: pulse 2, (scalar) peak field of the envelope (V/m)
# @param wavelength1: pulse 1, wavelength (meters)
# @param wavelength2: pulse 2, wavelength (meters)
# @param CEP1: pulse 1, carrier envelope phase
# @param CEP2: pulse 2, carrier envelope phase
# @param t1: instant for the peak field 1 (t1=0 by default)
# @param PulseDelay: delay between the amplitude maxima of pulse 1 and pulse 2 (seconds)
def PulseSquaredSinTemporalShapeDoublePulse(t, tau1, tau2, Efield1, Efield2, wavelength1, wavelength2, CEP1, CEP2, t1=0., PulseDelay=0.):
  omega1=2.*pi*c/wavelength1; omega2=2.*pi*c/wavelength2
  #sigmaTau1 = sigmaFWHM(tau1); sigmaTau2 = sigmaFWHM(tau2) #good for purely gaussian pulse, mmh? 
  t2 = t1 + PulseDelay
  H11        = step(t - t1 + tau1); H21 = step(t - t2 + tau2)
  H12        = step(t - t1 - tau1); H22 = step(t - t2 - tau2)
  FieldEnv1     = Efield1*np.sin(pi*(t-t1-tau1)/(2e0*tau1))**2 * H11 * (1.-H12) #could be bugged
  FieldEnv2     = Efield2*np.sin(pi*(t-t2-tau2)/(2e0*tau2))**2 * H21 * (1.-H22) #could be bugged
  Phase1 = np.exp(1e0j*(omega1*t+CEP1))
  Phase2 = np.exp(1e0j*(omega2*t+CEP2))
  #TotalEnvelope = np.sqrt( FieldEnv1*np.conj(FieldEnv1) + FieldEnv2*np.conj(FieldEnv2) + FieldEnv1*np.conj(FieldEnv2) * np.exp(1e0j*(omega1-omega2)*t) + np.conj(FieldEnv1)* FieldEnv2 * np.exp(1e0j*(omega2-omega1)*t) ) #complex square of the fields must provide the envelope
  #BUG: the phase does not work properly! Rederive the following formula! 
  TotalEnvelope = np.sqrt( FieldEnv1*np.conj(FieldEnv1) + FieldEnv2*np.conj(FieldEnv2) + FieldEnv1*np.conj(FieldEnv2) * np.exp(1e0j*((omega1-omega2)*t+CEP1-CEP2)) + np.conj(FieldEnv1)* FieldEnv2 * np.exp(1e0j*((omega2-omega1)*t+CEP2-CEP1)) ) #complex square of the fields must provide the envelope
  TotalField = FieldEnv1*Phase1 + FieldEnv2*Phase2
  return TotalEnvelope, TotalField

## Vectorial single pulse shape [table of peak_intensity(time)] evolution with time
# Defines the temporal shape of the laser pulse using a squared sinus law. 
# Polarization is taken into account. !  
# Outputs: <Array of real-valued field envelope, array of complex electric field>
## TODO: which formula to take to account for a complete vectorial pulse? 
## Single pulse shape [table of peak_intensity(time)] evolution with time
# Defines the temporal shape of the laser pulse using a squared sinus law. 
# Polarization is linear and constant. 
# Outputs: List of 3 <Temporal array of complex electric field>
# @param t: instants to output (can be a table)
# @param tau: pulse duration FWHM (s)
# @param PeakField: peak of the electric field envelope (V/m) (scalar only)
# @param wavelength: in meters (SI)
# @param PolarizationAngle (rad): projects the polarization in linear direction. 
# @param CEP: Carrier Envelope Phase inside the pulse envelope
# @param t0: central instant for the laser pulse (t0=0 by default)
# @param PulseDelay: temporal delay between 2 pulses (in seconds)
def PulseSquaredSinTemporalShape_vectorial_linear(t, tau, PeakField, wavelength, PolarizationAngle=0., CEP=0., t0=0., PulseDelay=0.):
  t1 = t0 + PulseDelay
  omega = 2e0*pi*c/wavelength
  H1 = step(t - t1 + tau) #! theer could be a mistake in pulse duration here!
  H2 = step(t - t1 - tau)
  EnvelopeX = PeakField*np.cos(PolarizationAngle)*np.sin(pi*(t-t1-tau)/(2e0*tau))**2 * H1 * (1.-H2)
  EnvelopeY = PeakField*np.sin(PolarizationAngle)*np.sin(pi*(t-t1-tau)/(2e0*tau))**2 * H1 * (1.-H2)
  Phase = np.exp(1e0j*(omega*t+CEP))
  FieldX = EnvelopeX * Phase
  FieldY = EnvelopeY * Phase
  return FieldX, FieldY, 0.

## TODO: when PolarizationAngle is time-dependent, we get elliptical polarization. 
## TODO: when phase velocity on various components are different, we get TWISTED LIGHT. This is same as bicolor interaction. 
  
step                                          = np.vectorize(step)
PulseGaussianTemporalShape                    = np.vectorize(PulseGaussianTemporalShape)
PulseSquaredSinTemporalShape                  = np.vectorize(PulseSquaredSinTemporalShape)
PulseSquaredSinTemporalShapeDoublePulse       = np.vectorize(PulseSquaredSinTemporalShapeDoublePulse)
PulseSquaredSinTemporalShape_vectorial_linear = np.vectorize(PulseSquaredSinTemporalShape_vectorial_linear)