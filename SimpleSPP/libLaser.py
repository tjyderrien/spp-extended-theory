#!/usr/bin/env python
#-*- coding: utf-8 -*-
## @package libLaser
# Functions to describe the laser pulse

import numpy as np
from scipy.constants import c, epsilon_0, mu_0, pi, e, m_e, h

## calculate laser frequency (Hz) from wavelength (m)
def omega(wavelength):#{{{
    return 2.0*pi*c/wavelength
#}}}

omega = np.vectorize(omega) 