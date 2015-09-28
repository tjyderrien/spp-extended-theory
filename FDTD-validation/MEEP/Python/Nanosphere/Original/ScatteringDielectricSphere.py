#!/usr/bin/env python
#-*- coding: utf-8 -*-
""" run_sim.py - an example python-meep simulation of a dielectric sphere scattering a broadband impulse, 
illustrating the use of the convenient functions provided by meep_utils.py 
(c) 2014 Filip Dominec, see http://fzu.cz/~dominecf/meep/ for more information """

# Initialisation, imports

import numpy as np
import time, sys, os
import meep_utils, meep_materials
from meep_utils import in_sphere, in_xcyl, in_ycyl, in_zcyl, in_xslab, in_yslab, in_zslab, in_ellipsoid
import meep_mpi as meep
#import meep
c = 2.997e8

sim_param, model_param = meep_utils.process_param(sys.argv[1:])


# Defining model

class SphereWire_model(meep_utils.AbstractMeepModel): #{{{
    def __init__(self, comment="", simtime=10e-12, resolution=6e-6, cells=1, monzc=0e-6, padding=50e-6,
            radius=25e-6, spacing=75e-6, wlth=10e-6, wtth=10e-6, Kx=0, Ky=0, dist=25e-6):
        meep_utils.AbstractMeepModel.__init__(self)        ## Base class initialisation
        self.simulation_name = "SphereWire"    
        monzd=spacing

        self.register_locals(locals())          ## Remember the parameters

        ## Constants for the simulation
        self.pml_thickness = 20e-6
        self.monitor_z1, self.monitor_z2 = (-(monzd*cells)/2+monzc-padding, (monzd*cells)/2+monzc+padding)  
        self.simtime = simtime      # [s]
        self.srcFreq, self.srcWidth = 1000e9, 2000e9     # [Hz] (note: gaussian source ends at t=10/srcWidth)
        self.interesting_frequencies = (0e9, 2000e9)     # Which frequencies will be saved to disk

        self.size_x = spacing 
        self.size_y = spacing
        self.size_z = 60e-6 + cells*monzd + 2*self.pml_thickness + 2*self.padding

        ## Define materials
        self.materials = [meep_materials.material_TiO2_THz(where = self.where_TiO2)]  
        if not 'NoMetal' in comment:
            self.materials += [meep_materials.material_Metal_THz(where = self.where_metal) ]
        self.TestMaterials()

        f_c = 2.997e8 / np.pi/self.resolution/meep_utils.meep.use_Courant()
        meep_utils.plot_eps(self.materials, mark_freq=[f_c])

    # each material has one callback, used by all its polarizabilities (thus materials should never overlap)
    def where_metal(self, r):
        if (in_yslab(r, cy=0e-6, d=self.wtth)) and in_zslab(r, cz=0, d=self.wlth):
            return self.return_value            # (do not change this line)
        return 0

    def where_TiO2(self, r):
        if (in_sphere(r, cx=0, cy=-self.size_x/2, cz=0, rad=self.radius) or
                in_sphere(r, cx=0, cy=self.size_x/2, cz=0, rad=self.radius)):
            return self.return_value             # (do not change this line)
        return 0
#}}}


# Initializing the instances model, voluem and structure
# Model selection
model = SphereWire_model(**model_param)
if sim_param['frequency_domain']: model.simulation_name += ("_frequency=%.4e" % sim_param['frequency'])

## Initialize volume and structure according to the model
vol = meep.vol3d(model.size_x, model.size_y, model.size_z, 1./model.resolution)
vol.center_origin()
s = meep_utils.init_structure(model=model, volume=vol, sim_param=sim_param, pml_axes=meep.Z)
