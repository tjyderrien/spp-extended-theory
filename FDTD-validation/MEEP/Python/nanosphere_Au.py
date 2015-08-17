#!/usr/bin/env python
#-*- coding: utf-8 -*-
""" run_sim.py - an example python-meep simulation of a dielectric sphere scattering a broadband impulse, 
illustrating the use of the convenient functions provided by meep_utils.py 
(c) 2014 Filip Dominec, see http://fzu.cz/~dominecf/meep/ for more information """

import numpy as np
import time, sys, os
import meep_utils, meep_materials
from meep_utils import in_sphere, in_xcyl, in_ycyl, in_zcyl, in_xslab, in_yslab, in_zslab, in_ellipsoid
import meep_mpi as meep
#import meep
c = 2.997e8

sim_param, model_param = meep_utils.process_param(sys.argv[1:])
class SphereWire_model(meep_utils.AbstractMeepModel): #{{{
    def __init__(self, comment="", simtime=100e-12, resolution=100e-9, cells=1, monzc=0e-6, padding=2e-6,
            radius=500e-9, spacing=1e-6, wlth=0.8e-6, wtth=1e-6, Kx=0, Ky=0, dist=2e-6):
        meep_utils.AbstractMeepModel.__init__(self)        ## Base class initialisation
        self.simulation_name = "SphereWire"    
        monzd=spacing

        self.register_locals(locals())          ## Remember the parameters

        ## Constants for the simulation
        self.pml_thickness = 2e-6
        self.monitor_z1, self.monitor_z2 = (-(monzd*cells)/2+monzc-padding, (monzd*cells)/2+monzc+padding)  
        self.simtime = simtime      # [s]
        self.srcFreq, self.srcWidth = 1E16, 1E16     # [Hz] (note: gaussian source ends at t=10/srcWidth)
        self.interesting_frequencies = (1E15, 1E16)     # Which frequencies will be saved to disk

        self.size_x = spacing 
        self.size_y = spacing
        self.size_z = 2e-6 + cells*monzd + 2*self.pml_thickness + 2*self.padding

        ## Define materials
        self.materials = [meep_materials.material_Au(where = self.where_TiO2)]  
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

# Model selection
model = SphereWire_model(**model_param)
if sim_param['frequency_domain']: model.simulation_name += ("_frequency=%.4e" % sim_param['frequency'])

## Initialize volume and structure according to the model
vol = meep.vol3d(model.size_x, model.size_y, model.size_z, 1./model.resolution)
vol.center_origin()
s = meep_utils.init_structure(model=model, volume=vol, sim_param=sim_param, pml_axes=meep.Z)

## Create fields with Bloch-periodic boundaries (any transversal component of k-vector 
## is allowed, but may not radiate)
f = meep.fields(s)
f.use_bloch(meep.X, -model.Kx/(2*np.pi))
f.use_bloch(meep.Y, -model.Ky/(2*np.pi))

## Add a source of the plane wave (see meep_utils for definition of arbitrary source shape)
if not sim_param['frequency_domain']:           ## Select the source dependence on time
    #src_time_type = meep.band_src_time(model.srcFreq/c, model.srcWidth/c, model.simtime*c/1.1)
    src_time_type = meep.gaussian_src_time(model.srcFreq/c, model.srcWidth/c)
else:
    src_time_type = meep.continuous_src_time(sim_param['frequency']/c)
srcvolume = meep.volume( 
        meep.vec(-model.size_x/2, -model.size_y/2, -model.size_z/2+model.pml_thickness),
        meep.vec( model.size_x/2,  model.size_y/2, -model.size_z/2+model.pml_thickness))
f.add_volume_source(meep.Ex, src_time_type, srcvolume)

## Define monitors planes and visualisation output
monitor_options = {'size_x':model.size_x, 'size_y':model.size_y, 'Kx':model.Kx, 'Ky':model.Ky}
monitor1_Ex = meep_utils.AmplitudeMonitorPlane(comp=meep.Ex, z_position=model.monitor_z1, **monitor_options)
monitor1_Hy = meep_utils.AmplitudeMonitorPlane(comp=meep.Hy, z_position=model.monitor_z1, **monitor_options)
monitor2_Ex = meep_utils.AmplitudeMonitorPlane(comp=meep.Ex, z_position=model.monitor_z2, **monitor_options)
monitor2_Hy = meep_utils.AmplitudeMonitorPlane(comp=meep.Hy, z_position=model.monitor_z2, **monitor_options)

slice_makers =  [meep_utils.Slice(model=model, field=f, 
		components=(meep.Dielectric), at_t=0, name='EPS')]
slice_makers += [meep_utils.Slice(model=model, field=f, 
		components=meep.Ex, at_x=0, min_timestep=.05e-12, outputgif=True)]
slice_makers += [meep_utils.Slice(model=model, field=f, 
		components=meep.Ex, at_t=2.5e-12)]

if not sim_param['frequency_domain']:       ## time-domain computation
    f.step()
    dt = (f.time()/c)
    meep_utils.lorentzian_unstable_check_new(model, dt)
    timer = meep_utils.Timer(simtime=model.simtime); meep.quiet(True) # use custom progress messages
    while (f.time()/c < model.simtime):                               # timestepping cycle
        f.step()
        timer.print_progress(f.time()/c)
        for monitor in (monitor1_Ex, monitor1_Hy, monitor2_Ex, monitor2_Hy): monitor.record(field=f)
        for slice_maker in slice_makers: slice_maker.poll(f.time()/c)
    for slice_maker in slice_makers: slice_maker.finalize()
    meep_utils.notify(model.simulation_name, run_time=timer.get_time())
else:                                       ## frequency-domain computation
    f.step()
    f.solve_cw(sim_param['MaxTol'], sim_param['MaxIter'], sim_param['BiCGStab']) 
    for monitor in (monitor1_Ex, monitor1_Hy, monitor2_Ex, monitor2_Hy): monitor.record(field=f)
    for slice_maker in slice_makers: slice_maker.finalize()
    meep_utils.notify(model.simulation_name)

## Get the reflection and transmission of the structure
if meep.my_rank() == 0:
    freq, s11, s12 = meep_utils.get_s_parameters(monitor1_Ex, monitor1_Hy, monitor2_Ex, monitor2_Hy, 
            frequency_domain=sim_param['frequency_domain'], frequency=sim_param['frequency'], 
            maxf=model.srcFreq+model.srcWidth, pad_zeros=1.0, Kx=model.Kx, Ky=model.Ky)
    meep_utils.savetxt(freq=freq, s11=s11, s12=s12, model=model)
    #import effparam        # process effective parameters for metamaterials

with open("./last_simulation_name.txt", "w") as outfile: outfile.write(model.simulation_name) 
meep.all_wait()         # Wait until all file operations are finished
