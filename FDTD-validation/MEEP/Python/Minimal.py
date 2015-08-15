#!/usr/bin/env python
import meep
x, y, z, voxelsize = 2e-6, 2e-6, 5e-6, 50e-9
vol = meep.vol3d(x, y, z, 1/voxelsize)

class Model(meep.Callback):
    def double_vec(self, r):
        if (r.x()-1e-6)**2+(r.y()-1e-6)**2+(r.z()-2.5e-6)**2 < .6e-6**2:
            return 10
        else: 
            return 1

model = Model()
meep.set_EPS_Callback(model.__disown__())
struct = meep.structure(vol, meep.EPS) 

f = meep.fields(struct)
f.add_volume_source(meep.Ex, 
        meep.gaussian_src_time(300e12/3e8, 500e12/3e8),    # (setting frequency = dividing by speed of light)
        meep.volume(meep.vec(0,0,.1e-6), meep.vec(2e-6,2e-6,.1e-6)))

while f.time()/3e8 < 30e-15:   # (retrieving time = dividing by speed of light)
    f.step()
    
print f.get_field(meep.Ex, meep.vec(.5e-6, .5e-6, 3e-6))

output = meep.prepareHDF5File("Ex.h5")
f.output_hdf5(meep.Ex, vol.surroundings(), output)
del(output)