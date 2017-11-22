import libKeldyshMcDonald
import numpy as np
import matplotlib.pyplot as plt

#Usage:
#libKeldyshMcDonald.solve(t, field_params, t_init, init_conds, K, T_2)
#---------------------------------------------------------------------
#t is the time parameter
#we are interested in n_v(t), n_c(t)
#---------------------------------------------------------------------
#field_params is a tuple:
#field_params = (peak, omega, FWHM, t0, phi)
#
#peak is the peak value of the field
#omega is the angular frequency
#FWHM is the pulse FWHM, total length of a squared cosine pulse is 2*FWHM
#t0 is the pulse center
#phi is the phase shift
#
#all according to this formula:
#Field(t) = peak * cos( pi*t/(2*FWHM) )^2 * cos(omega*t + phi)
#---------------------------------------------------------------------
#t_init is the initial time
#---------------------------------------------------------------------
#init_conds is a tuple:
#init_conds = (pi_init, n_v_init, n_c_init, S_init)
#
#pi_init is the initial value for pi in t=t_init
#etc.
#S_init it the initial value for the action (related to the divergence problem)
#---------------------------------------------------------------------
#K is a tuple:
#K = (K_x, K_y, K_z)
#---------------------------------------------------------------------
#T_2 is the damping time
#---------------------------------------------------------------------
#the function returns n_v(t), n_c(t)

#example:
times = np.linspace(0, 3, 200)
results = []
for t in times:
    results.append(libKeldyshMcDonald.solve(t,
                                            (1, 3, 1, 1, 0), #field_params
                                            0,               #initial time
                                            (1+1j, 2, 3, 0), #init_conds
                                            (0, 0, 0),       #K
                                            1                #T_2
                                            ))
plt.plot(times, results)
plt.grid(True)
plt.show()
