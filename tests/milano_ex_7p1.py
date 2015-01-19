# -*- coding: utf-8 -*-
"""
Created on Sun Dec  7 13:41:02 2014

@author: jmmauricio
"""

import sys, os

from pysimu import sim
import numpy as np
 
simu_1 = sim()

# parameters
p_m = 1.0
X = 0.5
e = 1.0
v = 1.0
H = 3.5
omega_s = 1.0
omega = omega_s
D = 1.0
Omega_b = 2.0*np.pi*50.0


# dynamic system
def f(t,x):
    
    delta = x[0]
    omega = x[1]
    
    p_e = e*v/X*np.sin(delta)
    
    ddelta = Omega_b*(omega - omega_s)
    domega = 1.0/(2*H)*(p_m - p_e - D*(omega - omega_s))
    
    return [ddelta, domega]


# outputs functions
def h(t,x):

    delta = x[0]
    omega = x[1]
    
    p_e = e*v/X*np.sin(delta)
    
    return np.array(p_e)


# initialization
p_e = p_m
delta_0 = np.arcsin(p_e*X/(e*v))  
omega_0 = omega_s
x_0 = np.array([delta_0, omega_0])

# system definition
simu_1.f = f
simu_1.x_0 = x_0
simu_1.h = h
simu_1.run(1.0)
v = 0.0

fault_ms = 200.0
simu_1.run(1.0+fault_ms/1000.0)
v = 1.0
simu_1.x_0 = simu_1.x
simu_1.run(5.0)

Delta = np.linspace(0.0, np.pi,100)
P_e = e*v/X*np.sin(Delta)


# plot results
import matplotlib.pyplot as plt
fig_1 = plt.figure( figsize=(14, 8))
    
ax_delta = fig_1.add_subplot(2,2,1)
ax_omega = fig_1.add_subplot(2,2,3)
ax_delta_omega = fig_1.add_subplot(2,2,(2,4))

ax_delta.plot(simu_1.T,simu_1.X[:,0], linewidth=2)
ax_omega.plot(simu_1.T,simu_1.X[:,1], linewidth=2)
ax_delta_omega.plot(Delta,P_e, label='$\sf p_e$')
ax_delta_omega.plot(Delta,P_e/P_e*p_m, label='$\sf p_m$', linewidth=2)
ax_delta_omega.plot(simu_1.X[:,0],simu_1.Y[:], 'r', linewidth=2)

ax_delta.set_ylabel('$\sf \delta $ (rad)')

ax_omega.set_xlabel('Time (s)')
ax_omega.set_ylabel('$\sf \omega $ (p.u.)')

ax_delta_omega.set_xlabel('$\sf \delta $ (rad)')
ax_delta_omega.set_ylabel('$\sf Power $ (p.u.)')

ax_delta.grid(True)
ax_omega.grid(True)
ax_delta_omega.grid(True)

publish = False
if publish: 
    import plotly.plotly as py

    py.sign_in("jmmauricio", "rwdnrmvuyg")
    plot_url = py.plot_mpl(fig_1, auto_open=True)
    
else:
    
    fig_1.savefig('machine_1_delta_omega_{:d}.png'.format(int(fault_ms)))