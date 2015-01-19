# -*- coding: utf-8 -*-
"""
Created on Sun Dec  7 13:41:02 2014

@author: jmmauricio
"""

import sys, os

from pysimu import sim
import numpy as np
from scipy.special import expit as sigmoid
import xlwt
import datetime
ezxf = xlwt.easyxf

sheet_name = "Sheet2"

book = xlwt.Workbook()
sheet = book.add_sheet(sheet_name)


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

p_max=8.0e6              # W

speed_max = 300.0/3.6    # m/s
accel_max = 0.3          # m/s2
# speed**2 = 2*accel*(pos)
# pos = 0.5*accel*t**2

Dpos = speed_max**2 /(2.0*accel_max)

lenght = 168.94e3
Position  = np.array([0.0,   Dpos, lenght-Dpos,    lenght, lenght+100.0])
Speed_ref = np.array([1.0, speed_max, speed_max,    0.0, -10])

f_max = p_max/speed_max*4 # N

# http://en.wikipedia.org/wiki/Rolling_resistance
# http://en.wikipedia.org/wiki/Drag_coefficient
# http://www.engineeringtoolbox.com/drag-coefficient-d_627.html

C_d = 1.2 
rho = 1.2    #
A = 3.0*4.0  # m2

K_p_pos = 100.0
K_p_speed = 5000.0
K_i_speed = 100.0
speed_ref = 300.0/3.6

C_rr = 0.0020
M = 300.0e3


# dynamic system
def f(t,x):
    
    pos = x[0]
    speed = x[1]
    int_speed = x[2]
    

    f_drag = 1.0/2*(rho*speed**2*C_d*A)
    f_roll = C_rr*M*9.81*sigmoid(100.0*speed)
        
    speed_ref = np.interp(pos,Position, Speed_ref )
    
    if speed_ref>speed_max:
        speed_ref = speed_max

    if speed_ref<(-speed_max):
        speed_ref = -speed_max   
        
    error_speed =(speed_ref - speed)
    f_accel = K_p_speed*error_speed + K_i_speed*int_speed
    p = f_accel*speed
    
    error_speed_sat = error_speed 
    if f_accel>f_max:
        f_accel = f_max
        error_speed_sat = 0.0

    if f_accel<(-f_max):
        f_accel = -f_max        
        error_speed_sat = 0.0

    if p>p_max*0.8:
        f_accel = p_max/speed  
        error_speed_sat = 0.0

    if p<(-p_max):
        f_accel = -p_max/speed  
        error_speed_sat = 0.0
        
    dpos   = speed
    dspeed = 1.0/(M)*(f_accel - f_drag - f_roll)
    dint_speed = error_speed_sat
    
    
    return [dpos,dspeed,dint_speed]


# outputs functions
def h(t,x):

    pos = x[0]
    speed = x[1]
    int_speed = x[2]
    speed_ref = np.interp(pos+1,Position, Speed_ref )
    error_speed =(speed_ref - speed)
    f_accel = K_p_speed*error_speed + K_i_speed*int_speed

    p = f_accel*speed
    

        

    f_drag = 1.0/2*(rho*speed**2*C_d*A)
    f_roll = C_rr*M*9.81*sigmoid(speed)
    
    return np.array([pos,speed,f_accel, speed_ref, f_drag, p, f_roll])


# initialization
p_e = p_m
delta_0 = np.arcsin(p_e*X/(e*v))  
omega_0 = omega_s
x_0 = np.array([0.0, 0.0,0.0])

# system definition
simu_1.f = f
simu_1.x_0 = x_0
simu_1.h = h
simu_1.max_step=10.0
simu_1.run(3100.0)

# plot results
import matplotlib.pyplot as plt
fig_1 = plt.figure( figsize=(14, 8))
    
ax_pos = fig_1.add_subplot(3,2,1)
ax_speed = fig_1.add_subplot(3,2,3)
ax_forces = fig_1.add_subplot(3,2,2)
ax_pos_speed = fig_1.add_subplot(3,2,4)
ax_powers = fig_1.add_subplot(3,2,5)

ax_pos.plot(simu_1.T,simu_1.X[:,0]/1.0e3, linewidth=2)
ax_pos.plot([simu_1.T[0],simu_1.T[-1]],[lenght/1000.0]*2, linewidth=2)

ax_speed.plot(simu_1.T,simu_1.X[:,1]*3.6, linewidth=2)
ax_speed.plot(simu_1.T,simu_1.Y[:,3]*3.6, linewidth=2)
ax_forces.plot(simu_1.T,simu_1.Y[:,2], label='$\sf f_{accel}$')
ax_forces.plot(simu_1.T,simu_1.Y[:,4], label='$\sf f_{drag}$')
ax_forces.plot(simu_1.T,simu_1.Y[:,6], label='$\sf f_{roll}$')
ax_pos_speed.plot(simu_1.X[:,0],simu_1.X[:,1], linewidth=2)
ax_pos_speed.plot(simu_1.X[:,0],simu_1.Y[:,3], linewidth=2)
ax_powers.plot(simu_1.X[:,0],simu_1.Y[:,5]/1.0e6, linewidth=2)
#ax_delta_omega.plot(Delta,P_e/P_e*p_m, label='$\sf p_m$', linewidth=2)
#ax_delta_omega.plot(simu_1.X[:,0],simu_1.Y[:], 'r', linewidth=2)

ax_pos.set_ylabel('$\sf Position $ (km)')

ax_speed.set_xlabel('Time (s)')
ax_speed.set_ylabel('Speed (km/h)')
#
#ax_delta_omega.set_xlabel('$\sf \delta $ (rad)')
#ax_delta_omega.set_ylabel('$\sf Power $ (p.u.)')

ax_pos.grid(True)
ax_speed.grid(True)
ax_forces.grid(True)
ax_powers.grid(True)
publish = False
if publish: 
    import plotly.plotly as py

    py.sign_in("jmmauricio", "rwdnrmvuyg")
    plot_url = py.plot_mpl(fig_1, auto_open=True)
    
else:
    
    fig_1.savefig('trajectory.png')
    
N = len(simu_1.X[:,0])    
for row in range(N):
    
    xf = xlwt.easyxf(num_format_str='MM:SS') 
    
    p = simu_1.Y[row,5]/1.0e3
    
    p_consumida = 0.0
    p_devuelta  = 0.0
    if p>0.0:
        p_consumida = p
        
    if p<0.0:
        p_devuelta  = -p
        
    sheet.write(row+4-1, 0, simu_1.X[row,0]/1000.0)
    sheet.write(row+4-1, 1, simu_1.T[row]/86400.0, xf)
    sheet.write(row+4-1, 2, simu_1.X[row,1]/3.6)
    sheet.write(row+4-1, 3, p_consumida)
    sheet.write(row+4-1, 4, p_devuelta)
    
    p = simu_1.Y[N-row-1,5]/1.0e3
    
    p_consumida = 0.0
    p_devuelta  = 0.0
    if p>0.0:
        p_consumida = p
        
    if p<0.0:
        p_devuelta  = -p    
    
    
    
    sheet.write(row+4-1, 6, simu_1.X[N-row-1,0]/1000.0)
    sheet.write(row+4-1, 7, simu_1.T[row]/86400.0, xf)
    sheet.write(row+4-1, 8, simu_1.X[N-row-1,1]/3.6)
    sheet.write(row+4-1, 9, p_consumida)
    sheet.write(row+4-1, 10, p_devuelta)
book.save('hola.xls')
