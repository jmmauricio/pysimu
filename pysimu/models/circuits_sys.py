# -*- coding: utf-8 -*-
"""
Created on Thu Aug 14 20:25:02 2014

@author: jmmauricio-m
"""
import numpy as np
import matplotlib.pyplot as plt
from pysimu.models import pwm
from scipy.integrate import odeint
dt_large = 0.9523809523809523e-3
dt_large = 100.0e-6
dt_short = 0.5e-6


class rl:
    
    def __init__(self):
        
        self.R = 0.5
        self.L = 30.0e-3

        
        self.chan_list = ['t','v','x', 'i_l']
        
        self.chan = {}
        for item in self.chan_list:
            self.chan.update({item:np.array([])})
            

    def setup(self):
        
        self.x = np.array([0.0])
        self.v = np.array([0.0])
        self.y = np.array([0.0])
        self.t = 0.0
        self.i_l = self.x
        
    def initialization(self):
        
        self.x = self.v
        self.y = self.x

        for item in self.chan_list:
            self.chan[item] = self.__dict__[item]
            
            
        return self.x
        
    def derivatives(self,t,x):

        v = self.v
        i_l = x 
        
        di_l = 1.0/self.L*(v - self.R*i_l)
        
        dx = di_l
        self.update(t,x)
        self.dx = dx
        return dx

    def update(self,t,x):
        
        y = x
        self.i_l = x 
        self.y = y     
        
    def output(self,t,x):
        
        self.t = t
        self.x = x

        for item in self.chan_list:
            self.chan[item] = np.hstack((self.chan[item],self.__dict__[item]))
#            
class rl_pwm:
    
    def __init__(self):
        
        self.R = 0.5
        self.L = 30.0e-3

        
        self.chan_list = ['t','v','x', 'i_l', 'eta']
        
        self.chan = {}
        for item in self.chan_list:
            self.chan.update({item:np.array([])})
            

    def setup(self):
        
        self.x = np.array([0.0])
        self.v = np.array([0.0])
        self.y = np.array([0.0])
        self.eta = np.array([0.0])
        self.t = 0.0
        self.i_l = self.x
        self.t_prev = 0.0
        
    def initialization(self):
        
        self.x = self.v
        self.y = self.x

        for item in self.chan_list:
            self.chan[item] = self.__dict__[item]
            
            
        return self.x
        
    def derivatives(self,t,x):
        
        
        v = self.v
        i_l = x 
        
        di_l = 1.0/self.L*(v - self.R*i_l)
        
        dx = di_l
        
        self.update(t,x)
        
        self.dx = dx
        return dx

    def update(self,t,x):
        eta = self.eta
        dt_large = t-self.t_prev
        dt_large = 10.0e-6
#        if dt_large<1e-6:
#            dt_large = 1e-6
        out_top, out_bottom, out_average, carrier_1, carrier_2 = pwm.pwm_1(t,5000.0,eta,dt_large,dt_short,dt_large)
        self.v = out_average*200.0 
        y = x
        self.i_l = x 
        self.y = y     
        self.t_prev = t
    def output(self,t,x):
        
        self.t = t
        self.x = x
        
        
        for item in self.chan_list:
            self.chan[item] = np.hstack((self.chan[item],self.__dict__[item]))

class csc_rl:
    
    def __init__(self):
        
        self.R_dc = 0.5
        self.L_dc = 30.0e-3
        self.R_l = 0.3 
        self.L_l = 30.0e-3
        self.C = 20.0e-6
        self.v = 100.0
        self.v_dc_p = self.v
        
#        self.chan_list = ['t','i_dc1','v_dc_p','i_inv1', 'v_c', 'eta', 'pwm_prom', 'i_load']
        self.chan_list = ['t','i_dc1','v_dc_p','i_inv1', 'v_c', 'eta', 'pwm_prom', 'i_load']
        self.chan = {}
        for item in self.chan_list:
            self.chan.update({item:np.array([])})
            

    def setup(self):
        
        self.x = np.zeros((3,1))
#        self.v = np.array([0.0])
        self.y = np.array([0.0])
        
        self.eta = np.array([0.0])
        
        self.t = 0.0
        self.i_l = self.x
        self.t_prev = 0.0
        
    def initialization(self):
        
        self.x = np.zeros((3,1))
        self.i_dc1  = self.x[0] 
        self.v_c    = self.x[1]
        self.i_load = self.x[2]
        self.i_inv1 = 0.0
        self.pwm_prom = 0.0 
        

        for item in self.chan_list:
            self.chan[item] = self.__dict__[item]
            
            
        return self.x
        
    def derivatives(self,t,x):
        
        
        v    = self.v        
        C    = self.C
        R_l  = self.R_l
        L_dc = self.L_dc
        R_dc = self.R_dc    
        L_l  = self.L_l 
        
        i_dc1  = x[0] 
        v_c    = x[1]
        i_load = x[2]
        
        eta = 0.5*np.sin(2.0*np.pi*50.0*t)
        self.eta =eta
        self.update(t,x)
        
        pwm_prom = self.pwm_prom
        i_inv1 = pwm_prom*i_dc1
        
        # p_dc = v_dc_p*i_dc1 = p_ac = v_c*i_inv1
        # p_dc = v_dc_p*i_dc1 = p_ac = v_c*eta*i_dc1
        
        
        v_dc_p = v_c*pwm_prom
        
        di_dc1  = 1.0/L_dc*(v - 1.0*R_dc*i_dc1 - v_dc_p)
        dv_c    = 1.0/C*(i_inv1 - i_load)
        di_load = 1.0/L_l*(v_c - R_l*i_load)
        
        dx = np.vstack((di_dc1,dv_c,di_load))
        
        self.v_dc_p = v_dc_p
        self.eta = eta
        self.pwm_signal = pwm_prom
        self.i_inv1 = i_inv1
        
        
        
        self.dx = dx
        return dx

    def update(self,t,x):
        
        self.i_dc1  = x[0] 
        self.v_c    = x[1]
        self.i_load = x[2]
        
        
        eta = self.eta
#        dt_large = t-self.t_prev
#        dt_large = 10.0e-6
##        if dt_large<1e-6:
##            dt_large = 1e-6
        out_top, out_bottom, out_average, carrier_1, carrier_2 = pwm.pwm_udc(t,21.0*50,eta,dt_large,dt_short,5.0e-6)
        self.pwm_prom = out_average 
#        y = x
#        self.i_l = x 
#        self.y = y     
#        self.t_prev = t
    def output(self,t,x):
        
        self.t = t
        self.x = x
        
        
        for item in self.chan_list:
            self.chan[item] = np.hstack((self.chan[item],self.__dict__[item]))
            
            
class solver(object):
    
    def __init__(self,derivatives):
        
        self.derivatives = derivatives
    
    def set_integrator(self, integrator, Dt=1.0e-6,max_step =100.0e-6, nsteps=20000):
        self.Dt = Dt
        self.integrator = integrator
        self.max_step = max_step
        self.nsteps = nsteps
        
    def set_solout(self,output):
        
        self.output = output
        
    def set_initial_value(self, initialize, x0):
        
        self.initialize = initialize
        self.x0 = x0
        self.x = self.x0
        self.t = 0.0
        
    def integrate(self,t_end):
        x_0 = self.x
        Dt = self.Dt
        t = self.t

        if self.integrator == 'euler':
            while t < t_end:  
                
                dx = self.derivatives(t,x_0)
                
                x = x_0 + Dt*dx
#                print t, x, dx, Dt
                self.t = t
                self.x = x
                self.output(t,x)                
                x_0 = x
                t += Dt

        if self.integrator == 'trapezoidal':
            while t < t_end:  
                
                dx_1 = self.derivatives(t,x_0)
                dx_2 = dx_1 
                for it_trap in range(2):
                    dx_2 = self.derivatives(t,x_0+Dt*dx_2)
                x = x_0 + Dt*0.5*(dx_1+dx_2)
#                print t, x, dx, Dt
                self.t = t
                self.x = x
                self.output(t,x)                
                x_0 = x
                t += Dt
                
        x_0 = x_0.reshape((3,))
        
        if self.integrator == 'odeint':
            while t < t_end:  
                
                t_array = np.arange(t,t+2*Dt,Dt)
 
                x = odeint(self.derivatives_2, x_0,t_array)
#                print x
                
#                print t, x, dx, Dt
                self.t = t
                self.x = x
                self.output(t,x)               
                x_0 = x[-1,:]
                t += Dt

    
    def derivatives_2(self,x,t):
#        print t,x
        dx = self.derivatives(t,x)
        return dx.reshape((3,))
            
def test_rl():
    
    from scipy.integrate import ode
    sys_rl = rl()

                  
    sys_rl.setup()                    
    r = solver(sys_rl.derivatives)
    r.set_integrator('euler', max_step =100.0e-6, nsteps=20000)
    r.set_solout(sys_rl.output)
    r.set_initial_value(sys_rl.initialization(),sys_rl.x)
    r.integrate(0.01)
    sys_rl.v=110.0
    r.integrate(0.8)    

    fig = plt.figure( figsize=(8, 8))
        
    ax_u  = fig.add_subplot(311)
    ax_y  = fig.add_subplot(312)
    ax_x = fig.add_subplot(313)

    ax_u.plot(sys_rl.chan['t'],sys_rl.chan['v'].T)
    ax_y.plot(sys_rl.chan['t'],sys_rl.chan['i_l'].T)
#    ax_powers.plot(sys_1.chan['t'],sys_1.chan['y']) 
    plt.show()    
    
    
def test_rl_pwm():
    
    from scipy.integrate import ode
    sys_rl = rl_pwm()

                  
    sys_rl.setup()                    
    r = solver(sys_rl.derivatives)
    r.set_integrator('euler', Dt=10.0e-6, max_step =100.0e-6, nsteps=2000000)
    r.set_solout(sys_rl.output)
    r.set_initial_value(sys_rl.initialization(),sys_rl.x)
    r.integrate(0.01)
    sys_rl.eta=0.8
    sys_rl.v_dc = 200.0
    r.integrate(0.8)    

    fig = plt.figure( figsize=(8, 8))
        
    ax_u  = fig.add_subplot(311)
    ax_y  = fig.add_subplot(312)
    ax_x = fig.add_subplot(313)

    ax_u.plot(sys_rl.chan['t'],sys_rl.chan['v'].T)
    ax_y.plot(sys_rl.chan['t'],sys_rl.chan['i_l'].T)
#    ax_powers.plot(sys_1.chan['t'],sys_1.chan['y']) 
    plt.show()        
    return sys_rl


def test_csc_rl():
    
    from scipy.integrate import ode
    sys_csc_rl = csc_rl()

                  
    sys_csc_rl.setup()                    
    r = solver(sys_csc_rl.derivatives)
#    r.set_integrator('odeint', Dt=5.0e-6, max_step =100.0e-6, nsteps=2000000)
    r.set_integrator('trapezoidal', Dt=dt_large, max_step =100.0e-6, nsteps=2000000)
    r.set_solout(sys_csc_rl.output)
    r.set_initial_value(sys_csc_rl.initialization(),sys_csc_rl.x)
    r.integrate(0.5)
#    sys_csc_rl.eta=0.8
##    sys_csc_rl.v_dc = 200.0
#    r.integrate(0.08)    
    data = np.load('/home/jmmauricio/Documents/private/prom/udec_pwm_rl_1.npz' ) 
    fig = plt.figure( figsize=(8, 8))
    fig_fourier   = plt.figure( figsize=(8, 8))
    
    ax_u  = fig.add_subplot(311)
    ax_y  = fig.add_subplot(312)
    ax_x = fig.add_subplot(313)
    
    ax_fourier = fig_fourier.add_subplot(111)

    ax_u.plot(sys_csc_rl.chan['t'],sys_csc_rl.chan['eta'].T)
    ax_u.plot(sys_csc_rl.chan['t'],sys_csc_rl.chan['pwm_prom'].T)
    ax_u.plot(data['Time'],data['tria1']) 
#    ax_y.plot(sys_csc_rl.chan['t'],sys_csc_rl.chan['v_c'].T)
    
    ax_y.plot(sys_csc_rl.chan['t'],sys_csc_rl.chan['i_inv1'].T)
    ax_y.plot(data['Time'],data['iinva1'])    
    ax_y.legend(['prom', 'psim'])
    

    ax_x.plot(sys_csc_rl.chan['t'],sys_csc_rl.chan['v_c'].T)
    ax_x.plot(data['Time'],data['vo']) 




    plt.show()        
    np.savez('/home/jmmauricio/Documents/private/prom/prom_sim.npz', t=sys_csc_rl.chan['t'],
                                                                     v_c=sys_csc_rl.chan['v_c'],
                                                                     i_inv1=sys_csc_rl.chan['i_inv1'],
                                                                      pwm_prom=sys_csc_rl.chan['pwm_prom'],
                                                                      eta=  sys_csc_rl.chan['eta']  )
    return sys_csc_rl
        
def test_fft():
    

    data = np.load('/home/jmmauricio/Documents/private/prom/udec_pwm_rl_1.npz' ) 
    data_chan = np.load('/home/jmmauricio/Documents/private/prom/prom_sim.npz')
    
    fig = plt.figure( figsize=(8, 8))
    fig_fourier   = plt.figure( figsize=(8, 8))
    
    ax_u  = fig.add_subplot(311)
    ax_y  = fig.add_subplot(312)
    ax_x = fig.add_subplot(313)
    
    ax_fourier = fig_fourier.add_subplot(111)


    ax_u.plot(data['Time'],data['tria1']) 

    ax_y.plot(data['Time'],data['iinva1'])    
    ax_y.legend(['prom', 'psim'])
    


    ax_x.plot(data['Time'],data['vo']) 



    n = len(data['Time']) # Number of data points
    dx = data['Time'][1]-data['Time'][0] # Sampling period (in meters)
    x = dx*np.arange(0,n) # x coordinates
    fx = data['vo']
    Fk = np.fft.fft(fx)
    nu = np.fft.fftfreq(n,dx) # Natural frequencies
    
    n_prom = len(data_chan['t']) # Number of data points
    dx_prom = data_chan['t'][2]-data_chan['t'][1] # Sampling period (in meters)
    x_prom = dx_prom*np.arange(0,n_prom) # x coordinates
    fx_prom = data_chan['v_c']
    Fk_prom = np.fft.fft(fx_prom)
    nu_prom = np.fft.fftfreq(n_prom,dx_prom) # Natural frequencies
    
#    Fk = fft.fftshift(Fk) # Shift zero freq to center
#    Fkabs=np.absolute(Fk)**2 ##Power spectrum
#    nu = fft.fftshift(nu) # Shift zero freq to center
    ax_fourier.bar(nu,np.abs(Fk)/n)
    ax_fourier.bar(nu_prom,np.abs(Fk_prom)/n_prom)
    ax_fourier.set_xlim((0,10000.0))
    plt.show()        

    
    
def read_psim():

    path_to_csv = '/home/jmmauricio/Documents/private/prom/udec_pwm_rl_1.csv'    
    data = np.genfromtxt(path_to_csv, dtype=float, delimiter=',', names=True) 
    
    Time = data['Time'][0:-1:10]
    iinva1 = data['iinva1'][0:-1:10]
    tria1 = data['tria1'][0:-1:10]
    m1 = data['m1'][0:-1:10]
    vo = data['vo'][0:-1:10]

    
    np.savez('/home/jmmauricio/Documents/private/prom/udec_pwm_rl_1.npz',Time = Time,
                                                                         iinva1 = iinva1,
                                                                         tria1 = tria1,
                                                                         m1 = m1,
                                                                         vo = vo )        

    
    data = np.load('/home/jmmauricio/Documents/private/prom/udec_pwm_rl_1.npz' )      
    
if __name__ == '__main__':   
        
#    test_rl()
#    sys_rl = test_rl_pwm()
#    read_psim()
#    sys_csc_rl = test_csc_rl()
    
#    test_fft()
    #test_observer()
    #test_first_order()
    test_csc_rl()
        
        
        
        
    
    
