# -*- coding: utf-8 -*-
"""
Created on Thu Aug 14 20:29:22 2014

@author: jmmauricio-m
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Aug  4 19:20:06 2014

@author: jmmauricio
"""
import numpy as np
from scipy.integrate import ode

class sys_freq_model_1:
    '''
    
    Power system simplified model for frequency studies
    The system is in pu with system base: S_sys
    
    
    '''
    
    def __init__(self):
        
        self.dx = np.zeros((2,1))        
        self.x = np.zeros((2,1))        

        self.n_x = len(self.x)
        
        # machine data
        self.S_g = 1000.0 # Machine based power MVA
        self.H = 10.0     # s (pu-m)
        self.D = 0.1      # (pu-m)   
        self.T_g = 2.0   # s
        self.R = 0.05     # droop (pu-m) 
        

        
        # power system data
        self.S_sys = 1.0  # System base power MVA
        self.H_sys = self.H*self.S_g/self.S_sys # s (pu-s)

        # load data        
        self.p_l_mw = 2000.0 # Load MW        
        self.p_l = self.p_l_mw/self.S_sys
        self.p_l_0 = self.p_l
        
        # initialization        
        self.p_c = self.p_l # Consigned power (pu-s)
        
        self.freq_0 = 1.0
        
        self.freq = self.freq_0     # frequency (pu-s)
        self.p_g = self.p_c        # Mechanical power (pu-s)
        
        # solver and output params
        self.ts_out = 0.1
        self.prev_t_out = 0.0
        
        self.p_nc = 0.0
        
        self.chan_list = ['t','freq','p_l','p_g']
        
        self.chan = {}
        for item in self.chan_list:
            self.chan.update({item:np.array([])})
        
    def ini(self):
        
        self.p_nc = self.p_g-self.p_l
        self.p_nc_0 = self.p_nc
        
        self.x = np.vstack((self.freq,self.p_g))
        
        return self.x
        
    def f(self,t,x):
        
        freq = x[0]   # pu     
        p_g = x[1]    # pu
          
        H_sys = self.H_sys
        D = self.D
        S_sys = self.S_sys
        S_g = self.S_g
        
        p_l = self.p_l
        T_g= self.T_g 
        p_c = self.p_c        
        p_nc = self.p_nc
        
        freq_0 = self.freq_0
        
        Dfreq = freq -freq_0
        
        p_ref = p_c - 1.0/self.R*Dfreq*S_g/S_sys 
      
        dfreq = 1.0/(2*H_sys)*(p_g + p_nc - p_l )
        dp_g = 1.0/T_g*(p_ref - p_g )

             
        dx = np.vstack((dfreq,dp_g))
        self.dx = dx
#        print dx   

        self.p_l = p_l
        self.freq = freq            
        self.prev_t_out = t
        
        return dx
        
    def update(self):
        self.freq = self.x[0]
        self.p_g = self.x[1] 
#        print self.freq 

        
    def out(self, t,x):
        self.t = t        
        for item in self.chan_list:
            self.chan[item] = np.hstack((self.chan[item],self.__dict__[item]))
            
    def perturbation(self,t):
        if t>0.1:
            self.p_l = self.p_l_0*1.1        


class gen_nc:
    
    def __init__(self):
        
        self.dx = np.zeros((1,1))        
        self.x = np.zeros((1,1))        
        self.n_x = len(self.x)

        # machine data
        self.K_f = 10.0 # Machine based power MVA
        self.p_nc = 0.0        
        self.freq = 1.0
        self.T_nc = 0.1
        
        self.chan_list = ['p_nc']
        
        self.chan = {}
        for item in self.chan_list:
            self.chan.update({item:np.array([])})
        
    def ini(self):
        
        self.x =self.p_nc  
        self.p_nc_0 = self.p_nc 
        return self.x
        
    def f(self,t,x):
        
        p_nc = x   # pu     
        p_nc_0 = self.p_nc_0
        T_nc = self.T_nc       
        freq = self.freq
        K_f = self.K_f
#        print K_f*(freq - 1) 
        p_f = - K_f*(freq - 1)
        if p_f>100.0:
            p_f = 100.0
        dp_nc = 1.0/(T_nc)*(p_nc_0  + p_f - p_nc )
               
        dx = dp_nc
        self.dx = dx
#        print self.dx
        return dx
     
    def update(self):
        self.p_nc = self.x


    def out(self, t,x):
        
        self.t = t
        self.p_nc = x[0]
        for item in self.chan_list:
            self.chan[item] = np.hstack((self.chan[item],self.__dict__[item]))

#        self.p_turbines=np.hstack((self.p_turbines, self.p_m))
#        self.power_edac=np.hstack((self.power_edac, self.p_edac))

#------------------------------------------------------------------


class gen_classic_inf:
    '''
    
    Model of classical generator connected to infinite bus
    
    
    '''
    
    def __init__(self):
        
        self.dx = np.zeros((2,1))        
        self.x = np.zeros((2,1))        

        self.n_x = len(self.x)
        
        # machine data
        self.S_g = 1000.0 # Machine based power MVA
        self.H = 10.0     # s (pu-m)
        self.D = 0.1      # (pu-m)   
      

        
        # power system data


        # load data        
        self.p_l_mw = 2000.0 # Load MW        
        self.p_l = self.p_l_mw/self.S_sys
        self.p_l_0 = self.p_l
        
        # initialization        
        self.p_c = self.p_l # Consigned power (pu-s)
        
        self.freq_0 = 1.0
        
        self.freq = self.freq_0     # frequency (pu-s)
        self.p_g = self.p_c        # Mechanical power (pu-s)
        
        # solver and output params
        self.ts_out = 0.1
        self.prev_t_out = 0.0
        
        self.p_nc = 0.0
        
        self.chan_list = ['t','freq','p_l','p_g']
        
        self.chan = {}
        for item in self.chan_list:
            self.chan.update({item:np.array([])})
        
    def ini(self):
        
        self.p_nc = self.p_g-self.p_l
        
        self.x = np.vstack((self.freq,self.p_g))
        
        return self.x
        
    def f(self,t,x):
        
        delta = x[0]   # pu     
        omega = x[1]    # pu

#        self.perturbation(t)
        
#        for item in self.__dict__:
#            to_exec = item + ' = self.' + item
#            print to_exec
#            exec(to_exec)    



               
        H = self.H
        D = self.D
        
        p_m = self.p_m

        v_inf = self.v_inf
        e = self.e
        x_pu = self.x_pu
        Omega_b = self.Omega_b
       
        p_e = e*v_inf/x_pu*np.sin(delta)
         
        ddelta = Omega_b*(omega-1)               
        domega = 1.0/(2.0*H)*(p_m - p_e )



             
        dx = np.vstack((ddelta,domega))
        
        self.dx = dx
#        print dx   
        self.delta = delta
        self.omega = omega

        
        return dx
        
    def out(self, t,x):
        self.t = t        
        for item in self.chan_list:
            self.chan[item] = np.hstack((self.chan[item],self.__dict__[item]))
            
    def perturbation(self,t):
        if t>0.1:
            self.p_l = self.p_l_0*1.1        


class gen_nc_2:
    
    def __init__(self):
        
        self.dx = np.zeros((1,1))        
        self.x = np.zeros((1,1))        
        self.n_x = len(self.x)

        # machine data
        self.K_f = 10.0 # Machine based power MVA
        self.p_nc = 0.0
        
        self.freq = 1.0
        
        self.T_nc = 0.5
        
        self.chan_list = ['p_nc']
        
        self.chan = {}
        for item in self.chan_list:
            self.chan.update({item:np.array([])})
        
    def ini(self):
        
        self.x =self.p_nc
        
        return self.x
        
    def f(self,t,x):
        
        p_nc = x[0]   # pu     
        
        T_nc = self.T_nc
        
        freq = self.freq
        K_f = self.K_f
        
        dp_nc = 1.0/(T_nc)*(K_f*(freq - 1) - p_nc )
        
        
        dx = dp_nc
        
#        if  (t-self.prev_t_out)>self.ts_out:
#            self.p_g = p_g
#            self.p_l = p_l
#            self.freq = freq            
#            self.out(t,x)
#            self.prev_t_out = t
            

#        self.p_m = p_m
        return dx
        
    def update(self):
        self.p_nc = self.x[0]

        
    def out(self, t,x):
        
        self.t = t
        self.p_nc = x[0]
        for item in self.chan_list:
            self.chan[item] = np.hstack((self.chan[item],self.__dict__[item]))

#        self.p_turbines=np.hstack((self.p_turbines, self.p_m))
#        self.power_edac=np.hstack((self.power_edac, self.p_edac))

#-----------------------------------------------------------------

    
class system:
    
    def __init__(self):
        
        self.dx = np.zeros((3,1))        
        self.x = np.zeros((3,1))        

        self.H_1 = 10.0     # s (pu-m)
        self.S_g1 = 1000.0 # MVA
        self.S_g2 = 1000.0 # MVA
        self.T_g1 = 10.0
        self.T_g2 = 10.0
        self.S_sys = 1.0 # MVA
        self.H_t = self.H_1*self.S_g1/self.S_sys
        self.p_l = 2000.0
        self.p_c_1 = self.p_l*self.S_sys/self.S_g1/2.0
        self.p_g1_pu = self.p_c_1
        self.p_c_2 = self.p_l*self.S_sys/self.S_g2/2.0
        self.p_g2_pu = self.p_c_2
        self.freq = np.array([])
        self.time = np.array([])
        self.p_turbines=np.array([])
        self.power_edac=np.array([])
        self.omega_pu_0 = 50.25/50.0
        self.p_m = 0.0
        self.max_step = 1.0
        self.t_pert = 2.0
        self.p_pert = -145.0
        self.R_1 = 0.15
        self.R_2 = 0.15
        self.p_edac = 0.0
        self.edac_f =  np.array([49.0, 48.9, 48.8, 48.7, 48.6, 48.5, 48.4, 48.3])
        self.edac_p_mw = np.array([45.3, 50.4, 89.8, 84.8, 87.0, 80.8, 123.4, 111.9])
    def ini(self):
        
        self.x = np.vstack((self.omega_pu_0,self.p_g1_pu ,self.p_g2_pu))
        
        return self.x
        
    def f(self,t,x):
        
        omega_pu = x[0]        
        p_g1_pu = x[1]
        p_g2_pu = x[2]
        
        H_t = self.H_t
        S_sys = self.S_sys
        S_g1 = self.S_g1
        S_g2 = self.S_g2
        p_l = self.p_l
        T_g1= self.T_g1 
        T_g2= self.T_g2 
        p_c_1 = self.p_c_1
        p_c_2 = self.p_c_2        
        Domega_pu = (omega_pu - self.omega_pu_0)
        
        p_1_ref = p_c_1 - 1.0/self.R_1*Domega_pu/S_sys
        p_2_ref = p_c_2 - 1.0/self.R_2*Domega_pu/S_sys   
        
        p_pert  = 0.0
        if t>self.t_pert:
            p_pert = self.p_pert
        
        f = omega_pu*50.0
        new_p_edac = np.sum(self.edac_p_mw[self.edac_f>f])
#        print new_p_edac,   self.p_edac      
        if self.p_edac<new_p_edac:
            self.p_edac = new_p_edac
            
            
            
        p_l = self.p_l - self.p_edac/S_sys + 000.0*(omega_pu -self.omega_pu_0)
        
        p_m = p_g1_pu*S_g1  + p_g2_pu*S_g2       
        
        domega_pu = 1.0/(2*H_t)*(p_m + p_pert - p_l )
        dp_g1_pu = 1.0/T_g1*(p_1_ref - p_g1_pu )
        dp_g2_pu = 1.0/T_g2*(p_2_ref - p_g2_pu )
        
        dx = np.vstack((domega_pu,dp_g1_pu,dp_g2_pu))
        self.p_m = p_m
        self.p_l = p_l
        return dx
        
    def out(self, t,x):
        
        self.freq=np.hstack((self.freq, 50.0*x[0]))
        self.time=np.hstack((self.time, t))
        self.p_turbines=np.hstack((self.p_turbines, self.p_m))
        self.power_edac=np.hstack((self.power_edac, self.p_edac))

def test_freq():
    sys_1 = sys_freq_model_1()
    
    r = ode(sys_1.f)
    r.set_integrator('dopri5', max_step =0.1, nsteps=20000)
    r.set_solout(sys_1.out)
    r.set_initial_value(sys_1.ini(),0.0)
    r.integrate(30.0)
    
    import matplotlib.pyplot as plt
    fig_freq = plt.figure( figsize=(8, 8))
        
    ax_freq = fig_freq.add_subplot(211)
    ax_freq.plot(sys_1.chan['t'],50.0*sys_1.chan['freq'])

    ax_powers = fig_freq.add_subplot(212)
    ax_powers.plot(sys_1.chan['t'],sys_1.chan['p_l'])
    ax_powers.plot(sys_1.chan['t'],sys_1.chan['p_g']) 
    plt.show()
    return sys_1
        
if __name__ == '__main__':
    
    sys_1 = test_freq()
    
    

#    ax_speed = fig_delta_speed.add_subplot(212)        
#        ax_V = fig_v.add_subplot(111)
#        ax_pg = fig_pg_qg.add_subplot(211)
#        ax_qg = fig_pg_qg.add_subplot(212)