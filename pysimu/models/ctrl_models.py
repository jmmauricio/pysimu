# -*- coding: utf-8 -*-
"""
Created on Thu Aug 14 20:25:02 2014

@author: jmmauricio-m
"""
import numpy as np
class first_order:
    
    def __init__(self):
        
        self.T = 0.001
        
        self.chan_list = ['t','u','y']
        
        self.chan = {}
        for item in self.chan_list:
            self.chan.update({item:np.array([])})
            

    def setup(self):
        
        self.x = np.array([0.0])
        self.u = np.array([0.0])
        self.y = np.array([0.0])
        self.t = 0.0
        
    def initialization(self):
        
        self.x = self.u
        self.y = self.x

        for item in self.chan_list:
            self.chan[item] = self.__dict__[item]
            
            
        return self.x
        
    def derivatives(self,t,x):
        T = self.T
        u = self.u
        
        dx = 1.0/T*(u - x)
        self.update(t,x)
        self.dx = dx
        return dx

    def update(self,t,x):
        
        y = x
        
        self.y = y     
        
    def output(self,t,x):
        
        self.t = t
        self.x = x
        
        
        for item in self.chan_list:
            self.chan[item] = np.hstack((self.chan[item],self.__dict__[item]))
#            
    
    
class lti:
    '''
    
    Implements lti model 
    
    
    '''
    
    def __init__(self):
        
         
        self.dx = np.zeros((2,1))        
        self.x = np.zeros((2,1))        

        self.n_x = len(self.x)

        self.A = np.array([]) 
        self.B = np.array([])
        self.C = np.array([])
        self.D = np.array([])
        
        
        self.chan_list = ['t','x','y','u']
        
        self.chan = {}
        for item in self.chan_list:
            self.chan.update({item:np.array([])})
    
    def setup(self,A,B,C,D):

        self.A = A
        self.B = B
        self.C = C
        self.D = D
        
        self.n_x = A.shape[1]
        self.n_u = B.shape[1]
        self.n_y = C.shape[0]
        
        self.x_0 = np.zeros((self.n_x,1))
        self.u_0 = np.zeros((self.n_u,1))
        
    def initialization(self):
        self.t = 0.0
        self.x = self.x_0
        self.u = self.u_0
        self.y = np.dot(self.C,self.x) + np.dot(self.D,self.u)
        
        for item in self.chan_list:
            self.chan[item] = self.__dict__[item]
            
            
        return self.x
        
    def derivatives(self,t,x):
        
#        self.perturbation(t)
        self.x=x.reshape((self.n_x,1))
        dx = np.dot(self.A,self.x) + np.dot(self.B,self.u)
#        self.perturbation(t)
        self.update(t,x)
        return dx

    def update(self, t,x):
        self.t = t 
        self.y = np.dot(self.C,self.x) + np.dot(self.D,self.u)

        
    def output(self, t,x):
        self.t = t 
        self.y = np.dot(self.C,self.x) + np.dot(self.D,self.u)
        
        for item in self.chan_list:
            self.chan[item] = np.hstack((self.chan[item],self.__dict__[item]))
#            
    def perturbation(self,t):
        self.u[0] = 0.0 
        if t>1.0:
            self.u[0] = 100.0 

 

class observer:
    '''
    
    Observer
    
    Imlements Luenberger observer, continuous-time case
    
    
    ..math::    
         {\dot \hat{x}} = A \hat x + B u + L \left( y - C \hat x \right)  
    
    http://en.wikipedia.org/wiki/State_observer
    
    '''
    
    def __init__(self):

        self.A = np.array([]) 
        self.B = np.array([])
        self.C = np.array([])
        self.D = np.array([])
        
         
        self.chan_list = ['t','x_est']
        
        self.chan = {}
        for item in self.chan_list:
            self.chan.update({item:np.array([])})
    
    def setup(self,A,B,C,D,K,L):

        self.A = A
        self.B = B
        self.C = C
        self.D = D
        self.K = K
        self.L = L
        
        self.n_x = A.shape[1]
        self.n_u = B.shape[1]
        self.n_y = C.shape[0]
        self.dx = np.zeros((self.n_x ,1))        
        self.x = np.zeros((self.n_x ,1))         
        self.x_0 = np.zeros((self.n_x,1))
        self.u_0 = np.zeros((self.n_u,1))
        
    def initialization(self):
        self.t = 0.0
        self.x_est = self.x_0
        self.u = self.u_0
        self.y = np.dot(self.C,self.x_est) + np.dot(self.D,self.u)
        
        for item in self.chan_list:
            print item
            self.chan[item] = self.__dict__[item]
            
            
        return self.x
        
    def derivatives(self,t,x):
        
        A,B,C,D,K,L = self.A,self.B,self.C,self.D,self.K,self.L       
        x_est=x.reshape((self.n_x,1))        
        u = self.u
        y = self.y
        
        dx_est = np.dot(A,x_est) + np.dot(B,u) + np.dot(L,(y-np.dot(C,x_est)))

        self.dx = dx_est
        return self.dx
        
    def output(self, t,x):
        self.t = t 
        self.x_est = x.reshape((self.n_x,1))
        
        for item in self.chan_list:
            self.chan[item] = np.hstack((self.chan[item],self.__dict__[item]))
            
            
            
            

def test_lti():
    from scipy.integrate import ode
    import matplotlib.pyplot as plt
    
    sys_1 = lti()
    
    X_pu = 0.3
    H=5.0
    Omega_b = 314.0
    D = 5.0
    A = np.array([[0 ,                     Omega_b],
                  [-1.0/(X_pu*2*H) ,     -D/(2*H)]])

    B = np.array([[0.0],
                  [1.0/(2*H)]])                  

    C = np.array([[1.0, 0.0],
                  [0.0, 1.0]])

    D = np.array([[0.0],
                  [0.0]])
                  
    sys_1.setup(A,B,C,D)                    
    r = ode(sys_1.derivatives)
    r.set_integrator('dopri5', max_step =0.1, nsteps=20000)
    r.set_solout(sys_1.output)
    r.set_initial_value(sys_1.initialization(),sys_1.x)
    r.integrate(30.0)
    
    import matplotlib.pyplot as plt
    fig = plt.figure( figsize=(8, 8))
        
    ax_delta_omega = fig.add_subplot(211)
    ax_powers = fig.add_subplot(212)

    ax_delta_omega.plot(sys_1.chan['t'],sys_1.chan['x'].T)
    ax_powers.plot(sys_1.chan['t'],sys_1.chan['u'].T)
#    ax_powers.plot(sys_1.chan['t'],sys_1.chan['y']) 
    plt.show() 

def test_first_order():
    
    from scipy.integrate import ode
    sys_lp = first_order()

                  
    sys_lp.setup()                    
    r = ode(sys_lp.derivatives)
    r.set_integrator('dopri5', max_step =0.1, nsteps=20000)
    r.set_solout(sys_lp.output)
    r.set_initial_value(sys_lp.initialization(),sys_lp.x)
    r.integrate(1.0)
    sys_lp.u=1.0
    r.integrate(2.0)    
    import matplotlib.pyplot as plt
    fig = plt.figure( figsize=(8, 8))
        
    ax_u  = fig.add_subplot(311)
    ax_y  = fig.add_subplot(312)
    ax_x = fig.add_subplot(313)

    ax_u.plot(sys_lp.chan['t'],sys_lp.chan['u'].T)
    ax_y.plot(sys_lp.chan['t'],sys_lp.chan['y'].T)
#    ax_powers.plot(sys_1.chan['t'],sys_1.chan['y']) 
    plt.show()    
    
    
    
def test_observer():
    from scipy.integrate import ode
    import matplotlib.pyplot as plt
    import control as ctrl
    

    X_pu = 0.3
    H=5.0
    Omega_b = 314.0
    D = 5.0
    A = np.array([[0 ,                     Omega_b],
                  [-1.0/(X_pu*2*H) ,     -D/(2*H)]])

    B = np.array([[0.0],
                  [1.0/(2*H)]])                  

    C = np.array([[1.0, 0.0],
                  [0.0, 1.0]])

    D = np.array([[0.0],
                  [0.0]])


    # L Observador
    #--------------------------------------------------------------------------
    R_o = 1.0 
    Q_o = 1.0
    Q_obs = Q_o*np.eye(A.shape[0])
    R_obs = R_o*np.eye(C.shape[0])
    
     
    K, S, w_ctrl = ctrl.lqr(A.T,C.T,Q_obs,R_obs)
    L_obs = K.T

    sys = lti()
    sys_obs = observer()
                  
    sys.setup(A,B,C,D)    
    sys_obs.setup(A,B,C,D,K,L_obs)       

    sys.x_0 = np.array([[ 0.0],
                      [ 0.1 ]])
                      
    r_sys = ode(sys.derivatives)
    r_sys.set_integrator('dopri5', max_step =0.05, nsteps=20000)
    r_sys.set_solout(sys.output)
    r_sys.set_initial_value(sys.initialization(),sys.x)

    
    r_obs = ode(sys_obs.derivatives)
    r_obs.set_integrator('dopri5', max_step =0.05, nsteps=20000)
    r_obs.set_solout(sys_obs.output)
    r_obs.set_initial_value(sys_obs.initialization(),sys_obs.x)
    


                   
    for t_end in np.arange(0.0,10.0,0.01):
        sys.u = np.array([[ 0.]])
        if t_end > 1.0:
            sys.u = np.array([[ 0.2]])
            
        
        sys_obs.y = sys.y
        sys_obs.u = sys.u
        
        r_sys.integrate(t_end) 
        r_obs.integrate(t_end) 
        
    
    
    import matplotlib.pyplot as plt
    fig = plt.figure( figsize=(8, 8))
        
    ax_delta  = fig.add_subplot(311)
    ax_omega  = fig.add_subplot(312)
    ax_powers = fig.add_subplot(313)

    ax_delta.plot(sys.chan['t'],sys.chan['x'][0,:].T)
    ax_delta.plot(sys_obs.chan['t'],sys_obs.chan['x_est'][0,:].T)
    ax_omega.plot(sys.chan['t'],sys.chan['x'][1,:].T)
    ax_omega.plot(sys_obs.chan['t'],sys_obs.chan['x_est'][1,:].T)
    ax_delta.legend(['$\sf \delta$','$\sf \hat \delta$'])
    ax_omega.legend(['$\sf \omega$','$\sf \hat \omega$'])
    ax_powers.plot(sys.chan['t'],sys.chan['u'].T)
#    ax_powers.plot(sys_1.chan['t'],sys_1.chan['y']) 
    plt.show()    
    
    
if __name__ == '__main__':   
        
#    test_lti()
    #test_observer()
    test_first_order()

    
    
