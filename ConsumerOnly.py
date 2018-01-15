#!/usr/bin/env python 

from pylab import * 
from sympy import Symbol
from sympy.matrices import Matrix, eye, zeros, ones 
from sympy.solvers import nsolve
from scipy.stats import norm as normdist 
import matplotlib.pyplot as plt 
from scipy.optimize import root, fsolve, fixed_point

RE      =   0.99
RX      =   0.8
SE      =   4.00
SX      =   3.50
B       =   0.99
H       =   0.30
PERIODS =   5
NSIMS   =   2



class Param:
    """
    Parameter collector
    """
    def __init__(self,re,rx,se,sx,b,H,vd=None):
        """
        Initialize by passing parameters (re,rx,ve,vx,b,H)   
        """
        self.re     = re
        self.rx     = rx
        self.se     = se
        self.sx     = sx
        self.b      = b
        self.H      = H
        self.theta      =   Matrix([-1,self.re,0])
        if vd is None:
            self.sd         =   Symbol('sd')
        else:
            self.sd         =   sd 
        self.m          =   'm'
        self.V          =   'V'
        self.mu0        =   zeros(3,3)
        self.Sigma0     =   Matrix([[self.se**2,0,0],[0,self.sd**2,0],[0,0,self.sx**2]])
        self.I          =   eye(3)
        self.Ii         =   Matrix([[0,0,0],[0,0,0],[0,0,1]])
        self.Ic         =   self.I - self.Ii 
        self.psi1       =   Matrix([1,self.rx-self.re,1])
        self.Psi1       =   self.bayes_coef(self.psi1,self.Sigma0)
        self.beliefs0   =   {self.m : self.mu0, self.V : self.Sigma0}
        self.beliefs1   =   self.bayes_update(self.Psi1,self.beliefs0)
        self.phi        =   (self.b*self.re/(1 - self.b*self.re)*self.theta.transpose()*(self.I - self.beliefs1[self.m])).transpose() + self.Ii*ones(3,1)
        self.phic       =   self.Ic*self.phi
        self.phii       =   self.Ii*self.phi
        #self.bstar      =   normdist.ppf(1-self.H,(loc=0,scale=self.rx**2/(1-self.rx**2) + self.phi[2]**2))
        self.psi2       =   self.phii +  Matrix([0,self.rx,0])
        self.psi2       =   self.phic + Matrix([0,-self.rx])
        self.psi2alt    =   self.phic - 1/(1-self.b*self.re)*self.theta
        self.Psi2       =   self.bayes_coef(self.psi2, self.beliefs1[self.V])
        self.beliefs2   =   self.bayes_update(self.Psi2, self.beliefs1)
        self.beliefs2alt=   self.bayes_update(self.bayes_coef(self.psi2alt,self.beliefs1[self.V]),self.beliefs1)
        self.delxta      =   (self.theta.transpose()*(self.I - self.beliefs2[self.m])).transpose()
        self.deltac     =   self.Ic*self.delta
        self.deltai     =   self.Ii*self.delta 
        self.sdvalue    =   (self.deltac.transpose()*self.beliefs0[self.V]*self.deltac)[0]**0.5
        self.stvalue    =   (self.deltai.transpose()*self.beliefs0[self.V]*self.deltai)[0]**0.5
        self.sol        =   self.solve_root()
        self.sd_sol     =   self.sol.x[0]
        self.Sigma0_sol =   self.Sigma0.subs(self.sd,self.sd_sol)
        self.phi_sol    =   self.phi.subs(self.sd,self.sd_sol)
        self.delta_sol  =   self.delta.subs(self.sd,self.sd_sol)
        

    def solve_root(self,init=SE*10,method='hybr'):
        """ Solve problem using root """ 
        self.root_sol = root(self.sol_eval,array([init]),method=method)
        return self.root_sol 
    def solve_fsolve(self,init=SE*10):
        ''' solve problem '''
        self.fsolve_sol = fsolve(self.sol_eval,array([init]))
        return self.fsolve_sol 
    def solve_nsolve(self,init=SE*10):
        self.nsolve_sol = nsolve(self.sdvalue - self.sd,self.sd,init)
        return self.nsolve_sol 
        
    def sol_eval(self,x):
        return float((self.sdvalue - self.sd).subs(self.sd,x)) 
    def sols_eval(self,x):
        retlist = list()
        for xx in x:
            retlist.append(self.sol_eval(xx))
        return array(retlist).reshape(x.shape)   
    def bayes_coef(self,v=zeros(3,1),S=eye(3)):
        """
        Create bayesian updating coefficient from a vector and variance matrix 
        """   
        return (v.transpose()*S*v).inv()[0]*S*v*v.transpose()
    
    def bayes_update(self,currcoef,initbeliefs):
        """
        Returns updated mean/variance given initial and new information
        """
        mupdate = initbeliefs[self.m] + currcoef*(self.I - initbeliefs[self.m])
        vupdate = (self.I - currcoef)*initbeliefs[self.V]*(self.I - currcoef).transpose() 
        return {self.m : mupdate, self.V : vupdate} 
        

def fp_eval(s,re,rx,se,sx,b):
    re     = re
    rx     = rx
    se     = se
    sx     = sx
    b      = b

    theta      =   Matrix([-1, 0, re,re])
    sd         =   s[0]
    st         =   s[1]
    m          =   'm'
    V          =   'V'
    mu0        =   zeros(4,4)
    Sigma0     =   Matrix([[se**2,0,0, 0],[0,sx**2,0, 0],[0,0,sd**2,0],[0,0,0,st**2]])
    I          =   eye(4)
    Ii         =   Matrix([[0,0,0,0],[0,1,0,0],[0,0,0,0],[0,0,0,1]])
    Ic         =   I - Ii 
    psi1       =   Matrix([1,1,rx-re,rx-re])
    Psi1       =   bayes_coef(psi1,Sigma0)
    beliefs0   =   {m : mu0, V : Sigma0}
    beliefs1   =   bayes_update(Psi1,beliefs0)
    phi        =   (b*re/(1 - b*re)*theta.transpose()*(I - beliefs1[m])).transpose() + Matrix([0,1,0,0])
    phic       =   Ic*phi
    phii       =   Ii*phi
    psi2       =   phii +  Matrix([0,0,rx,rx])
    psi2alt    =   phic - 1/(1-b*re)*theta
    Psi2       =   bayes_coef(psi2, beliefs1[V])
    beliefs2   =   bayes_update(Psi2, beliefs1)
    beliefs2alt=   bayes_update(bayes_coef(psi2alt,beliefs1[V]),beliefs1)
    delta      =   (theta.transpose()*(I - beliefs2[m])).transpose()
    deltac     =   Ic*delta
    deltai     =   Ii*delta 
    
    sdvalue    =   (deltac.transpose()*beliefs0[V]*deltac)[0]**0.5
    stvalue    =   (deltai.transpose()*beliefs0[V]*deltai)[0]**0.5

class Param4:
    """
    Parameter collector
    """
    def __init__(self,re,rx,se,sx,b,H):
        """
        Initialize by passing parameters (re,rx,ve,vx,b,H)   order: eps xi delta tau 
        """
        self.re     = re
        self.rx     = rx
        self.se     = se
        self.sx     = sx
        self.b      = b
        self.H      = H
        self.theta      =   Matrix([-1, 0, self.re,self.re])
        self.sd         =   Symbol('sd')
        self.st         =   Symbol('st')
        self.m          =   'm'
        self.V          =   'V'
        self.mu0        =   zeros(4,4)
        self.Sigma0     =   Matrix([[self.se**2,0,0, 0],[0,self.sx**2,0, 0],[0,0,self.sd**2,0],[0,0,0,self.st**2]])
        self.I          =   eye(4)
        self.Ii         =   Matrix([[0,0,0,0],[0,1,0,0],[0,0,0,0],[0,0,0,1]])
        self.Ic         =   self.I - self.Ii 
        self.psi1       =   Matrix([1,1,self.rx-self.re,self.rx-self.re])
        self.Psi1       =   self.bayes_coef(self.psi1,self.Sigma0)
        self.beliefs0   =   {self.m : self.mu0, self.V : self.Sigma0}
        self.beliefs1   =   self.bayes_update(self.Psi1,self.beliefs0)
        self.phi        =   (self.b*self.re/(1 - self.b*self.re)*self.theta.transpose()*(self.I - self.beliefs1[self.m])).transpose() + Matrix([0,1,0,0])
        self.phic       =   self.Ic*self.phi
        self.phii       =   self.Ii*self.phi
        self.psi2       =   self.phii +  Matrix([0,0,self.rx,self.rx])
        self.psi2alt    =   self.phic - 1/(1-self.b*self.re)*self.theta
        self.Psi2       =   self.bayes_coef(self.psi2, self.beliefs1[self.V])
        self.beliefs2   =   self.bayes_update(self.Psi2, self.beliefs1)
        self.beliefs2alt=   self.bayes_update(self.bayes_coef(self.psi2alt,self.beliefs1[self.V]),self.beliefs1)
        self.delta      =   (self.theta.transpose()*(self.I - self.beliefs2[self.m])).transpose()
        self.deltac     =   self.Ic*self.delta
        self.deltai     =   self.Ii*self.delta 
        
        self.sdvalue    =   (self.deltac.transpose()*self.beliefs0[self.V]*self.deltac)[0]**0.5
        self.stvalue    =   (self.deltai.transpose()*self.beliefs0[self.V]*self.deltai)[0]**0.5

    def solve_root(self,init=array([SE*10,SX*10]),method='hybr'):
        """ Solve problem using root """ 
        self.root_sol = root(self.sol_eval,init,method=method)
        return self.root_sol 
    def solve_fsolve(self,init=SE*10):
        ''' solve problem '''
        self.fsolve_sol = fsolve(self.sol_eval,array([init]))
        return self.fsolve_sol 
    def solve_nsolve(self,init=SE*10):
        self.nsolve_sol = nsolve(self.sdvalue - self.sd,self.sd,init)
        return self.nsolve_sol 
        
    def sol_eval(self,x):
        return array(Matrix([self.sdvalue - self.sd,self.stvalue-self.st]).subs({self.sd:x[0],self.st:x[1]})).astype(float).reshape(x.shape)
    def sols_eval(self,x):
        retlist = list()
        for xx in x:
            retlist.append(self.sol_eval(xx))
        return array(retlist).reshape(x.shape)   
    def bayes_coef(self,v=zeros(3,1),S=eye(3)):
        """
        Create bayesian updating coefficient from a vector and variance matrix 
        """   
        return (v.transpose()*S*v).inv()[0]*S*v*v.transpose()
    
    def bayes_update(self,currcoef,initbeliefs):
        """
        Returns updated mean/variance given initial and new information
        """
        mupdate = initbeliefs[self.m] + currcoef*(self.I - initbeliefs[self.m])
        vupdate = (self.I - currcoef)*initbeliefs[self.V]*(self.I - currcoef).transpose() 
        return {self.m : mupdate, self.V : vupdate} 
        
class Simulations:
    def __init__(self,param,n=1,t=1):
        """ Initialize simulation object """
        self.param  =   param 
        self.nsim   =   n
        self.nper   =   t 
        self.time   =   arange(0,t,1)
        self.states =   list()
        self.pact   =   list()
        self.peff   =   list()
        self.fig    =   plt.figure()
        
        for i in range(self.nsim):
            print("Running simulation %d\n" % (i+1))
            seed(i)
            self.simulate(i)
        self.states = array(self.states)
        self.pact = array(self.pact)
        self.peff = array(self.peff)
    def simulate(self,i):
        """ Seed RNG, draw initial values and iterate forward """
        simstates   =   list() # [etat-1 mut deltat-1 
        simpact     =   list()
        simpeff     =   list()

        currstate   =   Matrix([[normal(0,1/(1-self.param.re**2)*self.param.ve,1)],[normal(0,self.param.ve)],[normal(0,self.param.vd)]])
        
        simstates.append(currstate.copy())
        prices          =   self.compute_prices(currstate)
        simpact.append(prices[0])
        simpeff.append(prices[1])
        for i in range(0,self.nper-1):
            laststate       =   currstate.copy()
            currstate       =   currstate.copy()
            currstate[1]    =   normal(0,self.param.ve)
            currstate[0]    =   self.param.re*laststate[0] + currstate[1]            
            m = Matrix([[currstate[1]],[laststate[2]],[0]])
            mm=self.param.theta.transpose()*(self.param.I - self.param.beliefs2['m'])*self.param.Ic
            currstate[2]    =   mm*m
            prices          =   self.compute_prices(currstate)
            simpact.append(prices[0])
            simpeff.append(prices[1])
            simstates.append(currstate)
        
        self.states.append(array(simstates).transpose())
        self.pact.append(array(simpact))
        self.peff.append(array(simpeff))
            
            
    def compute_prices(self,state):
        pact   =   self.param.bstar/(1-self.param.b) + state[0]/(1-self.param.b*self.param.re)+self.param.phi[0]*state[1] + self.param.phi[1]*state[2]
        peff   =   self.param.bstar/(1-self.param.b) + state[0]/(1-self.param.b*self.param.re)
        return pact, peff
    def plot_prices(self,state):
        """ args: state: a state scalar """ 
        #self.clear_plot()
        plt.clf() 

        plt.plot(self.time+1,self.pact[state],'r',label='Actual Prices')
        plt.plot(self.time+1,self.peff[state],'b',label='Efficient Prices')
        plt.gca().set_xlabel('Period')
        plt.gca().set_ylabel('Price')
        plt.gca().legend() 
    def plot_state(self,state):
        """ args: state, scalar \nPlot state eps/delta """ 
        plt.clf()
        
        plt.plot(self.time+1,self.states[state,1,],'r',label="\varepsilon")
        plt.plot(self.time+1,self.states[state,2,],'b',label="\Delta_{t-1}")
        plt.gca().set_xlabel('Period')
        plt.gca().set_ylabel('\varepsilon/\Delta value')
        plt.gca().legend() 
    def clear_plot(self):
        """ Clear current plot object     """
        self.ax.clear() 
            
        
# PARAMETERS # arrayp




# param = Param(re=RE,rx=RX,se=SE,sx=SX,b=B,H=H)
param = Param(re=RE,rx=RX,se=SE,sx=SX,b=B,H=H)
# print('vd solution: %5.3f' % param.vd)
# print('Running %d simulations for %d periods' % (NSIMS,PERIODS))

# sim = Simulations(param,n=NSIMS,t=PERIODS)
 
# Symbolic System # 

# class SymbolicSystem:
#     """
#     Solve symbolic model 
#     """
#     def __init__(self,param):
#         """
#         Build symbolic matrices/vectors/etc.  
#         """
#         self.vd         =   Symbol('vd')
#         self.m          =   'm'
#         self.V          =   'V'
#         self.param      =   param 
#         self.Sigma0     =   Matrix([[param.ve,0,0],[0,self.vd,0],[0,0,param.vx]])
#         self.I          =   eye(3)
#         self.Ii         =   Matrix([[0,0,0],[0,0,0],[0,0,1]])
#         self.Ic         =   self.I - self.Ii 
#         self.psi        =   Matrix([1,param.rx-param.re,1])
#         self.Psi        =   self.bayes_coef(self.psi,self.Sigma0)
#         self.mu0        =   zeros(3,3)
#         self.beliefs0   =   {self.m : self.mu0, self.V : self.Sigma0}
#         self.beliefs1   =   self.bayes_update(self.Psi,self.beliefs0)
#         self.theta      =   Matrix([-1,param.re,0])
#         self.phi        =   (param.b*param.re/(1 - param.b*param.re)*self.theta.transpose()*(self.I - self.Psi)).transpose() + self.Ii*ones(3,1)
#         self.pi         =   self.phi - 1/(1-param.b*param.re)*self.theta + Matrix([0,param.rx,0]) 
#         self.Pi         =   self.bayes_coef(self.pi, self.beliefs1[self.V])
#         self.beliefs2   =   self.bayes_update(self.Pi, self.beliefs1)
#         self.vdvalue    =   (self.theta.transpose()*self.beliefs2[self.V]*self.theta)[0]
#         self.vd_solve   =   nsolve(self.vdvalue - self.vd,self.vd,50.0)
#         
#         
#     def bayes_coef(self,v=zeros(3,1),S=eye(3)):
#         """
#         Create bayesian updating coefficient from a vector and variance matrix 
#         """   
#         return (v.transpose()*S*v).inv()[0]*S*v*v.transpose()
#     def bayes_update(self,currcoef,initbeliefs):
#         """
#         Returns updated mean/variance given initial and new information
#         """
#         mupdate = initbeliefs[self.m] + currcoef*(self.I - initbeliefs[self.m])
#         vupdate = (self.I - currcoef)*initbeliefs[self.V]
#         return {self.m : mupdate, self.V : vupdate} 
#          

