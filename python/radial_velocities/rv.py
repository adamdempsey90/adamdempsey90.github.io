import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import scipy.signal as sig
from scipy.optimize import curve_fit,fsolve
from copy import copy

G = 4*np.pi**2 /(356.25)**2 # G in AU^3/(day^2 Solar Mass)

class Star():
    def __init__(self,name,mass,rv_data,mass_err=0,jitter=0):
        self.name = name
        self.mass = mass
        self.t = rv_data[:,0]
        self.vr = rv_data[:,1]
        self.vr_err = rv_data[:,2]
        self.mass_err = mass_err
        self.vr_err = np.sqrt( self.vr_err**2 + jitter**2)
        self.jitter = jitter
        
    def plot(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.errorbar(self.t,self.vr,yerr=self.vr_err,fmt='o')
        ax.set_xlabel('Time  (Julian Days)',fontsize=15)
        ax.set_ylabel('Radial Velocity (m/s)',fontsize=20)
        ax.set_title(self.name + ',    M = %.2f Msun' % self.mass)
        ax.minorticks_on()
        ax.grid()
        plt.show()

    

    def get_best_n(self):
    
        periods = 10**np.linspace(-3,4,1000)
        omega = 2*np.pi/periods
        res=sig.lombscargle(self.t-self.t[0],self.vr,omega)
        max_pers = np.argsort(res)
        return omega[max_pers[-1]]
    def periodogram(self):
        periods = 10**np.linspace(-3,4,1000)
        omega = 2*np.pi/periods
        power=sig.lombscargle(self.t-self.t[0],self.vr,omega)
        self.periodogram_plot(periods,power)

    def periodogram_plot(self,periods,power,pmax=None,ax=None):
        if ax == None:
            fig=plt.figure()
            ax=fig.add_subplot(111)
        ax.semilogx(periods,power)
        ax.set_xlabel('Period (days)',fontsize=15)
        ax.set_ylabel('Power',fontsize=20)
        if pmax != None:
            ax.axvline(pmax,color='k',linestyle='--',linewidth=3)
    def recover_params(self,*popt):
        """Recover the mass of the planet"""

        n = popt[0]
        k = popt[2]
        w = popt[3]
        e = popt[4]
       
        
        if e<0:
            w += np.pi
            e *= -1
        if k<0:
            k*=-1
            w += np.pi
            
        k1 = (G*n)**(1./3) * (1 - e**2)**(-1./2)
        mass_f = k/k1
        
        # Now make assumption that most of the mass is in the star
        mp = mass_f * self.mass**(2./3) # Planet mass in solar masses
        mp *= 9.543e-4
        p = 2*np.pi/n
        a = (G*self.mass/n**2)**(1./3)
        return mp,e,p,w,a
        
    def get_uncertainties(self,pcov,*popt):
        mp,e,p,w,a = self.recover_params(*popt)
        
        ms = self.mass
        ms_err = self.mass_err
        n = popt[0]
        k = popt[2]
        
        n_err = np.sqrt(pcov[0,0])
        k_err = np.sqrt(pcov[2,2])
        w_err = np.sqrt(pcov[3,3])
        p_err = p*np.abs(n_err/n)
        
        e_err = np.sqrt(pcov[4,4])
        a_err = a*np.sqrt( (2*n_err/(3*n))**2 + (ms_err/(3*ms))**2 )
        
        mp_err = mp*np.sqrt( (2*ms_err/(3*ms))**2 + (k_err/k)**2 + (n_err/(3*n))**2 + (e*e/(1-e*e))**2 *(e_err/e)**2)
        
        return mp_err,e_err,p_err,w_err,a_err
    def output_results(self,mp,mp_err,e,e_err,p,p_err,a,a_err):
    
        outstr = 'Mp > {:.3f} +\- {:.4f} MJ\n'.format(mp,mp_err)
        outstr += 'e = {:.3f} +\- {:.4f}\n'.format(e,e_err)
        outstr += 'P =  {:.3f} +\- {:.4f} days\n'.format(p,p_err)
        outstr += 'a = {:.3f} +\- {:.4f} AU\n'.format(a,a_err)
    
        print self.name
        print outstr

def solve_kep_eqn(l,e):
    """ Solve Keplers equation x - e*sin(x) = l for x"""
    try:
        l[0]
        res = np.zeros(l.shape)
        for i,li in enumerate(l):
            tmp,= fsolve(lambda x: x-e*np.sin(x) - li,li)
            res[i] = tmp
    except IndexError:
        res, = fsolve(lambda x: x - e*np.sin(x)-l,l)

    return res
        
    
        
def load_single_star(fname):
    with open(fname,'r') as f:
        header = f.readline()
    
    rv_data = np.loadtxt(fname)
    line = header.strip().split('#')[-1].strip()
    line = line.split('Mass')
    name = line[0]
    line = line[-1].split('Jitter')
    mass = line[0]
    jitter = float(line[-1])
    mass = float(line[0].split('(')[0])
    mass_err = float(line[0].split('(')[-1].split(')')[0])
    return Star(name,mass,rv_data,mass_err=mass_err,jitter=jitter)

def plot_data(star,t_fit,vr_fit,residuals,popt,pcov):
    t = star.t
    vr = star.vr
    vr_err = star.vr_err
    
    mp,e,p,w,a = star.recover_params(*popt)
    mp_err,e_err,p_err,w_err,a_err = star.get_uncertainties(pcov,*popt)
      
    legstr = 'Mp>{:.2f} +\- {:.2f} Mj, e={:.3f} +\- {:.2f}, P={:.1f} +\- {:.2f} days'.format(mp,mp_err,e,e_err,p,p_err)
        
    fig = plt.figure(figsize=(20,10))
    gs = gridspec.GridSpec(3,4)
    ax = fig.add_subplot(gs[:2,:])
    axr = fig.add_subplot(gs[-1,:])
    plt.subplots_adjust(hspace=0)


    ax.errorbar(t,vr,yerr=vr_err,fmt='o')
    ax.plot(t_fit,vr_fit,label=legstr)
    ax.set_xticklabels([]) #Remove x-tic labels for the first frame
    ax.get_yticklabels()[0].set_visible(False)
    ax.grid()
    ax.minorticks_on()
    ax.legend(loc='upper left')
    ax.set_ylabel('Radial Velocity (m/s)', fontsize=15)
    ax.set_title(star.name + ', M = {:.1f} Msun'.format(star.mass),fontsize=15)
            
    #Residual plot
    axr.errorbar(t,residuals,yerr=vr_err,fmt='o')
    axr.axhline(0,color='k',linewidth=2)
    axr.get_yticklabels()[-1].set_visible(False)
    axr.grid()
    axr.minorticks_on()
    axr.set_ylabel('Residuals',fontsize=15)
    axr.set_xlabel('Time (Julian Days)',fontsize=15)
    
    plt.show()