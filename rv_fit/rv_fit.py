import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import scipy.signal as sig
from scipy.optimize import curve_fit,fsolve
from copy import copy

G = 4*np.pi**2 /(356.25)**2 # G in AU^3/(day^2 Solar Mass)

class Star():
    def __init__(self,name,mass,rv_data):
        self.name = name
        self.mass = mass
        self.t = rv_data[:,0]
        self.vr = rv_data[:,1]
        self.vr_err = rv_data[:,2]
        
    def plot(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.errorbar(self.t,self.vr,yerr=self.vr_err,fmt='o')
        ax.set_xlabel('Time  (Julian Days',fontsize=15)
        ax.set_ylabel('$v_r$',fontsize=20)
        ax.set_title(self.name + ',    M = %.2f Msun' % self.mass)
        plt.show()

    def solve_kep_eqn(self,l,e):
        """ Solve Keplers equation u - e*sin(u) = l for u"""
        try:
            l[0]
            res = np.zeros(l.shape)
            for i,li in enumerate(l):
                tmp,= fsolve(lambda x: x-e*np.sin(x) - li,li)
                res[i] = tmp
        except IndexError:
            res, = fsolve(lambda x: x - e*np.sin(x)-l,l)
   
        return res

    def get_best_n(self,t,vr):
        periods = 10**np.linspace(-3,4,1000)
        omega = 2*np.pi/periods
        res=sig.lombscargle(t-t[0],vr,omega)
        max_pers = np.argsort(res)
        return omega[max_pers[-1]],periods,res

    def vr_func(self,t,n,tau,k,w,d):
        """Obtain the radial velocity due to a single planet.
        t = time of measurement,
        n = period of planet,
        tau = time of pericenter passage,
            k = amplitude of radial velocity (depends on planet mass and eccentricity),
        w = related to the argument of pericenter by a shift of pi.
        d = D.C offset, 
        The radial velocity at time t is given by
        vr = k*(cos(f + w)+d), where f is related to the number of periods since pericenter passage, n*(t-tau)"""
        e = d/np.cos(w)
        u = self.solve_kep_eqn(n*(t-tau),e)
        f = 2*np.arctan2(np.sqrt(1+e)*np.sin(u*.5),np.sqrt(1-e)*np.cos(u*.5))
    
        return k*(np.cos(f+w)+d)
    
    def fitting_func(self,n):
        def vr_func_n(t,tau,k,w,d):
            """Obtain the radial velocity due to a single planet.
            t = time of measurement,
            n = period of planet,
            tau = time of pericenter passage,
            k = amplitude of radial velocity (depends on planet mass and eccentricity),
            w = related to the argument of pericenter by a shift of pi.
            d = D.C offset, 
            The radial velocity at time t is given by
            vr = k*(cos(f + w)+d), where f is related to the number of periods since pericenter passage, n*(t-tau)"""
            
            e = d/np.cos(w)
            u = self.solve_kep_eqn(n*(t-tau),e)
            f = 2*np.arctan2(np.sqrt(1+e)*np.sin(u*.5),np.sqrt(1-e)*np.cos(u*.5))
    
            return k*(np.cos(f+w)+d)
        return vr_func_n
    
    def recover_params(self,n,tau,k,w,d):
        """Recover the masses of the planets"""
        if k < 0:
            w += np.pi
            k *= -1
        e = d/np.cos(w)
        k1 = (G*n)**(1./3) * (1 - e**2)**(-1./2)
        mass_f = k/k1
        # Now make assumption that most of the mass is in the star
        mp = mass_f * self.mass**(2./3) # Planet mass in solar masses
        mp *= 9.543e-4
        p = 2*np.pi/n
        a = (G*self.mass/n**2)**(1./3)
        return mp,e,p,a
    def periodogram(self):
        n,pers,power=self.get_best_n(self.t,self.vr)
        self.periodogram_plot(pers,power)
        plt.show()
    def periodogram_plot(self,periods,power,pmax=None,ax=None):
        if ax == None:
            fig=plt.figure()
            ax=fig.add_subplot(111)
        ax.semilogx(periods,power)
        ax.set_xlabel('Period (days)',fontsize=15)
        ax.set_ylabel('Power',fontsize=20)
        if pmax != None:
            ax.axvline(pmax,color='k',linestyle='--',linewidth=3)
        
        
    def plot_data(self,n,popt,pcov,pers,periodogram_res):
        t = self.t
        vr = self.vr
        vr_err = self.vr_err
        t_fit = np.linspace(t[0],t[-1],1e3)
        vr_fit = np.array([self.vr_func(x,n,*popt) for x in t_fit])
        resids = vr - np.array([self.vr_func(x,n,*popt) for x in t])
            
        fig = plt.figure(figsize=(20,10))
        gs = gridspec.GridSpec(3,4)
        ax = fig.add_subplot(gs[:2,:2])
        axr = fig.add_subplot(gs[-1,:2])
        axp = fig.add_subplot(gs[:2,-2:])
        axpr = fig.add_subplot(gs[-1,-2:])
        plt.subplots_adjust(hspace=0)


        ax.errorbar(t,vr,yerr=vr_err,fmt='o')
        ax.plot(t_fit,vr_fit)
        ax.set_xticklabels([]) #Remove x-tic labels for the first frame
        ax.get_yticklabels()[0].set_visible(False)
        ax.grid()
        ax.minorticks_on()
        ax.set_ylabel('$v_r$', fontsize=20)
        
        
    #Residual plot
#        axr=fig.add_axes((.1,.1,.8,.2))        
        axr.errorbar(t,resids,yerr=vr_err,fmt='o')
        axr.axhline(0,color='k',linewidth=2)
        axr.get_yticklabels()[-1].set_visible(False)
    #axr.plot(t_lin,ff(t_lin,*popt2))
        axr.grid()
        axr.minorticks_on()
        axr.set_ylabel('Residuals',fontsize=15)
        axr.set_xlabel('Time (Julian Days)',fontsize=15)
#        axp=fig.add_axes(((.91,.3,.6,.6)))
        axp.semilogx(pers,periodogram_res)
        axp.axvline(2*np.pi/n,linestyle='--',color='k',linewidth=3)
        axp.grid()
        axp.set_xticklabels([])
        axp.get_yticklabels()[0].set_visible(False)
        axp.minorticks_on()
        axp.set_title('Periodogram',fontsize=15)
        
        n_r,pers_r,periodo_r = self.get_best_n(t,resids)
#        axpr=fig.add_axes(((.91,.1,.6,.2)))
        axpr.semilogx(pers_r,periodo_r)
        axpr.axvline(2*np.pi/n_r,linestyle='--',color='k',linewidth=3)
        axpr.get_yticklabels()[-1].set_visible(False)
        axpr.grid()
        axpr.minorticks_on()
        axpr.set_xlabel('Period (Days)',fontsize=15)
#     plt.errorbar(t,vr,yerr=vr_err,fmt='o')
#     plt.plot(t_fit,vr_fit)

    def fit_data(self):
        n,pers,periodogram_res = self.get_best_n(self.t,self.vr)
        popt,pcov = curve_fit(self.fitting_func(n),self.t,self.vr, 
                                sigma=self.vr_err,absolute_sigma=True, 
                                p0=(self.t[len(self.t)/2],self.vr.max(),0,.1))
        popterr = np.sqrt(np.diag(pcov))
        popt_str_vals = []
        for p1,p2 in zip(popt,popterr):
            popt_str_vals.append(p1)
            popt_str_vals.append(p2)
        print 'Error in cosine fit:'
        print 'tau=%.2e pm %.2e\nK = %.2e pm %.2e\nomega = %.2e pm %.2e\nd = %.2e pm %.2e'% tuple(popt_str_vals)
        mp,e,p,a = self.recover_params(n,*popt)
        print '\n\n',self.name
        print 'Mp: >%.3f Mj\ne: %.3f\nPeriod: %.3f days\na: %.3f AU'  %(mp,e,p,a)
        self.plot_data(n,popt,pcov,pers,periodogram_res)
        return
    
def load_data(fname):
    with open(fname,'r') as f:
        lines = [line for line in f.readlines() if '#' not in line]
    data_lines = [line.strip().split('\t') for line in lines if 'HD' in line]
    rv_data=np.zeros((1,3))
    rv_star=[]
    stellar_data={}
    for line in data_lines:
        if len(line) > 4:
            # Stellar params
            stellar_data[line[0].strip()] = float(line[2].strip())        
        else:
            rv_star.append(line[0].strip())
            temp = np.array([float(lin.strip()) for lin in line[1:]])

            rv_data = np.vstack((rv_data,temp))

    rv_data = rv_data[1:,:]
    unique_stars = np.unique(rv_star)
    rv_star = np.array(rv_star)
    print 'Found data for stars: '
    print unique_stars
    
    stars=[]
    for name in unique_stars:
        stars.append(Star(name,stellar_data[name],rv_data[rv_star == name,:]))
    
    return stars