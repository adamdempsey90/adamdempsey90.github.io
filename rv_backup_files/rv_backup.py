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
    #    ax.set_title(self.name + ',    M = %.2f Msun' % self.mass)
        ax.minorticks_on()
        ax.grid()
        

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
        return omega[max_pers[-1]],periods,res/res[max_pers[-1]]

    def vr_func(self,t,t0,n,tau,h,c,v0,d):
        """Obtain the radial velocity due to a single planet.
        t = time of measurement,
        n = period of planet,
        tau = time of pericenter passage,
            k = amplitude of radial velocity (depends on planet mass and eccentricity),
        w = related to the argument of pericenter by a shift of pi.
        d = D.C offset, 
        The radial velocity at time t is given by
        vr = k*(cos(f + w)+d), where f is related to the number of periods since pericenter passage, n*(t-tau)"""
        
        k = np.sqrt(h*h + c*c)
        w = np.arctan2(-c,h)
        e = np.abs(v0/(k*np.cos(w)))
        u = self.solve_kep_eqn(n*(t-tau),e)
        f = 2*np.arctan2(np.sqrt(1+e)*np.sin(u*.5),np.sqrt(1-e)*np.cos(u*.5))
    
        return h*np.cos(f) + c*np.sin(f) + v0 + d*(t-t0)
    
    #def double_fitting_func(self,t0,fit_trend=True):
    #        
    #        def vr_func2(t,n1,tau1,h1,c1
    def fitting_func(self,t0,fit_trend=False):
        if fit_trend:
            def vr_func_n(t,n,tau,h,c,e,d):
                """Obtain the radial velocity due to a single planet.
                t = time of measurement,
                n = period of planet,
                tau = time of pericenter passage,
                k = amplitude of radial velocity (depends on planet mass and eccentricity),
                w = related to the argument of pericenter by a shift of pi.
   	            d = D.C offset, 
                The radial velocity at time t is given by
                vr = k*(cos(f + w)+d), where f is related to the number of periods since pericenter passage, n*(t-tau)"""
               
                k = np.sqrt(h*h + c*c)
                w = np.arctan2(-c,h)
   #             e = np.abs(v0/(k*np.cos(w)))
                u = self.solve_kep_eqn(n*(t-tau),e)
                f = 2*np.arctan2(np.sqrt(1+e)*np.sin(u*.5),np.sqrt(1-e)*np.cos(u*.5))
        
                return h*np.cos(f) + c*np.sin(f) +  k*e*np.cos(w)  + d*(t-t0)
        else:
            def vr_func_n(t,n,tau,k,w,e):
                """Obtain the radial velocity due to a single planet.
                t = time of measurement,
                n = period of planet,
                tau = time of pericenter passage,
                k = amplitude of radial velocity (depends on planet mass and eccentricity),
                w = related to the argument of pericenter by a shift of pi.
   	            d = D.C offset, 
                The radial velocity at time t is given by
                vr = k*(cos(f + w)+d), where f is related to the number of periods since pericenter passage, n*(t-tau)"""
            
 #               k = np.sqrt(h*h + c*c)
  #              w = np.arctan2(-c,h)
            
 #               e = np.abs(v0/(k*np.cos(w)))
                norm = -1 if e<0 else 1
                wnorm = np.pi if e<0 else 0
                knorm = -1 if k<0 else 1
                wnorm += np.pi if k<0 else 0
                
                u = self.solve_kep_eqn(n*(t-tau),e*norm)

                f = 2*np.arctan2(np.sqrt(1+e*norm)*np.sin(u*.5),np.sqrt(1-norm*e)*np.cos(u*.5))
                return knorm*k*(np.cos(f + w+wnorm) + norm*e*np.cos(w+wnorm))      
 #               return h*np.cos(f) + c*np.sin(f) + k*e*np.cos(w)        
        return vr_func_n
    
    def recover_params(self,*popt):
        #"""Recover the masses of the planets"""
        #if k < 0:
        #    w += np.pi
        #    k *= -1
        n = popt[0]
        tau = popt[1]
        k = popt[2]
       # c = popt[3]
        e = popt[4]
        w = popt[3]
        
        if e<0:
            w += np.pi
            e *= -1
        if k<0:
            k*=-1
            w += np.pi
#        k = np.sqrt(h*h + c*c)
 #       w = np.arctan2(-c,h)
        if e<0:
            e *= -1
            w += np.pi
  #      if np.sign(np.sin(w)) != np.sign(-c):
        
#        e = np.abs(v0/(k*np.cos(w)))

        k1 = (G*n)**(1./3) * (1 - e**2)**(-1./2)
        mass_f = k/k1
        # Now make assumption that most of the mass is in the star
        mp = mass_f * self.mass**(2./3) # Planet mass in solar masses
        mp *= 9.543e-4
        p = 2*np.pi/n
        a = (G*self.mass/n**2)**(1./3)
        return mp,e,p,w,a
    def get_uncertainties(self,pcov,*popt):
        return 0,0,0,0,0
  #      mp,e,p,w,a = self.recover_params(*popt)
  #      
  #     
  #      n = popt[0]
  #      h = popt[2]
  #      c = popt[3]
  #      e = np.abs(popt[4])
  #      
  #      k = np.sqrt(h*h+c*c)
  #       
  #      snn = pcov[0,0]
  #      sn = np.sqrt(snn)
  #      
  #      shh = pcov[2,2]
  #      scc = pcov[3,3]
  #      svv = pcov[4,4]
  #      
  #      shc = pcov[2,3] + pcov[3,2]
  #      shv = pcov[2,4] + pcov[4,2]
  #      scv = pcov[3,4] + pcov[4,3]
  #      
  #      p_err = np.abs(p * sn/n)
  #      a_err = a * np.abs(2*sn/(3*n))
  #      
  #      k_err = (1/k)**2 *( h**2 * shh + c**2 * scc + c*h*shc)
  #      k_err =np.sqrt(k_err)
  #      w_err = np.cos(w)**4 *(scc/(h*h) + (c/h**2)**2 * shh - c*shc/(h**3))
  #      w_err = np.sqrt(w_err)
  #      
  ##      e_err = e*e *( svv/(v0*v0) + shh/(h*h) - shv/(v0*h))
  ##      e_err = np.sqrt(e_err)
  #      e_err = np.sqrt(svv)     
  #      mp_err = mp**2*((k_err/k)**2 + (e**2/(1-e**2))**2*(e_err/e)**2 + (snn/(9*n*n))
  #                      -(e**2/(1-e**2))*(k_err/k)*(e_err/e) - (sn*k_err/(3*n*k))
  #                      + (e**2/(1-e**2))*(sn*e_err)/(3*n*e))
  #      mp_err += (2*mp* self.mass_err/(3*self.mass))**2
  #      mp_err =np.sqrt(mp_err)
  #      
  #      #mp_err = (mp/(1-e**2))**2*(skk/k**2 - e*e*skd/(k*d)-e*e*np.tan(w)*swk/k \
  #      #           + e**4 * sdd/d**2 + e**4 * np.tan(w)*sdw/d + e**4 * np.tan(w)**2 * sww)
  #      #mp_err = np.sqrt(mp_err)
  #      #e_err = e**2 *(sdd/d*2 + skk/k**2 + np.tan(w)**2 * sww - skd/(d*k) + np.tan(w)*sdw/d - np.tan(w)*swk/k)
  #      #e_err = np.sqrt(e_err)
  #      
  #      return mp_err,e_err,p_err,w_err,a_err
        
    def chi2(self,popt,pcov):
        return np.sum(((self.vr - self.func(self.t,*popt))/self.vr_err)**2)
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
        
        
    def plot_data(self,popt,pcov,pers,periodogram_res,chi2):
        t = self.t
        vr = self.vr
        vr_err = self.vr_err
        t_fit = np.linspace(t[0],t[-1],1e3)
        vr_fit = np.array([self.func(x,*popt) for x in t_fit])
        resids = vr - np.array([self.func(x,*popt) for x in t])
        mp,e,p,w,a = self.recover_params(*popt)
        mp_err,e_err,p_err,w_err,a_err = self.get_uncertainties(pcov,*popt)
        legstr = 'Mp>{:.2f} +\- {:.2f} Mj, e={:.3f} +\- {:.2f}, P={:.1f} +\- {:.2f} days'.format(mp,mp_err,e,e_err,p,p_err)
        
        fig = plt.figure(figsize=(20,10))
        gs = gridspec.GridSpec(3,4)
        ax = fig.add_subplot(gs[:2,:2])
        axr = fig.add_subplot(gs[-1,:2])
        axp = fig.add_subplot(gs[:2,-2:])
        axpr = fig.add_subplot(gs[-1,-2:])
        plt.subplots_adjust(hspace=0)


        ax.errorbar(t,vr,yerr=vr_err,fmt='o')
        ax.plot(t_fit,vr_fit,label=legstr)
        ax.set_xticklabels([]) #Remove x-tic labels for the first frame
        ax.get_yticklabels()[0].set_visible(False)
        ax.grid()
        ax.minorticks_on()
        ax.legend(loc='upper left')
        ax.set_ylabel('Radial Velocity (m/s)', fontsize=15)
        ax.set_title(self.name + ', M = {:.1f} Msun, chi2 = {:.3f}'.format(self.mass,chi2),fontsize=15)
        
        
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
        axp.axvline(p,linestyle='--',color='k',linewidth=3)
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
        
    def fit_data(self,fit_trend=False):
        n,pers,periodogram_res = self.get_best_n(self.t,self.vr)
        
        func = self.fitting_func(self.t[0],fit_trend)
        if fit_trend:
            init_guess = (n,self.t[self.vr == self.vr.max()],self.vr.max(),.2,0,0)
        else:
            init_guess = (2*np.pi/1043,17062,31.5,198*np.pi/180,.11)
         #   init_guess = (n,self.t[self.vr == self.vr.max()],self.vr.max(),-self.vr.max(),.11)
        
        self.func = func
        popt,pcov = curve_fit(func,self.t,self.vr, 
                                sigma=self.vr_err,absolute_sigma=True, 
                                p0=init_guess)
                                
        print 'Results'
        print popt
        print np.sqrt(np.diag(popt))
        print init_guess
        chi2 = self.chi2(popt,pcov)
        
        popterr = np.sqrt(np.diag(pcov))
        popt_str_vals = []
        for p1,p2 in zip(popt,popterr):
            popt_str_vals.append(p1)
            popt_str_vals.append(p2)
        if not fit_trend:
            popt_str_vals.append(0)
            popt_str_vals.append(0)

      #  print 'Error in cosine fit:'
      #  print 'tau=%.2e pm %.2e\nK = %.2e pm %.2e\nomega = %.2e pm %.2e\nd = %.2e pm %.2e\tslope = %.2e pm %.2e'% tuple(popt_str_vals)
        #out_str2 = 'n = {:.3f} +\- {:.3f}\n'.format(popt[0],popterr[0])
        #out_str2 += 'tau = {:.3f} +\- {:.3f}\n'.format(popt[1],popterr[1])  
        #out_str2 += 'h = {:.3f} +\- {:.3f}\n'.format(popt[2],popterr[2])  
        #out_str2 += 'c = {:.3f} +\- {:.3f}\n'.format(popt[3],popterr[3])  
        #out_str2 += 'k = {:.3f} +\- {:.3f}\n'.format(np.sqrt(popt[2]**2 + popt[3]**2),(popt[3]*popterr[3]+popt[2]*popterr[2])/np.sqrt(popt[2]**2 + popt[3]**2))  
        #out_str2 += 'v0 = {:.3f} +\- {:.3f}\n'.format(popt[4],popterr[4]) 
        #if fit_trend:
        #    out_str2 += 'd = {:.3f} +\- {:.3f}\n'.format(popt[5],popterr[5]) 
        #print out_str2
        mp,e,p,w,a = self.recover_params(*popt)
        mp_err,e_err,p_err,w_err,a_err = self.get_uncertainties(pcov,*popt)        
        
        print '\n\n',self.name
        out_str = 'Mp: >{:.3f} +\- {:.3f} Mj\n'.format(mp,mp_err)
        out_str += 'e: {:.3f} +\- {:.3f}\n'.format(e,e_err)
        out_str += 'P: {:.3f} +\- {:.3f} days\n'.format(p,p_err)
        out_str += 'a: {:.3f} +\- {:.3f} AU\n'.format(a,a_err)
        out_str += 'w: {:.3f} +\-{:.3f} deg\n'.format(w*180/np.pi,w_err*180/np.pi)
        print out_str
        print 'chi^2 = {:.3f}'.format(chi2)
        
        self.plot_data(popt,pcov,pers,periodogram_res,chi2)
        
        
        return popt,init_guess
def load_single(fname):
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