import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as sig
from scipy.optimize import curve_fit,fsolve
from rv import *




def fitting_function(t,n,tau,k,w,e):
    """Obtain the radial velocity due to a single planet.
        t = time of measurement,
        n = angular frequency of planet,
        tau = time of pericenter passage,
        k = amplitude of radial velocity (depends on planet mass and eccentricity),
        w = related to the argument of pericenter by a shift of pi.
        e = eccentricity of orbit
        The radial velocity at time t is given by
        vr = k*(cos(f + w)+e*cos(w)), 
        where f is related to the number of periods since pericenter passage, n*(t-tau)
    """
    

                
    e_anom = solve_kep_eqn(n*(t-tau),e)

    f = 2*np.arctan2(np.sqrt(1+e)*np.sin(e_anom*.5),np.sqrt(1-e)*np.cos(e_anom*.5))
   
    return k*(np.cos(f + w) + e*np.cos(w)) 



def fit_data(star,fitting_function):
    """ Fit the RV data with the given fitting function """
    
    # Set initial guesses to the parameters
    n0 = star.get_best_n()
          
    # Set the inital guess for the semi-amplitude equal to the maximum stellar velocity
    # Hint: The velocity data are stored in star.vr
    k0 = # FIX ME!
    tau0 = star.t[ star.vr == k0]
    w0 = 0
    e0 = 0.5
    initial_guess = (n0,tau0,k0,w0,e0)
        
    

    # Now do the fit
    # Pass the time and radial velocities to curve_fit
    # Hint: curve_fit is called with the arguments, curv_fit(func,x,y)
    # where func is the fitting_function, and x,y are the data arrays
    
    popt,pcov = curve_fit(fitting_function, #FIX ME!, #FIX ME!, 
                            sigma=star.vr_err,absolute_sigma=True, 
                            p0=initial_guess)
                            
    mp,e,p,w,a= star.recover_params(*popt)
    mp_err,e_err,p_err,w_err,a_err = star.get_uncertainties(pcov,*popt)
    t_fit = np.linspace(star.t[0],star.t[-1],1e3)
    vr_fit = np.array([fitting_function(x,*popt) for x in t_fit])
    residuals = star.vr - np.array([fitting_function(x,*popt) for x in star.t])
  
    # Output the results by calling the star.output_results
    star.output_results(mp,mp_err,e,e_err,p,p_err,a,a_err)                      

        
    # Finally plot the results by calling the plot data function
    star.periodogram()
    plot_data(star,t_fit,vr_fit,residuals,popt,pcov)
        
        
    return 


plt.close('all')
# Load the data file in here
# Pass the file name of the data file to the load_single_star function
star= load_single_star( #FIX ME! )

# You now have a star object
# Print out the star's name and mass
# Hint: star is a class defined in rv.py with attributes name and mass


# Plot the data using the star.plot function
# FIX ME! 
plt.show() 


# Now fit the data with the fit_data function defined above

#FIX ME!
plt.show()


# Now repeat these steps for the other stars.