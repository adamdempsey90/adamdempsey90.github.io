import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt


speed_of_light = 3e5
def chi_squared(func,xdata,ydata,yerr,*parameters):
    return np.sum( ((func(xdata,*popt)-ydata)/yerr)**2)

plt.close('all')

"""
 The supernovae data are laid out in 3 columns in the sn_data.txt file.
 The first is the redshift to the supernova, z
 The second is the distance to the superova, d
 The third is the error in the distance, d_err
 The first few lines of the file look like this,
#z       d          d_err
0.01	51.5229    8.77905
0.01	43.6516	   7.63887
0.01	55.4626	   9.96117
0.013	63.6796	   9.09091
"""


# Import the data file with numpy's loadtxt function
# Replace the FIX ME with the filename of the data file.
# Remember to put it in quotes!

sn_data = np.loadtxt( #FIX ME!  )

# The variable sn_data is a Matrix with three columns
# Assign each column of the data file to an appropriate variable

z = sn_data[:,0]
d = sn_data[:,1]

# Set d_err to the third column of the sn_data matrix
# Look at how we set z and d in the lines above for a hint.

d_err = #FIX ME!!

# Linear fitting function
def linear_fit_function(x,H):
    """
    This is Hubble's law.
    x is the redshift.
    H is Hubbles constant in km/s/Mpc
    returns the distance
    """
    y = speed_of_light * x / H
    return y

#  Part 3: Fitting to a quadratic function

# Define our quadratic fitting function.

def quadratic_fit_function(x,H,q):
    """
    This is next order correction to Hubble's law that describes an
    accelerated expansion.
    x is the redshift.
    H is Hubbles constant in km/s/Mpc
    q is the deceleration parameters. q<0 implies that the expansion rate is
    accelerating, while q>0 implies that the expansion is decelerating.
    returns the distance
    """
    y = speed_of_light * ( x + 0.5*(1 - q) * x*x)/H
    return y


# Lets restrict ourselves to supernovae with z < 0.5, set zmax to reflect this
zmax = # FIX ME!

# Replot the data in the domain z < zmax
plt.figure()

# Plot distance (d) versus redshift (z) in the range z < zmax below
# Make sure you are only plotting supernovae with z < 0.5 (look back to the previous part for a hint)

plt.errorbar( #FIX ME!, #FIX ME!,yerr=d_err[z<=zmax],fmt='o')


plt.xlabel('Redshift') # Labels
plt.ylabel('Distance (Mpc)')
plt.title(' z < %.1f' % zmax)
plt.show()

# Fit the data to a quadratic function.
# Replace the fix me with the quadratic_fit_function name (look back to the previous part for a hint)

popt,pcov = opt.curve_fit( #FIX ME!, z[z<=zmax], d[z<=zmax],  sigma = d_err[z<=zmax],absolute_sigma = True)


# Calculate errors
perr = np.sqrt(np.diag(pcov))
yval = quadratic_fit_function(z[z<=zmax],*popt)
dydq =  speed_of_light*z[z<=zmax]/(2*popt[0])
dydh = yval/popt[0]
yerr = np.sqrt( dydh**2 * pcov[0,0] + dydq**2 * pcov[1,1] + 2*dydq*dydh*pcov[0,1])


# Print out the results
print 'Quadratic fit Hubble constant for z < %.2f = %.2f +- %.2f' % (zmax,popt[0],perr[0])
print 'Quadratic fit for the deceleration constant q = %.2f +- %.2f' % (popt[1],perr[1])
print 'The percentage of matter-energy in the Universe that is in Dark Energy is %.2f' % (0.5*0.3 - popt[1])
plt.plot(z[z<=zmax],quadratic_fit_function(z[z<=zmax],*popt),'-k',label='H=%.2f+-%.2f, q=%.2f+-%.2f'%(popt[0],perr[0],popt[1],perr[1]))
plt.fill_between(z[z<=zmax],quadratic_fit_function(z[z<=zmax],*popt)-yerr,quadratic_fit_function(z[z<=zmax],*popt)+yerr,alpha=1, edgecolor='#3F7F4C', facecolor='#7EFF99')
plt.legend(loc='upper left')
plt.draw()

popt,pcov = opt.curve_fit(linear_fit_function, z[z<=zmax], d[z<=zmax],  sigma = d_err[z<=zmax],absolute_sigma = True)
perr = np.sqrt(np.diag(pcov))
yerr = np.abs( linear_fit_function(z[z<=zmax],*popt)/popt[0] * perr[0])
plt.plot(z[z<=zmax],linear_fit_function(z[z<=zmax],*popt),'--k',label='H=%.2f+-%.2f'%(popt[0],perr[0]))
plt.fill_between(z[z<=zmax],linear_fit_function(z[z<=zmax],popt[0])-yerr,linear_fit_function(z[z<=zmax],popt[0])+yerr,alpha=1, edgecolor='#CC4F1B', facecolor='#FF9848')
plt.legend(loc='upper left')
plt.draw()

