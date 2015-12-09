import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt


speed_of_light = 3e5 
def chi_squared(func,xdata,ydata,yerr,*parameters):
    return np.sum( ((func(xdata,*popt)-ydata)/yerr)**2)
    
plt.close('all')


# The supernovae data are laid out in 3 columns
# The first is the redshift to the supernova, z
# The second is the distance to the superova, d
# The third is the error in the distance, d_err
# The first few lines of the file look like this,
##z       L          L_err 
#0.01	51.5229    8.77905
#0.01	43.6516	   7.63887
#0.01	55.4626	   9.96117
#0.013	63.6796	   9.09091


# Import the data file with numpy's loadtxt function

sn_data = np.loadtxt('sn_data.txt')

# The variable sn_data is a Matrix with three columns
# Assign each column of the data file to an appropriate variable
z = sn_data[:,0]
d = sn_data[:,1]
d_err = sn_data[:,2]


## Part 1: A simple linear fit

# Plot redshfit versus distance
# Use the plotting function errorbar(x,y,yerr=yerr,fmt='o') to plot the 
# errorbars in the distances.

plt.figure()
plt.errorbar(z,d,yerr=d_err,fmt='o')
plt.xlabel('Redshift')
plt.ylabel('Distance (Mpc)')


# Linear fitting function
def linear_fit_function(x,H):
    y = speed_of_light * x / H
    return y  

# Fit the z,d data with a linear fit

popt,pcov = opt.curve_fit(linear_fit_function, z, d,  sigma = d_err,absolute_sigma = True)
perr = np.sqrt(np.diag(pcov))
print 'Linear fit Hubble constant = %.1f +- %.1f' % (popt[0],perr[0])
plt.plot(z,linear_fit_function(z,*popt),'-k')
plt.draw()



# Part 2: limiting our fitting domain
#Lets fit only the data with redshift < 0.1 

# Lets replot the data in this range of redshifts

plt.figure()
plt.errorbar( z[z<=0.1], d[z<=0.1],yerr=d_err[z<=0.1],fmt='o')
plt.xlabel('Redshift')
plt.ylabel('Distance (Mpc)')
plt.title(' z < 0.1')
plt.show()

popt,pcov = opt.curve_fit(linear_fit_function, z[z<=0.1], d[z<=0.1],  sigma = d_err[z<=0.1],absolute_sigma = True)
perr = np.sqrt(np.diag(pcov))
print 'Linear fit Hubble constant for z < 0.1 = %.1f +- %.1f' % (popt[0],perr[0])
plt.plot(z[z<=0.1],linear_fit_function(z[z<=0.1],*popt),'-k')
plt.draw()


#  Part 3: Fitting to a quadratic function

# Define our quadratic fitting function.

def quadratic_fit_function(x,H,q):
    y = speed_of_light * ( x + 0.5*(1 - q) * x*x)/H
    return y
def cubic_fit_function(x,H,q,a):
    y = speed_of_light * ( x + 0.5*(1 - (q+a*x)) * x*x)/H
    return y


# Lets restrict ourselves to supernovae with z < 0.5
zmax = 12
# Replot the data in the domain z < zmax
plt.figure()
plt.errorbar( z[z<=zmax], d[z<=zmax],yerr=d_err[z<=zmax],fmt='o')
plt.xlabel('Redshift')
plt.ylabel('Distance (Mpc)')
plt.title(' z < %.1f' % zmax)
plt.show()

popt,pcov = opt.curve_fit(quadratic_fit_function, z[z<=zmax], d[z<=zmax],  sigma = d_err[z<=zmax],absolute_sigma = True)
perr = np.sqrt(np.diag(pcov))
print 'Quadratic fit Hubble constant for z < %.2f = %.2f +- %.2f' % (zmax,popt[0],perr[0])
print 'Quadratic fit for the deceleration constant q = %.2f +- %.2f' % (popt[1],perr[1])
print 'The percentage of matter-energy in the Universe that is in Dark Energy is %.2f' % (0.5*0.3 - popt[1])
plt.plot(z[z<=zmax],quadratic_fit_function(z[z<=zmax],*popt),'-k')
plt.draw()

popt,pcov = opt.curve_fit(linear_fit_function, z[z<=zmax], d[z<=zmax],  sigma = d_err[z<=zmax],absolute_sigma = True)
perr = np.sqrt(np.diag(pcov))
plt.plot(z[z<=zmax],linear_fit_function(z[z<=zmax],*popt),'--k')
plt.draw()

popt,pcov = opt.curve_fit(cubic_fit_function, z[z<=zmax], d[z<=zmax],  sigma = d_err[z<=zmax],absolute_sigma = True)
perr = np.sqrt(np.diag(pcov))
print 'Cubic fit Hubble constant for z < %.2f = %.2f +- %.2f' % (zmax,popt[0],perr[0])
print 'Cubic fit for the deceleration constant q = %.2f +- %.2f' % (popt[1],perr[1])
print 'The percentage of matter-energy in the Universe that is in Dark Energy is %.2f' % (0.5*0.3 - popt[1])
plt.plot(z[z<=zmax],cubic_fit_function(z[z<=zmax],*popt),'-r')
plt.draw()

# Let's find out subset of the data is best fit by our quadratic fitting function
# To do this try out different zmax values and make a plot of zmax versus chi_squared

zmax_list = [ 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 
                0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
chi2_list = []
popt_list = []
perr_list = []
for zmax in zmax_list:
    popt,pcov = opt.curve_fit(quadratic_fit_function, z[z<=zmax], d[z<=zmax],  sigma = d_err[z<=zmax],absolute_sigma = True)
    perr = np.sqrt(np.diag(pcov))
    chi2 = chi_squared(quadratic_fit_function, z[z<=zmax], d[z<=zmax],d_err[z<=zmax], *popt)
    chi2_list.append(chi2)
    perr_list.append(perr)
    popt_list.append(popt)

for zmax in zmax_list:
    popt,pcov = opt.curve_fit(linear_fit_function, z[z<=zmax], d[z<=zmax],  sigma = d_err[z<=zmax],absolute_sigma = True)
    perr = np.sqrt(np.diag(pcov))
    chi2 = chi_squared(linear_fit_function, z[z<=zmax], d[z<=zmax],d_err[z<=zmax], *popt)
    chi2_list.append(chi2)
    perr_list.append(perr)
    popt_list.append(popt)            
# Make a plot of chi_squared versus zmax

plt.figure()
plt.plot(zmax_list,chi2_list[:len(zmax_list)],'-x',label='Quadratic Fit')
plt.plot(zmax_list,chi2_list[-len(zmax_list):],'-s',label='Linear Fit')
plt.xlabel('zmax')
plt.ylabel('Chi^2')
plt.legend(loc='best')
plt.show()

# Make a plot of the deceleration parameter versus zmax
plt.figure()
plt.errorbar(zmax_list,[x[1] for x in popt_list[:len(zmax_list)]],yerr=[x[1] for x in perr_list[:len(zmax_list)]],fmt='o')
plt.xlabel('zmax')
plt.ylabel('q')
plt.ylim(-1,1)
plt.show()

