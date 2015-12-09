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

# Linear fitting function
def linear_fit_function(x,H):
    y = speed_of_light * x / H
    return y  

#  Part 3: Fitting to a quadratic function

# Define our quadratic fitting function.

def quadratic_fit_function(x,H,q):
    y = speed_of_light * ( x + 0.5*(1 - q) * x*x)/H
    return y
def cubic_fit_function(x,H,q,a):
    y = speed_of_light * ( x + 0.5*(1 - (q+a*x)) * x*x)/H
    return y



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


