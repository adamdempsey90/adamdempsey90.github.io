import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt


speed_of_light = 3e5

plt.close('all') # Make sure all open plots are closed


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


## Part 1: A simple linear fit

# Plot redshfit versus distance
# Use the plotting function errorbar(x,y,yerr=yerr,fmt='o') to plot the
# errorbars in the distances.

plt.figure() # Bring up a new figure window


# Plot the variable z on the x-axis and the variable d on the y-axis below.
# hint to plot y vs x using the errorbar function you would write plt.errorbar(x,y)

plt.errorbar(#FIX ME!, #FIX ME!, yerr=d_err,fmt='o')


plt.xlabel('Redshift')   # axes labels
plt.ylabel('Distance (Mpc)')
plt.title('Full Data')
plt.show()

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

# Fit the z,d data with a linear fit

# Call curve_fit to fit the data to linear_fit_function
# Replace the FIX MEs with the redshift (z) and distance (d) variables

popt,pcov = opt.curve_fit(linear_fit_function, #FIX ME!, #FIX ME!,  sigma = d_err,absolute_sigma = True)


perr = np.sqrt(np.diag(pcov)) # The error in the fitted parameters
yerr = np.abs(linear_fit_function(z,*popt)*perr[0]/popt[0])

# Print out the results
print 'Linear fit Hubble constant = %.1f +- %.1f' % (popt[0],perr[0])

plt.plot(z,linear_fit_function(z,*popt),'-k',label='H=%.2f +- %.2f' % (popt[0],perr[0]))

# Fancy shading showing the errors
plt.fill_between(z,linear_fit_function(z,popt[0])-yerr,linear_fit_function(z,popt[0])+yerr,alpha=1, edgecolor='#3F7F4C', facecolor='#7EFF99')
plt.legend(loc='upper left')
plt.draw()



# Part 2: limiting our fitting domain
#Lets fit only the data with redshift < 0.1


# UNCOMMENT ALL OF THE LINES BELOW THIS AFTER YOU FINISH PART 1

#zmax = #FIX ME!   # Set the maximum redshift zmax equal to 0.1
#plt.figure()

## Replace the FIX ME with the new redshift upper bound zmax.
#plt.errorbar( z[z <= #FIX ME!], d[z <= #FIX ME! ],yerr=d_err[z <= #FIX ME!],fmt='o')

#plt.xlabel('Redshift')
#plt.ylabel('Distance (Mpc)')
#plt.title(' z < 0.1')
#plt.show()

#popt,pcov = opt.curve_fit(linear_fit_function, z[z<=zmax], d[z<=zmax], sigma=d_err[z<=zmax],absolute_sigma = True)
#perr = np.sqrt(np.diag(pcov))
#yerr = np.abs(linear_fit_function(z[z<=zmax],*popt)*perr[0]/popt[0])
#print 'Linear fit Hubble constant for z < %.1f = %.1f +- %.1f' % (zmax,popt[0],perr[0])
#plt.plot(z[z<=zmax],linear_fit_function(z[z<=zmax],*popt),'-k',label='H=%.2f +- %.2f' % (popt[0],perr[0]))
#plt.fill_between(z[z<=zmax],linear_fit_function(z[z<=zmax],popt[0])-yerr,linear_fit_function(z[z<=zmax],popt[0])+yerr,alpha=1, edgecolor='#3F7F4C', facecolor='#7EFF99')
#plt.legend(loc='upper left')
#plt.draw()



