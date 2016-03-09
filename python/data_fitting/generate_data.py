import numpy as np
from scipy.optimize import curve_fit
def fitting_function(x,a,b,c):
    return np.exp(-x/a) *(1 - (x/b)**2 + (x/c)**3)

def generate_data(func,xvals,*popt,**kargs):
    noise = kargs.pop('noise',.3)
    err = kargs.pop('err',.6)
    print 'Generating %d values with noise=%.2f and error noise=%.2f'%(len(xvals),noise,err)
    ytrue = func(xvals,*popt)
    #y  = ytrue*(1 + np.random.poisson(1.,len(ytrue)))
    y = np.random.normal(ytrue,noise)
    yerr = np.random.normal(0,err,len(y))
    return y,yerr

def get_sample_data(n=50):
    x = np.linspace(0,10,n)
    popt=(1.,.5,.7)
    y,yerr = generate_data(fitting_function,x,*popt)
    return x,y,yerr,popt

def fit_data(x,y,yerr):
    popt,pcov = curve_fit(fitting_function,x,y,sigma=yerr,absolute_sigma=True)
    perr = np.sqrt(np.diag(pcov))
    for i,(v,e) in enumerate(zip(popt,perr)):
        print 'Param %d: %f +- %f' % (i,v,e)
    return popt,perr
