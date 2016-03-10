import numpy as np
from scipy.optimize import curve_fit
import data_fits as df

def fitting_function_1(x,a,b,c):
    return np.exp(-x/a) *(1 - (x/b)**2 + (x/c)**3)

def fitting_function_2(x,a,b,c,d):
    return  x**2 * np.sin(b*x + c) + d

def fitting_function_3(x,a,b,c,d):
    return a*np.sin(b*x + c) + d

def fitting_function_4(x,a,b,c,d):
    return a*np.sin(b*x) + c*np.cos(d*x)


def generate_data(func,xvals,*popt,**kargs):
    noise = kargs.pop('noise',.3)
    err = kargs.pop('err',.6)
    print 'Generating %d values with noise=%.2f and error noise=%.2f'%(len(xvals),noise,err)
    ytrue = func(xvals,*popt)

    y = np.random.normal(ytrue,noise)
    yerr = np.random.normal(0,err,len(y))
    xerr = np.random.normal(0,err/2,len(y))
    return y,yerr,xerr

def get_sample_data(n=50,num=1):
    if num == 1:
        x = np.linspace(0,10,n)
        popt=(1.,.5,.7)
        y,yerr,xerr = generate_data(fitting_function_1,x,*popt)
    if num == 2:
        x = np.linspace(0,2*np.pi,n)
        popt=(.4,2.1,.7,0)
        y,yerr,xerr = generate_data(fitting_function_2,x,*popt)
    if num == 3:
        x = np.linspace(0,2*np.pi,n)
        popt=(1.,.5,.7,3)
        y,yerr,xerr = generate_data(fitting_function_3,x,*popt)
    if num == 4:
        x = np.linspace(0,5*np.pi,n)
        popt=(.5,5,-1, 3)
        y,yerr,xerr = generate_data(fitting_function_4,x,*popt,err=.1,noise=.1)
    return x,y,yerr,xerr,popt

def write_data(fname="data.dat",num=1,n=50):
    x,y,yerr,xerr,_ = get_sample_data(num=num,n=n)
    
    dat = np.vstack((x,y,yerr)).transpose()
    if dat.shape[1] > 3:
        head = '#x\ty\txerr\tyerr\n'
    else:
        head = '#x\ty\tyerr\n'
        
    with open(fname,'w') as f:
        f.write(head)
        f.write('\n'.join(['\t'.join(line.astype(str)) for line in dat])+'\n')


def fit_data(x,y,yerr):
    popt,pcov = curve_fit(fitting_function,x,y,sigma=yerr,absolute_sigma=True)
    perr = np.sqrt(np.diag(pcov))
    for i,(v,e) in enumerate(zip(popt,perr)):
        print 'Param %d: %f +- %f' % (i,v,e)
    return popt,perr
