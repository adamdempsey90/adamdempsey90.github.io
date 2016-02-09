import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
class DataPoints():
    def __init__(self,val,label,err=None):
        self.val = val
        self.err = err
        self.label = label



class Data(DataPoints):
    def __init__(self,Xdata=None, Ydata=None,from_file=False,filename=None,fitting_function = None):
        if from_file:
            self.load_from_file(filename)
        else:
            self.xdata = Xdata
            self.ydata = Ydata
        self._fitting_function = fitting_function
        
    @property
    def fitting_function(self):
        return self._fitting_function
    @fitting_function.setter
    def fitting_function(self,name):
        self._fitting_function = name
    @fitting_function.deleter
    def fitting_function(self):
        self._fitting_function = None
   
    def plot(self,ax=None,xlims=None,ylims=None,logx=False,logy=False,**kwargs):
        if ax == None:
            fig=plt.figure()
            ax = fig.add_subplot(111)
        fmt = kwargs.pop('fmt','o')
        fontsize = kwargs.pop('fontsize',15)
        if self.ydata.err != None:
            ax.errorbar(self.xdata.val,self.ydata.val,yerr=self.ydata.err,fmt=fmt,**kwargs)
        else:
            ax.plot(self.xdata.val,self.ydata.val,fmt,**kwargs)
        
        ax.set_xlabel(self.xdata.label,fontsize=fontsize)
        ax.set_ylabel(self.ydata.label,fontsize=fontsize)
        if logx:
            ax.set_xscale('log')
        if logy:
            ax.set_yscale('log')
        
        if xlims != None:
            ax.set_xlim(xlims)
        if ylims != None:
            ax.set_ylim(ylims)
        ax.minorticks_on()
        
            
    def load_from_file(self,filename):
        """ Load data from filename.
        The file should look like,
        # x y xerr yerr 
          0 1 0   0 
        """
        dat = np.loadtxt(filename)
        with open(filename,'r'):
            for line in f.readlines():
                if '#' in line:
                    header = line.split('#')[-1].split()
        no_xerr= False
        no_yerr= False
        if dat.shape[1] < 4:
            # Assume no xerr
            no_xerr = True
            if dat.shape[1] < 3:
                # Assume no yerr
                no_yerr = True
                if dat.shape[1] < 2:
                    # Not enough data points
                    print 'Not enough data points in {}!'.format(filename)
                    raise
                
        self.xdata = DataPoints(dat[:,0],header[0],err= None if no_xerr else dat[:,2])
        self.ydata = DataPoints(dat[:,1],header[1],err = None if no_yerr else dat[:,2] if no_xerr else dat[:,3])
        
        