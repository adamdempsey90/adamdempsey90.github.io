import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

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

    def statistics(self):
        print 'Data Statistics:'
        print'\tDomain: ({:.3e}, {:.3e})'.format(self.xdata.val.min(),self.xdata.val.max())
        print'\tRange: ({:.3e}, {:.3e})'.format(self.ydata.val.min(),self.ydata.val.max())
        print'\tMean {}: {:.3e}'.format(self.ydata.label,self.ydata.val.mean())
        print'\tStandard deviation {}: {:.3e}'.format(self.ydata.label,self.ydata.val.std())
        print'\tVariance {}: {:.3e}'.format(self.ydata.label,self.ydata.val.var())
        print'\tMedian {}: {:.3e}'.format(self.ydata.label,np.median(self.ydata.val))

    def plot(self,ax=None,xlims=None,ylims=None,logx=False,logy=False,**kwargs):
        if ax == None:
            fig=plt.figure()
            ax = fig.add_subplot(111)
        fmt = kwargs.pop('fmt','o')
        fontsize = kwargs.pop('fontsize',15)
        grid_off = kwargs.pop('grid_off',False)
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
        if not grid_off:
            ax.grid(which='both')


    def load_from_file(self,filename):
        """ Load data from filename.
        The file should look like,
        # x y xerr yerr
          0 1 0   0
        """
        dat = np.loadtxt(filename)
        with open(filename,'r') as f:
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

    def fit_data(self,initial_guess=None,fitting_function=None,**kwargs):
        if self.fitting_function == None:
            if fitting_function == None:
                print 'Fitting function not set!'
                return
            else:
                self.fitting_function = fitting_function
        abssigma = kwargs.pop('absolute_sigma',True)
        try:
            if self.ydata.err == None:
                popt,pcov = opt.curve_fit(self.fitting_function,
                                    self.xdata.val,self.ydata.val)
            else:
                popt,pcov = opt.curve_fit(self.fitting_function,
                                    self.xdata.val,self.ydata.val,
                                    sigma=self.ydata.err,absolute_sigma=abssigma,
                                    p0 = initial_guess)
        except RuntimeError:
           print 'Fit not converging. Try different inital_guess.'
           return


        resids = self.ydata.val - self.fitting_function(self.xdata.val,*popt)
        self.plot_fit(popt,pcov,resids)
        self.print_fit(popt,pcov)
        return popt,pcov,resids
    def print_fit(self,popt,pcov):
        for i,(v,e) in enumerate(zip(popt,np.sqrt(np.diag(pcov)))):
            print 'Parameter {:d}: {:.3e} +- {:.3e}'.format(i,v,e)

    def plot_fit(self,popt,pcov,resids,**kwargs):
        if self.fitting_function == None:
            print "Please define a fitting function first!"
            return
        fig = plt.figure(figsize=(20,10))
        gs = gridspec.GridSpec(3,3)
        ax = fig.add_subplot(gs[:2,:])
        axr = fig.add_subplot(gs[-1,:])
        plt.subplots_adjust(hspace=0)

        fontsize = kwargs.pop('fontsize',20)
        logx = kwargs.pop('logx',False)
        logy = kwargs.pop('logy',False)
        fmt = kwargs.pop('fmt','o')
        linewidth = kwargs.pop('linewidth',3)
        linestyle = kwargs.pop('linestyle','-')
        color = kwargs.pop('color','g')
        grid_off = kwargs.pop('grid_off',False)

        x_fit = np.linspace(self.xdata.val[0],self.xdata.val[-1],1000)
        y_fit = self.fitting_function(x_fit,*popt)
       # self.plot(ax=ax)
        ax.errorbar(self.xdata.val,self.ydata.val,yerr=self.ydata.err,fmt=fmt)
        ax.plot(x_fit,y_fit,color=color,linestyle=linestyle,linewidth=linewidth,**kwargs)
        axr.errorbar(self.xdata.val,resids,yerr=self.ydata.err,fmt=fmt)
        axr.axhline(0,color='k')
        axr.set_xlabel(self.xdata.label,fontsize=fontsize)
        axr.set_ylabel('Residuals',fontsize=fontsize)
        ax.set_ylabel(self.ydata.label,fontsize=fontsize)

        ax.set_xticklabels([]) #Remove x-tic labels for the first frame
        ax.get_yticklabels()[0].set_visible(False)
        axr.get_yticklabels()[-1].set_visible(False)

        if logx:
            ax.set_xscale('log')
            axr.set_xscale('log')
        if logy:
            ax.set_yscale('log')

        axr.minorticks_on()
        ax.minorticks_on()
        if not grid_off:
            axr.grid(which='both')
            ax.grid(which='both')
        fig.canvas.draw()

