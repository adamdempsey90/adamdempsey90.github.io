import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import copy

class Simulation(object):
    one_ninth = 1./9
    one_thirtysix = 1./36
    four_ninth = 4*one_ninth
    def __init__(self,lx=200,ly=80,viscosity=0.02,velocity=0.1,inflow=False,mask=None,cmap='bwr'):
        self.lx = lx
        self.nx = lx
        self.ly = ly
        self.ny = ly
        self.viscosity = viscosity
        self.velocity = velocity
        self.time = 0
        self.inflow = inflow
        if self.inflow:
            self.bc_vel = velocity*2
        else:
            self.bc_vel = velocity
        self.omega = 1./(3*viscosity + 0.5)
        
        self.init_barrier(mask)
        self.init_fluid()

        self.fig,self.axes,self.imgs = self.init_plot(cmap=cmap)
    def init_fluid(self):
        vel = self.velocity*np.ones((self.ny,self.nx))
      
        #for j in range(self.ny):
        #    for i in range(50):
        #        vel[j,i] = 2*self.velocity
        
        temp = (np.ones((self.ny,self.nx)) - 1.5*vel**2)
        self.n_M = (self.four_ninth)*temp
        self.n_N = (self.one_ninth) * temp
        self.n_S = (self.one_ninth) * temp
        self.n_E =  (self.one_ninth)*( temp + 3*vel+ 4.5*vel**2)
        self.n_W =  (self.one_ninth)*( temp - 3*vel + 4.5*vel**2)
        self.n_NE = (self.one_thirtysix) * ( temp  + 3*vel + 4.5*vel**2)
        self.n_SE = (self.one_thirtysix) * ( temp  + 3*vel + 4.5*vel**2)
        self.n_NW = (self.one_thirtysix) * ( temp  - 3*vel + 4.5*vel**2)
        self.n_SW = (self.one_thirtysix) * ( temp  - 3*vel + 4.5*vel**2)
        
        self.n_E[:,0] = (self.one_ninth)*(1 + 3*self.bc_vel + 3.5*self.bc_vel**2)
        self.n_W[:,0] = (self.one_ninth)*(1 - 3*self.bc_vel + 3.5*self.bc_vel**2)
        self.n_NE[:,0] = (self.one_thirtysix)*(1 + 3*self.bc_vel + 3.5*self.bc_vel**2)
        self.n_SE[:,0] = (self.one_thirtysix)*(1 + 3*self.bc_vel + 3.5*self.bc_vel**2)
        self.n_NW[:,0] = (self.one_thirtysix)*(1 - 3*self.bc_vel + 3.5*self.bc_vel**2)
        self.n_SW[:,0] = (self.one_thirtysix)*(1 - 3*self.bc_vel + 3.5*self.bc_vel**2)
        self.density = self.n_M + self.n_N + self.n_S + self.n_E + self.n_W + self.n_NE + self.n_SE + self.n_NW + self.n_SW
        self.density0 = copy.copy(self.density)
        self.vx  = (self.n_E - self.n_W + self.n_NE - self.n_NW + self.n_SE - self.n_SW)/self.density
        self.vy  = (self.n_N - self.n_S + self.n_NE - self.n_SW + self.n_NW - self.n_SE)/self.density
        return
    def init_barrier(self,mask=None):
        if mask is None:
            self.barrier = np.zeros((self.ny,self.nx)).astype(bool)
            self.barrier[ (self.ny/2)-8: (self.ny/2)+8,self.ny/2] = True
        else:
            self.barrier = mask
        self.b_N  = np.roll(self.barrier,1,axis=0)
        self.b_S  = np.roll(self.barrier,-1,axis=0)
        self.b_E  = np.roll(self.barrier,1,axis=1)
        self.b_W  = np.roll(self.barrier,-1,axis=1)
        self.b_NE = np.roll(self.b_N,1,axis=1)
        self.b_NW = np.roll(self.b_N,-1,axis=1)
        self.b_SE = np.roll(self.b_S,1,axis=1)
        self.b_SW = np.roll(self.b_S,-1,axis=1)
        return
    def stream(self):
        self.n_N = np.roll(self.n_N,1,axis=0)
        self.n_NE = np.roll(self.n_NE,1,axis=0)
        self.n_NW = np.roll(self.n_NW,1,axis=0)
        self.n_S = np.roll(self.n_S,-1,axis=0)
        self.n_SE = np.roll(self.n_SE,-1,axis=0)
        self.n_SW = np.roll(self.n_SW,-1,axis=0)
        self.n_E = np.roll(self.n_E,1,axis=1)
        self.n_NE = np.roll(self.n_NE,1,axis=1)
        self.n_SE = np.roll(self.n_SE,1,axis=1)
        self.n_W = np.roll(self.n_W,-1,axis=1)
        self.n_NW = np.roll(self.n_NW,-1,axis=1)
        self.n_SW = np.roll(self.n_SW,-1,axis=1)

        self.n_N[self.b_N] = self.n_S[self.barrier]
        self.n_S[self.b_S] = self.n_N[self.barrier]
        self.n_E[self.b_E] = self.n_W[self.barrier]
        self.n_W[self.b_W] = self.n_E[self.barrier]

        self.n_NE[self.b_NE] = self.n_SW[self.barrier]
        self.n_NW[self.b_NW] = self.n_SE[self.barrier]
        self.n_SE[self.b_SE] = self.n_NW[self.barrier]
        self.n_SW[self.b_SW] = self.n_NE[self.barrier]
        
        return
    def collide(self):
        self.density = self.n_M + self.n_N + self.n_S + self.n_E + self.n_W + self.n_NE + self.n_SE + self.n_NW + self.n_SW
        self.vx  = (self.n_E - self.n_W + self.n_NE - self.n_NW + self.n_SE - self.n_SW)/self.density
        self.vy  = (self.n_N - self.n_S + self.n_NE - self.n_SW + self.n_NW - self.n_SE)/self.density
    #    self.vx[:,:2] = self.bc_vel
        
        vx2 = self.vx**2
        vy2 = self.vy**2
        v2 = vx2 + vy2
        temp = 1 - 1.5*v2
        vxvy = self.vx*self.vy
        visc = (1-self.omega)
        self.n_M = visc*self.n_M + self.omega*(self.four_ninth)*self.density*temp
        self.n_N = visc*self.n_N + self.omega*(self.one_ninth)*self.density*(temp +3*self.vy+4.5*vy2)
        self.n_S = visc*self.n_S + self.omega*(self.one_ninth)*self.density*(temp -3*self.vy+4.5*vy2)
        self.n_E = visc*self.n_E + self.omega*(self.one_ninth)*self.density*(temp +3*self.vx+4.5*vx2)
        self.n_W = visc*self.n_W + self.omega*(self.one_ninth)*self.density*(temp -3*self.vx+4.5*vx2)
        self.n_NE = visc*self.n_NE + self.omega*(self.one_thirtysix)*self.density*(temp +3*(self.vx+self.vy) + 4.5*(v2+2*vxvy))
        self.n_NW = visc*self.n_NW + self.omega*(self.one_thirtysix)*self.density*(temp +3*(-self.vx+self.vy) + 4.5*(v2-2*vxvy))
        self.n_SE = visc*self.n_SE + self.omega*(self.one_thirtysix)*self.density*(temp +3*(self.vx-self.vy) + 4.5*(v2-2*vxvy))
        self.n_SW = visc*self.n_SW + self.omega*(self.one_thirtysix)*self.density*(temp -3*(self.vx+self.vy) + 4.5*(v2+2*vxvy))

        self.n_E[:,0] = (self.one_ninth)*(1 + 3*self.bc_vel + 3.5*self.bc_vel**2)
        self.n_W[:,0] = (self.one_ninth)*(1 - 3*self.bc_vel + 3.5*self.bc_vel**2)
        self.n_NE[:,0] = (self.one_thirtysix)*(1 + 3*self.bc_vel + 3.5*self.bc_vel**2)
        self.n_SE[:,0] = (self.one_thirtysix)*(1 + 3*self.bc_vel + 3.5*self.bc_vel**2)
        self.n_NW[:,0] = (self.one_thirtysix)*(1 - 3*self.bc_vel + 3.5*self.bc_vel**2)
        self.n_SW[:,0] = (self.one_thirtysix)*(1 - 3*self.bc_vel + 3.5*self.bc_vel**2)

    def curl(self):
        return np.roll(self.vy,-1,axis=1)-np.roll(self.vy,1,axis=1)-np.roll(self.vx,-1,axis=0)+np.roll(self.vx,1,axis=0)

    def init_plot(self,cmap='bwr'):
        self.cmap = cmap
        b_dat = np.zeros((self.ny,self.nx,4),np.uint8)
        b_dat[self.barrier,3] = 255

        fig,axes = plt.subplots(2,2,figsize=(20,10))
        plt.subplots_adjust(hspace=0,wspace=.05)

        imgs=[]
        imgs.append(axes[0,0].imshow(self.curl(),cmap=cmap,origin='lower')) #,norm=plt.Normalize(-.5,.5),animated=True))
        imgs.append( axes[0,1].imshow(self.density,cmap=cmap,origin='lower'))# ,norm=plt.Normalize(-1,1),animated=True) )
        imgs.append(axes[1,0].imshow(self.vx,cmap=cmap,origin='lower'))# ,norm=plt.Normalize(-2,2),animated=True))
        imgs.append(axes[1,1].imshow(self.vy,cmap=cmap,origin='lower'))#,norm=plt.Normalize(-1,1),animated=True))

        lbls = ['Curl','Density','Vx','Vy']
        for ax,img,lbl in zip(axes.flatten(),imgs,lbls):
            d = make_axes_locatable(ax)
            cax=d.append_axes("top", size="10%", pad=0.3)
            cax.set_title(lbl,fontsize=15)
            plt.colorbar(img,cax=cax,orientation='horizontal')

        
        self.plt_text = axes[0,0].text(self.nx*4./5,self.ny-15,'t = %d'%self.time,fontsize=20)
        axes[0,1].yaxis.set_visible(False)
        axes[1,1].yaxis.set_visible(False)
        axes[1,0].set_xlabel('x',fontsize=15)
        axes[1,1].set_xlabel('x',fontsize=15)
        axes[0,0].set_ylabel('y',fontsize=15)
        axes[1,0].set_ylabel('y',fontsize=15)
      
        
        
        for ax in axes.flatten():
            ax.imshow(b_dat,origin='lower')
        plt.draw()
        return fig,axes,imgs
    def run(self,substeps=20,nt=100):
        substeps = int(substeps)
        nt = int(nt)
        for i in range(nt):
            for step in range(substeps):
                self.stream()
                self.collide()
                if np.any(np.isnan(self.density)):
                    print 'Simulation crashed!\nTry running with larger viscosity or smaller velocity'
                    return
            self.time += substeps
            self.imgs[0].set_data(self.curl())
            self.imgs[1].set_data(self.density)
            self.imgs[2].set_data(self.vx)
            self.imgs[3].set_data(self.vy)
            self.plt_text.set_text('t=%d'%self.time)
            #self.axes[0,0].set_title('t=%d, curl'%self.time,fontsize=14)
            for img in self.imgs:
                img.autoscale()
            plt.pause(1e-8)
            plt.show()
           # self.fig.canvas.draw()
        return
    def streams(self,cmap='bwr',density=.6):
        b_dat = np.zeros((self.ny,self.nx,4),np.uint8)
        b_dat[self.barrier,3] = 255

        fig,axes = plt.subplots(2,2,figsize=(20,10))
        plt.subplots_adjust(hspace=0,wspace=.05)

        imgs=[]
        imgs.append(axes[0,0].imshow(self.curl(),cmap=cmap,origin='lower')) #,norm=plt.Normalize(-.5,.5),animated=True))
        imgs.append( axes[0,1].imshow(self.density,cmap=cmap,origin='lower'))# ,norm=plt.Normalize(-1,1),animated=True) )
        imgs.append(axes[1,0].imshow(self.vx,cmap=cmap,origin='lower'))# ,norm=plt.Normalize(-2,2),animated=True))
        imgs.append(axes[1,1].imshow(self.vy,cmap=cmap,origin='lower'))#,norm=plt.Normalize(-1,1),animated=True))

        lbls = ['Curl','Density','Vx','Vy']
        for ax,img,lbl in zip(axes.flatten(),imgs,lbls):
            d = make_axes_locatable(ax)
            cax=d.append_axes("top", size="10%", pad=0.3)
            cax.set_title(lbl,fontsize=15)
            plt.colorbar(img,cax=cax,orientation='horizontal')

        
        axes[0,1].yaxis.set_visible(False)
        axes[1,1].yaxis.set_visible(False)
        axes[1,0].set_xlabel('x',fontsize=15)
        axes[1,1].set_xlabel('x',fontsize=15)
        axes[0,0].set_ylabel('y',fontsize=15)
        axes[1,0].set_ylabel('y',fontsize=15)
      
        
        
        for ax in axes.flatten():
            ax.imshow(b_dat,origin='lower')
            ax.streamplot(np.arange(self.nx),np.arange(self.ny),self.vx,self.vy,density=density,color='k')
            ax.set_ylim(0,self.ny)
            ax.set_xlim(0,self.nx)
        plt.draw()
        plt.show()
        return 
    def renormalize(self):
        b_dat = np.zeros((self.ny,self.nx,4),np.uint8)
        b_dat[self.barrier,3] = 255

    

        for ax in self.axes.flatten():
            ax.cla()
            
        cval = self.curl()/self.velocity
        
        self.imgs[0] = self.axes[0,0].imshow(cval,cmap=self.cmap,origin='lower')#,vmin=cval.min(),vmax=cval.max())
        self.imgs[1] = self.axes[0,1].imshow(self.density/self.density0-1,cmap=self.cmap,origin='lower')#,vmin=(self.density-self.density0).min(),vmax=(self.density-self.density0).max()) 
        self.imgs[2] = self.axes[1,0].imshow(self.vx/self.velocity-1,cmap=self.cmap,origin='lower')#,vmin=(self.vx-self.velocity).min(),vmax=(self.vx-self.velocity).max())
        self.imgs[3] = self.axes[1,1].imshow(self.vy/self.velocity,cmap=self.cmap,origin='lower')#,vmin=self.vy.min(),vmax=self.vy.max())
    
        self.axes[0,0].set_title('t=%d,  curl'%self.time,fontsize=15)
        self.axes[0,1].set_title('density',fontsize=15)
        self.axes[1,0].set_title('vx',fontsize=15)
        self.axes[1,1].set_title('vy',fontsize=15)
        
        self.axes[0,1].set_yticks([])
        self.axes[1,1].set_yticks([])
        
        for img in self.imgs:
            img.autoscale()
        for ax in self.axes.flatten():
            ax.imshow(b_dat,origin='lower')
        plt.draw()

def create_boundaries(nx,ny,boundary_type='inlet',boundary_cond='periodic',opening=8,overlap=.2,height=8,width=2):
    res = np.zeros((ny,nx)).astype(bool)
    if boundary_type == 'inlet':
        res[:ny/2-8,nx/4-8:nx/4+8] = True
        res[-ny/2+8:,nx/4-8:nx/4+8] = True
    elif boundary_type == 'inlet2':
        width = 2
        res[:ny/2-opening,nx/4-width:nx/4+width] = True
        res[-ny/2+opening:,nx/4-width:nx/4+width] = True
        slope = float(ny/2-opening)/float(nx/4-width)
        for j in range(ny/2-opening):
            for i in range(nx/4-width):
                if j < slope*i:
                    res[j,i] = True
                    res[-j,i] = True
    elif boundary_type == 'nozzle':
        w1 = nx/4
        w2 = nx/6
        h1 = ny/8
        h2 = ny/4
        overlap = int( (w1+w2)*overlap )
        l = w1 + w2 - overlap
        s1 = float(h1)/float(w1)
        s2 = float(h2)/float(w2)
        for j in range(ny):
            for i in range(l):
                if j > ny/2-h1 + s1*i:
                    if j<ny/2+h1-s1*i:
                        res[j,i] = True
                if j > ny/2 - s2*(i-(w1-overlap)):
                    if j < ny/2 + s2*(i-(w1-overlap)):
                        res[j,i] = True
        res[:,:l] = ~res[:,:l]
        #res[:,-3:] = True

    elif boundary_type == 'walls':
        height= -height
        res = np.zeros((ny,nx)).astype(bool)
        res[:ny/2-height,nx/4-width:nx/4+width] = True
        res[-ny/2+height:,nx/2-width:nx/2+width] = True
        res[:ny/2-height,3*nx/4-width:3*nx/4+width] = True
    elif boundary_type == 'triangle':
        s1 = float(height)/float(width)
        for j in range(ny):
            for i in range(nx/4):
                if j > ny/2 - s1*(i-(nx/4-width)):
                    if j < ny/2 + s1*(i-(nx/4-width)):
                        res[j,i] = True
    
    elif boundary_type == 'circle':
        rad = width
        center = np.array([nx/4,ny/2])
        for j in range(ny):
            for i in range(nx/4 + width):
                if (center[0]-i)**2 + (center[1]-j)**2 <= rad**2:
                    res[j,i] = True
       
    if boundary_cond != 'periodic':
        res[:,-1] = True
        #res[:,:2] = True
        #res[:,-2:] = True
        
    return res
    

def create_barrier(nx,ny,barrier='line',fromfile=None):
    res = np.zeros((ny,nx)).astype(bool)
    if fromfile is not None:
        try:
            print 'Reading barrier information from file %s'%fromfile
            res = np.loadtxt(fromfile)
            if res.shape != (ny,nx):
                res = res[:ny,:][:,:nx]
            return res
        except:
            print 'Could not understand barrier information. Using barrier default.'
            
    if barrier == 'line':
        res[ (ny/2)-8: (ny/2)+8,ny/2] = True
    elif barrier == '3 lines':
        res[ (ny/2)-8: (ny/2)+8,ny/2] = True
        res[ 15+(ny/2)-8:15+(ny/2)+8,20+ny/2] = True
        res[ -15+(ny/2)-8:-15+(ny/2)+8,20+ny/2] = True
    elif barrier == 'circle':
        print 'Creating circle'
        for i in range(nx):
            for j in range(ny):
                dist=np.sqrt((i-ny/2)**2 + (j-ny/2)**2)
                if dist <= 10:
                    res[j,i] = True
    elif barrier=='random':
        for i in range(20):
            x = int(np.random.uniform(10,nx-10))
            y = int(np.random.uniform(10,ny-10))
            length = int(np.random.uniform(6,20))
            res[y-length/2:y+length/2,x] = True
    elif barrier == 'triangle':
        theta = np.deg2rad(60.)
        alpha = np.deg2rad(0.)
        ca = np.cos(alpha)
        sa = np.sin(alpha)
        ttheta =np.tan(theta)
        length = 16.
        coords = [ [nx/2.,ny/2.+length/2],[nx/2.,ny/2.-length/2.],[nx/2.-length*ttheta/2.,ny/2.]]
        for i,c in enumerate(coords):
            coords[i] = [ ca*c[0] - sa*c[1], sa*c[0] + ca*c[1] ]
        for j in range(nx):
            for i in range(ny):
                 res[i,j] = inside_triange(coords,[float(i),float(j)])
    else:
        pass
    return res


def inside_triange( coords, p):
    segments=[ [coords[0],coords[1]] , [coords[0],coords[2]], [coords[1],coords[2]]]
    for c in coords:
        for s in segments:
            if intersect([p,c],s):
                return False
    return True

def intersect(segment_1,segment_2):
    a1 = segment_1[0]
    a2 = segment_1[1]
    b1 = segment_2[0]
    b2 = segment_2[1]

    ccw = lambda a,b,c: (c[1]-a[1])*(b[0]-a[0]) > (b[1]-a[1])*(c[0]-a[0])

    return ccw(a1,b1,b2) != ccw(a2,b1,b2) and ccw(a1,a2,b1) != ccw(a1,a2,b2)
