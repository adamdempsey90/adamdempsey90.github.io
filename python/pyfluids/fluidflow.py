import numpy as np
import matplotlib.pyplot as plt

one_ninth = 1./9
one_thirtysix = 1./36
four_ninth = 4*one_ninth

class Simulation(object):
    def __init__(self,lx=200,ly=80,viscosity=0.02,velocity=0.1,mask=None):
        self.lx = lx
        self.nx = lx
        self.ly = ly
        self.ny = ly
        self.viscosity = viscosity
        self.velocity = velocity

        self.omega = 1./(3*viscosity + 0.5)
        self.init_fluid()
        self.init_barrier(mask)
        self.fig,self.axes,self.imgs = self.init_plot()
    def init_fluid(self):
        temp = (np.ones((self.ny,self.nx)) - 1.5*self.velocity**2)
        self.n_M = (four_ninth)*temp
        self.n_N = (one_ninth) * temp
        self.n_S = (one_ninth) * temp
        self.n_E =  (one_ninth)*( temp + 3*self.velocity + 4.5*self.velocity**2)
        self.n_W =  (one_ninth)*( temp - 3*self.velocity + 4.5*self.velocity**2)
        self.n_NE = (one_thirtysix) * ( temp  + 3*self.velocity + 4.5*self.velocity**2)
        self.n_SE = (one_thirtysix) * ( temp  + 3*self.velocity + 4.5*self.velocity**2)
        self.n_NW = (one_thirtysix) * ( temp  - 3*self.velocity + 4.5*self.velocity**2)
        self.n_SW = (one_thirtysix) * ( temp  - 3*self.velocity + 4.5*self.velocity**2)
        self.density = self.n_M + self.n_N + self.n_S + self.n_E + self.n_W + self.n_NE + self.n_SE + self.n_NW + self.n_SW
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

        vx2 = self.vx**2
        vy2 = self.vy**2
        v2 = vx2 + vy2
        temp = 1 - 1.5*v2
        vxvy = self.vx*self.vy
        visc = (1-self.omega)
        self.n_M = visc*self.n_M + self.omega*(four_ninth)*self.density*temp
        self.n_N = visc*self.n_N + self.omega*(one_ninth)*self.density*(temp +3*self.vy+4.5*vy2)
        self.n_S = visc*self.n_S + self.omega*(one_ninth)*self.density*(temp -3*self.vy+4.5*vy2)
        self.n_E = visc*self.n_E + self.omega*(one_ninth)*self.density*(temp +3*self.vx+4.5*vx2)
        self.n_W = visc*self.n_W + self.omega*(one_ninth)*self.density*(temp -3*self.vx+4.5*vx2)
        self.n_NE = visc*self.n_NE + self.omega*(one_thirtysix)*self.density*(temp +3*(self.vx+self.vy) + 4.5*(v2+2*vxvy))
        self.n_NW = visc*self.n_NW + self.omega*(one_thirtysix)*self.density*(temp +3*(-self.vx+self.vy) + 4.5*(v2-2*vxvy))
        self.n_SE = visc*self.n_SE + self.omega*(one_thirtysix)*self.density*(temp +3*(self.vx-self.vy) + 4.5*(v2-2*vxvy))
        self.n_SW = visc*self.n_SW + self.omega*(one_thirtysix)*self.density*(temp -3*(self.vx+self.vy) + 4.5*(v2+2*vxvy))

        self.n_E[:,0] = (one_ninth)*(1 + 3*self.velocity + 3.5*self.velocity**2)
        self.n_W[:,0] = (one_ninth)*(1 - 3*self.velocity + 3.5*self.velocity**2)
        self.n_NE[:,0] = (one_thirtysix)*(1 + 3*self.velocity + 3.5*self.velocity**2)
        self.n_SE[:,0] = (one_thirtysix)*(1 + 3*self.velocity + 3.5*self.velocity**2)
        self.n_NW[:,0] = (one_thirtysix)*(1 - 3*self.velocity + 3.5*self.velocity**2)
        self.n_SW[:,0] = (one_thirtysix)*(1 - 3*self.velocity + 3.5*self.velocity**2)

    def curl(self):
        return np.roll(self.vy,-1,axis=1)-np.roll(self.vy,1,axis=1)-np.roll(self.vx,-1,axis=0)+np.roll(self.vx,1,axis=0)

    def init_plot(self):
        #fig = plt.figure(figsize=(8,3))
        #ax = fig.add_subplot(111)
        b_dat = np.zeros((self.ny,self.nx,4),np.uint8)
        b_dat[self.barrier,3] = 255

        fig,axes = plt.subplots(2,2,figsize=(8,4))
        plt.subplots_adjust(hspace=0,wspace=0)

        imgs=[]
        imgs.append(axes[0,0].imshow(self.curl(),origin='lower',norm=plt.Normalize(-self.velocity,self.velocity)))
        imgs.append( axes[0,1].imshow(self.density,origin='lower',norm=plt.Normalize(self.density.min()*.7,self.density.max()*1.2)) )
        imgs.append(axes[1,0].imshow(self.vx,origin='lower',norm=plt.Normalize(-self.velocity*.1,self.velocity*3)))
        imgs.append(axes[1,1].imshow(self.vy,origin='lower',norm=plt.Normalize(-self.velocity,self.velocity)))

        axes[0,0].set_title('t=0,  curl',fontsize=15)
        axes[0,1].set_title('density',fontsize=15)
        axes[1,0].set_title('vx',fontsize=15)
        axes[1,1].set_title('vy',fontsize=15)
        for ax in axes.flatten():
            ax.imshow(b_dat,origin='lower')
        fig.canvas.draw()
        return fig,axes,imgs
    def run(self,substeps=20,nt=100):

        for i in range(nt):
            for step in range(substeps):
                self.stream()
                self.collide()

            self.imgs[0].set_array(self.curl())
            self.imgs[1].set_array(self.density)
            self.imgs[2].set_array(self.vx)
            self.imgs[3].set_array(self.vy)
            self.axes[0,0].set_title('t=%d, curl'%(i*substeps),fontsize=14)
            self.fig.canvas.draw()
        return
def create_barrier(nx,ny,barrier='line'):
    res = np.zeros((ny,nx)).astype(bool)
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
        #for i,c in enumerate(coords):
        #    coords[i] = [ ca*c[0] - sa*c[1], sa*c[0] + ca*c[1] ]
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
