import fluidflow as ff
import numpy as np



def nozzle_model(lx=100,ly=50,opening=.1,velocity=.1,viscosity=.07,cmap='bwr',inflow=True):
    wall = ff.create_boundaries(lx,ly,boundary_type='nozzle',overlap=opening)
    return ff.Simulation(velocity=velocity,viscosity=viscosity,mask=wall,lx=lx,ly=ly,cmap=cmap)

def inlet_model(lx=100,ly=50,velocity=.1,viscosity=.07,cmap='bwr',opening=5,inflow=True):
    wall = ff.create_boundaries(lx,ly,boundary_type='inlet2',opening=opening)
    return ff.Simulation(velocity=velocity,viscosity=viscosity,mask=wall,lx=lx,ly=ly,cmap=cmap)

def pipe_model(lx=100,ly=50,velocity=.1,viscosity=.07,cmap='bwr',width=2,height=8,inflow=True):
    wall = ff.create_boundaries(lx,ly,boundary_type='walls',width=width,height=height)
    return ff.Simulation(velocity=velocity,viscosity=viscosity,mask=wall,lx=lx,ly=ly,cmap=cmap)

def wing_model(lx=100,ly=50,velocity=.1,viscosity=.07,cmap='bwr',width=2,height=8,inflow=True):
    wall = ff.create_boundaries(lx,ly,boundary_type='triangle',width=width,height=height)
    return ff.Simulation(velocity=velocity,viscosity=viscosity,mask=wall,lx=lx,ly=ly,cmap=cmap)
    
def cylinder_model(lx=100,ly=50,velocity=.1,viscosity=.07,cmap='bwr',radius=8,inflow=True):
    wall = ff.create_boundaries(lx,ly,boundary_type='circle',width=radius)
    return ff.Simulation(velocity=velocity,viscosity=viscosity,mask=wall,lx=lx,ly=ly,cmap=cmap)