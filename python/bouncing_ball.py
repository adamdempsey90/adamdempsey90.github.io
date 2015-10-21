import matplotlib.pyplot as plt
import numpy as np



g = 9.8

def evolve(vi,theta,h,t_end,dt):
    """ Evolution function given the initial velocity (speed and angle),
    initial height, ending time, and time step."""

# Initialize everything
    x,y,vx,vy = initialize(vi,theta,h)

    fig,axes,data = initialize_plot(h,x,y,vx,vy)

# Start moving the ball from t=0 to t=t_end in steps of dt
    t = 0
    while t <= t_end:
        x,y,vx,vy,t = update_simple(x,y,vx,vy,t,dt)
        update_plot(fig,axes,data,t,h,x,y,vx,vy)

    return

def initialize(vi,theta,h):
    """ Initialize the position and velocity of the ball.
    theta is in degrees.  """
    vx = vi*np.cos(theta * np.pi/180)
    vy = vi*np.sin(theta * np.pi/180)
    x = 0
    y = h
    return x,y,vx,vy


def update_simple(x,y,vx,vy,t,dt):
    """ Simple update scheme that advnaces the ball forward in
    time by dt. """

# Set the variables to their new values using v = dx/dt and a = dv/dt
    vx_new =            # FIX ME!
    vy_new =            # FIX ME!
    x_new =             # FIX ME!
    y_new =            # FIX ME!


# What happens if the ball hits the wall.
# What are these lines doing?
    if y_new < 0:
        y_new = 0   # Set the y coordinate to zero
        dt_new = -y/vy  # Calculate the time it takes to hit the floor
        vy_new = -(vy-g*dt_new)  # Calculate the velocity that the ball hits the floor and reverse it.
        x_new = x + vx*dt_new
# Advance the simulation time t by dt  on the next 2 lines
        t =  t     # FIX ME!
    else:
        t =  t      # FIX ME!

    return x_new,y_new,vx_new,vy_new,t


def update_leapfrog(x,y,vx,vy,t,dt):
    """ Update scheme based on a drift-kick-drift method. """

    # Set the variables to their new values using the drift-kick-drift scheme

    # First Drift

#    vx_new =        # FIX ME!
#    vy_new =        # FIX ME!

    # Now Kick

#    x_new =         # FIX ME!
#    y_new =         # FIX ME!

    # Second Drift

#    vx_new =        # FIX ME!
#    vy_new =        # FIX ME!


# Special updates for when the ball hits the floor
    if y_new < 0:
        y_new = 0
        dt_new = (vy + np.sqrt(vy*vy +2*g*y) )/(g)
        vy_new = - (vy - g*dt_new)
        x_new = x + vx*dt_new

# Advance the simulation time t by dt  on the next 2 lines
        t = t    # FIX ME!
    else:
        t = t    # FIX ME!

    return     # FIX ME!

def update_dissipative(x,y,vx,vy,t,dt):
    """ Update the position and velocity of the bouncing ball over a time dt,
    taking into account a loss of kinetic energy as the ball bounces. """

    return



# Begin plotting functions. Nothing to be done here.

def initialize_plot(h,x,y,vx,vy):
    """ Plot is set up as 2x2 square"""
    fig,axes = plt.subplots(2,2,figsize=(15,13))
    data = []
# First do the x,y ball plot
    axes[0,0].set_xlabel('x')
    axes[0,0].set_ylabel('y')
    axes[0,0].set_title('t = 0')
    axes[0,0].plot([0,0],[0,h],'-k',linewidth=3)
    axes[0,0].plot([-3,0],[h,h],'-k',linewidth=3)
    line, = axes[0,0].plot(x,y,'ob',markersize=20)
    axes[0,0].set_xlim(-3,4*h)
    axes[0,0].set_ylim(0,3*h)
    data.append(line)

# Next do x,y versus time,

    axes[1,0].set_xlabel('t')
    axes[1,0].set_ylabel('x,y')
    axes[1,0].set_title('x,y position of ball')
    linex,=axes[1,0].plot(0,x,'-r',label=r'x-position')
    liney,=axes[1,0].plot(0,y,'-b',label=r'y-position')
    data.append( (linex,liney))
    axes[1,0].set_xlim(0,1)
    axes[1,0].legend(loc='upper right')
# Next do vx,vy versus time

    axes[0,1].set_xlabel('t')
    axes[0,1].set_ylabel('vx,vy')
    axes[0,1].set_title('velocities of ball')
    linex,=axes[0,1].plot(0,vx,'-r',label=r'x-velocity')
    liney,=axes[0,1].plot(0,vy,'-b',label=r'y-velocity')
    data.append( (linex,liney))
    axes[0,1].set_xlim(0,1)
    axes[0,1].legend(loc='upper right')

# Finally add energy plot
    ke = .5*(vx**2 + vy**2)
    pe = g*y
    axes[1,1].set_xlabel('t')
    axes[1,1].set_ylabel('E')
    axes[1,1].set_title('Kinetic and potential energies of ball')
    linex,=axes[1,1].plot(0,ke,'-r',label=r'Kinetic Energy')
    liney,=axes[1,1].plot(0,pe,'-b',label=r'Potential Energy')
    data.append( (linex,liney))
    axes[1,1].set_xlim(0,1)
    axes[1,1].legend(loc='upper right')

    plt.show()

    return fig,axes,data

def update_plot(fig,axes,data,t,h,x,y,vx,vy):
# update x,y plot
    data[0].set_xdata(x)
    data[0].set_ydata(y)
    axes[0,0].set_title('t = %.2f' % t)
    if x > axes[0,0].get_xlim()[1]:
        axes[0,0].set_xlim(-3,2*x)
    if y > axes[0,0].get_ylim()[1]:
        axes[0,0].set_ylim(0,2*y)

# update position plot
    linex,liney = data[1]
    linex.set_xdata( np.append(linex.get_xdata(), t))
    linex.set_ydata(np.append(linex.get_ydata(),x))
    liney.set_xdata( np.append(liney.get_xdata(), t))
    liney.set_ydata(np.append(liney.get_ydata(),y))
    if t > axes[1,0].get_xlim()[1]:
        axes[1,0].set_xlim(0,2*t)
    if y > axes[1,0].get_ylim()[1] or x > axes[1,0].get_ylim()[1]:
        axes[1,0].set_ylim(0,2*max(x,y))

# update velocity plot
    linex,liney = data[2]
    linex.set_xdata( np.append(linex.get_xdata(), t))
    linex.set_ydata(np.append(linex.get_ydata(),vx))
    liney.set_xdata( np.append(liney.get_xdata(), t))
    liney.set_ydata(np.append(liney.get_ydata(),vy))
    if t > axes[0,1].get_xlim()[1]:
        axes[0,1].set_xlim(0,2*t)
    if vy > axes[0,1].get_ylim()[1] or vx > axes[0,1].get_ylim()[1]:
        axes[0,1].set_ylim(-2*max(abs(vx),abs(vy)),2*max(abs(vx),abs(vy)))

# update energy plot
    ke = .5*(vx**2 + vy**2)
    pe = g*y
    linex,liney = data[3]
    linex.set_xdata( np.append(linex.get_xdata(), t))
    linex.set_ydata(np.append(linex.get_ydata(),ke))
    liney.set_xdata( np.append(liney.get_xdata(), t))
    liney.set_ydata(np.append(liney.get_ydata(),pe))
    if t > axes[1,1].get_xlim()[1]:
        axes[1,1].set_xlim(0,2*t)
    if  ke > axes[1,1].get_ylim()[1] or pe > axes[1,1].get_ylim()[1]:
        new_upper = 2*max(ke,pe)
    else:
        new_upper = axes[1,1].get_ylim()[1]
    if ke < axes[1,1].get_ylim()[0] or pe < axes[1,1].get_ylim()[0]:
        new_lower = min(ke,pe)
    else:
        new_lower = axes[1,1].get_ylim()[0]
    axes[1,1].set_ylim(new_lower,new_upper)


    plt.pause(.0001)
    plt.show()
    return
