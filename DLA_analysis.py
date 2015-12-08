#!/usr/bin/env python

#==============================================================================
#  Max Byers, Jon Ueki, Jonas Kaufman
#  December 7 2015
#  Program to analyze and plot DLA data.
#==============================================================================

from pylab import * # this includes numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d
import matplotlib.colors as clrs
import matplotlib.cm as cm

def drawSphere(xCenter, yCenter, zCenter, r):
    """" Sets up a sphere drawing for plotting """
    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    x=np.cos(u)*np.sin(v)
    y=np.sin(u)*np.sin(v)
    z=np.cos(v)
    # shift and scale sphere
    x = r*x + xCenter
    y = r*y + yCenter
    z = r*z + zCenter
    return (x,y,z)

def plotParticles(positions,times):
    """ Creates 3 different perspective plots of the simulation results """
    # Set up figure
    fig = plt.figure(facecolor='w')
    fig.set_figwidth(8)
    ax0 = fig.add_subplot(131, projection='3d',aspect='equal')
    ax1 = fig.add_subplot(132, projection='3d',aspect='equal')
    ax2 = fig.add_subplot(133, projection='3d',aspect='equal')
    (bx,by,bz) = (length/2,length/2,length/2) # box.bounds
    subplots = [ax0,ax1,ax2]
    fSize=10
    for i in range(len(subplots)):
        ax = subplots[i]
        ax.set_xlim3d([-bx, bx])
        ax.set_ylim3d([-by, by])
        ax.set_zlim3d([-bz, bz])
        ax.view_init(elev=10.0, azim=i*120) #
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks([])
    plt.subplots_adjust(wspace=0.01,hspace=0.01)

    # Organize simulation results for plotting
    totalTime = times[-1]
    xs = [p[0] for p in positions]
    ys = [p[1] for p in positions]
    zs = [p[2] for p in positions]
    colors = [float(t)/totalTime for t in times]

    # Plot it
    for ax in subplots:
        for (xi,yi,zi,ci) in zip(xs,ys,zs,colors):
            (x,y,z) = drawSphere(xi,yi,zi,colRad)
            ax.plot_surface(x, y, z,color=cm.rainbow(ci),
                rstride=4, cstride=4,linewidth=0,shade=False)
    
    # Set up colorbar ###### STILL WORKING ON THIS
    """
    norm = clrs.Normalize(vmin=0,vmax=totalTime)
    sm = cm.ScalarMappable(cmap=cm.rainbow, norm=norm)
    sm.set_array([])
    fig.colorbar(sm)
    """
    fileName = simName+'_pp.png'
    plt.savefig(fileName, dpi=300, bbox_inches='tight')

def analyzeData(positions,times):
    """ looks at particles stuck and max radius vs. time and plots them"""
    nStuck = range(len(times)) # number of particles stuck at each time
    # Find max radius at each time
    maxRadii = []
    stickRadii = [np.linalg.norm(pos) for pos in positions]
    maxR = stickRadii[0]
    for r in stickRadii:
        if r > maxR:
            maxR = r
        maxRadii += [maxR]
    print 'Maximum radius = %f'%maxR
    gyrX = radiusOfGyration(positions, 0)
    gyrY = radiusOfGyration(positions, 1)
    gyrZ = radiusOfGyration(positions, 2)
    print ("Radius of Gyration x: %f"%gyrX)
    print ("Radius of Gyration y: %f"%gyrY)
    print ("Radius of Gyration z: %f"%gyrZ)
    # Set up plot
    fig = plt.figure(facecolor='w')
    plt.subplot(2, 1, 1)
    plt.plot(times, nStuck,'b-')
    plt.ylabel('Particles Stuck')
    plt.subplot(2, 1, 2)
    plt.plot(times, maxRadii, 'b-')
    plt.xlabel('Time Steps')
    plt.ylabel('Maximum Radius')
    # Save figure
    fileName = simName+'_anl.png'
    plt.savefig(fileName, dpi=300, bbox_inches='tight')
    
def radiusOfGyration(positions, axis):
    """ Calculate the radius of gyration for a collection of particles
    about a given axis: 0=x, 1=y, 2=z"""
    inertia = 0
    numParticles = len(positions)
    a = (axis+1)%3
    b = (axis+2)%3
    for pos in positions:
        inertia += (pos[a]**2 + pos[b]**2)
    radius = inertia/numParticles
    return radius
 

#==============================================================================
#  Main Program
#==============================================================================
# Put the data here
#length = 50
simName = 'plttest2'

# Plot the data
#plotParticles(positions,times)

# Analyze data and plot
analyzeData(positions,times)