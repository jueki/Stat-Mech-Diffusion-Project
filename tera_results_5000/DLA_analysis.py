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
from tera_data import pList, tList, names

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
        ax.w_xaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
        ax.w_yaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
        ax.w_zaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
    plt.subplots_adjust(wspace=0.01,hspace=0.01)


    # Organize simulation results for plotting
    totalTime = times[-1]
    xs = [p[0] for p in positions]
    ys = [p[1] for p in positions]
    zs = [p[2] for p in positions]
    colors = [float(t)/totalTime for t in times]

    # Plot it
    for ax in subplots:
        ax.plot([bx,bx],[-bx,-bx],[-bx,bx],'-',linewidth=1,color='darkblue')
        ax.plot([bx,bx],[-bx,bx],[-bx,-bx],'-',linewidth=1,color='darkgreen')
        ax.plot([-bx,bx],[bx,bx],[-bx,-bx],'-',linewidth=1,color='darkred')
        for (xi,yi,zi,ci) in zip(xs,ys,zs,colors):
            (x,y,z) = drawSphere(xi,yi,zi,colRad)
            ax.plot_surface(x, y, z,color=cm.rainbow(ci),
                rstride=4, cstride=4,linewidth=0,shade=False)

  # Set up colorbar ###### STILL WORKING ON THIS
    norm = clrs.Normalize(vmin=0,vmax=totalTime)
    sm = cm.ScalarMappable(cmap=cm.rainbow, norm=norm)
    sm.set_array([])
    cax = fig.add_axes([0.92, 0.375, 0.01, 0.25])
    cb = fig.colorbar(sm, cax=cax)     
    cb.outline.remove() 
    cb.ax.tick_params(labelsize=5) 
    cb.ax.yaxis.offsetText.set(size=5)
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
    rms = rmsRadius(positions)
    print ("RMS radius: %f"%rms)
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

    
def rmsRadius(positions):
    """ Calculate the radius of gyration for a collection of particles
    about center of mass"""
    inertia = 0
    N = len(positions)
    (rx,ry,rz) = (0,0,0)
    for pos in positions:
        rx += pos[0]
        ry += pos[1]
        rz += pos[2]
    rx = rx/N
    ry = ry/N
    rz = rz/N
    for pos in positions:
        inertia += ((pos[0]-rx)**2 + (pos[1]-ry)**2 + (pos[2]-rz)**2)
    radius = np.sqrt(inertia/N)
    return radius
 

 

#==============================================================================
#  Main Program
#==============================================================================
# Put the data here
colRad = 1
length = 120

simName = '5000_120len'

for i in range(len(tList)):
    simName = '5000_120len_'+names[i]
    positions = pList[i]
    print 'Prob = ' + names[i]
    stickRadii = [np.linalg.norm(pos) for pos in positions]
    maxR = stickRadii[0]
    for r in stickRadii:
        if r > maxR:
            maxR = r
    print 'Maximum radius = %f'%maxR
    rms = rmsRadius(positions)
    print ("RMS radius = %f"%rms)
    print
    plotParticles(pList[i],tList[i])

# Analyze data and plot

"""
nStuck = range(len(times02)) # number of particles stuck at each time
# Find max radius at each time
maxList = []
for positions in pList:
    maxRadii = []
    stickRadii = [np.linalg.norm(pos) for pos in positions]
    maxR = stickRadii[0]
    for r in stickRadii:
        if r > maxR:
            maxR = r
        maxRadii += [maxR]
    maxList += [maxRadii]
# Set up plot
colors = ['r-','y-','g-','b-','m-']
fig = plt.figure(facecolor='w')
ax1 = plt.subplot(2, 1, 1)
for i in range(len(tList)):
    plt.plot(tList[i], nStuck,colors[i],label=labels[i])
plt.ylabel('Particles Stuck')
plt.legend(loc=4,bbox_to_anchor=(0.95, 0.1),fontsize = 10,frameon = False)
plt.yticks(fontsize = 10)
plt.setp(ax1.get_xticklabels(), visible=False)
plt.subplot(2, 1, 2)
for i in range(len(tList)):
    plt.plot(tList[i], maxList[i],colors[i])
plt.xlabel('Time Steps')
plt.xticks(fontsize = 10)
plt.yticks(fontsize = 10)
plt.ylabel('Maximum Radius')
# Save figure
fileName = simName+'_anl.png'
plt.savefig(fileName, dpi=300, bbox_inches='tight')
"""