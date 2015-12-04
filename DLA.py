#!/usr/bin/env python

#==============================================================================
#  Max Byers, Jon Ueki, Jonas Kaufman
#  November 27 2015
#  Program to simulate diffusion-limited aggregation in three dimensions.
#==============================================================================

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d
import matplotlib.colors as clrs
import matplotlib.cm as cm
import random
import time

def fromSphericalCoordinates(r,theta,phi):
	""" converts from spherical to cartesian coordinates """
	x = r*np.sin(theta)*np.cos(phi)
	y = r*np.sin(theta)*np.sin(phi)
	z = r*np.cos(theta)
	return np.array([x,y,z])

def randomPosition(r):
	""" gives a random position on a sphere of radius r """
	theta = np.random.uniform(0,np.pi,1)[0]
	phi = np.random.uniform(0,2*np.pi,1)[0]
	return fromSphericalCoordinates(r,theta,phi)

class Particle:
	""" holds information for a single particle """

	def __init__(self,position=[0,0,0],radius=1,time=0, stickProb=1):
		self.position = np.array(position)
		self.radius = colRad
		self.time = time
		self.stickProb = stickProb

	def move(self,newPosition):
		self.position = np.array(newPosition)
  
	def decision(self):
		""" checks if particle has collided with another particle """
		return random.random() < self.stickProb

	def collided(self,other):
		""" checks if particle has collided with another particle """
		collisionRadius = self.radius + other.radius
		distance = np.linalg.norm(np.array(self.position) - np.array(other.position))
		if distance > collisionRadius:
			return False
		elif self.decision(): #Sticking prob turns to be true
			return True
		else: #Sticking prob turn out not to be true
			return False

	def step(self,dr):
		""" moves the particle by dr in random direction """
		newPosition = self.position + randomPosition(dr)
		self.move(newPosition)
  
	def getPosition(self):
		""" Returns the (x,y,z) position of the particle """
		return self.position
  
class Grid:
    """The 3d grid that holds the position of each partice"""
    
    def __init__(self, length):
        """A 3D cube with each length^3 many cells"""
        self.cell = [[[[]*length for x in range(length)] for y in range(length)]
                    for z in range(length)]
        self.length = length
        
    def inGrid(self, x, y, z):
        """Checks if this continous position is a valid position inside the grid"""
        gridPos = self.CTG([x,y,z])
        x = gridPos[0]
        y = gridPos[1]
        z = gridPos[2]
        
        if (x < 0 or x >= self.length):
            return False
        
        elif (y < 0 or y >= self.length):
            return False
        
        elif (z < 0 or z >= self.length):
            return False
            
        else:  
            return True
        
    def CTG(self, position):
        """Takes the continous position and transforms it into
            a position in the grid"""
        xGrid = int((position[0]) + (self.length/2))
        yGrid = int((position[1]) + (self.length/2))
        zGrid = int((position[2]) + (self.length/2))
        
        return [xGrid, yGrid, zGrid]
        
    
    def getParticles(self, x, y, z):
        """Gets all of the particles in the cell at position (x, y, z)"""
        return self.cell[x][y][z]
    
    def push(self, thisParticle):
        """Push a particle into a cell in the grid"""
        gridPos = self.CTG(thisParticle.getPosition())
        xPos = gridPos[0]
        yPos = gridPos[1]
        zPos = gridPos[2]
        
        self.cell[xPos][yPos][zPos].append(thisParticle)
        
class Simulation:
    def __init__(self, length, colRad, spawnRad, stepSize,time=0, stickProb=1):

        self.grid = Grid(length)
        self.allParticles = [] #A list of all of the particles existing
        self.stepSize = stepSize
        self.length = length
        self.time = time
        self.colRad = colRad
        self.spawnRad = spawnRad
        self.prob = stickProb #the stickign probability
        seed = Particle(radius=self.colRad)
        
        #Create a seed at (0,0,0)
        self.grid.push(seed)
        self.allParticles.append(seed)
        
        
    def collision(self, current):
        """Checks if a particle has collided with any particle in all 26
            neighboring cells"""
        pos = current.getPosition()
        pos = self.grid.CTG(pos)
        x = pos[0]
        y = pos[1]
        z = pos[2]
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    currentCell = self.grid.cell[x-1+i][y-1+j][z-1+k]
                    for l in range(len(currentCell)):
                        if(currentCell[l].collided(current)):
                            return True
                        else:
                            return False
            
    def inBounds(self, particle):
        """Checks if the particle is within the boundry, where the boundry
            is when we decided the particle wandered off too far"""
        pos = particle.getPosition()
        x = self.grid.CTG(pos)[0]
        y = self.grid.CTG(pos)[1]
        z = self.grid.CTG(pos)[2]
        
        # Bounds are one less than the grid size in order to prevent out of bounds
        # error during collision check
        if x < 0 or x >= self.length-1:
            return False
            
        elif y < 0 or y >= self.length-1:
            return False
            
        elif z < 0 or z >= self.length-1:
            return False
            
        else:
            return True
    
    def newParticle(self):
        """Creates a new particle and evolves it until it collides with seed"""
        particle = Particle(randomPosition(self.spawnRad),radius=self.colRad, stickProb=self.prob)

        while(self.inBounds(particle) and not self.collision(particle)):
            particle.step(self.stepSize)
            self.time += dt
        
        # If the particle is still inbounds, that means it collided with the seed
        # So fix it into the grid and append it to the list of existing particles
        if (self.inBounds(particle)):
            particle.time = self.time
            self.allParticles.append(particle)
            self.grid.push(particle)
            print 'Particle %d stuck!'%(len(self.allParticles) - 1)
            
    def run(self, numParticles):
        """Runs the simulation until we have a certain amount of particles stuck
            to the seed"""
            
        while(len(self.allParticles) <= numParticles):
            self.newParticle()
            
    def getParticles(self):
        """Returns all of the particles in the simulation"""
        return self.allParticles
        
    def printParticlePositions(self):
        print("Positions: ")
        for x in range(len(self.allParticles)):
            print(self.allParticles[x].getPosition())

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

def saveData(positions,times,filename='data'):
    """ Saves simulation data to a file for later use"""
    positions = [p.tolist() for p in positions]
    s = 'length = %d\nrSpawn = %d\ndr = %f\n'%(length,rSpawn,dr)
    s += 'colRad = %f\nstickProb = %f\nsimName = %s\n'%(colRad,stickProb,simName)
    s += 'positions = '+ str(positions)
    s += '\ntimes = '+ str(times)
    f = open(filename, 'w')
    f.write(s)
    f.close()

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
    # Set up plot
    fig = plt.figure(facecolor='w')
    plt.subplot(2, 1, 1)
    plt.plot(times, nStuck,'-')
    plt.ylabel('Particles Stuck')
    plt.subplot(2, 1, 2)
    plt.plot(times, maxRadii, '-')
    plt.xlabel('Time Steps')
    plt.ylabel('Maximum Radius')
    # Save figure
    fileName = simName+'_anl.png'
    plt.savefig(fileName, dpi=300, bbox_inches='tight')
    
#==============================================================================
#  Main Program
#==============================================================================
# Get inputs from user (with defaults)
length = raw_input('Box side length: ')
if not length: length = 100
else: length = int(length)
print length

while True: 
    rSpawn = raw_input('Spawn radius: ')
    if not rSpawn:
        rSpawn = 45
        break
    elif int(rSpawn) > length/2:
        print "Spawn radius too large."
    else:
        rSpawn = int(rSpawn)
        break
print rSpawn

dt = 1 # time step (arbitrary units)

dr = raw_input('Distance step: ')
if not dr: dr = 1
else: dr = float(dr)
print dr

colRad = raw_input('Particle radius: ')
if not colRad: colRad = 1
else: colRad = float(colRad)
print colRad

stickProb = raw_input('Sticking probability (0-1): ')
if not stickProb: stickProb = 1
else: stickProb = float(stickProb)
print stickProb

numParticles = raw_input('Number of particles: ')
if not numParticles: numParticles = 100
else: numParticles = int(numParticles)
print numParticles

simName = raw_input('Simulation name: ')
if not simName: simName = 'test'
print simName

# Set up the simulation and run it
print 'Running simulation...'
start_time = time.time()
sim = Simulation(length, colRad, rSpawn, dr, stickProb)
sim.run(numParticles)
runTime = time.time() - start_time # the runtime in seconds

# Obtain Particle Information
particles = sim.getParticles()
positions = [p.position for p in particles]
radii = [p.radius for p in particles]
times = [p.time for p in particles]

# Save the data to file
saveData(positions,times,filename=simName+'_data')

# Plot the data
print 'Simulation done in %d seconds. Plotting...'%runTime
plotParticles(positions,times)
print 'Done.'

# Analyze data and plot
analyzeData(positions,times)