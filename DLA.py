#!/usr/bin/env python

#==============================================================================
#  Max Byers, Jon Ueki, Jonas Kaufman
#  November 27 2015
#  Program to simulate diffusion-limited aggregation in three dimensions.
#==============================================================================

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation


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

	def __init__(self,position=[0,0,0],radius=1,free=False,color=(0,0,1)):
		self.position = np.array(position)
		self.radius = radius
		self.free = free
		self.color = color

	def fix(self):
		self.free = False

	def move(self,newPosition):
		self.position = np.array(newPosition)

	def collided(self,other):
		""" checks if particle has collided with another particle """
		collisionRadius = self.radius + other.radius
		distance = np.linalg.norm(np.array(self.position) - np.array(other.position))
		if distance > collisionRadius:
			return False
		else:
			return True

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
    def __init__(self, length, colRad, spawnRad, stepSize):
        self.grid = Grid(length)
        self.allParticles = [] #A list of all of the particles existing
        self.stepSize = stepSize
        self.length = length
        self.colRad = colRad
        self.spawnRad = spawnRad
        seed = Particle()
        
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
        particle = Particle(randomPosition(self.spawnRad),radius=self.colRad,free=True)
        
        while(self.inBounds(particle) and not self.collision(particle)):
            particle.step(self.stepSize)
        
        # If the particle is still inbounds, that means it collided with the seed
        # So fix it into the grid and append it to the list of existing particles
        if (self.inBounds(particle)):
            self.allParticles.append(particle)
            self.grid.push(particle)
            particle.fix()
            
            
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
            

# setting up 3d figure
fig = plt.figure()
ax = p3.Axes3D(fig)
(bx,by,bz) = (10,10,10) # box.bounds
ax.set_xlim3d([-bx, bx])
ax.set_xlabel('x')

ax.set_ylim3d([-by, by])
ax.set_ylabel('y')

ax.set_zlim3d([-bz, bz])
ax.set_zlabel('z')

ax.set_title('DLA Test')


length = 20 #length of grid
rSpawn = 5 #spawn radius
dt = 50 # time step in milliseconds
dr = 0.25 # distance step
colRad = .4 #Collision radius
time = 0 #starting time
numParticles = 500

#Set Up the simulation and run it
sim = Simulation(length, colRad, rSpawn, dr)
sim.run(numParticles)

#Obtain Particle Information
particles = sim.getParticles()
positions = [p.position for p in particles]
radii = [p.radius for p in particles]
colors = [p.color for p in particles]

xs = [p[0] for p in positions]
ys = [p[1] for p in positions]
zs = [p[2] for p in positions]
ax.scatter(xs,ys,zs,s=radii) # plot initial configuration

pFree = Particle(randomPosition(rSpawn),radius=1,free=True) # spawn particle 
pos = pFree.position
scat = ax.scatter([pos[0]],[pos[1]],[pos[2]],s=100*radii,c=colors) #plot it

def animate(i):
	global time
	time += dt
	global scat
	global pFree
	global particles
	collide = False
	for p in particles: # loop through to check for collisions (change this)
		if pFree.collided(p):
			collide = True
	if collide:
		pFree.fix()
		particles += [pFree]
		pFree = Particle(randomPosition(rSpawn),radius=1,free=True) # spawn new particle
		pos = pFree.position
		r = pFree.radius
		xs = [pos[0]]
		ys = [pos[1]]
		zs = [pos[2]]
		scat = ax.scatter(xs,ys,zs,s=r) # plot new particle
	else:
		scat.remove() # clear previous plot
		pFree.step(dr)
		pos = pFree.position
		r = pFree.radius
		xs = [pos[0]]
		ys = [pos[1]]
		zs = [pos[2]]
		scat = ax.scatter(xs,ys,zs,s=r) # replot

# Creating the Animation object
ani = animation.FuncAnimation(fig, animate, interval=dt, blit=False)
plt.show()