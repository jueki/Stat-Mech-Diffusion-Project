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
  
  
    def position(self):
        """Returns the (x,y,z) position of the particle"""
        return self.position

class Grid:
    """The 3d grid that holds the position of each partice"""
    def __init__(self, length):
        """A 3D cube with each length^3 many cells"""
        self.cell = [[[]*length for x in range(length)] for y in range(length)]
        self.length = length
        
        
    def CTG(position):
        """Takes the continous position and transforms it into
            a position in the grid"""
        xGrid = int(position[0]) + (self.length/2)
        yGrid = int(position[1]) + (self.length/2)
        zGrid = int(position[2]) + (self.length/2)
        
        return [xGrid, yGrid, zGrid]
        
        
        
    
    def getParticles(self, x, y, z):
        """Gets all of the particles in the cell at position (x, y, z)"""
        return self.cell[x][y][z]
    
    def push(thisParticle):
        """Push a particle into a cell in the grid"""
        gridPos = CTG(thisParticle.position)
        xPos = gridPos[0]
        yPos = gridPos[1]
        zPos = gridPos[2]
        
        self.cell[xPos][yPos][zPos].append(thisParticle)
        


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


rSpawn = 3 #spawn radius
dt = 50 # time step in milliseconds
dr = 0.25 # distance step

time = 0 #starting time

particles = [Particle()] #initial config, one particle fixed at origin
positions = [p.position for p in particles]
radii = [p.radius for p in particles]
colors = [p.color for p in particles]

xs = [p[0] for p in positions]
ys = [p[1] for p in positions]
zs = [p[2] for p in positions]
ax.scatter(xs,ys,zs,s=radii) # plot initial configuration

pFree = Particle(randomPosition(rSpawn),radius=1,free=True) # spawn particle 
pos = pFree.position
scat = ax.scatter([pos[0]],[pos[1]],[pos[2]],s=10*radii,c=colors) #plot it

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