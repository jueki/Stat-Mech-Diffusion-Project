#!/usr/bin/env python

#==============================================================================
#  Max Byers, Jon Ueki, Jonas Kaufman
#  November 27 2015
#  Program to simulate diffusion-limited aggregation in three dimensions.
#==============================================================================

import numpy as np
from matplotlib import pyplot as plt
#from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D

def fromSphericalCoordinates(r,theta,phi):
	""" converts from spherical to cartesian coordinates """
	x = r*np.sin(theta)*np.cos(phi)
	y = r*np.sin(theta)*np.sin(phi)
	z = r*np.cos(theta)
	return [x,y,z]

def randomPosition(r):
	""" gives a random position on a sphere of radius r """
	theta = np.random.uniform(0,np.pi,1)[0]
	phi = np.random.uniform(0,2*np.pi,1)[0]
	return fromSphericalCoordinates [r,theta,phi]

class Particle:
	""" holds information for a single particle """

	def __init__(self,position=[0,0,0],radius=1.0,free=True):
		self.position = position
		self.radius = radius
		self.free = free
		#self.color = color

	def __repr__(self):
		s = "%s %s %s"%(self.position,self.radius,self.free)
		return s

	def fix(self):
		self.free = False

	def move(self,newPosition):
		self.position = newPosition

	def collided(self,other):
		""" checks if particle has collided with another particle """
		collisionRadius = self.radius + other.radius
		distance = np.linalg.norm(np.array(self.position) - np.array(other.position))
		if distance > collisionRadius:
			return False
		else:
			return True




