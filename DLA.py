#!/usr/bin/env python

#==============================================================================
#  Max Byers
#  Jon Ueki
#  Jonas Kaufman
#  November 27 2015
#  Program to simulate diffusion-limited aggregation in three dimensions.
#==============================================================================

import numpy as np

def fromSphericalCoordinates(r,theta,phi):
	""" converts from spherical to cartesian coordinates """
	x = r*np.sin(theta)*np.cos(phi)
	y = r*np.sin(theta)*np.sin(phi)
	z = r*np.cos(theta)
	return (x,y,z)

def randomPosition(r):
	""" gives a random position on a sphere of radius r """
	theta = np.random.uniform(0,np.pi,1)[0]
	phi = np.random.uniform(0,2*np.pi,1)[0]
	return fromSphericalCoordinates(r,theta,phi)

print randomPosition(1)