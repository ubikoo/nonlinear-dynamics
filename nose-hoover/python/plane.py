#!/usr/bin/env ipython

#
# untitled.py
#
# Copyright (c) 2020 Carlos Braga
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the MIT License.
#
# See accompanying LICENSE.md or https://opensource.org/licenses/MIT.
#

import sys
import atexit
import logging

import numpy as np
from matplotlib import pyplot, cm
from mpl_toolkits.mplot3d import axes3d
from scipy.integrate import odeint, ode
from PIL import Image

#sys.path.append('/Users/carlos/usr/local/example/lib')
#import example

## ----------------------------------------------------------------------------
# Plane
class Plane:
    """ Class doc """
    ##
    # Constructor/destructor
    def __init__(self, centre, normal):
        self.centre = centre
        self.normal = self.normalize(normal)

        # Compute the orthonormal basis set
        self.o_u = np.array([1.0, 0.0, 0.0])
        self.o_v = np.array([0.0, 1.0, 0.0])
        self.o_w = np.array([0.0, 0.0, 1.0])

        eps = 1.0E-6
        c = self.normal.dot(self.o_w)
        if (c > -1.0):
            n = np.cross(self.o_w, self.normal)
            K = np.array([[  0.0, -n[2],  n[1]],
                          [ n[2],   0.0, -n[0]],
                          [-n[1],  n[0],   0.0]])
            R = np.eye(3) + K + K.dot(K) / (1.0 + c)
            self.o_u = R.dot(self.o_u)
            self.o_v = R.dot(self.o_v)
            self.o_w = R.dot(self.o_w)
        else:
            self.o_u = -self.o_u
            self.o_v =  self.o_v
            self.o_w = -self.o_w

    def __del__(self):
        pass

    ##
    # Normalize a vector
    def normalize(self, v):
        norm = np.linalg.norm(v)
        if norm == 0:
            return v
        return v / norm

    ##
    # Return the intersection point between a ray and the plane
    def intersect(self, ray_orig, ray_dir):
        alpha = self.normal.dot(self.centre - ray_orig) / self.normal.dot(ray_dir)
        r_isect = ray_orig + alpha * ray_dir
        if alpha < 0.0:
            print(alpha, ray_orig, ray_dir, r_isect)
        return r_isect

    ##
    # Return side indicator of a given point wrt plane normal:
    #   s = -1 if pos is on the interior side
    #   s =  1 if pos is on the exterior side
    #   s =  0 if pos is on the plane
    def side(self, pos):
        s = self.normal.dot(pos - self.centre)
        return -1 if (s < 0.0) else 1 if (s > 0.0) else 0

    ##
    # generate
    # Generate n points in a plane region of size width and height
    def generate(self, n_points, width, height):
        """ Function doc """
        half_width = 0.5 * width
        half_height = 0.5 * height
        step_width = width / float(n_points - 1)
        step_height = height / float(n_points - 1)
        return np.array([
            self.centre +
            self.o_u * (-half_width + i * step_width) +
            self.o_v * (-half_height + j * step_height)
            for i in range(n_points)
            for j in range(n_points)])

## ----------------------------------------------------------------------------
# draw function
def draw(r_isect):
    """ Function doc """
    fig = pyplot.figure(figsize=(7,6))
    ax = fig.add_subplot(111, projection='3d')

    # set axes limits
    rmin = np.array([-2, -2, -2])
    rmax = np.array([ 2,  2,  2])
    # rmin = np.min(r_isect,0)
    # rmax = np.max(r_isect,0)
    xmin, xmax = rmin[0], rmax[0]
    ymin, ymax = rmin[1], rmax[1]
    zmin, zmax = rmin[2], rmax[2]
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_zlim(zmin, zmax)

    tics_maj = 1
    tics_min = 1
    ax.set_xticks(np.linspace(xmin, xmax, tics_maj+1))
    ax.set_xticks(np.linspace(xmin, xmax, tics_min+1), minor=True)
    ax.set_yticks(np.linspace(ymin, ymax, tics_maj+1))
    ax.set_yticks(np.linspace(ymin, ymax, tics_min+1), minor=True)
    ax.set_zticks(np.linspace(zmin, zmax, tics_maj+1))
    ax.set_zticks(np.linspace(zmin, zmax, tics_min+1), minor=True)

    # ax.grid(which='both')
    # ax.grid(which='minor', alpha=0.2)
    # ax.grid(which='major', alpha=0.8)
    ax.grid(None)

    # set axes labels
    # ax.legend(frameon=False, fontsize="small")
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    # set axes view angle
    #ax.view_init(elev=20., azim=-35)

    # plot data
    #ax.scatter(x, y, zs=0, zdir='z', label='(x,y)')
    #ax.scatter(x, y, zs=0, zdir='y', label='(x,z)')
    #ax.scatter(x, y, zs=0, zdir='x', label='(y,z)')
    ax.scatter(r_isect[:,0], r_isect[:,1], r_isect[:,2], s=1.0, marker='.', label='(x,y,z)')

    pyplot.show()
    # pyplot.savefig('fig.pdf', bbox_inches='tight')

## ----------------------------------------------------------------------------
# main function
def main (argv):
    """ Function doc """
    ##
    # Create a plane
    centre = np.array([0.0, 0.0, 0.0])
    normal = np.array([0.0, 0.0, -1.0])
    plane = Plane(centre, normal)

    ##
    # Intersect the plane with a set of rays
    n_rays = 1000;
    ray_orig = centre + 2.0 * normal
    ray_dir = np.array([-normal + 0.1*np.random.randn(3) for i in range(n_rays)])
    for i,item in enumerate(ray_dir):
        ray_dir[i] = item / np.linalg.norm(item)

    r_isect = np.zeros(np.shape(ray_dir))
    for i,dir in enumerate(ray_dir):
        r_isect[i] = plane.intersect(ray_orig, dir)

    ##
    # plot the interesection points
    draw(r_isect)
    draw(plane.generate(100, 10, 10))

    # Done
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv))
