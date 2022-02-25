#!/usr/bin/env ipython
#
# nosehoover.py
#
# Copyright (c) 2020 Carlos Braga
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the MIT License.
#
# See accompanying LICENSE.md or https://opensource.org/licenses/MIT.
#

import sys
import numpy
import matplotlib
from matplotlib import pyplot
from mpl_toolkits.mplot3d import axes3d
from scipy.integrate import odeint, ode
from PIL import Image

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
        self.ortho_u = numpy.array([1.0, 0.0, 0.0])
        self.ortho_v = numpy.array([0.0, 1.0, 0.0])
        self.ortho_w = numpy.array([0.0, 0.0, 1.0])

        eps = 1.0E-6
        c = self.normal.dot(self.ortho_w)
        if (c > -1.0):
            n = numpy.cross(self.ortho_w, self.normal)
            K = numpy.array([[  0.0, -n[2],  n[1]],
                          [ n[2],   0.0, -n[0]],
                          [-n[1],  n[0],   0.0]])
            R = numpy.eye(3) + K + K.dot(K) / (1.0 + c)
            self.ortho_u = R.dot(self.ortho_u)
            self.ortho_v = R.dot(self.ortho_v)
            self.ortho_w = R.dot(self.ortho_w)
        else:
            self.ortho_u = -self.ortho_u
            self.ortho_v =  self.ortho_v
            self.ortho_w = -self.ortho_w

    def __del__(self):
        pass

    ##
    # Normalize a vector
    def normalize(self, v):
        norm = numpy.linalg.norm(v)
        if norm == 0:
            return v
        return v / norm

    ##
    # Return the intersection point between a ray and the plane
    def intersect(self, ray_orig, ray_dir):
        alpha = self.normal.dot(self.centre-ray_orig) / self.normal.dot(ray_dir)
        r_isect = ray_orig + alpha * ray_dir
        if alpha < 0.0:
            print(alpha, ray_orig, ray_dir, r_isect)
        return r_isect

    ##
    # Compute the vector projection onto the plane normal
    def project(self, pos):
        return self.normal.dot(pos)

    ##
    # Return side indicator of a given point wrt plane normal:
    #   s = -1 if pos is on the interior side
    #   s =  1 if pos is on the exterior side
    #   s =  0 if pos is on the plane
    def side(self, pos):
        s = self.project(pos - self.centre)
        return -1 if (s < 0.0) else 1 if (s > 0.0) else 0

    ##
    # generate
    # Generate n points in a plane region of size width and height
    def generate(self, n_points, width, height):
        half_width = 0.5 * width
        half_height = 0.5 * height
        step_width = width / float(n_points - 1)
        step_height = height / float(n_points - 1)
        return numpy.array([
            self.centre +
            self.ortho_u * (-half_width + i * step_width) +
            self.ortho_v * (-half_height + j * step_height)
            for i in range(n_points)
            for j in range(n_points)])

    ##
    # coord
    # Return the
    def coord(self, pos, scale, width, height):
        return numpy.array([
            0.5 + 0.5 * self.ortho_u.dot(scale * pos - self.centre) / width,
            0.5 + 0.5 * self.ortho_v.dot(scale * pos - self.centre) / height])



## ----------------------------------------------------------------------------
# NoseHoover
class NoseHoover:
    """ Class doc """
    ##
    # Constructor/destructor
    def __init__(self, gamma = 0.0):
        self.gamma = gamma

    def __del__(self):
        pass

    ##
    # Compute the Nose Hooover EoM
    def deriv(self, x):
        dxdt = numpy.zeros(numpy.size(x))
        dxdt[0] =  x[1]
        dxdt[1] = -x[0] - (x[2] * x[1])
        dxdt[2] =  (x[1] * x[1]) - (1.0 + self.gamma * numpy.tanh(x[0]))
        return dxdt

    ##
    # Perform an integration step
    def step(self, x_0, t_step, max_err = 1.0E-12, max_iter = 12):
        err = numpy.finfo(float).max
        z_val = 0.0
        n_iter = 0

        while (err > max_err and n_iter < max_iter):
            dx_dt = self.deriv(x_0 + 0.5 * z_val)
            z_new = t_step * dx_dt

            err = numpy.linalg.norm(z_new - z_val)
            z_val = z_new
            n_iter += 1

        x_1 = x_0 + z_val
        return x_1

## ----------------------------------------------------------------------------
# Poincare
class Poincare:
    """ Class doc """
    ##
    # Constructor/destructor
    def __init__(self, centre, normal, width = 1.0, height = 1.0):
        self.plane = Plane(centre, normal)
        self.width = width
        self.height = height
        self.x_points = None
        self.x_isect = []
    def __del__(self):
        pass

    ##
    # generate
    # Generate a set of initial points on the a region in the Poincare plane
    def generate(self, n_points):
        self.x_points = self.plane.generate(n_points, self.width, self.height)

    ##
    # intersect
    # Return 1 if both points are on the same side of the plane.
    # Return -1 if they are on opposite sides.
    # Return 0 if at least one point is on the plane
    def intersect(self, x_0, x_1):
        s_0 = self.plane.side(x_0)
        s_1 = self.plane.side(x_1)
        return s_0 * s_1

    ##
    # bisect
    # Compute the intersection point of a trajectory using a bisection method.
    def bisect(self, model, x_0, t_step, max_err = 1.0E-3, max_iter = 5):
        s_0 = self.plane.side(x_0)

        t_lo = 0.0
        t_hi = t_step
        n_iter = 0

        while n_iter < max_iter:
            t = 0.5 * (t_lo + t_hi)
            x_1 = model.step(x_0, t)

            if abs(self.plane.project(x_1)) < max_err:
                break

            s_1 = self.plane.side(x_1)
            if s_0 * s_1 < 0:
                t_hi = t
            else:
                t_lo = t

            n_iter += 1

        self.x_isect.append(x_1)

    ##
    # Generate an image from the Poincare intersection points.
    def image(self, scale = 1.0, img_width = 512, img_height = 512):
        # Scan the array and compute the intersection points.
        palette = self.denormalize(matplotlib.cm.get_cmap("turbo").colors)

        image = Image.new(
            mode="RGB", size=(img_width, img_height), color=(24,24,24))
        for ix,x in enumerate(self.x_isect):
            coord = self.plane.coord(x, scale, self.width, self.height)

            x = int(img_width * coord[0])
            y = int(img_height * coord[1])
            x = max(0, min(x, img_width - 1))
            y = max(0, min(y, img_height - 1))

            color = palette[ix % len(palette)]
            # r = (255,   0, 0)
            # g = (  0, 255, 0)
            # image.putpixel((x,y), r if self.plane.side(x) < 0 else g)
            image.putpixel((x,y), color)
        return image

    def denormalize(self, palette):
        return [
            tuple(int(channel * 255) for channel in color)
            for color in palette]

## ----------------------------------------------------------------------------
# main function
def main (argv):
    """ Function doc """
    ##
    # Create a poincare map and generate a collection of initial points.
    centre = numpy.array([0.0, 0.0, 0.0])
    normal = numpy.array([1.0, 0.0, 0.0])
    poincare = Poincare(centre, normal, width = 1.0, height = 1.0)
    poincare.generate(n_points = 10)

    ##
    # Compute the Poincare section of a Nose Hoover trajectory starting at
    # each initial point.
    t_step = 0.01
    n_steps = 1000
    nosehoover = NoseHoover(gamma = 0.25)
    for i,x_0 in enumerate(poincare.x_points):
        for n in range(n_steps):
            x_1 = nosehoover.step(x_0, t_step);
            if poincare.intersect(x_0, x_1):
                poincare.bisect(nosehoover, x_0, t_step)
            x_0 = x_1
        poincare.x_points[i] = x_0

        print(i, ": ", numpy.size(poincare.x_isect))
        image = poincare.image(
            scale = 0.2, img_width = 2048, img_height = 2048)
        image.save("out.pdf")

    # Done
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv))
