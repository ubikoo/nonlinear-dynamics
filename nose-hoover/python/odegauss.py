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
import numpy
import matplotlib
from matplotlib import pyplot
from mpl_toolkits.mplot3d import axes3d
from scipy.integrate import odeint, ode
from PIL import Image

## ----------------------------------------------------------------------------
# Kepler
class Kepler:
    """ Class doc """
    def __init__(self, delta):
        self.delta = delta

    def __del__(self):
        pass

    ##
    # Compute the Nose Hooover EoM
    def deriv(self, x):
        r = numpy.sqrt(x[0]*x[0] + x[1]*x[1])
        inv_r3 = 1.0 / (r*r*r)
        inv_r5 = 1.0 / (r*r*r*r*r)

        dxdt = numpy.zeros(numpy.size(x))
        dxdt[0] = x[2]
        dxdt[1] = x[3]
        dxdt[2] = -1.0*(inv_r3 + 3.0 * self.delta * inv_r5)*x[0]
        dxdt[3] = -1.0*(inv_r3 + 3.0 * self.delta * inv_r5)*x[1]
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

        print("n_iter: ", n_iter, " err: ", err)

        x_1 = x_0 + z_val
        return x_1

    ##
    # plot function
    @staticmethod
    def plot (x, y, title):
        """ Function doc """
        # plot data
        fig = pyplot.figure()
        ax = fig.add_subplot(1,1,1)

        xmin, xmax = numpy.min(x), numpy.max(x)
        ymin, ymax = numpy.min(y), numpy.max(y)

        tics_maj = 5
        tics_min = 5

        ax.set_xticks(numpy.linspace(xmin, xmax, tics_maj+1))
        ax.set_xticks(numpy.linspace(xmin, xmax, tics_min+1), minor=True)

        ax.set_yticks(numpy.linspace(ymin, ymax, tics_maj+1))
        ax.set_yticks(numpy.linspace(ymin, ymax, tics_min+1), minor=True)

        # ax.grid(which='both')
        ax.grid(which='minor', alpha=0.2)
        ax.grid(which='major', alpha=0.8)

        pyplot.axis([xmin, xmax, ymin, ymax])
        #pyplot.legend(frameon=False, fontsize="small")

        n_points = numpy.size(x)
        c = numpy.arange(n_points)
        pyplot.plot(x, y, '-', linewidth=0.1)
        ax.scatter(x, y, marker='.', s=1, c=c, cmap='RdBu')
        ax.scatter(x[0], y[0], marker='*', s=16, c='red')
        ax.scatter(x[-1:], y[-1:], marker='o', s=16, c='green')

        pyplot.title(title)
        pyplot.show()

# -----------------------------------------------------------------------------
# main function
#
def main (argv):
    """ Function doc """
    ##
    # Model data
    delta = 0.001
    e = 0.6
    x0 = numpy.array([1.0 - e, 0.0, 0.0, numpy.sqrt((1.0+e)/(1.0-e))])
    t_max = 100.0
    n_steps = int(abs(t_max)/0.01)
    t_val = numpy.concatenate([numpy.linspace(0, t_max, n_steps),
                            numpy.linspace(t_max, 0, n_steps)])
    x_val = numpy.zeros((numpy.size(t_val), numpy.size(x0)))
    x_val[0,:] = x0

    ##
    # Kepler model
    kepler = Kepler(delta)
    for i in range(numpy.size(t_val)-1):
        t_step = t_val[i+1] - t_val[i]
        x_val[i+1] = kepler.step(x_val[i], t_step)

    pyplot.close('all')
    kepler.plot(x_val[:,0], x_val[:,1], 'kepler')


if __name__ == '__main__':
    sys.exit(main(sys.argv))
