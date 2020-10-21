#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 11:14:04 2019

Produce surface meshes for any parametric spherical function.
These are useful for modeling, because thickness between
mathematical curves can be determined theoretically.

Resulting meshes have perfect vertex-to-vertex correspondence
and are perfectly smooth (so long as the given function is differentiable).

These objects can be manipulated with param2xfm or
converted to MINC volumes with surface_mask2.

Uses physics (ISO) convention for spherical coordinates
https://en.wikipedia.org/wiki/Spherical_coordinate_system#/media/File:3D_Spherical_2.svg

@author: Jennings Zhang <jenni_zh@protonmail.com>
"""

import sys
import numpy as np
from numpy.linalg import norm
from math import *
import argparse
from pybicpl import MniObj

def euclidean_distance(*args):
    """
    Wrapper for numpy.linalg.norm
    """
    return norm(args)


def cart2sphere(x, y, z):
    """
    non-vectorized adaptation of dipy.core.geometry.cart2sphere
    """
    # doesn't matter for sinusoidal functions, but it would be nice to
    # change these angles to be in the range of
    # 0 <= theta < pi
    # 0 <= phi < 2*pi
    r = euclidean_distance(x, y, z)
    # center will be missing if it exists at whole numbers
    theta = acos(z / r)
    phi = atan2(y, x)
    return r, theta, phi


def sphere2cart(r, polar, azimuth):
    x = r * sin(polar) * cos(azimuth)
    y = r * sin(polar) * sin(azimuth)
    z = r * cos(polar)
    return x, y, z


def project_mesh(function, sampling_tensor):
    def proj(point):
        r, theta, phi = cart2sphere(*point)
        r = function(theta, phi)
        return sphere2cart(r, theta, phi)
    return np.array([proj(p) for p in sampling_tensor], dtype=np.float32)



def project_mesh_xyz(f, sampling_tensor):
    def proj(point):
        r, theta, phi = cart2sphere(*point)
        x = f[0](theta, phi)
        y = f[1](theta, phi)
        z = f[2](theta, phi)
        return x, y, z
    return np.array([proj(p) for p in sampling_tensor], dtype=np.float32)


def func(equation_string):
    return eval('lambda theta, phi: ' + equation_string)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Produce a surface from a spherical function.')
# dependency on create_tetra instead of requiring an input data file
#     ap.add_argument('-sphere', metavar='template.obj', default=unit_sphere,
#                     help='Pre-existing unit sphere to reference for which '
#                          'directional angles to sample. '
#                          '(Use inflate_to_sphere on the extracted surface.)')
# option deprecated because it's better to just use create_tetra and param2xfm
#    ap.add_argument('-ellipsoid', action='store_true',
#                    help='Interpret the equations as directional scalings '
#                         'and draw an ellipse (or sphere)')
    ap.add_argument('equations', nargs='+', type=func,
                    help='Supply one string in the form of r(polar, azimuth)=?'
                         ' OR three to represent x=?, y=?, z=?')
    ap.add_argument('output_filename', metavar='surface.obj')
    args = ap.parse_args()

    n_eq = len(args.equations)

#    if args.ellipsoid:
#        a = args.equations[0](None, None)
#        if n_eq > 2:
#            b = args.equations[1](None, None)
#            c = args.equations[2](None, None)
#        else:
#            b = a
#            c = a
#        args.equations = [
#            lambda theta, phi: a * sin(theta) * cos(phi),
#            lambda theta, phi: b * sin(theta) * sin(phi),
#            lambda theta, phi: c * cos(theta)
#        ]

    obj = MniObj()

    if n_eq == 1:
        obj.points = project_mesh(args.equations[0], obj.points)
    else:
        if n_eq != 3:
            print('error: expected 3 parametric functions, got {n_eq}')
            sys.exit(1)
        obj.points = project_mesh_xyz(args.equations, obj.points)

    obj.recompute_normals()
    obj.save(args.output_filename)
