#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 16:10:18 2019

@author: Jennings Zhang <jenni_zh@protonmail.com>
"""

from os.path import isfile, join
import numpy as np
from math import acos
import subprocess
from tempfile import TemporaryDirectory
import argparse


def run(command_array):
    return subprocess.run(command_array, stdout=subprocess.DEVNULL, check=True)


def depth_potential_normals(input_filename, buffer, dtype=np.float32):
    run(['depth_potential', '-normals', input_filename, buffer])
    return np.loadtxt(buffer, dtype=dtype)


def angle_between(vectors):
    return np.arccos(np.clip(np.dot(*vectors), -1.0, 1.0))

def ideal_angles(inner_filename, outer_filename, tmpdir):
    inner = depth_potential_normals(inner_filename, join(tmpdir, 'inner.txt'))
    outer = depth_potential_normals(outer_filename, join(tmpdir, 'outer.txt'))
    # angle between two normal vectors
    return [angle_between(vectors) for vectors in zip(inner, outer)]


if __name__ == '__main__':
    # perfectionist: wasted temporary directory if user is just running --help
    with TemporaryDirectory() as tmpdir:
        ap = argparse.ArgumentParser(description='Runner for surface_angles '
                                     + 'which then recalculates values to '
                                     + 'represent the values as radians.')
        def input_file(filename):
            if isfile(filename):
                return filename
            else:
                msg = 'Required input file "' + filename + '" does not exist.'
                ap.error(msg)
        ap.add_argument('-mid', metavar='mid.obj',
                        default=join(tmpdir, 'middle_surface.obj'),
                        help='save the mid surface created by average_objects')
        ap.add_argument('-ideal', action='store_true',
                        help='Calculate angles between corresponding normals.')
        ap.add_argument('-error', action='store_true',
                        help='Subtract distortion from the angle between '
                        + 'corresponding normal vectors.')
        ap.add_argument('-regular', metavar='t.txt',
                        help='multiply output with normalized thickness.'
                        + '(Distortion of points which don\'t move far '
                        + 'is not important to consider, this produces '
                        + 'a more reasonable map for visualization.)')
        ap.add_argument('-view', action='store_true',
                        help='open brain-view')
        ap.add_argument('inner', metavar='inner.obj', type=input_file)
        ap.add_argument('outer', metavar='outer.obj', type=input_file)
        ap.add_argument('output', metavar='angles.txt', type=str)
        args = ap.parse_args()


        run(['average_objects', args.mid, args.inner, args.outer])
        run(['adapt_object_mesh', args.mid, args.mid, '0', '10', '0', '0'])
        if args.ideal:
            angles = ideal_angles(args.inner, args.outer, tmpdir)
        else:
            run(['surface_angles', args.inner, args.mid, args.outer, args.output])

            invert = np.vectorize(lambda a: acos(1.0 / a))

            angles = np.loadtxt(args.output, dtype=np.float64)
            angles = invert(angles)

            if args.error:
                #angles *= ideal_angles(args.inner, args.outer, tmpdir)
                raise NotImplementedError('gotta think about this')

        if args.regular:
            thickness = np.loadtxt(args.regular, dtype=np.float32)
            # normalize values between 0 and 1
            thickness -= thickness.min()
            thickness /= thickness.max()
            angles *= thickness

    np.savetxt(args.output, angles, fmt='%f')

    if args.view:
        obj = args.mid if isfile(args.mid) else args.outer
        command = 'brain-view {} {}'.format(obj, args.output)
        print(command)
        subprocess.run(command + '&', shell=True)
