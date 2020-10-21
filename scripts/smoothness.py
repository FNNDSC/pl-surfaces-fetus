#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 15:19:52 2019

@author: Jennings Zhang <jenni_zh@protonmail.com>
"""

import sys
import pybicpl as bicpl
import numpy as np

def smoothness(filename):
    data = bicpl.depth_potential(filename, '-mean_curvature')
    obj = bicpl.MniObj(filename)
    return bicpl.difference_average(obj.neighbor_graph(), data)


if __name__ == '__main__':

    usage = 'usage: ' + sys.argv[0] + ' surface.obj smoothness.txt'

    if len(sys.argv) > 1 and '-h' in sys.argv[1]:
        print(usage)
        print()
        print('Local smoothness can be defined as having a small change in')
        print('curvature between any vertex and its neighbors.')
        print('This program quantifies smoothness errors (or bumpiness?)')
        print('by calculating the average change in mean curvature')
        print('for every vertex.')
        print('A surface can be said to be "smooth" when its smoothness')
        print('values are close to 0.')
        sys.exit(0)

    elif len(sys.argv) < 2:
        print(usage)
        print('error: missing filenames')
        sys.exit(1)

    result = smoothness(sys.argv[1])
    if len(sys.argv) > 2:
        bicpl.write_file(sys.argv[2], result)
    else:
        print(round(np.mean(np.fromiter(result, float)), 3))