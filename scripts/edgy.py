#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 13:51:24 2019

@author: Jennings Zhang <jenni_zh@protonmail.com>
"""

import sys
import pybicpl as bicpl
import numpy as np


def edgy(filename):
    obj = bicpl.MniObj(filename)
    # can't figure out the output order of surface-stats
    #data = check_output(['surface-stats', '-edge_length', filename])
    #data = np.float32(data.split())
    def average_length(point):
        index, neighbors = point
        coord = obj.points[index]
        l = [np.linalg.norm(coord - obj.points[n]) for n in neighbors]
        return np.mean(l)
    return [average_length(p) for p in enumerate(obj.neighbor_graph())]


if __name__ == '__main__':
    
    usage = 'usage: ' + sys.argv[0] + ' surface.obj edges.txt'
    
    if len(sys.argv) > 1 and '-h' in sys.argv[1]:
        print(usage)
        print()
        print('Average edge length at every vertex.')

        sys.exit(0)
    
    elif len(sys.argv) < 3:
        print(usage)
        print('error: missing filenames')
        sys.exit(1)

    bicpl.write_file(sys.argv[2], edgy(sys.argv[1]))
