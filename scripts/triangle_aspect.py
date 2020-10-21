#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 15:58:24 2019

@author: Jennings Zhang <jenni_zh@protonmail.com>
"""


import sys
import pybicpl as bicpl
import numpy as np


def aspect_ratio(t):
    a = np.linalg.norm(t[0] - t[1])
    b = np.linalg.norm(t[1] - t[2])
    c = np.linalg.norm(t[2] - t[0])
    s = (a + b + c) / 2
    return a * b * c / (8 * (s - a) * (s - b) * (s - c))

def triangle_aspect(filename):
    obj = bicpl.MniObj(filename)
    # can't figure out the output order of surface-stats
    #data = check_output(['surface-stats', '-edge_length', filename])
    #data = np.float32(data.split())
    triangles = np.reshape(obj.indices, (obj.n_items, 3))
    return [aspect_ratio([obj.points[i] for i in t]) for t in triangles]


if __name__ == '__main__':
    
    usage = 'usage: ' + sys.argv[0] + ' surface.obj [aspect_ratios.txt]'
    
    if len(sys.argv) > 1 and '-h' in sys.argv[1]:
        print(usage)
        print()
        print('Aspect ratio of triangles given by the')
        print('ratio of the circumradius to twice its inradius')
        print('=a*b*c/(8*(s-a)*(s-b)*(s-c)) where s=(a+b+c)/2')

        sys.exit(0)
    
    elif len(sys.argv) < 2:
        print(usage)
        print('error: missing filenames')
        sys.exit(1)

    result = triangle_aspect(sys.argv[1])
    
    q1 = np.quantile(result, 0.25)
    q3 = np.quantile(result, 0.75)
    iqr = q3 - q1
    mask = [a < q1 - 1.5 * iqr or a > q3 + 1.5 * iqr for a in result]
    m = np.ma.array(result, mask=mask)
    count = np.ma.count_masked(m)
    cent = count / len(result) * 100
    nums = (count, len(result), cent)
    
    print('all triangles, n={:d}'.format(len(result)))
    print('mean={:.3f}  std={:.3f}'.format(np.mean(result), np.std(result)))
    print('{} out of {} ({:.1f}%) outliers masked'.format(*nums))
    print('mean={:.3f}  std={:.3f}'.format(np.mean(m), np.std(m)))

    if len(sys.argv) > 2:    
        bicpl.write_file(sys.argv[2], result)
