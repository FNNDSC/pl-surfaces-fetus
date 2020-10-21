#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 16:14:27 2019

diencephalon is found where the intermediate zone (label=4) is in direct
contact with the CSF or background (label<1.5)
If the input MINC file has been artificially corrected to have a one-voxel
layer of label 3 covering the intermediate zone at this border,
we would have to use a different strategy.
Either the number of dilations of the intermediate zone should be increased,
or the binary mask of the intermediate zone should be transformed linearly
depending on side to exceed the boundary.

MINC tools must come from the CIVET-2.1.0 quarantine.
mincinfo -version
program: 2.1.02
libminc: 2.1.02
netcdf : 3.6.1
HDF5   : 1.8.8

Vlad's minc-toolkit (minc-toolkit-1.9.16-20180117-Ubuntu_16.04-x86_64) is buggy
and dilate_volume will not work.

@author: Jennings Zhang <jenni_zh@protonmail.com>
"""

import numpy as np
import argparse
from os.path import join, isfile
from tempfile import TemporaryDirectory
import subprocess


def run(command_array):
    return subprocess.run(command_array, stdout=subprocess.DEVNULL)


minccalc = ['minccalc', '-quiet', '-clobber', '-byte', '-unsigned', '-expr']


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Marks points on the given '
                                 + 'mesh to exclude regions where the '
                                 + 'subplate (SP) is discontinuous.')
    def input_file(filename):
        if isfile(filename):
            return filename
        else:
            ap.error(f'Required input file "{filename}" does not exist.')
    ap.add_argument('-keep', metavar='die.mnc', type=str,
                    help='Keep the intermediate file which highlights '
                    + 'the diencephalonic surface in voxel space.')
    ap.add_argument('-view', action='store_true',
                    help='Open the result in Display and brain-view')
    ap.add_argument('-boundary', metavar='label', type=int, default=1,
                    help='Possible values: 1 or 2 (default). Represents the '
                    'outer label which neighbors the IZ at label=4 over '
                    'regions that should be excluded. label=1 (CSF) is the '
                    'more strict option which usually covers the diencephalon.'
                    ' This can be interpreted as where the SP is '
                    'discontinuous. Setting label=2 (CP) will select a larger '
                    'area to be masked out, where the SP thickness is '
                    'expected to be 0 or close to 0.')
    ap.add_argument('-invert', action='store_true',
                    help='Flip 0/1 to create a negative mask so that '
                    'the lateral surface is painted 0.')
    ap.add_argument('labels', metavar='labels.mnc', type=input_file,
                    help='Painted segmentation volume with labels 1-6')
    ap.add_argument('obj', metavar='surface.obj', type=input_file,
                    help='Extracted mesh for subplate surface.')
    ap.add_argument('output', metavar='sp_mask.txt', type=str,
                    help='Mask in vertex space to be created for the mesh.')
    
    args = ap.parse_args()
    
    with TemporaryDirectory() as tmpdir:
        labels = args.labels
        iz = join(tmpdir, 'iz.mnc')
        pre = join(tmpdir, 'die.mnc')
        die = args.keep if args.keep else join(tmpdir, 'die_final.mnc')
        run(minccalc + ['A[0]>3.5',  labels, iz])
        run(['dilate_volume', iz, iz, '1', '6', '1'])

        boundary = str(args.boundary + 0.5)
        run(minccalc + ['A[0]>0.5&&A[1]<' + boundary, iz, labels, pre])
        # max_connect is passed as the third argument to mincdefrag.
        # larger values can compensate detections that happen on the 
        # lateral surface that result from messy segmentation
        # "garbage in, garbage out." --Claude
        # if these values don't work for a brain, just do this all by hand
        run(['mincdefrag', pre, pre, '1', '1', '14'])
        run(['dilate_volume', pre, pre, '1', '26', '1'])
        run(['dilate_volume', pre, pre, '1', '6', '1'])
        # just in case, this shouldn't do anything but it would remove
        # parts that cover the subplate zone
        run(minccalc + ['A[0]>0.5&&(A[1]<2.5||A[1]>3.5)', pre, labels, die])
        
        run(['volume_object_evaluate', '-linear', die, args.obj, args.output])
        
    # volume_object_evaluate produces inexact values around the border
    # list comprehension fixes values as 0 or 1 depending on a threshold
    array = np.loadtxt(args.output, dtype=np.float16)
    k, b = (0, 1,) if args.invert else (1, 0,)
    array = [k if a < 0.6 else b for a in array]
    np.savetxt(args.output, array, fmt='%i')

    if args.view:
        command = 'brain-view {} {}'.format(args.obj, args.output)
        print(command)
        subprocess.run(command + '&', shell=True)
        if args.keep:  # die.mnc shouldn't exist anymore unless it's kept
            command = 'Display -spectral ' + args.labels
            command += ' -label ' + args.keep
            command += ' ' + args.obj
            print(command)
            subprocess.run(command + '&', shell=True)

