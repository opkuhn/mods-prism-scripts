#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description = '')

parser.add_argument('fileroot', type=str, help = 'fileroot')
parser.add_argument('ext', type=int, help = 'extension: 1=obj,2=sig,3=sky')
parser.add_argument('scl', type=float, help = 'scale factor')


args = parser.parse_args()

fileroot = args.fileroot
ext = args.ext
scl = args.scl

infile = fileroot + "_wav.dat"

if ext==1:
   ord = 5
if ext==2:
   ord = 8
if ext==3:
   ord = 11

if fileroot.find("m1b") != -1:
   x1=3200
   x2=6000
   clr = 'cornflowerblue'
elif fileroot.find("m1r") != -1:
   x1=5000
   x2=10500
   clr = 'darkred'

a = np.genfromtxt(infile)

ind = [i for i in range(len(a)) if a[i,1]>x1 and a[i,1]<x2]

buffer = 1.0 * np.std(scl*a[ind,5])

if ext==1 or ext==2:
   y1 = np.min(scl*a[ind,5]) - buffer
   y2 = np.max(scl*a[ind,5]) + buffer
elif ext==3:
   y1 = 0
   y2 = 50

plt.axis([x1,x2,y1,y2])

plt.plot(a[ind,1],scl*a[ind,ord],ls="solid",ds="steps-mid",color=clr)
