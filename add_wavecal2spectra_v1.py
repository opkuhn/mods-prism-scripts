#!/usr/bin/env python
#

import numpy as np
import scipy as sp
from scipy.interpolate import interp1d
import matplotlib.pyplot as pl
import argparse

# add_wavecal2spectra.py will manipulate the output of listpix and add to it 
# the wav(pix) solution, shifted by the appropriate number of pixels. 
#
# listpix.cl: for each normalized spectrum, list the extracted object spectrum and associated sigma spectrum
# listpix ("myspec[*,*,1]", wcs="logical", formats="", verbose=no, >"myspec_obj.list")
# listpix ("myspec[*,*,4]", wcs="logical", formats="", verbose=no, >"myspec_sig.list")
# listpix ("myspec[*,*,3]", wcs="logical", formats="", verbose=no, >"myspec_sky.list")
# and these have 3 columns each: 
#
#  object: pixel, 1, flux1s(ADU)     
#  and
#  sigma: pixel, 1, sig_flux1s(ADU)
#
# DualBlue_wav.txt and DualRed_wav.txt each contain 3 columns: pixel, wavelength, delta_wav
      
# correct for flexure by using np.roll to shift object, sigma and sky spectra 
# by an integer number of pixels 
#
# np.roll(array,d,axis=0) where d<0 will shift the rows up, adding the first rows to the end
#                           and d>0 will shift the rows down, brigning the last rows up to the top

# inputs: (1) root name of object and sigma spectra 
#         (2) shift (pix)

#
# Main Program starts here...
#

# Parse the command-line arguments 

parser = argparse.ArgumentParser(description = 'wavelength calibrate prism spectra')
parser.add_argument('fileroot', type=str, help = 'root name of object, sky and sigma spectra')
parser.add_argument('flexshift', type=int, help = 'flexure shift (positive or negative pixels)')
parser.add_argument('--nosig', action="store_true", help = 'indicates no sigma spectrum is available')

args = parser.parse_args()

fileroot = args.fileroot
shift = args.flexshift
nosig = False

if args.nosig:
   # there is no sigma spectrum, create it by sqrt(obj)
   nosig = True

if fileroot.find('m1b') != -1:
   wavcal = "DualBlue_wav.txt"
if fileroot.find('m1r') != -1:
   wavcal = "DualRed_wav.txt"

objdat = fileroot + "_obj.list"
skydat = fileroot + "_sky.list"
outdat = fileroot + "_wav.dat"

if nosig:
   objarr = np.genfromtxt(objdat)
   skyarr = np.genfromtxt(skydat)
   wavarr = np.genfromtxt(wavcal)
if not nosig:
   objarr = np.genfromtxt(objdat)
   skyarr = np.genfromtxt(skydat)
   wavarr = np.genfromtxt(wavcal)

if not nosig:
   sigdat = fileroot + "_sig.list"
   sigarr = np.genfromtxt(sigdat)
else:
   sigarr = np.zeros_like(objarr)
   sigarr[:,0:2] = objarr[:,0:2]
   for i in range(len(objarr)):
      if objarr[i,2] > 0:
         sigarr[i,2] = np.sqrt(objarr[i,2])
      elif objarr[i,2] <= 0:
         sigarr[i,2] = 0

#newarr = np.hstack((np.roll(wavarr,-shift,axis=0),objarr,sigarr,skyarr))
newarr = np.hstack((wavarr,np.roll(objarr,shift,axis=0),np.roll(sigarr,shift,axis=0),np.roll(skyarr,shift,axis=0)))

np.savetxt(outdat,newarr,fmt="%d %10.5f %6.4f %d %d %.5f %d %d %.5f %d %d %.5f")
