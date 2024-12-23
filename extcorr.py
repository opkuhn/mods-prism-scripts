#!/usr/bin/env python

import argparse
import numpy as np
import scipy as sp
from scipy.interpolate import interp1d


parser = argparse.ArgumentParser(description = "Prism Extinction Correction ")
parser.add_argument("infile", type=str, help = "Name of input file - observations ")
parser.add_argument("airmass", type=str, help = "airmass")
args = parser.parse_args()

infile = args.infile # infile is wavelength-calibrated, flexure-corrected data
# 
# infile name must start with "m1b" or "m1r", it will be of the form m1r_2023YO1_A_cert_iavg_1s_wav.dat
#
#  (0)     (1)         (2)       (3)  (4)  (5)      (6 ) (7)   (8)        (9)  (10)  (11) 
# pixel wavelength disp[A/pix] Pix_obs 1 Cnts_obs Pix_obs 1 SigCnts_obs Pix_obs 1 SkyCnts_obs
#...
#2000 4444.46352 6.8999 2001 1 0.36532 2001 1 0.04948 2001 1 3.83307
#2001 4451.37919 6.9315 2002 1 0.27496 2002 1 0.04889 2002 1 3.96984
#2002 4458.32659 6.9634 2003 1 0.34175 2003 1 0.04970 2003 1 4.00730
#2003 4465.30594 6.9954 2004 1 0.30518 2004 1 0.04983 2004 1 4.10196
#2004 4472.31746 7.0277 2005 1 0.27718 2005 1 0.04950 2005 1 4.11982
#2005 4479.36138 7.0602 2006 1 0.26521 2006 1 0.04980 2006 1 4.19614
#2006 4486.43792 7.0929 2007 1 0.34195 2007 1 0.04993 2007 1 4.08614
#2007 4493.54732 7.1259 2008 1 0.36480 2008 1 0.05037 2008 1 4.14512
#2008 4500.68980 7.1591 2009 1 0.36093 2009 1 0.05018 2009 1 4.11015
#2009 4507.86562 7.1926 2010 1 0.31341 2010 1 0.04985 2010 1 4.16603
#2010 4515.07500 7.2263 2011 1 0.38784 2011 1 0.04999 2011 1 4.06039
#2011 4522.31819 7.2602 2012 1 0.34013 2012 1 0.04987 2012 1 4.12991
#2012 4529.59544 7.2944 2013 1 0.33333 2013 1 0.04968 2013 1 4.11050
#2013 4536.90699 7.3288 2014 1 0.36000 2014 1 0.04995 2014 1 4.11665
#...
#
airmass = float(args.airmass)

outfile = infile[:infile.find('.dat')] + '_ext.dat'

blue_wlim1=3550
blue_wlim2=5650
red_wlim1=5650
red_wlim2=10000

if infile.find('m1b') != -1:
   wlim1 = blue_wlim1
   wlim2 = blue_wlim2
if infile.find('m1r') != -1:
   wlim1 = red_wlim1
   wlim2 = red_wlim2

# Extinction Curve
#
#extinct = np.array([[3200,0.866],[3500,0.511],[4000,0.311],[4500,0.207],[5000,0.153],[5500,0.128],[6000,0.113],[6450,0.088],[6500,0.085],[7000,0.063],[7500,0.053],[8000,0.044],[8210,0.043],[8260,0.042],[8370,0.041],[8708,0.026],[10256,0.020]],float)
extinct = np.genfromtxt("/Users/olga/Dropbox/MODS/CalTables/fluxcal/NewTables/LBTO_AtmosExtModel_2024.txt",skip_header=15)

iext = sp.interpolate.interp1d(extinct[:,0],extinct[:,1],kind='linear')

spec = np.genfromtxt(infile)
speccorr = np.zeros_like(spec)

# iw covers only the pixels of interest. This is hardcoded and was chosen to be
# where the per pixel dispersion is increasing, the limits defined as blue,red_wlim1,2. 

iw = [i for i in range(len(spec)) if spec[i,1]<wlim2 and spec[i,1]>wlim1]

speccorr[:,0:5] = spec[:,0:5]
speccorr[:,6:8] = spec[:,6:8]
speccorr[:,9:11] = spec[:,9:11]
for i in iw:
   extmag = iext(spec[i,1])*(airmass-1)
   extflux = 10**(0.4*extmag)   # positive because we are removing the effect of extinction
   print ("%d %f %f" % (i,extmag,extflux))
   speccorr[i,5] = spec[i,5]*extflux
   speccorr[i,8] = spec[i,8]*extflux
   speccorr[i,11] = spec[i,11]*extflux

np.savetxt(outfile,speccorr,fmt="%d %.5f %.4f %d %d %.5f %d %d %.5f %d %d %.5f")
