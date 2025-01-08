#!/usr/bin/env python
# /opt/anaconda/bin/python
#
# This program takes the output (database/id*) of IRAF's identify
# task, which has 5 columns: pixel (selected), pixel (fit) and
# wavelength (also fwidth and 1/0 indicating, when a fit has been made
# in identify, whether the point has been used or not), 
# and it fits a polynomial to the data, 
# first to the wavelength as a function of pixel and 
# second to the pixel as a function of wavelength.
#
# The plots output show (at top) the datapoints, the fit and 
# (at bottom) the fit residuals. 
#
# The number of points used is indicated as well as the fit rms.
#
# update from prismwavecal.py:  
#
# update plot style: 2 panels share an axis 
#  and to use Tex style fonts 
# fit only one channel at a time 
# input the fit order 
# use argparse 
#
# 2020-05-19 opk

import os
from sys import argv, exit
import getopt
import numpy as np
import scipy as sp
from scipy.interpolate import interp1d
import math
import matplotlib.pyplot as pl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator) 
import argparse

# Version number and date

versNum = "0.0.1"
versDate = "2020-05-18"

#----------------------------------------------------------------
#
# Main Program starts here...
#
# Parse the command-line arguments

parser = argparse.ArgumentParser(description = "Prism wavelength fits")
parser.add_argument("infile", type=str, help = "Name of input file")
parser.add_argument("ord", type=int, help = "order of polynomial") 
parser.add_argument("inv", type=str, help = "no=lambda(pix) or yes=pix(lambda)?") 
parser.add_argument("niter", type=int, help = "# of iterations with sigma-clipping?") 
parser.add_argument("sigclip", type=float, help = "sigma for sigma-clipping (2.5)?") 
args = parser.parse_args()

f_chanid = args.infile
ord = args.ord
inv = args.inv
niter = args.niter
sig = args.sigclip

if f_chanid[:2] == "bd":
	config= "Blue Direct Prism"
	if inv == "yes":
		plfile = "BlueDirRev.png"
	if inv == "no":
		plfile = "BlueDir.png"
		pixwav = "BlueDir_wav.txt"
if f_chanid[:2] == "db":
	config= "Dual Blue Prism"
	if inv == "yes":
		plfile = "DualBluRev.png"
	if inv == "no":
		plfile = "DualBlue.png"
		pixwav = "DualBlue_wav.txt"
if f_chanid[:2] == "rd":
	config= "Red Direct Prism"
	if inv == "yes":
		plfile = "RedDirRev.png"
	if inv == "no":
		plfile = "RedDir.png"
		pixwav = "RedDir_wav.txt"
if f_chanid[:2] == "dr":
	config= "Dual Red Prism"
	if inv == "yes":
		plfile = "DualRedRev.png"
	if inv == "no":
		plfile = "DualRed.png"
		pixwav = "DualRed_wav.txt"

# Read the input file into chanid
# the first 3 columns of this file will have 
# pixel(selected), pixel(fit) and wavelength ID.
# use pixel(fit) instead of pixel(selected)
# multiple measures of the same line have been combined by combid.py

chanid = np.genfromtxt(f_chanid)


################## Polynomial fits #####################

###### lam(pix) for inv=no. Writes out a text file also
###### pix(lam) for inv=yes

#sig= 2.5

if (inv == "no"): # wavelength(pix)
  keep = range(len(chanid))
  for k in range(niter):
   cfit,cvar,crank,csv,crcond = np.polyfit(chanid[keep,1],chanid[keep,2],ord,full=True)
   xcp = np.poly1d(cfit)
   cres = chanid[keep,2] - xcp(chanid[keep,1])
   rms = np.sqrt(cvar/len(keep))
   if k < niter-1:
      keep = [i for i in range(len(cres)) if np.fabs(cres[i]/rms) < sig]

elif (inv == "yes"): # pix(wavlength)
  keep = range(len(chanid))
  for k in range(niter):
   cfit,cvar,crank,csv,crcond = np.polyfit(chanid[keep,2],chanid[keep,1],ord,full=True)
   xcp = np.poly1d(cfit)
   cres = chanid[keep,1] - xcp(chanid[keep,2])
   rms = np.sqrt(cvar/len(keep))
   if k < niter-1:
      keep = [i for i in range(len(cres)) if np.fabs(cres[i]/rms) < sig]

# Write the pix-to-wavelength files which will be used to 
# wavelength calibrate the data.
#
if (inv == "no"):
	xpix = np.linspace(0,4095,4096,endpoint=True)
	pcwav = np.zeros((len(xpix),3))
	pcwav[:,0] = xpix+1
	pcwav[:,1] = xcp(xpix+1)
	#pcwav[:,2] = xcp(xpix+0.5) - xcp(xpix-0.5)
	pcwav[:,2] = xcp(xpix+1.5) - xcp(xpix+0.5)
	np.savetxt(pixwav,pcwav,fmt="%.1f %.5f %.5f")

################## Make plots #####################

### Generate a smooth polynomial curve

# Set the limits over which to plot the best-fit polynomials 
#
c1 = int(np.rint(np.min(chanid[:,1]))-1)
c2 = int(np.rint(np.max(chanid[:,1]))+1)
nc = c2-c1
w1 = int(np.rint(np.min(chanid[:,2]))-1)
w2 = int(np.rint(np.max(chanid[:,2]))+1)
nw=w2-w1

if (inv == "no"):
	cx = np.linspace(c1,c2,nc,endpoint=True)

if (inv == "yes"):
	cx = np.linspace(w1,w2,nw,endpoint=True)

cmod=  xcp(cx)

### Set up the plot parameters

buf = 50

fig,ax = pl.subplots()
divider = make_axes_locatable(ax)
ax2 = divider.append_axes("bottom",size="20%",pad=0)
ax.figure.add_axes(ax2)

if (inv == "no"):
	ax.set_xlim([c1-buf,c2+buf])
	ax.set_ylim([w1-2*buf,w2+2*buf])
	ax.set_ylabel("Wavelength [Ang]")
	ax.set_xticklabels([])
	#ax.xaxis.set_major_locator(MultipleLocator(50))
	#ax.xaxis.set_minor_locator(MultipleLocator(10))
	ax.plot(chanid[keep,1],chanid[keep,2],marker='o',color='k',linestyle="none",fillstyle='none')
	ax.plot(chanid[keep,1],chanid[keep,2],marker='x',color='k',linestyle="none",fillstyle='none')
	ax.plot(cx,cmod,color='k',linestyle=":",label='%s: order = %d N = %d rms = %.2f Ang' % (config,ord, len(chanid[keep,:]),rms))
	ax.legend()

	ax2.xaxis.set_major_locator(MultipleLocator(100))
	ax2.xaxis.set_minor_locator(MultipleLocator(100))
	ax2.set_xlim([c1-buf,c2+buf])
	ax2.hlines(0,c1-2*buf,c2+2*buf,linestyle="dotted")
	lim = np.rint(max(np.fabs(min(cres)),np.fabs(max(cres)))) +.5
	ax2.set_ylim([-lim,lim])
	ax2.set_ylabel("res[Ang]")
	ax2.plot(chanid[keep,1],cres,marker='x',linestyle="none",color='k')
	ax2.set_xlabel("X [pix]")

if (inv == "yes"):
	ax.set_xlim([w1-2*buf,w2+2*buf])
	ax.set_ylim([c1-buf,c2+buf])
	ax.set_ylabel("Pix")
	ax.set_xticklabels([])
	ax.plot(chanid[keep,2],chanid[keep,1],marker='o',color='k',linestyle="none",fillstyle='none')
	ax.plot(chanid[keep,2],chanid[keep,1],marker='x',color='k',linestyle="none",fillstyle='none')
	ax.plot(cx,cmod,color='k',linestyle=":",label='%s: order = %d N = %d rms = %.2f pix' % (config,ord,len(chanid[keep,:]),rms))
	ax.legend()

	ax2.xaxis.set_major_locator(MultipleLocator(500))
	ax2.xaxis.set_minor_locator(MultipleLocator(100))
	ax2.set_xlim([w1-2*buf,w2+2*buf])
	ax2.hlines(0,w1-2*buf,w2+2*buf,linestyle="dotted")
	lim = np.rint(max(np.fabs(min(cres)),np.fabs(max(cres)))) +.5
	ax2.set_ylim([-lim,lim])
	ax2.set_ylabel("res[pix]")
	ax2.plot(chanid[keep,2],cres,marker='x',linestyle="none",color='k')
	ax2.set_xlabel("Wavelength [Ang]")
	
pl.savefig(plfile)

#pl.show()
