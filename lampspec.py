#!/usr/bin/env python

import numpy as np
import argparse
import matplotlib.pyplot as plt
import matplotlib
import scipy as sp
from scipy.interpolate import interp1d

parser = argparse.ArgumentParser(description = 'create line plots')

parser.add_argument('infile', type=str, help = 'ascii file (m1b_dual_lamps_otf.dat)')

args = parser.parse_args()

infile = args.infile

pngfile = infile[:infile.find(".dat")] + ".png"

lineids = "database/id" + infile[:infile.find(".dat")]

lamproot = infile.split("_")[2]
lamps = lamproot[:lamproot.find(".dat")]
instr = infile.split("_")[0]
dich = infile.split("_")[1]

if str.lower(dich) == "dual":
   mode = "Dual"
if str.lower(dich) == "red":
   mode = "Red"
if str.lower(dich) == "blue":
   mode = "Blue"
if instr == "m1r":
   inst = "MODS1R"
if instr == "m2r":
   inst = "MODS2R"
if instr == "m1b":
   inst = "MODS1B"
if instr == "m2b":
   inst = "MODS2B"
if str.lower(lamps) == "hgar":
   lamptitle = "Hg(Ar) lamp"
if str.lower(lamps) == "ArHg": 
   lamptitle = "Ar + Hg(Ar) lamps"
if lamps == "NeHg":
   lamptitle = "Ne + Hg(Ar) lamps"
if str.lower(lamps) == "ne":
   lamptitle = "Ne lamp"
if str.lower(lamps) == "kr":
   lamptitle = "Kr lamp"
if str.lower(lamps) == "xe":
   lamptitle = "Xe lamp"

a = np.genfromtxt(infile)
ltab = np.genfromtxt(lineids,skip_header=10)
#a = np.genfromtxt("m1r_dual_Xe.dat")
#ltab = np.genfromtxt("database/idm1r_dual_Xe",skip_header=10)

ip = [i for i in range(len(a)) if a[i,0] > 1750 and a[i,0] < 2450]

# interpolate line spectrum for plotting labels only
ia = sp.interpolate.interp1d(a[ip,0],a[ip,1],kind='linear')
ilog10a = sp.interpolate.interp1d(a[ip,0],np.log10(a[ip,1]),kind='linear')

log10a = np.log10(a[ip,1])

fig,ax = plt.subplots(1,1)

lwidth = 1.0
matplotlib.rcParams.update({'font.size':10})
matplotlib.rcParams.update({'axes.titlesize':'large'})
matplotlib.rcParams.update({'xtick.labelsize':'large'})
matplotlib.rcParams.update({'ytick.labelsize':'large'})
matplotlib.rc('axes',linewidth=lwidth)

#plt.rc('text', usetex=True)
#plt.rc('font', **{'family':'serif','serif':['Times'],'weight':'bold','size':'14'})
plt.rcParams['xtick.major.pad']='10'
plt.rcParams['ytick.major.pad']='10'
plt.rcParams['axes.labelpad'] = '10'


plotHeight = 3000 # plot width and height in pixels
plotWidth = 4000
dpi = 300 # dpi = 300 gives finer resolution and nicer plots.
fig.set_dpi(dpi)
wDisp = plotWidth
hDisp = plotHeight
wInches = float(wDisp)/float(dpi)
hInches = float(hDisp)/float(dpi)
fig.set_size_inches(wInches,hInches,forward=True)


frac  = 0.05
#off= frac*(np.max(a[ip,1]) - np.min(a[ip,1]))
off=  frac*(np.max(log10a) - np.min(log10a))
y1 =  np.min(log10a) - off
y2 =  np.max(log10a) + 4*off

ax.plot(a[ip,0],log10a,linestyle="solid",color="k")
#ax.axis([1750,2400,2,8])
ax.axis([1750,2400,y1,y2])

ax.set_xlabel("Pixel")
ax.set_ylabel("log10(Counts[ADU])")


for l,f in zip(ltab[:,1],ltab[:,2]):
    ax.annotate('%s' % f,xy=(l,ilog10a(l)+off),xycoords='data',xytext=(0,5),textcoords='offset points', ha='center',va='bottom',rotation=90,color="black",size='medium')

plt.title("%s %s %s" % (inst, mode, lamptitle))

plt.savefig(pngfile)
