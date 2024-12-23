import numpy as np
import argparse
import matplotlib.pyplot as plt
import scipy as sp
from scipy.interpolate import interp1d

parser = argparse.ArgumentParser(description = 'create line plots')

parser.add_argument('infile', type=str, help = 'infile')

args = parser.parse_args()

infile = args.infile

lineids = "database/id" + infile[:infile.find(".dat")]

lamps = infile.split("_")[2]
instr = infile.split("_")[0]
dich = infile.split("_")[1]
if dich == "dual":
   mode = "Dual"
if dich == "red":
   mode = "Red"
if dich == "blue":
   mode = "Blue"
if instr == "m1r":
   inst = "MODS1R"
if instr == "m2r":
   inst = "MODS2R"
if instr == "m1b":
   inst = "MODS1B"
if instr == "m2b":
   inst = "MODS2B"
if lamps == "ArHg":
   lamptitle = "Ar + Hg(Ar) lamps"
if lamps == "NeHg":
   lamptitle = "Ne + Hg(Ar) lamps"
if lamps == "Kr":
   lamptitle = "Kr lamp"
if lamps == "Xe":
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
plt.plot(a[ip,0],log10a,linestyle="solid",color="k")

plt.xlabel("Pixel")
plt.ylabel("log10(Counts[ADU])")

frac  = 0.05
#off= frac*(np.max(a[ip,1]) - np.min(a[ip,1]))
off=  frac*(np.max(log10a) - np.min(log10a))

for l,f in zip(ltab[:,1],ltab[:,2]):
    plt.annotate('%s' % f,xy=(l,ilog10a(l)+off),xycoords='data',xytext=(0,5),textcoords='offset points', ha='center',va='bottom',rotation=90,color="black",size='small')

plt.title("%s %s %s" % (inst, mode, lamptitle))
