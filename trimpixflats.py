#!/usr/bin/env python


from astropy.io import fits
import matplotlib.pyplot as plt
import astropy.stats as stats
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='trim 8Kx3K pixel flats to 4Kx3K')

parser.add_argument('inFile',type=str,help='input file 1')

args = parser.parse_args()

inFile  = args.inFile

outSuffix = "_tr4K"
outFile = inFile[:inFile.find(".fits")] + outSuffix + ".fits"

f1 = fits.open(inFile)
d1 = f1[0].data
h1 = f1[0].header


t1 = d1[:,2048:6144]
h1.add_history('Overscan-subtracted 8Kx3K trimmed to 4K, region [:,2048:6144]')

hdu = fits.PrimaryHDU(data=t1,header=h1)
hdu.writeto(outFile)
