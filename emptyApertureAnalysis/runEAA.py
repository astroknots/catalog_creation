import astropy.io.fits as pyfits
import astropy.table
table = astropy.table.Table()
import os, random, time
import numpy as np
from photutils import CircularAperture
from photutils import aperture_photometry
import math_fits
import math_models
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import matplotlib as mpl
mpl.rcParams['axes.linewidth'] = 1.5
import emptyApertureAnalysis as eaa
import sys

#name and path of image to input
nn=sys.argv[1]
#base name for output files (there will be one per aperture, with the aperture size appended to it)
of=sys.argv[2]
#pixel scale, from the header
ps=sys.argv[3]

apertures=[0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0]

ref_aper=0.7

sigmas,alpha,beta=eaa.run(imname=nn,outfile=of,pixscl=float(ps),apertures=apertures,ref_aper=ref_aper,naper=1000,figures=False,readfromfile=False)

print('successful run, look for outfiles')
