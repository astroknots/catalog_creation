#latest updates July 2021
#########the part of the program that does the calculations
import numpy as np
import astropy.table
import matplotlib as mpl
import matplotlib.pyplot as plt
import os,random,time
import matplotlib.backends.backend_pdf
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import fits
import scipy.stats as ss
from astropy.modeling import models, fitting
import sys
from array import array
import pylab as pl
from IPython.display import IFrame



def depth(totfile,photvals):
    valsin=np.loadtxt(totfile,dtype=float,delimiter='\t',usecols=2,comments='#')
    #########calculate average flux and stdev
    stdev_empty=np.std(valsin)
    avg_emptyflux=np.average(valsin)

    #########set up parameters
    temp_photflam=photvals[0]
    temp_photplam=photvals[1]

    #########calculate abmag_depth
    abmagzpt=-2.5*np.log10(temp_photflam) - 21.10 - 5.0*np.log10(temp_photplam) + 18.692
    abmag_depth=-2.5*np.log10(5.0*stdev_empty)+abmagzpt
    
    return abmag_depth,stdev_empty

###########the part of the program that calculates the magnitude given a flux
def magcalc(fluxfile,photvals):
    fluxvalsin=np.loadtxt(fluxfile,dtype=float,delimiter='\t',usecols=0,comments='#')
#that assumes that the flux values will be put in via a single column txt file
    ####set up parameters
    temp_photflam=photvals[0]
    temp_photplam=photvals[1]

    #####calculate abmag_depth
    abmagzpt=-2.5*np.log10(temp_photflam) - 21.10 - 5.0*np.log10(temp_photplam) + 18.692
    abmag=-2.5*np.log10(fluxvalsin)+abmagzpt
    return abmag

#########gaussian fit to histogram
def fithists(totfile,binnum):
    valsin=np.loadtxt(totfile,dtype=float,delimiter='\t',usecols=2,comments='#')
    binnys=binnum
    #####set up the name to reference the dictionary with
    tmpfile=totfile
    tmplist=tmpfile.split('/')
    tmpname=tmplist[len(tmplist)-1]
    outname=tmpname.replace('.dat','')

    #####make a histogram out of the empty aperture data
    hist1,bins1=np.histogram(valsin,binnys)
    
    ######fit a guassian to the histogram
    holdhist={}
    model=models.Gaussian1D()
    fitter=fitting.LevMarLSQFitter()
    bestfit=fitter(model,bins1[0:(len(bins1)-1)],hist1)
    holdhist={'hist1':valsin,'bins':bins1,'bestfit':bestfit}
    return holdhist,outname

#############calculate RMS for Labbeplot
def rmscalc(fits_sci,fits_mask,fits_seg,apertures,pixscale=0.1):

####Labbe03 plot setup
#take full image (hdusci) and select only those pixels that are:
#data (from hdumask) and that are not sources (from segdata)
	####### pull out the data
	hdusci=fits.open(fits_sci)
	hdumask=fits.open(fits_mask)
	hduseg=fits.open(fits_seg)
	apers=np.array(apertures)
	scidata=hdusci[0].data
	maskdata=hdumask[0].data
	segdata=hduseg[0].data
	
    ####### select 'valid' points that are not parts of objects or off-image
	validpts=np.where((segdata<1) & (maskdata==1))
	scivalid=scidata[validpts]

    
    #####actual calculation
	rmss=np.std(scivalid)
	rmsactual=np.sqrt(np.mean(scivalid**2.0))
	print('STDEV of all "valid" pixels in scidata:',rmss)
	print('rms actual (square root of the mean of the squared values of all valid data points):',rmsactual)
	linearsize=0.5*np.array(apers)/pixscale 
	#times 1/2 for 'radius' rather than 'diameter' linear size
	#is my linear size correct? N = sqrt(Area) ~ radius. above is radius in pixels, so, should be.
	return rmss,linearsize


