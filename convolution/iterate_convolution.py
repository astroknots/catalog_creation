#a program to test successive convolution kernels to find the best FWHM-image matching kernel
import numpy as np
from astropy.wcs import WCS
from astropy.io import fits
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
import sewpy
import astropy.table
import matplotlib as mpl
import matplotlib.pyplot as plt
import os,random,time
import os.path
from os import path
from astropy.io import fits
import scipy.stats as ss
from astropy.modeling import models, fitting
import sys
from array import array
import pathlib


# okay what do I want this to do?

def conv2Dgauss(fitsfilein,fwhm,i=0,ident='sept21_'):
    fits1=fits.open(fitsfilein)
    y=fwhm/(2*np.sqrt(2*np.log(2)))
    kern=Gaussian2DKernel(y)
    img=fits1[0].data
    tmpfile=fitsfilein
    tmpname=tmpfile.split('/')
    tmp2=tmpname[len(tmpname)-1]
    filename=tmp2.replace('.fits','')
    outfile='/home/k689j329/HSTdata/'+filename+'_conv_'+str(ident)+str(i)+'.fits'
    filenameout=outfile
    fil=pathlib.Path(outfile)
    if fil.exists():
        print('already convolved')
        nameout=outfile
        return nameout
    else:
        print('convolving...')
        astrconv=convolve(img,kern)
        fits1[0].data=astrconv
        filenameout=outfile
        fits1.writeto(filenameout)
        return filenameout
    
# another thing to run an input image through sewpy and produce a catalog
def makecat(fits2,see=0.195,ps=0):
    sew=sewpy.SEW(params=["ALPHA_J2000","DELTA_J2000","FWHM_IMAGE","FLUX_APER","X_IMAGE","Y_IMAGE","FLAGS","CLASS_STAR"],config={"DETECT_MINAREA":5, "DETECT_THRESH":1.5,"ANALYSIS_THRESH":1.5,"PHOT_APERTURES":5,"SEEING_FWHM":see,"PIXEL SCALE":ps},sexpath="sextractor")
    outcat=sew(fits2)
    return outcat["table"]

def calcfwhm(table1,star_limit=0.8,maglow=-2.0,maghigh=-0.5,fwhmlim=4.5,plottin=False):
    #define 'stellar sequence' limits from visually assessed convolved image. 
#these  limits were defined by convolution with y=0.15/cdarcblue, and I plotted its FWHM_IMAGE vs logflux to do it
    ylimit=table1['FWHM_IMAGE'] < fwhmlim
    xlimit=np.logical_and(table1['MAG'] > maglow, table1['MAG'] < maghigh)
    goodstar=table1['CLASS_STAR'] > star_limit
    noflags=table1['FLAGS'] < 0.9
    gooddata=np.logical_and(goodstar,noflags)
    limitdata=np.logical_and(xlimit,ylimit)
    limittable=table1[np.where(np.logical_and(gooddata,limitdata))]
    print('filtered data is size:',len(limittable))
    if plottin == True:
        plt.plot(limittable['MAG'],limittable['FWHM_IMAGE'],'ro')
        plt.ylim(0.0,fwhmlim)
        plt.xlim(maglow,maghigh)
        plt.show()
    print('Median FWHM_IMAGE',np.median(limittable['FWHM_IMAGE']),'Mean FWHM_IMAGE',np.mean(limittable['FWHM_IMAGE']))
    return limittable
# a third thing to calculate the fwhm_image average/median in that catalog (for only the good shit)
