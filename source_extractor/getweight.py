#latest updates Sept 2021
#########program to pull out the weight values
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
#from IPython.display import IFrame



############calculate weight value at a given (x,y)
def wi(fits_wht,x,y): #,apertures,pixscale=0.1):

####Labbe03 plot setup
#take weight image and produce a wi value for the given x,y

    hduwht=fits.open(fits_wht)
    whtdata=hduwht[0].data
    wi=whtdata[x-1,y-1]
    #need the -1 bc fits images are gridded from pixel 1,1 and the array of data in whtdata is from 0,0
    return wi
    
    
def getmax(fits_wht):
    hduwht=fits.open(fits_wht)
    whtdata=hduwht[0].data
    validpts= (whtdata >= 1)
    dat=whtdata[validpts]
    #print(np.shape(dat))
    #print(np.shape(whtdata))
    sortwht=np.sort(dat)
    pt9perct=round(0.01*len(sortwht))
    pt1perct=len(sortwht)-pt9perct
    testmax=sortwht[-1]
    #print(testmax)
    
    wimax=sortwht[int(pt1perct)]
    
    return wimax,testmax



