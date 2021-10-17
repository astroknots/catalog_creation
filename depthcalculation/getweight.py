#latest updates Sept 2021
#########program to pull out the weight values
import numpy as np
import astropy.table
import matplotlib as mpl
import matplotlib.pyplot as plt
import os,random,time
from astropy.io import fits
import scipy.stats as ss
from astropy.modeling import models, fitting
import sys
from array import array
import pylab as pl
from IPython.display import IFrame



############calculate weight value at a given (x,y)
def wi(fits_wht,x,y): #,apertures,pixscale=0.1):

####Labbe03 plot setup
#take weight image and produce a wi value for the given x,y

    hduwht=fits.open(fits_wht)
    whtdata=hduwht[0].data
    wi=whtdata[x,y]
    
    ##### get wimax too but need to adjust for how relates to invs variance vals
    #### to check, wimax should be the same for all (x,y) values)
   
    validpts= (whtdata >= 2)
    dat=whtdata[validpts]
    sortwht=np.sort(dat)
    pt9perct=round(0.01*len(sortwht))
    pt1perct=len(sortwht)-pt9perct
    
    wimax=sortwht[int(pt1perct)]
    
    return wi,wimax


