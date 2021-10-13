import numpy as np
from reproject import reproject_exact
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.io import fits
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve

fileblue=sys.argv[1]
filered=sys.argv[2]
fileout=sys.argv[3]
hrp=fits.open(filered)
hop=fits.open(fileblue)
array2,footprint1= reproject_exact(hop,hrp[0].header)

fits.writeto(fileout, array2, hrp[0].header, overwrite=True)
print('done reprojectin sci')


tmpblue=fileblue
tmplist=tmpblue.split('/')
tmpname=tmplist[len(tmplist)-1]
whtfileblue=tmpname.replace('_sci.fits','_wht.fits')

tmpred=filered
tmprlist=tmpred.split('/')
tmprname=tmprlist[len(tmprlist)-1]
whtfilered=tmprname.replace('_sci.fits','_wht.fits')


tmpo=fileout
tmpolist=tmpo.split('/')
tmponame=tmpolist[len(tmpolist)-1]
outwht=tmponame.replace('.fits','_wht.fits')


wb=fits.open(whtfileblue)
wr=fits.open(whtfilered)
whtarray2,footprintwht=reproject_exact(wb,wr[0].header)
fits.writeto(outwht,whtarray2,wr[0].header,overwrite=True)
print('done reprojectin wht')
