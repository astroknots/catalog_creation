import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import astropy
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
#from astropy.wcs.utils import skycoord_to_pixel
from astropy import units as u

filename=sys.argv[1]
hdu1=astropy.io.fits.open(filename)
wcoords=WCS(hdu1[0].header)

dat1=hdu1[0].data
#pos1=(34.58875,-5.174166667)
pos2=SkyCoord('02h18m21.300s -05d10m27.000s',frame='fk5')
siz1=u.Quantity((878,547),u.arcsec)


cutout=Cutout2D(dat1,position=pos2,size=siz1,wcs=wcoords,mode='trim')
print(cutout.shape)

print((cutout.position_original, cutout.position_cutout))



outfile=sys.argv[2]
hdu1[0].data=cutout.data
hdu1[0].header.update(cutout.wcs.to_header())
hdu1.writeto(outfile,overwrite=False)


