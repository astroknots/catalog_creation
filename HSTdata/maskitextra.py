import os
import sys
import numpy as np
import scipy as sp
import astropy
from astropy.io import fits
import matplotlib.pyplot as plt

##########read in the trimmed science image and the trimmed weight image
filename=sys.argv[1]
whtfile=sys.argv[2]
hdu1=astropy.io.fits.open(filename)
hdu2=astropy.io.fits.open(whtfile)
data1=hdu1[0].data
data2=hdu2[0].data



#########construct the normalization factor fn
#first, select only the parts of the weight image where data exist
dat=data2[np.where(data2 >=2)]
#next, flatten that to an array of all values in the data, and sort it
dat2=dat.flatten()
sortwht=np.sort(dat2)
#find the position of the 1% highest value within the sorted array; it should be at the 1% position in the array
pt9perct=round(0.01*len(sortwht))
pt1perct=len(sortwht)-pt9perct
print('length and 1percent point',len(sortwht),pt1perct)
print('former position and value',pt9perct,sortwht[int(pt9perct)])
val1perct=sortwht[int(pt1perct)]
#compare that value (val1perct) with the value if I'd rounded up or down instead--this proves the values are virtually identical, so it doesn't matter
print('1percent points?',val1perct,sortwht[int(pt1perct)-1])
#so then the normalization factor fn is:
fn=np.sqrt(data2)/np.sqrt(val1perct)


#######perform the normalization and output the data
outdata=data1*fn
hdout=fits.PrimaryHDU(outdata,hdu1[0].header)
tmpfile=filename
tmplist=tmpfile.split('/')
tmpname=tmplist[len(tmplist)-1]
outname=tmpname.replace('.fits','')
newname=outname+'_newnorm3.fits'
print(newname)
hdout.writeto(newname)


#####below here, create a mask to use for removing non-zero values outside the data image created by convolution in sextractor
data3=hdu2[0].data

i1=np.where(data3>2)
i0=np.where(data3<2)

data3[i1]=1
data3[i0]=0
#######and output it to a new file
hdo2=fits.PrimaryHDU(data3,hdu2[0].header)
newname2=outname+'_mask.fits'
hdo2.writeto(newname2)
