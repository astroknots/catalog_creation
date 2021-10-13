import numpy as np
import pyraf
from pyraf import iraf
from astropy import table
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.io import fits
import sewpy
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.table import Table,join
import astropy.table
import catalog_match_im_trans as cmit
import iterate_convolution as ic

file1='/home/k689j329/HSTdata/j210604m5845/j210604m5845-f606w_drc_sci.fits'
filered='/home/k689j329/HSTdata/j210604m5845/j210604m5845-f105w_drz_sci.fits'
#filereprjd='/home/k689j329/HSTdata/j2106_blue2red_exact.fits'
filereprjd='/home/k689j329/HSTdata/j2106_blue2red_exact_conv_jan217.fits'

hop=fits.open(file1)
hrp=fits.open(filered)
hprpjd=fits.open(filereprjd)

hdr_reprj=hprpjd[0].header
fwhmpix_reprj=1.2
cd1_reprj=hdr_reprj['CD1_1']
cdarc_reprj=abs(cd1_reprj)*3600.0
fwhmarc_reprj=fwhmpix_reprj*cdarc_reprj
table_reprj=ic.makecat(filereprjd,see=fwhmarc_reprj,ps=cdarc_reprj)
mags2=(-1.0)*np.log10(table_reprj['FLUX_APER'])
table_reprj['MAG']=mags2
limtab_reprj=table_reprj[np.where(table_reprj['FLAGS'] < 0.9)]
#limflags_reprj=table_reprj[np.where(table_reprj['FLAGS'] < 0.9)]
#limtab_reprj=limflags_reprj[np.where(limflags_reprj['FWHM_IMAGE'] < 5)]
#limtab_reprj=lim2[np.where(lim2['MAG'] < 0.5)]

hdr_red=hrp[0].header
fwhmpix_red=1.95
cd1_red=hdr_red['CD1_1']
cdarc_red=abs(cd1_red)*3600.0
fwhmarc_red=fwhmpix_red*cdarc_red
table_red=ic.makecat(filered,see=fwhmarc_red,ps=cdarc_red)
magsr=(-1.0)*np.log10(table_red['FLUX_APER'])
table_red['MAG']=magsr
limflags_red=table_red[np.where(table_red['FLAGS'] < 0.9)]
limtab_red=limflags_red[np.where(limflags_red['MAG'] < 0.0)]


xref=limtab_red['X_IMAGE']
yref=limtab_red['Y_IMAGE']
xin=limtab_reprj['X_IMAGE']
yin=limtab_reprj['Y_IMAGE']

#septol needs to be in pixels but I want to initially define it in arcsec
sarc=2.0
septol=sarc/cdarc_red
xoldmatch,yoldmatch,xreprjmatch,yreprjmatch,lims=cmit.cat_im_match(xref,yref,xin,yin,septol)

difx=xoldmatch-xreprjmatch
dify=yoldmatch-yreprjmatch



#plotting
fig,ax=plt.subplots(2,2)
fig.suptitle('X, Y vs DelX, DelY for J2106')
ax[0,0].scatter(xoldmatch,difx)
ax[0,0].set_title('X vs deltax')
ax[0,1].scatter(xoldmatch,dify)
ax[0,1].set_title('X vs deltay')
ax[1,0].scatter(yoldmatch,dify)
ax[1,0].set_title('Y vs deltay')
ax[1,1].scatter(yoldmatch,difx)
ax[1,1].set_title('Y vs deltax')
plt.show()


#cmit.match_diff_im_plot(xoldmatch,yoldmatch,xreprjmatch,yreprjmatch)
print('guessthatworked')
