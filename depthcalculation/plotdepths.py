import numpy as np
import scipy
import astropy
import matplotlib.pyplot as plt
import math
#from matplotlib.backends.backend_pdf import Pdfpages
import sys

file1='output_irc0218_f160w_v1.txt'
data160=np.genfromtxt(file1)
#data814=np.genfromtxt('data_j16_f814w_adjustexptime2.txt')


apertures=np.array([0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0])
pixs=  0.0999999999999972
pixs2= 0.0499999999999968
ap1=0.5*apertures/pixs
ap2=0.5*apertures/pixs2
#data1=np.array(data160)/(ap1*ap1)
#data8=np.array(data814)/(ap2*ap2)

#fig=plt.figure(0)
#plt.plot(apertures)
#plt.plot(ap1,data1,'r*')
#plt.plot(ap2,data8,'b^')
#plt.title('sigma/area as a function of linear size')
#plt.ylabel('sigma per unit area in e-/s/pix^2')
#plt.xlabel('pixels')
#plt.savefig('j160116m0029_sigmaperarea_linear.pdf')
plt.plot(apertures,data160,'r*')
#plt.plot(apertures,data814,'b^')
plt.ylabel('depth in AB magnitude')
plt.xlabel('aperture diameter in arcsec')
plt.title('5 sigma empty aper fluxes')
#plt.savefig('j160116m0029_depths_v7.pdf')
plt.show()
plt.close()
