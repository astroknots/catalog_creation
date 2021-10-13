import numpy as np
import astropy.io.fits
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
import os, sys

#open up data
infile=sys.argv[1]
valsin=np.loadtxt(infile,dtype=float,delimiter='\t',usecols=2,comments='#')
tmpfile=infile
tmplist=tmpfile.split('/')
tmpname=tmplist[len(tmplist)-1]
outname=tmpname.replace('.dat','')


binnys=100
n,bins,patchs=plt.hist(valsin,binnys)
#print(np.size(bins))


#fit a gaussian to the resultant data
model=models.Gaussian1D()
fitter=fitting.LevMarLSQFitter()
bestfit=fitter(model,bins[0:np.size(bins)-1],n)
print(bestfit)
print(bestfit.stddev.value)

print(np.where(n > 5.0*bestfit.stddev.value))
#plot histograms of vals#
outfile=outname+str('.pdf')
#plt.figure()
#plt.hist(valsin,binnys)
#plt.plot(bins[0:np.size(bins)-1],bestfit(bins[0:np.size(bins)-1]))
#plt.title(outname)
#plt.xlabel('electrons/s (intrinsic HST brightness in the filter)')
#plt.ylabel('counts')
#plt.savefig(outfile)
#plt.close()

