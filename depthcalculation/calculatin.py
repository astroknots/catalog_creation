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



def calc_reduced_chi_square(fit, x, y, yerr, N, n_free):
#    fit (array) values for the fit
#    x,y,yerr (arrays) data
#    N total number of points
#    n_free number of parameters we are fitting
    return 1.0/(N-n_free)*sum(((fit - y)/yerr)**2)
