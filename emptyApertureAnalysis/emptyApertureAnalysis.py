import astropy.io.fits as pyfits
import astropy.table
table = astropy.table.Table()
import os, random, time
import numpy as np
from photutils import CircularAperture
from photutils import aperture_photometry
import math_fits
import math_models
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import matplotlib as mpl
mpl.rcParams['axes.linewidth'] = 1.5

def emptyApertureFluxes(im, outfile='j160116m0029_f160w', pixscl=0.128, aper=3.00, naper=1000, sigclip=5., regfile=False):
    """This function finds the x, y coordinates of apertures encircling only nonzero pixels. It then writes the x, y coordinates and the flux within each
    aperture to an output file.

    Keyword arguments:
    im -- input image (for empty aperture analysis, this image should be noise normalized and object subtracted with masked pixels set to zero)
    outfile -- string added to beginning of output file names (e.g. 'macs0416clu_f160w')
    pixscl -- pixel scale of input image (arcsec/pixel)
    aper -- aperture diameter to be used for analysis (arcsec)
    naper -- target number of aperture positions
    sigclip -- apertures prevented from containing pixels with values +/- (sigma * sigclip). Used to prevent selection of apertures containing cosmic rays
    regfile -- creates a DS9 region file of circles at x, y positions (used to visually inspect apertures in DS9)

    Returns: list of measured fluxes.
    """

    # make directory for output files
    try:
        os.stat('aperflux')
    except:
        os.mkdir('aperflux')
    # add an underscore to outfile
    if outfile != '':
        if outfile[-1] != '_':
            outfile = outfile+'_'

    # check image for NaN values and set to zero
    im[np.where(im != im)] = 0

    # calculate sigma for image
    im_mean = np.mean(im[im != 0])
    im_std = np.std(im[im != 0])
    sig_upper = im_mean + (im_std * sigclip)
    sig_lower = im_mean - (im_std * sigclip)

    # aperture radius in pixels, rounded up
    aper_r = aper/2./pixscl
    r = int(np.ceil(aper_r))

    # initial x and y coordinate lists (pixels that are nonzero, between the sigma limits, and more than an aperture radius away from the image edge)
    coordinates = np.where((im != 0) & (im > sig_lower) & (im < sig_upper))
    xlist = coordinates[1][np.where((coordinates[1] > r) & (coordinates[1] < im.shape[1]-r) & (coordinates[0] > r) & (coordinates[0] < im.shape[0]-r))]
    ylist = coordinates[0][np.where((coordinates[1] > r) & (coordinates[1] < im.shape[1]-r) & (coordinates[0] > r) & (coordinates[0] < im.shape[0]-r))]

    # number of good/bad apertures found
    good_aper_count = 0
    # list of good coordinates
    good_x = []
    good_y = []

    while (good_aper_count < naper) & (len(xlist) > 0):
        # randomly select from the initial x, y lists
        i = random.randint(0, len(xlist)-1)
        x = xlist[i]
        y = ylist[i]
        # check if pixels within aperture radius are nonzero and between the sigma limits
        if sum(sum(im[y-r:y+r, x-r:x+r] == 0)) + sum(sum(im[y-r:y+r, x-r:x+r] < sig_lower)) + sum(sum(im[y-r:y+r, x-r:x+r] > sig_upper)) == 0:
            good_aper_count += 1
            good_x.append(xlist[i])
            good_y.append(ylist[i])
            xlist = np.delete(xlist, i)
            ylist = np.delete(ylist, i)
        else:
            xlist = np.delete(xlist, i)
            ylist = np.delete(ylist, i)

    # open output file where the measured fluxes will be saved
    file = open('aperflux/%semptyaperflux_%s.dat' %(outfile, ('%.2f' %(aper)).replace('.','')), 'w')
    file.write('# aperflux run %s, aper=%.2f, sigclip=%.1f\n' %(time.ctime(), aper, sigclip))
    file.write('# x y flux\n')
    if regfile:
        rfile = open('aperflux/%semptyaperflux_%s.reg' %(outfile, ('%.2f' %(aper)).replace('.','')), 'w')
        rfile.write('# aperflux run region file %s, aper=%.2f, sigclip=%.1f\n' %(time.ctime(), aper, sigclip))
        rfile.write('image\n')
    fluxes = []
    # measure flux within aperture at each position and write to save file
    for i in range(len(good_x)):
        x = good_x[i]
        y = good_y[i]
        a = CircularAperture((x, y), r=aper_r)
        flux = (aperture_photometry(im, a))['aperture_sum'][0]
        if flux != 0:
            fluxes.append(flux)
            file.write('%d\t%d\t%f\n' %(x, y, flux))
        if regfile:
            rfile.write('circle(%d, %d, %.2f)\n' %(x, y, aper_r))
    file.close()
    if regfile:
        rfile.close()

    return fluxes


def run(imname='', outfile='', pixscl=0.128, apertures=[], ref_aper=0.7, naper=1000, figures=True, title='', nbins=35, readfromfile=False):
    """This function calculates the image's empty aperture noise properties and returns a list of sigma values and power-law fit parameters.

    Keyword arguments:
    imname -- name of the input image (for empty aperture analysis, this image should be noise normalized and object subtracted with masked pixels set to zero)
    outfile -- string added to beginning of output file names (e.g. 'macs0416clu_f160w')
    pixscl -- pixel scale of input image (arcsec/pixel)
    apertures -- list of aperture diameters to be used for analysis (arcsec) default set to list of apertures used in Shipley et al. empty aperture analysis.
    ref_aper -- aperture used for normalization, must be included in list of apertures
    naper -- target number of aperture positions
    figures -- if True, makes Figure 9 from Shipley et al. 2018
    title -- title in figures (e.g. 'M0416-clu F160W')
    nbins -- number of bins used in histogram figure
    readfromfile -- if True, reads flux from previously saved file

    Returns: list of measured sigma values for each aperture, as well as, the power-law fit parameters, alpha and beta.
    """

    # make directory for output files
    if figures:
        try:
            os.stat('figures')
        except:
            os.mkdir('figures')

    # add an underscore to outfile
    if outfile != '':
        if outfile[-1] != '_':
            outfile = outfile+'_'

    try:
        im = pyfits.getdata(imname)
    except:
        print("Error loading image data")
        return []

    if (apertures == []) | (apertures == [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.75, 2.0]):
        apertures = [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.75, 2.0]
        excluded_apertures = [0.1, 0.15, 1.75, 2.0]
        included_apertures = [aper for aper in apertures if aper not in excluded_apertures]
    else:
        excluded_apertures = []

    if figures:
        figure_apertures = [aper for aper in [0.2, 0.5, 0.7, 1.0, 1.2, 1.5] if aper in apertures]

    sigmas = []
    excluded_sigmas = []
    included_sigmas = []
    NN = []
    included_NN = []
    for aper in apertures:
        if readfromfile:
            try:
                fluxfile = table.read('aperflux/%semptyaperflux_%s.dat' %(outfile, ('%.2f' %(aper)).replace('.','')), format='ascii')
                try:
                    f = fluxfile['aper_fluxes']
                except:
                    f = fluxfile['col3']
                print("%.2f aperture read from file" %(aper))
            except:
                print("working on %.2f aperture" %(aper))
                f = emptyApertureFluxes(im, outfile=outfile, pixscl=pixscl, aper=aper, naper=naper, sigclip=5., regfile=True)
        else:
            print("working on %.2f aperture" %(aper))
            f = emptyApertureFluxes(im, outfile=outfile, pixscl=pixscl, aper=aper, naper=naper, sigclip=5., regfile=True)
        sigma = math_fits.Smad(f)
        n, bins, patches = plt.hist(f, nbins, histtype='stepfilled')
        bsize = bins[1] - bins[0]
        gbins = bins[:-1]+bsize/2
        gfit = list(math_fits.gaussfit(gbins, n, peak=np.max(n), std=sigma))
        gmu = gfit[1]
        gsigma = abs(gfit[2])
        sigmas.append(gsigma)
        r_aper = aper/2./pixscl
        area = np.sqrt(np.pi * r_aper**2)
        NN.append(area)
        if aper in excluded_apertures:
            excluded_sigmas.append(gsigma)
        else:
            included_sigmas.append(gsigma)
            included_NN.append(area)
        if aper == ref_aper:
            ref_sigma = gsigma

    xx = np.log(included_NN)
    yy = np.log(included_sigmas/ref_sigma)
    res = math_fits.robust_linefit(xx, yy)
    alpha = np.exp(res[2])
    beta = res[0]

    if figures:
        emptyAperGaussHist(outfile=outfile, title=title, figure_apertures=figure_apertures, nbins=nbins)
        emptyAperSigmaLinefitFigure(alpha, beta, outfile=outfile, title=title, included_apertures=included_apertures, excluded_apertures=excluded_apertures, included_sigmas=included_sigmas, excluded_sigmas=excluded_sigmas, NN=NN)

    return sigmas, alpha, beta


def emptyAperGaussHist(outfile='', title='', figure_apertures=[0.2, 0.5, 0.7, 1.0, 1.2, 1.5], nbins=35):
    # add an underscore to outfile
    if outfile != '':
        if outfile[-1] != '_':
            outfile = outfile+'_'

    fig = plt.figure(figsize=(8,7))
    plt.xlabel('counts [ADU]', fontsize=16)
    plt.ylabel('N apertures', fontsize=16)

    # hard-coded colors
    colors = ['red','orange','green','blue','cyan','purple']

    handles = []
    legend_names = []
    ymax = 0
    for i in range(len(figure_apertures)):
        aper = figure_apertures[i]
        fluxfile = table.read('aperflux/%semptyaperflux_%s.dat' %(outfile, ('%.2f' %(aper)).replace('.','')), format='ascii')
        try:
            fluxes = fluxfile['aper_fluxes']
        except:
            fluxes = fluxfile['col3']
        sigma = math_fits.Smad(fluxes)
        nbins = nbins
        n, bins, patches = plt.hist(fluxes, nbins, histtype='step', stacked =True, fill=False, color=colors[i%len(colors)], linewidth=2)
        bsize = bins[1] - bins[0]
        gbins = bins[:-1]+bsize/2
        gfit = list(math_fits.gaussfit(gbins, n, peak=np.max(n), std=sigma))
        gmu = gfit[1]
        gsigma = abs(gfit[2])
        gy = math_models.gauss(gbins, peak=gfit[0], center=gfit[1], std=gfit[2])
        gplot, = plt.plot(gbins, gy, '-', color=colors[i], linewidth=2)
        if np.max(n) > ymax:
            ymax = np.max(n)
        handles.append(gplot)
        legend_names.append(r'%.1f$^{\prime\prime}$' %(aper))
        if aper == figure_apertures[-1]:
            ylimit = ymax*1.1
            plt.ylim(0, ylimit)
            xlimit = np.max((np.ceil(abs(bins[0]*7.5))/10., np.ceil(abs(bins[-1]*7.5))/10.))
            plt.xlim(xlimit*-1., xlimit)
            plt.text(xlimit*0.225, ylimit*0.925, title, fontsize=16, path_effects=[pe.withStroke(linewidth=2,foreground='w')])

    plt.legend(handles, legend_names, loc='upper left', frameon=False, title='aperture diameters')
    fig.savefig('figures/%ssigma_aper_hist.eps' %(outfile), bbox_inches='tight')


def emptyAperSigmaLinefitFigure(alpha, beta, outfile='', title='', included_apertures=[], excluded_apertures=[], ref_aper=0.7, included_sigmas=[], excluded_sigmas=[], NN=[]):

    # check lengths of aperture, sigma
    if (len(included_apertures) != len(included_sigmas)) | (len(excluded_apertures) != len(excluded_sigmas)):
        print('Error: all "included" lists must be the same size and all "excluded" lists must be the same size')
        return

    # add an underscore to outfile
    if outfile != '':
        if outfile[-1] != '_':
            outfile = outfile+'_'

    apertures = included_apertures[:]
    apertures.extend(excluded_apertures)
    apertures.sort()
    sigmas = included_sigmas[:]
    sigmas.extend(excluded_sigmas)
    sigmas.sort()

    ref_sigma = included_sigmas[np.where(np.array(included_apertures) == ref_aper)[0][0]]
    sigmafit = ref_sigma * alpha * (np.array(NN) ** (beta))
    sigfitup = ref_sigma * alpha * (np.array(NN) ** 1.0)
    sigfitlow = ref_sigma * alpha * (np.array(NN) ** 2.0)
    fig = plt.figure(figsize=(8, 7))
    plt.xlabel('aperture diameter [arcsec]', fontsize=16)
    plt.ylabel(r'$\sigma$ [ADU]', fontsize=16)
    plt.plot(included_apertures, included_sigmas, 'k^', markersize=12)
    plt.plot(excluded_apertures, excluded_sigmas, 'k^', markersize=12, markerfacecolor="white", markeredgecolor="black")
    plt.plot(apertures, sigmafit, 'r', linewidth=2)
    plt.plot(apertures, sigfitup, 'k--', linewidth=2)
    plt.plot(apertures, sigfitlow, 'k--', linewidth=2)
    plt.xlim(0, max(apertures))
    plt.ylim(0, (max(sigmas)*11/10))

    if title == '':
        plt.text(min(apertures), np.max(sigmas)*1.02, r'$\alpha$ = %.3f' %(alpha), fontsize=14, path_effects=[pe.withStroke(linewidth=8,foreground='w')])
        plt.text(min(apertures), np.max(sigmas)*0.95, r'$\beta$ = %.3f' %(beta), fontsize=14, path_effects=[pe.withStroke(linewidth=8,foreground='w')])
        plt.text(min(apertures), np.max(sigmas)*0.89, r'$\sigma$ (D = %.1f$^{\prime\prime}$) = %.3f' %(ref_aper, ref_sigma), fontsize=14, path_effects=[pe.withStroke(linewidth=8,foreground='w')])
    else:
        plt.text(min(apertures), np.max(sigmas)*1.02, title, fontsize=16, path_effects=[pe.withStroke(linewidth=8,foreground='w')])
        plt.text(min(apertures), np.max(sigmas)*0.95, r'$\alpha$ = %.3f' %(alpha), fontsize=14, path_effects=[pe.withStroke(linewidth=8,foreground='w')])
        plt.text(min(apertures), np.max(sigmas)*0.89, r'$\beta$ = %.3f' %(beta), fontsize=14, path_effects=[pe.withStroke(linewidth=8,foreground='w')])
        plt.text(min(apertures), np.max(sigmas)*0.83, r'$\sigma$ (D = %.1f$^{\prime\prime}$) = %.3f' %(ref_aper, ref_sigma), fontsize=14, path_effects=[pe.withStroke(linewidth=8,foreground='w')])

    fig.savefig('figures/%ssigma_aper_figure.eps' %(outfile), bbox_inches='tight')
