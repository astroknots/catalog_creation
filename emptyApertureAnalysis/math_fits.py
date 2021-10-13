import numpy as N
import math_models as mm
from scipy.optimize import leastsq
import scipy.stats as stats
import bootstrap as bs
from statsmodels.formula.api import rlm

def robust_linefit(x, y):
    """Calculates and returns linear fit to given (x,y) pairs using a robust linear fit routine.
    Return order is slope, standard deviation of slope, intercept and standard deviation of intercept."""
    regression = rlm("y ~ x", data=dict(y=y, x=x)).fit()
    slope = regression.params[1]; yint = regression.params[0] # slope and intercept fit values
    cov = regression.bcov_scaled.x[0]
    
    sig_slope = regression.bse[1]; sig_yint = regression.bse[0]/N.sqrt(regression.nobs)*.5
    # normalized standard deviation errors of the central values for the given sample.
    
    s_yint = N.sqrt(N.sum(regression.resid**2)/regression.nobs)*.6745#/1.105
    # s_yint and the 1.105 is an empirical value taken from robust_linefit.pro to match IDL's routine
    # for the standard deviations of the given sample.

    return slope, sig_slope, yint, sig_yint, s_yint, cov

def linearfit(x, y, runs=1e4):
    """Calculates and returns linear fit to given (x,y) pairs.
    Return order is slope, standard deviation of slope, intercept and standard deviation of intercept."""
    x *= N.ones(1); y *= N.ones(1)
    slope, yint, r, prob2, see = stats.linregress(x, y)
    scale = N.mean(y) - N.mean(x)
    nsig = Smad(N.log10(10**y/(10**x*10**scale)))
    bs_slope = N.zeros(runs); bs_yint = N.zeros(runs)
    for i in range(bs_slope.size):
        bsrun = bs.bootstrap_resample(N.arange(x.size))
        bsuse = N.where((10**y[bsrun]/(10**x[bsrun]*10**scale)>=10**(-1.*nsig)) & (10**y[bsrun]/(10**x[bsrun]*10**scale)<=10.**nsig))
        bs_slope[i], bs_yint[i] = stats.linregress(x[bsuse], y[bsuse])[:2]
    sd_slope = bs_slope.std(); sd_yint = N.sqrt(sd_slope**2/bs_yint.size*N.sum(bs_yint**2))
    
    return slope, sd_slope, yint, sd_yint

def linearfit_s1(x, y):
    """Calculates and returns linear fit to given (x,y) pairs for slope equals 1.
    Return order is slope, standard deviation of slope, intercept and standard deviation of intercept.
    Uses robust_linefit to calculate linear parameters."""
    x *= N.ones(1); y *= N.ones(1)
    mx = N.mean(x); my = N.mean(y)
    slope = 1.; yint = my-mx
    rl_slope, sd_slope, rl_yint, sd_yint, s_yint, cov = robust_linefit(x, y)
    
    return slope, sd_slope, yint, sd_yint, s_yint, cov

def Smad(x):
    """Computes the median absolute deviation (see Beers et al. 1990) from a sample median."""
    x *= N.ones(1)
    MAD = N.median(N.abs( x - N.median(x) ))
    Smad = MAD/0.6745

    return Smad

def gaussfit(x, y, peak=1., center=0., std=.1):
    """Calculates and returns a gaussian fit to given (x,y) pairs with intial guesses.
    Default initial guesses are used unless supplied by user."""
    def res(p, y, x):
        top1, m1, std1 = p
        y_fit = mm.gauss(x, top1, m1, std1)
        err = y - y_fit
        return err
    p = [peak, center, std]  # Initial guesses for leastsq
    plsq = leastsq(res, p, args = (y, x), maxfev=2000)
    
    return plsq[0]

def gaussfit2(x, y, peak=1., center=0., std=.1):
    """Calculates and returns a gaussian fit to given (x,y) pairs with intial guesses.
    Default initial guesses are used unless supplied by user."""
    def res(p, y, x):
        top1, top2, m1, m2, std1, std2 = p
        y_fit = mm.gauss(x, top1, m1, std1) + mm.gauss(x, top2, m2, std2)
        err = y - y_fit
        return err
    p = [peak[0], peak[1], center[0], center[1], std[0], std[1]]  # Initial guesses for leastsq
    plsq = leastsq(res, p, args = (y, x), maxfev=2000)
    
    return plsq[0]

def gaussfit3(x, y, peak=1., center=0., std=.1):
    """Calculates and returns a gaussian fit to given (x,y) pairs with intial guesses.
    Default initial guesses are used unless supplied by user."""
    def res(p, y, x):
        top1, top2, top3, m1, m2, m3, std1, std2, std3 = p
        y_fit = mm.gauss(x, top1, m1, std1) + mm.gauss(x, top2, m2, std2) + mm.gauss(x, top3, m3, std3)
        err = y - y_fit
        return err
    p = [peak[0], peak[1], peak[2], center[0], center[1], center[2], std[0], std[1], std[2]]  # Initial guesses for leastsq
    plsq = leastsq(res, p, args = (y, x), maxfev=2000)
    
    return plsq[0]

def lorentzfit(x, y, peak=1., center=0., std=.1):
    """Calculates and returns a lorentzian fit to given (x,y) pairs with intial guesses.
    Default initial guesses are used unless supplied by user."""
    def res(p, y, x):
        top1, m1, std1 = p
        y_fit = mm.lorentz(x, top1, m1, std1)
        err = y - y_fit
        return err
    p = [peak, center, std]  # Initial guesses for leastsq
    plsq = leastsq(res, p, args = (y, x), maxfev=2000)
    
    return plsq[0]
