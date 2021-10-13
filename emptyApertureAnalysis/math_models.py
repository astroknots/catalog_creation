import numpy as N

#######################################
## Gaussian
def gauss(x, peak=1., center=0., std=.1):
    norm = []
    for i in range(x.size):
        #norm += [peak/(std*N.sqrt(2*N.pi))*N.exp(-(x[i] - center)**2/(2*std**2))]
        norm.append(peak*N.exp(-(x[i] - center)**2/(2*std**2)))
    return N.array(norm)

## Lorentzian
def lorentz(x, peak=1., center=0., std=.1):
    norm = []
    for i in range(x.size):
        norm += [peak/N.pi*( (N.sqrt(2*N.log(2)*std))/((x[i] - center)**2 + 2*N.log(2)*std**2) )]
    return N.array(norm)
