import sys

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

from scipy import integrate
from math import sin, log, pi

def lobachevsky(theta):
    if theta == 0.0:
        return(0.0)
    
    result = integrate.quad(lambda u: log(abs(2*sin(u))), 0, theta)
    
    return(-result[0])