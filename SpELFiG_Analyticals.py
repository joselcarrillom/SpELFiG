import numpy as np
from scipy.integrate import quad
from scipy.optimize import fsolve

c = 299792.458  # Speed of light in km/s (vacuum)
inf = np.inf

# ANALYTICAL FUNCTIONS FOR THE SPECTRAL MODELS. ººººººººººººººººººººººººººººººººººººººººººººººººººº

def gauss(loc, a0, sd, x):
    diff = x-loc
    return  a0 * np.exp(-0.5*(((diff) / sd)**2))

def composed_gaussians(Ncomp, x, *params):
    pline = np.zeros_like(x)
    start = 0
    # print('%%% HERE composed_gaussians')
    # print(f'Ncomp: {Ncomp}')
    # print('params')
    # print(params)
    # print('len params', len(params))
    _params = np.array(*params[start:start+3])
    # print('_params')
    # print(_params)
    for i in range(Ncomp):
        sub_params = _params[start:start+3]
        loc_i = sub_params[0]
        amp_i = sub_params[1]
        sigma_i = sub_params[2]
        pline += gauss(loc_i, amp_i, sigma_i, x)
        start += 3
    return pline

def continuum_function(x, *params):
    a, b, c = params[0], params[1], params[2]
    if b == 0:
        a == 0
    return a*(x/b)**(-c)

def intflux_gauss(line, ncomp):
    '''
    Function to calculate the integrated fluxes of a given line, fitted with a gaussian (or linear
    sum of gaussians) profile
    '''
    A, SD = [], []
    for i in range(ncomp):
        #print('This A inside loop', line[(3*i)+1])
        A.append(line[(3*i)+1])
        SD.append(line[(3*i)+2])
    #print('ºººººº', A)
    A = abs(np.array(A))
    SD = abs(np.array(SD))
    fluxint = np.sum(A*SD)*np.sqrt(2*np.pi)
    return fluxint

def mean_loc(line, ncomp):
    '''
    Calculate the mean centroid of multicomponent gaussian line
    '''
    LOC, A, SD = [], [], []
    for i in range(ncomp):
        LOC.append(line[3*i])
        A.append(line[(3*i)+1])
        SD.append(line[(3*i)+2])
    LOC = np.array(LOC)
    A = np.array(A)
    SD = np.array(SD)
    loc_c = np.sum(LOC*A*SD)/np.sum(A*SD)
    return loc_c


# VELOCITY TRANSFORMATIONS AND FLUX FRACTION VELOCITY CALCULATION ºººººººººººººººººººººººººººººººº

def vel_correct(l0, l):
    '''
    Transform wavelength to units of velocity (km/s) given a reference wavelength l0 (in case of a single gaussian-fitted-line analysis, it is the centroid).
    '''
    vel = c*((l-l0)/l0)
    return vel

def vel_correct_sigma(sigma, l0):
    '''
    Transform sigma (D(lambda)) to velocity units
    '''
    vel = c*(sigma/l0)
    return vel

def V_flux(NCOMP, linepars, r, LOC):
    '''
    Function to solve for flux fraction velocity of a given fitted line.
    The radial velocity corresponds to the line-of-sight component of the velocity
    # vector of the object. In this case, it is represented by the median velocity of the line, i. e.; the value that bisects the total integrated flux of the line.
    '''
    def eq_int_vel(x, Ncomp, linepars, r):
        g = lambda x: composed_gaussians(Ncomp, x, linepars)
        F = intflux_gauss(linepars, Ncomp)
        y, err = quad(g,-inf, x)
        return y - r*F

    sol = fsolve(eq_int_vel, LOC, args=(NCOMP, linepars, r))
    vel_rad = vel_correct(LOC, sol)
    return vel_rad
