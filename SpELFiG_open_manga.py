'''
This is the module written with the purpose of defining the opening and extracting of spectra from the fits files. In this case, manga DAP cubes and integrated spectra are considered. Feel free to add new functions in order to open and extracting the spectra from your dataset according to its propper structure. DO NOT FORGET TO MODIFY THE CORRECT FUNCTION CALL IN THE Spelfic_Config_Temp file
'''

import os
import numpy as np
from astropy.io import fits
from scipy.ndimage import gaussian_filter1d

smooth_sigma = 0.4 # To smooth the data with a gaussian kernel.

# Lets open the dap table to extract the dap table:
basepath = os.getcwd()
path_to_dapall = basepath+'/Auxiliary/dapall-v3_1_1-3.1.0.fits'
dap_all = fits.open(path_to_dapall)
daptable = dap_all[1].data


# Redshift correction for DAP spectra:
def redshift_DAP_correct(plate, ifu, wavelenght):
    # Obtain the redshift from DAP catalogue:

    get = np.where(daptable['PLATEIFU']==str(plate)+'-'+str(ifu))
    redshift = daptable['STELLAR_Z'][get][0]

    if redshift<0:
        redshift = dap['Z'][get][0]

    wave0 = wavelenght/(1.0+redshift)

    return wave0

def openDAP_single_spec(inputfile, plate, ifu):

    # Open data:
    data = fits.open(inputfile)
    data = data[1].data

    # Redshift correction:
    wavelenght = redshift_DAP_correct(plate, ifu, data['wavelength'])
    flux = data['emlines']

    spec = np.array([wavelenght, flux, data['flux_error']]).T
    spec = gaussian_filter1d(spec, smooth_sigma, axis=0)
    return spec

def openDAP_cube_spectra(inputfile, plate, ifu):
    '''
    Open DAP MaNGA cube spectra, specifically the VOR10 cubes
    '''
    # First open the file:
    cube = fits.open(inputfile)
    # Obtaining the list of bins by number:
    bin_list = list(np.unique(cube['BINID'].data)[1:])
    # Extract the wavelength vector:
    wavelength = cube['WAVE'].data
    wavelength = redshift_DAP_correct(plate, ifu, wavelength)
    extracted_expectra = []

    # Extract the spectra in each bin:
    for bin in bin_list:
        # Obtaining the spaxels inside the bin:
        where = np.where(bin == cube['BINID'].data[0,:,:])
        x_pos, y_pos = where[0][0], where[1][0]

        # Obtaining the S/N ratio if later useful:
        # S_N = cube['BIN_SNR'].data[x_pos, y_pos]
        # Obtaining the AR and DEC position of the bin, if later useful:
        RA, DEC = cube[0].header['OBJRA'], cube[0].header['OBJDEC']

        flux = cube['EMLINE'].data[:, x_pos, y_pos]
        ivar = cube['IVAR'].data[:, x_pos, y_pos]

        flux_err = np.sqrt(1.0/(ivar))
        spec = np.array([wavelength, flux, flux_err]).T

        # Now obtaining redshift and make redshit correction:
        spec_i = {'spectrum': spec, 'bin number': bin, 'RA': RA, 'DEC': DEC}
        extracted_expectra.append(spec_i)

    return extracted_expectra
