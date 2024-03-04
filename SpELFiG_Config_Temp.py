#!/usr/bin/env python3

'''
MODULE OF CONFIGURATION FOR THE SPELFIC CODE EXECUTION
Please follow up the instructions in order to correctly configurate Spelfic for its execution.
'''


# These variables are going to be modified by the RUN script.
ID_OBJ = '[ID_OBJ]'
INPUT_FILE = '[INPUT_FILE]'
spectype = '[SPECTYPE_I]'

# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                      I.    ::::::::::: LISTS OF EMISSION LINES :::::::::
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The following lists in form of dictionaries contain the emission lines to be fitted to the spectra
# Please review them as you can add more lines to be fitted at the spectral range you are going to
# provide. In such case, please add the line with same structure (one item per line, along with its
# restframe wavelenght, we recommend to mimic the existing IDs for other lines, using ,1 ,2 for
# doublets). The restframe wavelengths used were taken from FIREFLY (Wilkinson et al., 2017)
# archive; you may modify and update as you wish. ALL_ELM contains all lines, BLR_EML and NLR_EML
# contains the lines from the Broad Line Region and Narrow Line Region respectively, and OUT_EML
# contain the lines associated also with kinematic outflows.
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


# All emission lines: ------------------------------------------------------------------------------
ALL_EML = {
	'He-II,1':  [3202.15],
	'He-II,2': [4685.74],
	'Ne-V,1':   [3345.81],
	'Ne-V,2':   [3425.81],
	'O-II,1':   [3726.03],
	'O-II,2':  [3728.73],
	'Ne-III,1': [3868.69],
	'Ne-III,2':  [3967.40],
	'H-ζ':     [3889.05],
	'H-ε':     [3970.07],
	'H-δ':    [4101.73],
	'H-γ':     [4340.46],
	'O-III,0':   [4363.15],
	'O-III,1':   [4958.83],
	'O-III,2':   [5006.77],
	'Ar-IV,1':  [4711.30],
	'Ar-IV,2':  [4740.10],
	'H-β':     [4861.32],
	'N-I,1':    [5197.90],
	'N-I,2':    [5200.39],
	'He-I':   [5875.60],
	'O-I,1':   [6300.20],
	'O-I,2':   [6363.67],
	'N-II,1':   [6547.96],
	'N-II,2':   [6583.34],
	'H-α':     [6562.80],
	'S-II,1':   [6716.31],
	'S-II,2':   [6730.68],
	'Ar-III': [7135.67],
	}


# Broad Line Region Lines (FOR AGNs): --------------------------------------------------------------
BLR_EML = {
	'He-II,1': [3202.15],
	'He-II,2': [4685.74],
	'H-ζ': [3889.05],
	'H-ε': [3970.07],
	'H-δ': [4101.73],
	'H-γ': [4340.46],
	'H-β': [4861.32],
	'He-I': [5875.60],
	'H-α': [6562.80]
}

# Narrow Line Region Lines (FOR AGNs): -------------------------------------------------------------
NLR_EML = {
	'Ne-V,1': [3345.81],
	'Ne-V,2': [3425.81],
	'O-II,1': [3726.03],
	'O-II,2': [3728.73],
	'Ne-III,1': [3868.69],
	'Ne-III,2': [3967.40],
	'O-III,0': [4363.15],
	'O-III,1': [4958.83],
	'O-III,2': [5006.77],
	'Ar-IV,1': [4711.30],
	'Ar-IV,2': [4740.10],
	'N-I,1': [5197.90],
	'N-I,2': [5200.39],
	'O-I,1': [6300.20],
	'O-I,2': [6363.67],
	'N-II,1': [6547.96],
	'N-II,2': [6583.34],
	'S-II,1': [6716.31],
	'S-II,2': [6730.68],
	'Ar-III': [7135.67]
}

# Outlfow Associated Lines (FOR AGNs): -------------------------------------------------------------
OUT_EML = {
	'O-III,1': [4958.83],
	'O-III,2': [5006.77]
}


# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                  II.    ::::::::::: RECIPROCAL RATIOS BETWEEN DOUBLETS :::::::::
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Now here the doublets tie conditions dictionary:
# BOTH LINES of the doublet must be present in the dictionary with the correspondant values of the
# amplitude and sigma ratios (reciprocal between them). With the free entrance, set True or False
# depending on which line you want to set with free parameters for the fits to run. Add ratios for
# other doublets as you wish, in order for be considered by Spelfic.
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

doublets_info = {
	'O-III,1': {'free': True, 'amplitude_ratio': 1./3., 'sigma_ratio': 1.},
	'O-III,2': {'free': False, 'amplitude_ratio': 3., 'sigma_ratio': 1.},
	'N-II,1': {'free': True, 'amplitude_ratio': 1./3., 'sigma_ratio': 1.},
	'N-II,2': {'free': False, 'amplitude_ratio': 3., 'sigma_ratio': 1.}
}

# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                  III.    ::::::::::: LOCAL PATHS AND FILE WRITTING FLAGS :::::::::
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# In the bellow variables, please fill as indicated the local paths of the directories as needed
# Also, ensure to save in the GALS variable an array, list-like object with the list of manga
# galaxies to fit. Modify and/or add your own lines of code if you wish to do so.
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

import os
basepath = os.getcwd()

if spectype == 'S':
	respath = basepath+f'/Outputs/Spectra/{ID_OBJ}'

elif spectype == 'C':
	respath = basepath+f'/Outputs/Cubes/{ID_OBJ}'

inputfile = f'{INPUT_FILE}'
outputname = respath +  f'/SpL-{ID_OBJ}'
prior_file = respath + f'/priorsMCMC-{ID_OBJ}.csv'

# ··································      FLAGS    ·················································
redshift_flag = True # Set this to true or false depending if you want to correct the input spectra
# for redshit. At this version, DAP redshifts from the dapfile given above are used.

# Modify this to True if you want to automatically estimate parameters and
# produce the priors and ranges file for the mcmc iteration regardless if there is one already. This
# will cause that modifying by hand this file will be useless.
rewrite_priors_flag = True


# ****************  PRODUCTION OF THE OUTPUT FILE FLAGS: *******************************************
# If the priors_plots_flag is set as True, the execution of the code will automatically produce the plots
# of the fits, one per each spectral zoom specified in specranges:
priors_plots_flag = False
samples_plots_flag = False
fit_plots_flag = True


# If the samples_flag is set as True, the samples file will be written. This file occupies BIG
# memory space (~ 150 MB with the chain production associated with the default configuration). if
# you do not have a lot of storing space we recommend to set it as False.
samples_file_flag = False

# If the properties_flag is set as True, the derivated lines properties file will be written.
properties_file_flag = True

Lines_for_global_tables = ['H-β', 'O-III,2', 'O-I,1', 'H-α', 'N-II,2', 'S-II,2']

# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#               IV.    ::::::::::: SPECTRAL GENERAL FEATURES FOR THE FITS :::::::::
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# In the next lines, some important configuration assets for the fits are settled. In the variable
# specrange0, specify the spectral range (in Agnstrom units) to be fitted. In the tuple zooms you
# can modify/add some spectral ranges in order to produce plots with zoom in those regions (one
# plot per each one). The variable niter is for the number of chains in the emcee cycles. Then, the
# maxncomp variables set the maximun number of components per each subgroup of lines (the adding of
# components will stop either when the BIC stops improving or when maximun number of components is
# reached). Finally, the Lines_for_global_tables list must contain the IDs of the lines (the above
# lists same ones) for which a global properties of the line table (for the entire objects sample)
# is desired at the end of the execution.
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

# JEFFREY SCALE FOR THE BIC:
# Choose the value of the BIC cut the difference in BIC

DBIC_cut = 10. # Weak
# DBIC_cut = 10**(1.5) # Definite
# DBIC_cut = 10**(2) # Strong
# DBIC_cut = 10**(3) # Very Strong


# Define the spectral range to fit:
specrange0 = [4800, 6800]


# Define a range of plane continuum INSIDE specrange0:
continuum_data0 = [5050, 5150]

# Define the spectral range to fit:
plotzooms = ([4800.0, 5100.0],
	[6200.0, 6800.0],
	[5050, 5150])

specranges = (specrange0, plotzooms)

# IMPORTANT: Number of samples to be calculated in the emcee cycles: n parameters * niter
'''
RECOMENDATION: Adjust between 1000 and 500, depending the number of free parameters. If HIGH number
of free parameters, a lower niter is better for the best performance of the code.
'''

niter = 1000

smooth_sigma = 0.4
# Factors for the estimations. The following values will serve for the estimation of new components
# as indicated:  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Maximal number of components for the Narrow Line Region:
maxncomp_nlr = 1
# IMPORTANT: This is a first estimation of the maximal fwhm of the narrow first component of the
# lines. This is in Angstroms.
max_width_narrow = 10.0

# Maximal number of components corresponding to emission from the Broad Line Region:
maxncomp_blr = 2
# IMPORTANT: This is a first estimation of the maximal fwhm of the broad emission component(s) of
# BLR lines. This is in Angstroms.
max_width_blr = 50.0
f_sigma_blr = 2. # The factor of sigma for making the initial guesses in terms of the narrow
# component, i.e., sigma_blr ~ f_sigma_blr * sigma_nlr

# Maximal number of components for Outflow Associated Lines:
maxncomp_outflow = 2
# IMPORTANT: This is a first estimation of the maximal fwhm of the outflow associated component(s)
# of the outflow lines. This is in Angstroms.
max_width_outflow = 30.0
f_sigma_out = 1.2 # The factor of sigma for making the initial guesses in terms of the narrow
# component, i.e., sigma_out ~ f_sigma_out * sigma_nlr


delta_loc = 0.0 # The centroid of the new component(s) will be shifted by this value (use a
# negative sign for a blueshift). For symetric gaussian profiles, the code works well with 0.0. You
# can, for example, modify this for a spectra where a systematic and noticeable shift of the extra
# component(s) is present.

# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#      V.    ::::::::::: LAST SECTION: SPECTRUM READING AND OTHER CONFIGURATION FLAGS :::::::::
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FROM HERE, modify things carefully. In the first section the spectrum must be read from the input
# file. There is a redshift correction, in case it is necessary if the fits file comes from the DAP
# pipeline. Modify this lines as needed, but the final SPEC object must have the exact same form.
#
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

# FIRST: Read the spectra: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import SpELFiG_open_manga as splo

# Extract plate & ifu FOR A MaNGA ID:
group = ID_OBJ.split('-')
plate, ifu = group[1], group[2]

if spectype == 'S':
	spec = splo.openDAP_single_spec(inputfile, plate, ifu)
	spec = spec[(spec[:,0]>=specrange0[0]) & (spec[:,0]<=specrange0[1])]

elif spectype == 'C': spec = splo.openDAP_cube_spectra(inputfile, plate, ifu)
