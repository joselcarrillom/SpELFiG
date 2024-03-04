#!/usr/bin/env python3
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from astropy.cosmology import WMAP9 as cosmo
from scipy.interpolate import CubicSpline

import os
import shutil
import tempfile
import re
import math
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

import warnings
warnings.filterwarnings("ignore")


# Enable LaTeX rendering
plt.rc('text', usetex=True)
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.size'] = 16

def read_and_extract(file_path, lines_of_interest):
    # Read the CSV file into a DataFrame
    galaxy_data = pd.read_csv(file_path)
    match = re.search(r'manga-(\d+)-(\d+)', file_path)

    if match:
        plate, ifu = match.groups()
        galaxy_id = f'{plate}-{ifu}'
    else:
        # Default value if no match is found
        galaxy_id = 'Unknown'

    # Extract rows for the lines of interest
    lines_data = galaxy_data[galaxy_data['Line Name'].isin(lines_of_interest)]
    lines_data.loc[:, 'Galaxy ID'] = galaxy_id
    # Print the resulting lines_data

    return lines_data

def get_files(src_directory, template):
    matching_files = []

    # Iterate through the source directory and its subdirectories
    for root, _, files in os.walk(src_directory):
        # Filter files based on the template
        matching_files.extend([os.path.join(root, file) for file in files if template in file])

    return matching_files

'''
def properties_extraction(num_threads, csv_files, lines_of_interest):
    aggregated_data_dict = {}
    # Use ThreadPoolExecutor to parallelize the reading process
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        # Submit the reading and extracting task for each CSV file
        futures = [executor.submit(read_and_extract, csv_file, lines_of_interest) for csv_file in csv_files]
        # Retrieve the results as completed
        for future in futures:
            try:
                lines_data = future.result()

                # Iterate through each line of interest
                for line in lines_of_interest:
                    # If the line is not already a key in the dictionary, create it
                    if line not in aggregated_data_dict:
                        aggregated_data_dict[line] = pd.DataFrame()

                    # Filter the lines_data for the current line
                    line_data = lines_data[lines_data['Line Name'] == line]

                    # Append the line data to the corresponding DataFrame in the dictionary
                    aggregated_data_dict[line] = aggregated_data_dict[line].append(line_data, ignore_index=True)
            except Exception as e:
                print(f"Error processing file: {e}")
    return aggregated_data_dict
'''

def properties_extraction(num_threads, csv_files, lines_of_interest):
    aggregated_data_dict = {}

    # Use ThreadPoolExecutor to parallelize the reading process
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        # Submit the reading and extracting task for each CSV file
        futures = [executor.submit(read_and_extract, csv_file, lines_of_interest) for csv_file in csv_files]
        # Retrieve the results as completed
        for future in futures:
            try:
                lines_data = future.result()

                # Iterate through each line of interest
                for line in lines_of_interest:
                    # If the line is not already a key in the dictionary, create it
                    if line not in aggregated_data_dict:
                        aggregated_data_dict[line] = pd.DataFrame()

                    # Create a DataFrame with all galaxies and the line of interest
                    galaxies = lines_data['Galaxy ID'].unique()
                    all_galaxies_data = pd.DataFrame({'Galaxy ID': galaxies, 'Line Name': line})

                    # Merge with the lines_data for the specific line and galaxy
                    line_data = pd.merge(all_galaxies_data, lines_data, how='left', on=['Galaxy ID', 'Line Name'])

                    # Append the line data to the corresponding DataFrame in the dictionary
                    if not aggregated_data_dict[line].empty:
                        aggregated_data_dict[line] = pd.concat([aggregated_data_dict[line], line_data], ignore_index=True)
                    else:
                        aggregated_data_dict[line] = line_data
            except Exception as e:
                print(f"Error processing file: {e}")

    return aggregated_data_dict

# For extracting the maximal sigmas of the lines:

def extract_max_min_values(dataframe, line_name, parameter):
    # Filter the DataFrame for the specified line name
    line_data = dataframe[dataframe['Line Name'] == line_name]

    if line_data.empty:
        return [np.nan, np.nan]

    # Find the maximum and minimum values for the specified parameter
    if len(line_data[parameter]) == 1:
        return line_data[parameter].values
    else:
        max_value = line_data[parameter].max()
        min_value = line_data[parameter].min()
        return [min_value, max_value]

# PLOTTING FUNCTION : ··············································································
# MOSAIC FORM ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

def plot_comparison_mosaic(data_tuples, titles):
    num_plots = len(data_tuples)
    rows = (num_plots + 2) // 3  # Calculate the number of rows with a maximum of 3 subplots per row

    fig, axes = plt.subplots(rows, 3, figsize=(15, 5 * rows))

    for i, (data, title) in enumerate(zip(data_tuples, titles), 1):
        ax = axes.flatten()[i - 1]

        # Scatter plot
        ax.scatter(data[0], data[1], alpha=0.7)

        # Identity line
        min_val = min(np.min(data[0]), np.min(data[1]))
        max_val = max(np.max(data[0]), np.max(data[1]))
        # if np.mean(data[0]) == 0. & np.mean(data[1]) == 0.:
        #     min_val = (-0.5, -0.5)
        #     max_val = (0.5, 0.5)
        ax.plot([min_val, max_val], [min_val, max_val], linestyle='--', color='red', label='Identity Line')

        ax.set_title(title)
        ax.legend(frameon=False)

    # Adjust layout and show the plot
    plt.tight_layout()
    return fig

# ··································································································
# Defining the output
spectype = 'S'
if spectype == 'S': specdir = 'Spectra'
elif spectype == 'C': specdir = 'Cubes'

basepath = os.getcwd()
outputs_directory = basepath + f'/Outputs/{specdir}'

# outputs_directory = f'/home/joselcarrillom/SpELFiG/Outputs/{specdir}'
template_name = '_line_properties.csv'

# Define the list of lines for the subgroup (replace with the lines you're interested in)
lines_of_interest = ['O-III,2', 'H-β', 'O-I,1', 'N-II,2', 'H-α', 'S-II,2']

# Get a list of CSV file paths
csv_files = get_files(outputs_directory, template_name)

# Number of threads to use (you can adjust this based on your system's capabilities)
num_threads = 2

PROPTAB = properties_extraction(num_threads, csv_files, lines_of_interest)

for tab in PROPTAB:
    df = PROPTAB[tab]
    df = df.sort_values(by='Galaxy ID').reset_index(drop=True)

# Loading DAP - NSA table for obtaining NSA properties.
nsafile = basepath+'/Auxiliary/manga_DR17_nsa.txt'
ALL = pd.read_csv(nsafile, delim_whitespace=True)

# Depurate to the fitted galaxies and sort by galaxy ID:
GALAXIES = PROPTAB['H-β']['Galaxy ID']
ALL = ALL[ALL['plateifu'].isin(GALAXIES)]
ALL = ALL.sort_values(by='plateifu').reset_index(drop=True)

# Extracting the max and min sigmas:
# List of files:
par_files = get_files(outputs_directory, '_results.csv')

# Creating the array for the sigmas:
SIGMAS = pd.DataFrame()
SIGMAS['plateifu'] = None
for line in lines_of_interest:
    key1, key2 = line+' s_min', line+' s_max'
    SIGMAS[key1], SIGMAS[key2] = None, None

# Now extracting and storing them:

for i, file in enumerate(par_files):
    match = re.search(r'manga-(\d+)-(\d+)', file)
    if match:
        plate, ifu = match.groups()
        galaxy_id = f'{plate}-{ifu}'
        SIGMAS.loc[galaxy_id, 'plateifu'] = galaxy_id  # Set the index value
    par_df = pd.read_csv(file)
    for line in lines_of_interest:
        sigmas = extract_max_min_values(par_df, line, 'Sigma (km/s)')
        key1, key2 =  line+' s_max', line+' s_min'
        SIGMAS.at[galaxy_id, key1] = sigmas[0]
        SIGMAS.at[galaxy_id, key2] = np.nan if len(sigmas) == 1 else sigmas[1]

SIGMAS.reset_index(drop=True, inplace=True)
SIGMAS = SIGMAS.apply(pd.to_numeric, errors='coerce')



# Importing the stellar mass:
STELLAR_MASS = ALL['nsa_sersic_mass']
Z = ALL['nsa_z']
LOG_STELLAR_MASS = np.log10(STELLAR_MASS)

ALL['OIII/Hb'] = np.log10(PROPTAB['O-III,2']['Flux']/PROPTAB['H-β']['Flux'])
ALL['NII/Ha'] = np.log10(PROPTAB['N-II,2']['Flux']/PROPTAB['H-α']['Flux'])
ALL['SII/Ha'] = np.log10(PROPTAB['S-II,2']['Flux']/PROPTAB['H-α']['Flux'])
ALL['OI/Ha'] = np.log10(PROPTAB['O-I,1']['Flux']/PROPTAB['H-α']['Flux'])

# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# ..................................................................................................
# ::::::::::::::::::::::::::::::::::  ASTROPHYSICAL MODELS  ::::::::::::::::::::::::::::::::::::::::
# ··································································································
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

# COSMOLOGICAL PARAMETERS: _________________________________________________________________________
Om_mat_cero = 0.307
Om_lambda = 0.693
Om_baryons = 0.048
sigma8 = 0.829
h = 0.678
delta_c = 1.686

Cosmology = [0, Om_mat_cero, Om_lambda, Om_baryons, sigma8, h, delta_c]

# :::::::::::::::::::::::::::  ASTROPHYSICAL CALCULATIONS ::::::::::::::::::::::::::::::::::::::::::

def log_Mvir(Ms, logM1, logMs0, b, d, g):
    '''
    Calculates the mean logMvir (log of Halo Mass) in function of the stellar mass
    '''
    Ms0 = 10**(logMs0) # solar masses
    MS = Ms/Ms0
    B = b*np.log10(MS)
    C = ((MS)**d)/(1+MS**(-g))
    return logM1 + B + C - (1/2)

def Vmax(log10Mh, z):
    '''
    Return the maximal velocity of equilibrium Vmax inside de galaxy-halo potential.
    '''
    log10Mh = log10Mh + np.log10(Cosmology[5]);
    afac = 1 / (1 + z);
    beta = 2.209059 + 0.059565 * afac -0.021009 * afac * afac;
    alpha = 0.346032 -0.058864 * afac + 0.024920 * afac * afac;

    log10vmax = beta + (log10Mh-12) * alpha + 0.5 * np.log10(Cosmology[2] + Cosmology[1] \
        * pow(1+z,3)) * alpha

    return pow(10, log10vmax)

def Mvir(Ms):
    '''
    Returns directly the Virial Mass of the Halo in function of the stellar mass, according to the
    z = 0.10 values of the table 6 from Rodríguez-Puebla et al., 2017
    '''
    logMvir = log_Mvir(Ms, 12.58, 10.90, 0.48, 0.29, 1.52)
    Mvir = 10**logMvir
    return Mvir

def a0r(Z):
    '''
    Returns the a0r factor given the redshift (z << zeq)
    '''
    c_Ho = 4420 #Mpc
    fz = lambda z: 1/(Cosmology[2]+Cosmology[1]*((1+z)**3))**(1/2)
    a0r = c_Ho*int.quad(fz, 0, Z)
    return a0r

def luminosity_distance(Z):
    '''
    Returns the Luminosity distance given a redshift.
    '''
    DL = cosmo.luminosity_distance(Z)
    DL = np.array([DL[i].value for i in range(len(Z))])
    DL = DL*3.086e+24 #Mpc to cm
    return DL

def logOH(NII_HA):
    '''
    Return the log([O/H]) abundance given the log([NII]6584/Halpha) ratio, given by two calibrators
    '''
    logOH1 = 8.90 + 0.57*NII_HA - 12
    return logOH1

def Luminosity(Flux, DL):
    '''
    Returns the Luminosity given the flux and the luminosity distance
    '''
    Lux = (4*np.pi*(DL**2)*Flux*1e-17)
    return Lux

def Moutflow(FOIII, MET, DL, ne3):
    '''
    Returns an estimate of the outflow mass according to the equation of Cano-Díaz et al, 2012,
    given the flux of the OIII line.
    '''
    LOIII = (4*np.pi*(DL**2)*FOIII*1e-17)/(1e44)
    MOUT = 5.33e7*(LOIII/ne3*(10**MET)) # Solar Masses
    return MOUT


def Mout_dot(FOIII, MET, DL, ne3, Vout, Rout):
    '''
    Returns an estimate of the outflow mass ejection rate
    '''
    LOIII = (4*np.pi*(DL**2)*FOIII*1e-17)/(1e44)
    #eLOIII = (4*np.pi*(DL**2)*eFOIII*1e-17)/(1e44)
    MOUTDOT = 164*((LOIII*Vout*1e-3)/(ne3*(10**MET)*Rout))
    #eMOUTDOT = 164*(((eLOIII*Vout+LOIII*eVout)*1e-3)/(ne3*(10**MET)*Rout)) #Solar Masses/yr

    return MOUTDOT

def SFR_Ha(FHA, DL):
    '''
    Returns an estimate of the SFR (Star Formation Rate) by Kennicut et al., 1994
    '''
    LHA = (4*np.pi*(DL**2)*FHA*1e-17)
    #eLHA = (4*np.pi*(DL**2)*eFHA*1e-17)
    SFR = 7.9e-42*LHA
    #eSFR = 7.9e-42*eLHA
    return SFR

def MBH_Ha(FHA, FWHM_HA, DL):
    LHA = (4*np.pi*(DL**2)*FHA*1e-17)/(1e42)
    FWHM_HA = FWHM_HA/1e3
    MBH = (2.0)*1e6*(LHA**(0.55))*(FWHM_HA**(2.06))
    return MBH

def MBH_Hb(FHB, FWHM_HB, DL):
    LHB = (4*np.pi*(DL**2)*FHB*1e-17)/(1e42)
    FWHM_HB = FWHM_HB/1e3
    MBH = (3.6)*1e6*(LHB**(0.56))*(FWHM_HB**(2.0))
    return MBH


# RELATIONS FROM LITERATURE: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Star formation Main Sequence relation: -----------------------------------------------------------
def SFR_M(logMstar, ncurve):
    '''
    The fit for the Star Formation Main Sequence, taken from Cano Díaz et al., 2019
    '''
    if ncurve == 1:
        da, db = 0.0, 0.0
    elif ncurve == 2:
        da, db = 0.15, 0.01
    elif ncurve == 3:
        da, db = -0.15, -0.01
    logSFR = -(7.64+da) + (0.74+db)*logMstar
    return logSFR

def sSFR_M(logMstar, ncurve):
    '''
    The fit for the Star Formation Main Sequence, taken from Cano Díaz et al., 2019
    '''
    if ncurve == 1:
        da, db = 0.0, 0.0
    elif ncurve == 2:
        da, db = 0.15, 0.01
    elif ncurve == 3:
        da, db = -0.15, -0.01
    logsSFR = -(7.64+da) - (0.26+db)*logMstar
    return logsSFR

# Black Hole Mass Scaling relation: ----------------------------------------------------------------

def logMBH(logMStar, ncurve, alpha, dalpha, beta, dbeta):
    '''
    The fit for Black Hole Mass - Stellar Mass Relation, taken from Reines & Volonteri, 2016
    '''
    if ncurve == 1:
        dalpha, dbeta = 0.0, 0.0
    elif ncurve == 2:
        dalpha, dbeta = dalpha, dbeta
    elif ncurve == 3:
        dalpha, dbeta = -dalpha, -dbeta
    logMBH = (alpha + dalpha) + (beta+dbeta)*logMStar
    return logMBH

# Black Hole Mass vs Stellar Sigma relation: -------------------------------------------------------

def MBH_sigma(sigma):
    MBH = (10**(8.12))*((sigma/200)**(4.24))
    return MBH


# CALCULATIONS WITH ALDO'S PROCEDURE (Mgas COOLING) ------------------------------------------------
# Functions to calculate M_cooling (from Aldo Rodriguez): ------------------------------------------

# Constants: ---------------------------------------------------------------------------------------
v_c = 299792 #km/s
T_sfr = 150E6
epsilon_rad = 0.1
G = 4.299E-6 #kpc/Msol*(km/s)^2
H0 = 100        #Today's Hubble constant km/s/Mpc
f_bar = 0.16
Om_mat = 0.3
Om_lam = 1. - Om_mat
h = 0.7
Ks_cte = -4.40909
n_KS = 1.4
mu = 0.59 #number weighted proton mass
mp = 1.673E-24 # g proton mass
kB = 1.381E-16 # erg / K #Boltzmann constant
E_ers = 1.989E43 # 1 Msol km^2 /s^2 = 1.989E^43 erg
Lsolar = 3.82710E33 #ergs/s
yr = 31556926


# --------------------------------------------------------------------------------------------------
def Delta_vir(Om_mat_cero,Om_lambda_cero,z):

    x = Om_mat_cero*(1+z)**3 / ( Om_mat_cero* (1+z)**3 + Om_lambda_cero ) - 1.

    return ( 18 * math.pi * math.pi + 82*x - 39*x*x) / (1+x)

def Om_m(Om_mat,Om_lambda,z):        #Omega Matter
    return Om_mat * ( 1. + z )**3/( Om_lambda + Om_mat * ( 1. + z )**3 )

def Om_l(Om_mat,Om_lambda,z):        #Omega lambda
    return Om_lambda / ( Om_lambda + Om_mat * ( 1. + z )**3 )

def H(Om_mat,Om_lambda,z):
    return H0 * np.sqrt( Om_lambda + Om_mat * ( 1. + z )**3 )    #Hubble constant km/s/Mpc

def rho_crit(Om_mat,Om_lambda,z):
    return 3 * H(Om_mat,Om_lambda,z)**2 / 8 / math.pi / ( G / 1E3 )    #Critical density M_sun / Mpc**3

def rho_m(Om_mat,Om_lambda,z):
    return Om_m(Om_mat,Om_lambda,z)*rho_crit(Om_mat,Om_lambda,z)


def Rvir(Mvir,z):   #virial radius in kpc
    log10Mvir = np.log10(Mvir)
    rv = log10Mvir - np.log10( 4. * math.pi * rho_m(Om_mat,Om_lam,z) * Delta_vir(Om_mat,Om_lam,z) / 3.)
    rv = rv / 3. + 3;
    return 10**rv / h

def t_dyn(z): #dynamical time in seconds
    time = G * 1E-3 * (3.2407792896664E-20)**2 * rho_m(Om_mat,Om_lam,z) * Delta_vir(Om_mat,Om_lam,z)
    return 1./np.sqrt(time) / 31556926 / 1E9

def dMvir_dt(Mh,z): #Dynamical time averaged halo accreation rate RP16 in Msol/yr
    log10Mh = np.log10(Mh)
    log10Mh = log10Mh + np.log10(h);
    afac = 1 / (1 + z);
    beta = 2.730 - 1.828 * afac + 0.654 * afac * afac
    alpha = 1. + 0.329 * afac + 0.206 * afac * afac
    logdMvir_dt = beta + (log10Mh-12) * alpha + 0.5 * np.log10(Om_lam + Om_mat * (1+z)**3);
    return 10**logdMvir_dt

################################Gaseous halo##########################################

def facc_hot(Mvir,x_rvir): #fraction of baryons that were accreated and are hot. Correa, Schaye, van de Voort et al. 2018
    M_12 = 1E12*10**(-0.15)
    x = Mvir/M_12
    frac = 1 + x**(-1.86)
    return x_rvir/frac

def facc_cold(Mvir,x_rvir): #fraction of baryons that were accreated and are cold. #Correa, Schaye, van de Voort et al. 2018
    return 1.- facc_hot(Mvir,x_rvir)

def epsilon_gas(Mvir,x_rvir):
    return facc_hot(Mvir,x_rvir) + facc_cold(Mvir,1.)

def Mhot(Mvir): #Hot mass for a given halo mass.
    x = np.log10(Mvir/1E12)
    fhot = -0.79 + 0.52 * x - 0.05 * x * x
    fhot = 10**fhot
    fhot = np.where(fhot >= 1, 1, fhot)
    return fhot * f_bar * Mvir

def Mcold_flow(Mvir): #Cold mass for a given halo mass.
    fcold = 1 - Mhot(Mvir) / f_bar / Mvir
    return fcold * f_bar * Mvir

def Thalo(logVmax): #Gas temperature in hydrostatic equilibrium in K
    return 35.9 * 10**( 2.* logVmax)

def t_cool(Th,cooling_rate,nH): #cooling_rate is in units of erg cm^-3 * s^-1 and nH of cm^-3
    cooling_rate = cooling_rate / nH / nH
    return 1.5 * (Th/1E6) * (1E-3/nH) * (1E-23/cooling_rate)


def x_cool(Mvir,Th,cooling_rate,nH,z): #cooling radius in terms of the virial radius
    cooling_rate = cooling_rate / nH / nH
    td = t_dyn(z) * 31556926 * 1e9 #sec
    M_hot = Mhot(Mvir) * 1.99E33 #g
    rv = Rvir(Mvir,z) * (1E3 * 3.086E18) #cm^3
    rv = rv**3
    x_up = td * M_hot * cooling_rate # s * g * erg cm^3 * s^-1
    x_down = 6 * math.pi * mu * mp * kB * Th * rv  # g * erg / K * K * cm^3
    x = np.sqrt(x_up / x_down)
    x = np.where(x > 1, 1, x)
    return x

def Mcooling(Mvir,x_cooling,z): #Cooling mass
    return Mhot(Mvir) * x_cooling

def Mcold_gas(Mvir,x_cooling,z): #Total cold gas mass
    return Mcold_flow(Mvir) + Mhot(Mvir) * x_cooling

# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# ··································································································
# ::::::::::::::::::::::::::::::: CALCULATIONS WITH THE DATA  ::::::::::::::::::::::::::::::::::::::
# ..................................................................................................
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

# lines_of_interest = ['O-III,2', 'H-β', 'O-I,1', 'N-II,2', 'H-α', 'S-II,2']
# Luminosity distance:

LD = luminosity_distance(Z)
ALL['Luminosity distance'] = LD

# Virial mass:
VIRIAL_MASS = Mvir(STELLAR_MASS)
LOG_VIRIAL_MASS = np.log10(VIRIAL_MASS)
ALL['Mvir'] = VIRIAL_MASS

# Vmax for dataset and Vmax curve:
# Dataset:
VMAX = Vmax(LOG_VIRIAL_MASS, Z)
ALL['Vmax'] = VMAX
# Curve:
sm_mock = np.linspace(8.0, 12.0, 10000)
virmock = Mvir(10**sm_mock)
Z_0 = np.random.uniform(min(Z), max(Z), 10000)
VMAX_CURVE =(sm_mock, np.log10(Vmax(np.log10(virmock), Z_0)))

# [OIII] outflow mass:
ne3 = 1.0
OH = logOH(ALL['NII/Ha'])
ALL['log OH'] = OH
OH = 10**OH
MOUT_OIII = Moutflow(PROPTAB['O-III,2']['Flux'], OH, LD, ne3)
ALL['Mout ([OIII])'] = MOUT_OIII
# [OIII] mass-out rate:
f_reff = 0.5
Vout = SIGMAS['O-III,2 s_max']
Rout = f_reff * ALL['nsa_sersic_th50']
ALL['Mout_dot ([OIII])'] = Mout_dot(PROPTAB['O-III,2']['Flux'], OH, LD, ne3, Vout, Rout)

# OIII Luminosity:
ALL['L ([OIII])'] = Luminosity(PROPTAB['O-III,2']['Flux'], LD)

# Star Formation Rate:
SFR = SFR_Ha(PROPTAB['H-α']['Flux'], LD)
ALL['SFR'] = SFR
ALL['sSFR'] = SFR/STELLAR_MASS

# Black hole mass estimation:
ALL['M_BH'] = MBH_Hb(PROPTAB['H-β']['Flux'], PROPTAB['H-β']['W80 (km/s)'], LD)

# One cubic interpolation: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
file_path = basepath+'/Auxiliary/RadLoss_-1.00_0.47_0.00_1.00e-04.rates'
# Define the column names
column_names = ["temperature", "cooling_rate", "heating_rate"]
# Read the table into a pandas DataFrame
df = pd.read_csv(file_path, sep="\t", names=column_names, skiprows=2)
cooling_rate = CubicSpline(df.temperature, df.cooling_rate - df.heating_rate)

# Cool gas analysisis: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ALL['Rvir'] = Rvir(ALL['Mvir'], 0.)
ALL['T_Halo'] = Thalo(ALL['Vmax'])
ALL['lambda_cool'] = cooling_rate(ALL['T_Halo'])
ALL['x_cool'] = x_cool(ALL['Mvir'], ALL['T_Halo'], ALL['lambda_cool'], 1E-4, 0.)
ALL['x_cool'] = np.where(ALL['x_cool'] > 10.*ALL['nsa_sersic_th50']/ALL['Rvir'],
    ALL['x_cool'], 10.*ALL['nsa_sersic_th50']/ALL['Rvir'])
ALL['Mcold_flow'] = Mcold_flow(ALL['Mvir'])
ALL['Mcooling'] = Mcooling(ALL['Mvir'],ALL['x_cool'], 0.)
ALL['Mcold'] = Mcold_gas(ALL['Mvir'], ALL['x_cool'], 0.) - ALL['nsa_sersic_mass']
ALL['epsilon_gas'] = epsilon_gas(ALL['Mvir'], ALL['x_cool'])
ALL['Mcold_rate'] = ALL['epsilon_gas'] * f_bar * dMvir_dt(ALL['Mvir'], 0.)

# ··································································································
# ::::::::::::::::::::::::::::::::::::::::  PLOTS  :::::::::::::::::::::::::::::::::::::::::::::::::
# ··································································································

# Some generic function for produce scatter plots:

def horizontal_scatter_plots(entire_sample, y_label, x_labels,
        n_subplots_per_row=3, x0=None, y0=None, extra_subsets=None, y_range=None, x_ranges=None,
        title=None):

    sns.set(style="ticks")
    sns.set_context("paper", font_scale=1.2)

    n_plots = len(x_labels)
    n_rows = (n_plots - 1) // n_subplots_per_row + 1
    fig, axs = plt.subplots(n_rows, n_subplots_per_row, sharey=True, gridspec_kw={'wspace': 0.0},
        figsize=(13, 5 * n_rows))

    # Flatten the axs array if it has multiple rows
    axs = axs.flatten()
    tot_plots = len(axs)

    row_counter = 0
    for i in range(n_plots):
        ax = axs[i]
        # If this subplot is not used, make it invisible
        # Scatter plot for entire sample
        y_entire = entire_sample[0]
        x_entire = entire_sample[1][i]
        ax.scatter(x_entire, y_entire, label='Sample', color='deepskyblue', marker='H', alpha=0.7)

        # Plot straight line
        if x0 is not None:
            ax.axvline(x=x0, color='red', linestyle='--', label=f'x0 = {x0}')
        if y0 is not None:
            ax.axhline(y=y0, color='green', linestyle='--', label=f'y0 = {y0}')

        if y_range:
            ax.set_ylim(y_range)
        if x_ranges and i < len(x_ranges):
            ax.set_xlim(x_ranges[i])

        # Customize plot
        ax.set_xlabel(x_labels[i])
        if i == row_counter:
            ax.set_ylabel(y_label)
            row_counter += n_subplots_per_row
        if i == 0: ax.legend(frameon=False)

    # Plot extra subsets
    if extra_subsets is not None:
        for l in range(len(extra_subsets)):
            subset_l = extra_subsets[l]
            for i in range(len(subset_l)):
                x = subset_l[i][0]
                y = subset_l[i][1]
                label = subset_l[i][2]
                color = subset_l[i][3]
                linestyle = subset_l[i][4][0]
                marker = subset_l[i][4][1]
                axs[l].plot(x, y, label=label, color=color, marker=marker, linestyle=linestyle,
                    linewidth=2)

    for i in range(tot_plots):
        if i >= n_plots:
            axs[i].set_visible(False)

    if title:
        fig.suptitle(title)

    # Adjust layout
    plt.tight_layout()

    # Show the plots
    return fig

def vertical_scatter_plots(entire_sample, x_label, y_labels,
                            n_subplots_per_col=3, x0=None, y0=None, extra_subsets=None,
                            x_range=None, y_ranges=None, title=None):

    sns.set(style="ticks")
    sns.set_context("paper", font_scale=1.2)

    n_plots = len(y_labels)
    if n_plots < n_subplots_per_col: hplots = n_plots
    else: hplots = n_subplots_per_col
    n_cols = (n_plots - 1) // n_subplots_per_col + 1
    fig, axs = plt.subplots(n_subplots_per_col, n_cols, sharex=True,
                            gridspec_kw={'hspace': 0.0}, figsize=(5.5 * n_cols, 4*hplots))

    # Flatten the axs array if it has multiple columns
    axs = axs.flatten()
    tot_plots = len(axs)

    for i in range(n_plots):
        ax = axs[i]

        # Scatter plot for entire sample
        y_entire = entire_sample[1][i]
        x_entire = entire_sample[0]
        ax.scatter(x_entire, y_entire, label='Sample', color='deepskyblue', marker='H', alpha=0.7)

        # Plot straight line
        if x0 is not None:
            ax.axvline(x=x0, color='red', linestyle='--', label=f'x0 = {x0}')
        if y0 is not None:
            ax.axhline(y=y0, color='green', linestyle='--', label=f'y0 = {y0}')

        if x_range:
            ax.set_xlim(x_range)
        if y_ranges and i < len(y_ranges):
            ax.set_ylim(y_ranges[i])

    # Plot extra subsets
    if extra_subsets is not None:
        for l in range(len(extra_subsets)):
            subset_l = extra_subsets[l]
            for i in range(len(subset_l)):
                x = subset_l[i][0]
                y = subset_l[i][1]
                label = subset_l[i][2]
                color = subset_l[i][3]
                linestyle = subset_l[i][4][0]
                marker = subset_l[i][4][1]
                if len(subset_l) > 5:
                    alpha = subset_l[i][5]
                else:
                    alpha = 1.0
                axs[l].plot(x, y, label=label, color=color, marker=marker, linestyle=linestyle,
                            linewidth=2, alpha=alpha)
                axs[l].legend(frameon=False)

    # Customize plot
    col_count = n_plots - n_cols
    for i in range(n_plots):
        if i >= col_count:
            axs[i].set_xlabel(x_label)
        axs[i].set_ylabel(y_labels[i])
        if i == 0: axs[i].legend(frameon=False)
        if i >= n_plots:
            axs[i].set_visible(False)

    if title:
        fig.suptitle(title)

    # Adjust layout
    plt.tight_layout()

    # Show the plots
    return fig


def simple_scatter_plot(main_set, y_label, x_label, extra_subsets=None, y_range=None, x_range=None,
    x0=None, y0=None, title=None):
    sns.set(style="ticks")
    sns.set_context("paper", font_scale=1.2)

    fig, ax = plt.subplots(figsize=(5, 4))

    # Scatter plot for main set
    y_main = main_set[1]
    x_main = main_set[0]
    ax.scatter(x_main, y_main, label='Main Set', color='deepskyblue', marker='H', alpha=0.7)

    # Plot extra subsets
    if extra_subsets is not None:
        for subset in extra_subsets:
            x_extra = subset[0]
            y_extra = subset[1]
            label_extra = subset[2]
            color_extra = subset[3]
            linestyle_extra = subset[4][0]
            marker_extra = subset[4][1]
            ax.plot(x_extra, y_extra, label=label_extra, color=color_extra, marker=marker_extra,
                    linestyle=linestyle_extra, linewidth=2)

    # Plot straight lines
    if x0 is not None:
        ax.axvline(x=x0, color='gray', linestyle='--', label=f'x0 = {x0}')
    if y0 is not None:
        ax.axhline(y=y0, color='black', linestyle='--', label=f'y0 = {y0}')

    # Customize plot
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)

    if y_range:
        ax.set_ylim(y_range)
    if x_range:
        ax.set_xlim(x_range)

    if title:
        fig.suptitle(title)

    ax.legend(frameon=False)

    # Adjust layout
    plt.tight_layout()

    # Show the plot
    return fig

# Define the analysis results path:
respath = basepath+f'/Analysis/{specdir}'
# Create directory if this not exists:
if not os.path.exists(respath):
    os.makedirs(respath)

# Todas las líneas:
# sigma vs Mstar, sigma / vmax vs Mstar plots: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
X0 = np.log10(ALL['nsa_sersic_mass'])
x_label = r'log $M_\star$ ($M_\odot$)'

# Comparison sigma stellar curve, from Zahid et al., 2018:

def logSigma_Zahid(logMs, logSb, logMb, a1, a2):
    A = np.where(logMs > logMb, a2, a1)
    result = A * logMs + (logSb - A * logMb)
    return result

logSb = 2.073
logMb = 10.26
a1 = 0.403
a2 = 0.293

log_Sigma_Stellar = logSigma_Zahid(sm_mock, logSb, logMb, a1, a2)

EXP_CURVE = (sm_mock, log_Sigma_Stellar-VMAX_CURVE[1]) # log (sigma_stellar / Vmax)

# Vmax curve and stellar sigma:

vmax_curve = {
    0:[[VMAX_CURVE[0], VMAX_CURVE[1], r'log $V_{max}$ (M$_\star$, z); Rodriguez-Puebla et al. (2017)', 'darkorange', ('-', ''), 0.4]],
    1:[[EXP_CURVE[0], EXP_CURVE[1], r'log ($\sigma_\star$ / V$_{max}$) (from Zahid et al, 2018)', 'crimson', ('-', ''), 0.4]]}


for line in lines_of_interest:
    # Convert the column to numeric, replacing non-numeric values with NaN
    column_data_numeric = pd.to_numeric(SIGMAS[f'{line} s_max'], errors='coerce')
    # Apply log10 only to non-NaN values
    Y1 = np.log10(column_data_numeric)
    Y2 = np.log10(column_data_numeric / ALL['Vmax'])
    entire_sample = [X0, (Y1, Y2)]
    title = f'{line}'
    y_labels = [r'log $\sigma_{max}$ (kms$^{-1}$)', r'log ($\sigma_{max}$/$V_{max}$)']

    # Rest of your plotting code
    fig = vertical_scatter_plots(entire_sample, x_label, y_labels, extra_subsets=vmax_curve,
        x_range=[8.5, 12.], n_subplots_per_col=2)

    fig.savefig(respath + f'/sigma_velmax_{line}.png', format='png')


# A & K plots for OIII, Hb, and Ha lines: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
lines_toplot = ['O-III,2', 'H-β', 'H-α']

x_labels = [r'log W$_{80}$ (kms$^{-1}$)', r'log $V_{med}$ (kms$^{-1}$)',
    r'log $\sigma_{max}$ (kms$^{-1}$)']

for line in lines_toplot:
    X1 = np.log10(PROPTAB[line]['W80 (km/s)'])
    X2 = np.log10(PROPTAB[line]['V_med (km/s)'])
    column_data_numeric = pd.to_numeric(SIGMAS[f'{line} s_max'], errors='coerce')
    X3 = np.log10(column_data_numeric)
    # A figure:
    Y0 = PROPTAB[line]['A']
    entire_sample = [Y0, (X1, X2, X3)]
    y_label = 'A'
    fig = horizontal_scatter_plots(entire_sample, y_label, x_labels)
    fig.savefig(respath+f'/A_vs_vel_{line}.png', format='png')
    # K figure:
    Y0 = PROPTAB[line]['K']
    entire_sample = [Y0, (X1, X2, X3)]
    y_label = 'K'
    fig = horizontal_scatter_plots(entire_sample, y_label, x_labels)
    fig.savefig(respath+f'/K_vs_vel_{line}.png', format='png')

# sSFR vs OIII parameters: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Y0 = np.log10(ALL['sSFR'])
X1 = np.log10(SIGMAS['O-III,2 s_max'])
X2 = np.log10(ALL['L ([OIII])'])
X3 = np.log10(ALL['Mout ([OIII])'])
entire_sample = [Y0, (X1, X2, X3)]

# Ranges:
y_range = [-13., -8.]
x_ranges = ([1., 4.], [36., 44.], [0., 8.])

y_label = r'log sSFR (yr$^{-1}$)'
x_labels = [r'log $\sigma_{[OIII]}$ (kms$^{-1}$)', r'log L$_{[OIII]}$ (erg/s)',
    r'log M$_{[OIII]}$ ($M_\odot$)']
fig = horizontal_scatter_plots(entire_sample, y_label, x_labels, y0=-9.8, y_range=y_range,
    x_ranges=x_ranges)
fig.savefig(respath+f'/sSFR_vs_OIII.png', format='png')

# sigma & MOIII vs z: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

X0 = ALL['nsa_z']
Y1 = np.log10(ALL['Mout ([OIII])'])
Y2 = np.log10(SIGMAS['O-III,2 s_max']/ALL['Vmax'])
entire_sample = [X0, (Y1, Y2)]

x_label = r'z'
y_labels = [r'log M$_{[OIII]}$ ($M_\odot$)', r'log ($\sigma_{[OIII]}$/$V_{max}$)']
fig = vertical_scatter_plots(entire_sample, x_label, y_labels, n_subplots_per_col=2)
fig.savefig(respath+f'/MOII_sigma_vs_z.png', format='png')

# SFR vs Ms: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
X0 = np.log10(ALL['nsa_sersic_mass'])
Y0 = np.log10(ALL['SFR'])

x_label = r'log M$_\star$ (M$_\odot$)'
y_label = r'log SFR (M$_\odot$/yr)'
main_set = [X0, Y0]

x_range = [8.0, 12.0]
# SFR main sequence:
logsfr_mock = SFR_M(sm_mock, 1)
sfr_sequence = [[sm_mock, logsfr_mock, r'Star Formation Main Sequence (Cano-Diaz et al, 2019)',
    'seagreen', ('--', ' ')]]

fig = simple_scatter_plot(main_set, y_label, x_label, extra_subsets=sfr_sequence, x_range=x_range)
fig.savefig(respath+f'/SFR.png', format='png')

# sigma OIII vs L OIII: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
X0 = np.log10(ALL['L ([OIII])'])
Y0 = np.log10(SIGMAS['O-III,2 s_max'])
main_set = [X0, Y0]

x_label = r'log L$_{[OIII]}$ (erg/s)'
y_label = r'log $\sigma_{OIII}$ (kms$^{-1}$)'

fig = simple_scatter_plot(main_set, y_label, x_label)
fig.savefig(respath+f'/sigmaOIII_vs_LOIII.png', format='png')

# sigma OIII vs L OIII: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
X0 = np.log10(ALL['L ([OIII])'])
Y0 = np.log10(SIGMAS['O-III,2 s_max']/ALL['Vmax'])
main_set = [X0, Y0]

x_label = r'log L$_{[OIII]}$ (erg/s)'
y_label = r'log ($\sigma_{OIII}$/$V_{max}$)'

fig = simple_scatter_plot(main_set, y_label, x_label)
fig.savefig(respath+f'/sigmaOIIIVmax_vs_LOIII.png', format='png')

# Mstar vs M OIII: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
X0 = np.log10(ALL['Mout ([OIII])'])
Y0 = np.log10(ALL['nsa_sersic_mass'])
main_set = [X0, Y0]

x_label = r'log M$_{[OIII]}$ ($M_\odot$)'
y_label = r'log M$_\star$ (M$_\odot$)'

fig = simple_scatter_plot(main_set, y_label, x_label)
fig.savefig(respath+f'/Ms_vs_MOIII.png', format='png')

# BPT: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# BPT limits:
def BPT_lim_N(log_n2):
    log_o3 = ((0.61)/(log_n2 - 0.47)) + 1.19
    return log_o3

def BPT_lim_S(log_s2):
    log_o3 = ((0.72)/(log_s2 - 0.32)) + 1.30
    return log_o3

def BPT_lim_O(log_oI):
    log_o3 = ((0.73)/(log_oI + 0.59)) + 1.33
    return log_o3


# BPT curves:
x01 = np.linspace(-3.0, 0.3, 10000)
x02 = np.linspace(-3.0, -0.7, 10000)
y01 = BPT_lim_N(x01)
y02 = BPT_lim_S(x01)
y03 = BPT_lim_O(x02)

# Data:
Y0 = ALL['OIII/Hb']
X1 = ALL['NII/Ha']
X2 = ALL['SII/Ha']
X3 = ALL['OI/Ha']
entire_sample = [Y0, (X1, X2, X3)]
y_label = r'log ([OIII]/H$\beta$)'
x_labels = [r'log ([NII]/H$\alpha$)', r'log ([SII]/H$\alpha$)', r'log ([OI]/H$\alpha$)']

# Ranges:
y_range = [-2.5, 2.5]
x_ranges = ([-2.5, 2.5], [-2.5, 2.5], [-3., 1.5])

bpt_curves = {
    0:[[x01, y01, '', 'tomato', ('-', '')]],
    1:[[x01, y02, '', 'tomato', ('-', '')]],
    2:[[x02, y03, r'Keweley et al. (2001)', 'tomato', ('-', '')]]}

fig = horizontal_scatter_plots(entire_sample, y_label, x_labels, extra_subsets=bpt_curves,
    y_range=y_range, x_ranges=x_ranges)
fig.savefig(respath+f'/BPT.png', format='png')

# Figuras Aldo:
############################### M cool plot ########################################################

fig, axs = plt.subplots(1,figsize=(10,8))
# Set tick font size
for label in (axs.get_xticklabels() + axs.get_yticklabels()):
    label.set_fontsize(20)

vmin = -3.5
vmax = -2.5
clrmap = 'coolwarm'
edg_col = 'gray'

sc_plt = axs.scatter(ALL['nsa_sersic_mass'], ALL['Mcold_rate'] / ALL['SFR'], c= np.log10( ALL['Mout_dot ([OIII])']), cmap=clrmap, marker='o', linewidth=3, s=50,vmin=vmin, vmax=vmax)

cbar = plt.colorbar(sc_plt,ax=axs, extend='max')

#cbar.set_label(r'$\log (M_{\rm BH} / M_{\ast})$',fontsize=20)
cbar.set_label(r'log($\dot{M}_{\rm [OIII]}$)',fontsize=20)

cbar.ax.tick_params(labelsize=15)  # Set the label font size

#axs.scatter(ALL['Ms'], ALL['Mcold_rate'] / ALL['SFR'], c=np.log10(ALL['M_bh'] / ALL['Ms'] ), cmap=clrmap, marker='H', linewidth=3, s=50, vmin=vmin, vmax=vmax)

#axs.scatter(OUTFLOWS['Ms'], OUTFLOWS['Mcold_rate'] / OUTFLOWS['SFR'], c=np.log10( OUTFLOWS['Mbh'] / OUTFLOWS['Ms']), cmap=clrmap, marker='D', s=50, linewidth=3, vmin=vmin, vmax=vmax)

axs.set_ylabel(r'$\epsilon_{\rm cold,gas}f_{\rm bar} \dot{M}_{\rm vir,dyn}/{\rm SFR}$',fontsize=20)
axs.set_xlabel(r'$M_{\rm \ast} \; [{\rm M}_\odot]$',fontsize=20)
axs.legend(loc='lower left',fontsize=20, frameon=False)
axs.set_yscale("log")
axs.set_xscale("log")
axs.axis([1E8,1E12,1E-2,1E5])

plt.savefig(respath+'/Mcool_ratio.png', format='png', bbox_inches='tight')


## plot 2: ---------------------------------------------------------------------------------------
fig, axs = plt.subplots(1,figsize=(10,8))
# Set tick font size
for label in (axs.get_xticklabels() + axs.get_yticklabels()):
    label.set_fontsize(20)

sc_plt = axs.scatter(ALL['nsa_sersic_mass'],
    (ALL['Mcold_rate']-ALL['sSFR']) / ALL['Mout_dot ([OIII])'],
    c=np.log10( ALL['Mout_dot ([OIII])']), cmap=clrmap, marker='o', linewidth=3,
    s=50,vmin=vmin, vmax=vmax)

cbar = plt.colorbar(sc_plt,ax=axs, extend='max')

#cbar.set_label(r'$\log (M_{\rm BH} / M_{\ast})$',fontsize=20)
cbar.set_label(r'log($\dot{M}_{\rm [OIII]}$)',fontsize=20)
cbar.ax.tick_params(labelsize=15)  # Set the label font size

#axs.scatter(ALL['Ms'], ALL['Mcold_rate'] / ALL['SFR'], c=np.log10(ALL['M_bh'] / ALL['Ms'] ), cmap=clrmap, marker='H', linewidth=3, s=50, vmin=vmin, vmax=vmax)

# axs.scatter(OUTFLOWS['Ms'], (OUTFLOWS['Mcold_rate']-OUTFLOWS['sSFR']) / OUTFLOWS['Mout_dot'],
#    c=np.log10(OUTFLOWS['Mbh'] / OUTFLOWS['Ms']), cmap=clrmap, marker='D', s=50, linewidth=3,
#    vmin=vmin, vmax=vmax)

axs.set_ylabel(r'$\epsilon_{\rm cold,gas}f_{\rm bar} \dot{M}_{\rm vir,dyn} - SFR/{\dot{M}_{\rm OIII}}$',fontsize=20)
axs.set_xlabel(r'$M_{\rm \ast} \; [{\rm M}_\odot]$',fontsize=20)
axs.legend(loc='lower left',fontsize=20, frameon=False)
axs.set_yscale("log")
axs.set_xscale("log")
axs.axis([1E8,1E12,1E-2,1E5])


plt.savefig(respath+'/Mcool_Mout_{}.png', format='png', bbox_inches='tight')


## plot 3: MOIII vs MBH ------------------------------------------------------------------------

vmin = -12.
vmax = -9.

fig, axs = plt.subplots(1,figsize=(10,8))
# Set tick font size
for label in (axs.get_xticklabels() + axs.get_yticklabels()):
    label.set_fontsize(20)

sc_plt = axs.scatter(ALL['M_BH'], ALL['Mout_dot ([OIII])'], c=np.log10(ALL['sSFR']), cmap=clrmap,
    marker='o', linewidth=3, s=50,vmin=vmin, vmax=vmax)

cbar = plt.colorbar(sc_plt,ax=axs, extend='max')

cbar.set_label(r'$\log ({\rm sSFR})$',fontsize=20)
cbar.set_label(r'log(SFR)',fontsize=20)

cbar.ax.tick_params(labelsize=15)  # Set the label font size

#axs.scatter(OUTFLOWS['M_BH'], OUTFLOWS['Mout_dot'], c=np.log10(OUTFLOWS['sSFR']), cmap=clrmap,
#    marker='D', s=50, linewidth=3, vmin=vmin, vmax=vmax)

axs.set_ylabel(r'$\dot{M}_{\rm OIII}$',fontsize=20)
axs.set_xlabel(r'$M_{\rm BH} \; [{\rm M}_\odot]$',fontsize=20)
axs.legend(loc='lower left',fontsize=20, frameon=False)
axs.set_yscale("log")
axs.set_xscale("log")
axs.axis([1E6,1E10,1E-5,1E+00])


plt.savefig(respath+'/Mbh_MOIIIdot.png', format='png', bbox_inches='tight')


print('FINISHED!')
