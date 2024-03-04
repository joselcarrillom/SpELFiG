#!/usr/bin/env python3
# Record the start time: ········································································
import time
start_time = time.time()

# º Manage data libraries: ·········································································
import os
import sys
import importlib.util
import pandas as pd
import numpy as np

# º Plots libraries : ·············································································
import matplotlib.pyplot as plt
import seaborn as sns
import random

# º Font configuration for the plots: ······························································
from matplotlib import rc
rc('text', usetex=True)
# rc('text.latex', preamble=r'\usepackage{amsmath} \usepackage{amssymb} \usepackage{amsfonts} \usepackage[scaled=1]{helvet} \usepackage[helvet]{sfmath} \everymath={\sf}')

# Spelfic classes to perform the fits: ·····························································
from SpELFiG_fit import SpELFiG_emcee as spl_emcee
from SpELFiG_fit import SpELFiG_Scipy as spl_scipy

# To supress warnings: ·············································································
import warnings
warnings.filterwarnings("ignore")

# Import the initial configuration objects: ························································

# Get the script path from the command-line argument
config_path = sys.argv[1]

# Import the config script of the galaxy as a module
spec = importlib.util.spec_from_file_location("config", config_path)
config = importlib.util.module_from_spec(spec)
spec.loader.exec_module(config)


# Import some variables from config: ·······························································
# Spectra and paths:
SPEC_I = config.spec
specrange0 = config.specranges[0]
respath0 = config.respath
outputname = config.outputname

# Lines dictionaries:
ALL = config.ALL_EML
BLR =   config.BLR_EML
OUTFLOW = config.OUT_EML

# Lines max component numbers:
NCOMP0 = config.maxncomp_nlr
N_EXTRA_BLR = config.maxncomp_blr
N_EXTRA_OUT = config.maxncomp_outflow
NITER = config.niter

width_NLR = config.max_width_narrow
width_BLR = config.max_width_blr
width_OUT = config.max_width_outflow

# Create the outputpath
if not os.path.exists(respath0):
	os.makedirs(respath0)

# Functions to produce the plots: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def specfitplot(specfit, zoom_range=None, goodness_marks=True, scipy=False):
	# Extract data from the DataFrame
	obs_data = specfit.data
	chi2 = specfit.goodness['Reduced Chi-Squared']
	BIC =  specfit.goodness['BIC']
	x_fit = np.linspace(min(obs_data[0]) - 10, max(obs_data[0]) + 10, 10000)
	y_fit = np.zeros_like(x_fit)

	# Create a figure with two subplots (upper and lower panels)
	fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 5), sharex=True, gridspec_kw={'hspace': 0.0, 'height_ratios': [5, 2]})

	# Upper panel: Observed spectrum, model, and individual Gaussian components
	for _, row in specfit.lines_parameters.iterrows():
		line_name = row['Line Name']
		component = int(row['Component'])
		amplitude = row['Amplitude']
		centroid = row['Centroid']
		sigma = row['Sigma']
		print()

		# Generate Gaussian component for the current row
		component_y = amplitude * np.exp(-((x_fit - centroid) / sigma)**2 / 2)
		# y_fit += component_y

		color = (random.random(), random.random(), random.random())
		ax1.plot(x_fit, component_y, linestyle='--', linewidth=0.6, color=color)

	# Add the total fitted spectrum
	params = specfit.best_fit.flatten()
	if scipy is True:
		y_fit = specfit.spectral_model(x_fit, params)
		y_fit_obs = specfit.spectral_model(obs_data[0], params)
	else:
		y_fit = specfit.spectral_model(params, x_fit)
		y_fit_obs = specfit.spectral_model(params, obs_data[0])
	ax1.plot(x_fit, y_fit, color='crimson', linewidth=2.0, label='Total Fitted Spectrum')
	ax1.errorbar(obs_data[0], obs_data[1], yerr=obs_data[2], color='black', linestyle='-',
		marker='.', alpha=0.6, markersize=2, linewidth=0.6, label='Observed Spectrum')

	# Set the y-axis label and legend for the upper panel
	ax1.set_ylabel(r'I [$10^{-17}$ erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]')
	ax1.legend(loc='upper right', frameon=False)

	# Add Chi-squared and BIC as labels with transparency:
	if goodness_marks is True:
		ax1.text(0.05, 0.85, f'Chi-squared: {chi2:.2f}', transform=ax1.transAxes, fontsize=12,
			color='gray', alpha=0.8)
		ax1.text(0.05, 0.78, f'BIC: {BIC:.2f}', transform=ax1.transAxes, fontsize=12,
			color='gray', alpha=0.8)


	residuals = (obs_data[1] - y_fit_obs) / y_fit_obs # Compute percentage residuals
	# ax2.axhline(0.0, color='k', linestyle='--', label='Zero Line')
	e_range = [0, 1]
	ax2.plot(obs_data[0], abs(residuals), marker='.', color='steelblue', linestyle='--', alpha=0.5,
		markersize=2, label='Percentage error')
	ax2.set_ylim(e_range)
	# ax2.axhline(np.mean(np.abs(residuals)), color='m', linestyle='--', label='Mean Residuals')

	# Set the x-axis and y-axis labels for the lower panel
	ax2.set_xlabel(r'$\lambda$ [$\AA$]')
	ax2.set_ylabel('Residuals')
	ax2.legend(loc='upper right', frameon=False)

	# Set the zoom range if provided
	if zoom_range:
		ax1.set_xlim(zoom_range)
		ax2.set_xlim(zoom_range)

	# Fine-tune the plot layout and remove top and right spines
	plt.tight_layout()
	ax1.spines['top'].set_visible(False)
	ax1.spines['right'].set_visible(False)
	ax2.spines['top'].set_visible(False)
	ax2.spines['right'].set_visible(False)

	return fig


def plot_mcmc_chains(specfit, zoom_range=None):

	x, y, dy = specfit.data[0], specfit.data[1], specfit.data[2]
	xfit = np.linspace(min(x), max(x), 10000)
	samples = specfit.samples

	fig, ax = plt.subplots(figsize=(8, 5))

	ax.errorbar(x, y, dy, label='Observed flux')
	for theta in samples[np.random.randint(len(samples), size=100)]:
		ax.plot(xfit, specfit.spectral_model(theta, xfit), color="purple", alpha=0.1)

	ax.set_xlabel(r'$\lambda$')
	ax.set_ylabel(r'Flux')
	ax.legend(loc='upper right', frameon=False)

	# Set the zoom range if provided
	if zoom_range:
		ax.set_xlim(zoom_range)

	plt.tight_layout()
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)

	return fig


# Create the linelist object: ······································································
'''
In this instance we create a list of dictionaries, one per emission line to be fitted in the
spectra. Each dictionary will store the necessary information for the fits in each line.
'''

def filter_lines_by_presence(spec, linelist, std_cont, mean_cont, window_width=20,
	min_snr_threshold=5.0):
	filtered_linelist = []

	for line in linelist:
		line_wavelength = line['wavelength']
		name = line['name']

		# Extract the data within the specified window
		window_data = spec[
			(spec[:,0] >= line_wavelength - window_width / 2) &
			(spec[:,0] <= line_wavelength + window_width / 2)]

		# Calculate the signal-to-noise ratio (SNR) based on the maximum flux and standard deviation
		max_flux = window_data[:,1].max() - mean_cont
		noise_std = std_cont

		snr = max_flux / noise_std if noise_std > 0 else 0

		# Check if the SNR is above the threshold
		if snr >= min_snr_threshold:
			filtered_linelist.append(line)

	return filtered_linelist

NCOMP0 = config.maxncomp_nlr
max_width_narrow = config.max_width_narrow


def create_linelist(eml_dict, wavelength_range, max_width, maxflux, ncomp):
	sigma_max = max_width / (2.0 * np.sqrt(2.0 * np.log(2.0)))
	min_wavelength, max_wavelength = wavelength_range
	linelist = []
	for name, wavelength in eml_dict.items():
		if min_wavelength <= wavelength[0] <= max_wavelength:
			lineloc_min = wavelength[0] - 2 * sigma_max
			lineloc_max =  wavelength[0] + 2 * sigma_max
			line = {'name': name, 'wavelength': wavelength[0],
				'min_loc': lineloc_min, 'max_loc': lineloc_max,
				'min_sd': 1.0, 'max_sd': sigma_max, 'max_flux': maxflux, 'Ncomp': ncomp}
			linelist.append(line)
	return linelist

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CYCLE OF FITS FUNCTIONS:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
delta_loc = config.delta_loc
f_sigma_blr = config.f_sigma_blr
f_sigma_out = config.f_sigma_out

def add_component(linelist, dfparams, linenames, max_width, f_sigma):
	'''
	Function that returns the updated parameter dataframe (to serve as initial guesses for a next
	fit) and the linelist updated object.
	'''
	sigma_max = max_width / (2.0 * np.sqrt(2.0 * np.log(2.0)))
	new_rows = []
	linelist_ = linelist[:]

	for i, line, in enumerate(linelist_):
		# Create the new component in the dataframe:
		line_name = line['name']
		if line_name in linenames:
			Ncomp = line['Ncomp']
			index = dfparams.index[(dfparams['Line Name'] == line['name']) &
				(dfparams['Component'] == line['Ncomp'])].tolist()[0]
			line_row = dfparams.loc[index]
			new_row = line_row.copy()

			# Update the necessary linelist elements:
			lineloc_min = line['wavelength'] - 2 * sigma_max
			lineloc_max =  line['wavelength'] + 2 * sigma_max
			line['Ncomp'] += 1
			if line['max_sd'] < sigma_max:
				line['max_sd'] = sigma_max
			if line['max_loc'] < lineloc_max:
				line['max_loc'] = lineloc_max
			if line['min_loc'] > lineloc_min:
				line['min_loc'] = lineloc_min

			# Update values in dataframe:
			amplitude_value = line_row['Amplitude'] / line['Ncomp']
			dfparams.at[index, 'Amplitude'] = amplitude_value
			new_row['Component'] = line['Ncomp']
			new_row['Centroid'] = new_row['Centroid'] + delta_loc
			new_row['Amplitude'] = amplitude_value
			new_row['Sigma'] = sigma_max * random.uniform(0.5, 1.0)
			new_row['min_loc'] = lineloc_min
			new_row['max_loc'] = lineloc_max
			new_row['minflux'] = 0.0
			new_row['maxflux'] = line['max_flux']
			new_row['min_sd'] = 1.0
			new_row['max_sd'] = sigma_max

			dfparams = pd.concat([dfparams.loc[:index], pd.DataFrame([new_row]),
				dfparams.loc[index + 1:]]).reset_index(drop=True)

	return dfparams, linelist_


def add_components_bybic(spec, linelist0, dfparams, continuum, n_extra_max, LINES_EML, max_width,
	f_sigma):
	# Updating the linelist object and the dfparams object:
	n_extra = 1
	extra_comp_lines = list(LINES_EML.keys())
	previousBIC = np.inf
	linelist = linelist0
	dfparams0 = dfparams

	while n_extra <= n_extra_max:
		# Add the component:
		dfparams, linelist_ = add_component(linelist, dfparams0, extra_comp_lines, max_width,
			f_sigma)
		# Perform the new fit:
		specfit = spl_scipy(spec, linelist_, prev_dataframe=dfparams, continuum0=continuum,
			line_first_estimations=False)
		# Evaluate goodness:
		goodness = specfit.goodness
		BIC = abs(goodness['BIC'])
		chi2 = goodness['Reduced Chi-Squared']
		DeltaBIC = BIC - previousBIC
		# Store the current fit in case the max number of components is reached:
		if (n_extra == n_extra_max and abs(DeltaBIC) >= config.DBIC_cut and DeltaBIC < 0.):
			specfit0 = specfit
		# Stop the cicle and return the fit according to the cut criteria:
		if (DeltaBIC >= 0. or abs(DeltaBIC) < config.DBIC_cut) or n_extra == n_extra_max:
			# Update correctly number of components.
			df = specfit0.lines_parameters
			for line in specfit0.linelist:
				line_name = line['name']
				component = df[df['Line Name'] == line_name]['Component'].max()
				line['Ncomp'] = component
			return specfit0

		# Update for next iteration:
		specfit0 = specfit
		dfparams0 = dfparams
		linelist = linelist_
		n_extra += 1
		previousBIC = BIC

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# ····························  &&& SPELFIC GLOBAL FLOW EXECUTION &&&  ·····························
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
'''
Now, from this instance we will execute the entire Spelfic Flow Algorithm: We will create the
consecutive adding of components to the different groups of lines with the scipy class; then, from
those results, generate the file dataframe that stores THE FINAL initial guesses, for the final
iteration: the MCMC iteration performed with the emcee class. If the file of the final initial
guesses dataframe ALREADY EXISTS by the execution of this script, the code goes directly to the
MCMC iteration depending on the REWRITE flag.
'''

# PRIOR_FILE = config.prior_file
# REWRITE_FLAG = config.rewrite_priors_flag

# spec_i = {'spectrum': spec, 'bin number': bin, 'S/N': S_N, 'RA': RA, 'DEC': DEC}
#        extracted_expectra.append(spec_i)

def run_cube_spelfig(SPEC):
	#  Number of bins:
	N_spec = len(SPEC)
	empty_tables = [[] for _ in range(len(config.Lines_for_global_tables))]
	bins_list = []
	for i in range(N_spec):
		spec = SPEC[i]['spectrum']
		n_bin = SPEC[i]['bin number']
		bins_list.append(n_bin)
		#sn = SPEC[i]['S/N']
		respath = respath0+f'/bin{n_bin}'
		if not os.path.exists(respath):
			os.makedirs(respath)
		outputname =  respath+f'/SpL-{config.ID_OBJ}_bin{n_bin}'
		PRIOR_FILE = respath+f'/priors-{config.ID_OBJ}_bin{n_bin}.csv'
		print(f'>>> Processing bin: {n_bin}')
		maxflux = max(spec[:,1])
		continuum_data0 = spec[(spec[:,0] >= config.continuum_data0[0]) &
			(spec[:,0] <= config.continuum_data0[1])]
		std_cont = continuum_data0[:,1].std()
		mean_cont = continuum_data0[:,1].mean()


		linelist_0 = create_linelist(ALL, specrange0, width_NLR, maxflux, NCOMP0)

		linelist_1 = filter_lines_by_presence(spec, linelist_0, std_cont, mean_cont)

		SPECFIT1 = spl_scipy(spec, linelist_1, line_first_estimations=True)
		# Extract the new parameters and linelist object:
		dfresults_1 = SPECFIT1.lines_parameters
		continuum = SPECFIT1.continuum
		# Now adding the components to the BLR and outflow lines:
		# print(' >>>>>>>   STEP TWO: Adding BLR components   <<<<<<<')
		SPECFIT2_BLR = add_components_bybic(spec, linelist_1, dfresults_1, continuum, N_EXTRA_BLR,
			BLR, width_BLR, f_sigma_blr)

		linelist_2, dfresults_2 = SPECFIT2_BLR.linelist, SPECFIT2_BLR.lines_parameters
		# print(' >>>>>>>   STEP THREE: Adding OUTFLOW components   <<<<<<<')
		SPECFIT2_OUT = add_components_bybic(spec, linelist_2, dfresults_2, continuum, N_EXTRA_OUT,
			OUTFLOW, width_OUT, f_sigma_out)

	# Finally extract the parameter dataframe, include the ranges, and save it in the PRIOR_FILE:
		RESULTSDF = SPECFIT2_OUT.lines_parameters
		LINELISTF = SPECFIT2_OUT.linelist
		CONTINUUMF = SPECFIT2_OUT.continuum

		if len(config.specranges[1]) > 0 and config.priors_plots_flag is True:
			for i, zoom in enumerate(config.specranges[1]):
				galfigi = specfitplot(SPECFIT2_OUT, zoom_range=zoom, scipy=True)
				galfigi.savefig(outputname+f'PRIORS_zoom_{i}.jpg', format='jpg')

		# STORE THE LINE RANGES IN THE DATAFRAME:
		# Create new columns
		RESULTSDF['min_sd'] = 0.0
		RESULTSDF['max_sd'] = 0.0
		RESULTSDF['max_flux'] = 0.0
		RESULTSDF['min_loc'] = 0.0
		RESULTSDF['max_loc'] = 0.0

		for line in LINELISTF:
			line_name = line['name']
			line_index = RESULTSDF.index[RESULTSDF['Line Name'] == line_name].tolist()

			# Check if the line is present in the dataframe:
			if line_index:
				for l in range(len(line_index)):
					l_index = line_index[l]
					# Update the values in the dataframe based on the linelist line
					RESULTSDF.at[l_index, 'min_sd'] = line['min_sd']
					RESULTSDF.at[l_index, 'max_sd'] = line['max_sd']
					RESULTSDF.at[l_index, 'max_flux'] = line['max_flux']
					RESULTSDF.at[l_index, 'min_loc'] = line['min_loc']
					RESULTSDF.at[l_index, 'max_loc'] = line['max_loc']

		# STORE THE CONTINUUM:
		RESULTSDF['Continuum'] = np.nan
		RESULTSDF.loc[:2, 'Continuum'] = CONTINUUMF
		dfparams_mcmc = RESULTSDF
		linelist_mcmc = LINELISTF
		# FINALLY STORE THE FILE
		RESULTSDF.to_csv(PRIOR_FILE, index=False)
		continuum_mcmc = dfparams_mcmc.loc[:2, 'Continuum']

		# Finally, the MCMC iteration: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		SPECFIT_FINAL = spl_emcee(spec, linelist_mcmc, dfparams_mcmc, continuum_mcmc, NITER)

		# Producing and saving the plots according to the flags:

		if len(config.specranges[1]) > 0:
			for i, zoom in enumerate(config.specranges[1]):
				if config.fit_plots_flag is True:
					galfigi = specfitplot(SPECFIT_FINAL, zoom_range=zoom)
					galfigi.savefig(outputname+f'_zoom_{i}.png', format='jpg')
				if config.samples_plots_flag is True:
					samplesfigi = plot_mcmc_chains(SPECFIT_FINAL, zoom_range=zoom)
					samplesfigi.savefig(outputname+f'samples_zoom_{i}.jpg', format='jpg')

		# Producing and saving final output files: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		SPECFIT_FINAL.results_df.to_csv(outputname + '_results.csv', index=False)
		if config.samples_file_flag is True:
			SPECFIT_FINAL.CHAINS.to_csv(outputname + '_samples.csv', index=False)

		if config.properties_file_flag is True:
			galfit_prop = SPECFIT_FINAL.lines_properties_df
			for i, line in enumerate(config.Lines_for_global_tables):
				# Appending global properties:
				gal_row = galfit_prop[galfit_prop['Line Name'] == line]
				empty_tables[i].append(gal_row)

	global_tables = [pd.concat(table) for table in empty_tables]
	# Save these tables:
	for i, table in enumerate(global_tables):
		line = config.Lines_for_global_tables[i]
		table['bin number'] = bins_list
		cols = ['bin number'] + [col for col in table.columns if col != 'bin number']
		table = table[cols]
		table.to_csv(respath+f'/{line}_global_properties.csv',
			sep=',', index=False, encoding='utf-8')

# Execute the fits for the entire cube:
run_cube_spelfig(SPEC_I)


'''
elif (os.path.exists(PRIOR_FILE)) & (REWRITE_FLAG is False):
	# This case implies the file already existed and the user wish not to rewrite (repeat the
	# previous fits). This allows the user to modify BY HAND the priors and ranges for the mcmc
	# iteration.
	dfparams_mcmc = pd.read_csv(PRIOR_FILE)
	linelist_mcmc = linelist_1
	max_component_rows = dfparams_mcmc.groupby('Line Name')['Component'].idxmax()
	df_filtered = dfparams_mcmc.loc[max_component_rows]
	for i, line in enumerate(linelist_mcmc):
		dfline = df_filtered.loc[df_filtered['Line Name'] == line['name']]
		line['min_loc'] = dfline['min_loc']
		line['max_loc'] = dfline['max_loc']
		line['min_sd'] = dfline['min_sd']
		line['max_sd'] = dfline['max_sd']
		line['max_flux'] = dfline['max_flux']
		line['Ncomp'] = dfline['Component']
'''

# Record the end time:
end_time = time.time()

# Calculate the elapsed time in seconds
elapsed_time_seconds = end_time - start_time
elapsed_time_minutes = elapsed_time_seconds / 60  # 60 seconds in a minute
elapsed_time_hours = elapsed_time_minutes / 60  # 60 minutes in an hour
print(f"Execution Time (minutes): {elapsed_time_minutes:.2f}")
print(f"Execution Time (hours): {elapsed_time_hours:.2f}")
