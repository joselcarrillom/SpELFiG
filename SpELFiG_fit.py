# LIBRARIES:
import random
import numpy as np
import pandas as pd
import emcee
import SpELFiG_Analyticals as spl
import scipy.integrate as spi
import matplotlib.pyplot as plt

# from astropy.modeling.models import Voigt1D
from astropy.io import fits
from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import curve_fit
from scipy.stats import truncnorm

nan = np.nan
inf = np.inf

class SpELFiG_emcee():
	'''
	This class performs only a final mcmc round, taking a dataframe of parameters with a fixed
	number of components per line, and iterates over the associated spectral model in order to
	achieve the final model
	'''

	def __init__(self, spec, linelist, priors_dataframe, continuum0, niter):

		self.spec = spec
		self.nlambda = len(spec[:,0])
		self.linelist = linelist
		self.dfparams = priors_dataframe
		self.ndim = len(priors_dataframe)*3 + 3
		self.nwalkers = 2*self.ndim
		self.continuum0 = continuum0
		self.niter = niter
		self.data = (spec[:,0], spec[:,1], spec[:,2])
		self.minspec = min(self.data[0])
		self.maxspec = max(self.data[0])
		self.maxflux = max(self.data[1])

		# Execute the fit:
		self.dataframe_to_theta()
		self.best_fit, self.errors = self.mcmc_chains()

		# Create the atributes with the dataframes of results:
		self.goodness = self.goodness()
		# self.line_errors = self.calculate_line_errors()

		self.results_df = self.dataframe_results()
		self.lines_properties_df = self.calculate_line_properties()

	def dataframe_to_theta(self):
		# Extract all line parameters:
		theta_init = self.dfparams[['Centroid', 'Amplitude', 'Sigma']].values.flatten()
		self.lines_priors_min = self.dfparams[['min_loc', 'minflux', 'min_sd']].values.flatten()
		self.lines_priors_max = self.dfparams[['max_loc', 'maxflux', 'max_sd']].values.flatten()

		# Add the continuum:
		theta_init = np.concatenate((theta_init, self.continuum0))
		# components = self.prev_results.groupby('Line Name')['Component'].max()
		# self.components = components.tolist()

		# Generate the theta0 array:
		self.theta0 = np.array([theta_init + 1e-3 * np.random.randn(self.ndim) for _ in range(self.nwalkers)])


	def spectral_model(self, theta, x):
		# Initialize the model with an empty array
		model = np.zeros_like(x)

		# Start index for unpacking theta
		param_start = 0
		for i, line in enumerate(self.linelist):
			Ncomp = line['Ncomp']
			# So, the counter of parameters for this line:
			params_per_line = Ncomp*3
			# Extract parameters for the current component
			component_params = theta[param_start:param_start + params_per_line]

			model += spl.composed_gaussians(Ncomp, x, component_params)

			# Move the start index to the next component
			param_start += params_per_line

		# Add the continuum to the model
		continuum = theta[param_start:]
		model += spl.continuum_function(x, *continuum)

		self.model = model

		return self.model

	# EMCEE METHODS : %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	def lnprob(self, theta, x, y, dy):
		g = self.spectral_model(theta, x)
		return -0.5 * np.sum(((g - y) / dy)**2)

	def chi2(self, theta, x, y, dy):
		'''
		Calculates the chi squared of a model
		'''
		g = self.spectral_model(theta, x)
		return np.sum(((y - g) / dy)**2)

	def bic(self, xi_s, k, n):
		'''
		Calculates the bayesian information criterion
		'''
		BIC = xi_s + k*np.log(n)
		if n < 50:
			BIC = BIC + (2*k*(k+1))/(n-k-1)
		return BIC

	def residuals(self, theta, data):
		model = self.spectral_model(theta, data[0])
		return data[1] - model

	def lnprior(self, theta):
		# if self.prior_list_df(theta) is True:
		#	return 0.0
		if self.prior_list_df(theta) is True:
			return 0.0
		return -inf

	def log_probability(self, theta, x, y, dy):
		lp = self.lnprior(theta)
		if not np.isfinite(lp):
			return -inf
		g = self.lnprob(theta, x, y, dy)
		return lp + g

	def prior_list_linelist(self, theta):
		for i, line in enumerate(self.linelist):
			a_sum = 0
			Ncomp = line['Ncomp']
			for j in range(Ncomp):
				loc, a0, sd = theta[i * Ncomp + j * 3: i * Ncomp + (j + 1) * 3]
				a_sum += a0
				if sd < 0:
					sd = -sd
					theta[i * Ncomp + j * 3 + 2] = sd  # Update theta with the positive sigma
					# Check other conditions (min_loc, max_loc, maxflux, etc.)
				if not (line['min_loc'] <= loc <= line['max_loc'] and
					line['min_sd'] <= sd <= line['max_sd'] and a0 >= 0.0 and
					a_sum <= line['max_flux']):
					return False

				# Check the ratio conditions in case of the doublets:

				#if line['name'] in self.doublets_info:

					#ratio_info = self.doublets_info[line['name']]
					#amplitude_ratio = ratio_info['amplitude_ratio']
					#sigma_ratio = ratio_info['sigma_ratio']

					#if abs(a0 / sd - amplitude_ratio / sigma_ratio) > self.noise_scale:
						##print('This the invalid parameter combnination ratio check')
						##print(abs(a0 / sd - amplitude_ratio / sigma_ratio))
						#return False  # Invalid parameter combinatio

		return True

	def prior_list_df(self, theta):
		lines_params = theta[:-3]
		for i, param in enumerate(lines_params):
			min = self.lines_priors_min[i]
			max = self.lines_priors_max[i]
			if not (min <= param <= max):
				return False
		continuum = theta[-3:]
		cont_mins = [0.0, self.minspec, -2.0]
		cont_maxs = [self.maxflux, self.maxspec, 2.0]
		for i, param in enumerate(continuum):
			min, max = cont_mins[i], cont_maxs[i]
			if not (min <= param <= max):
				return False

				# Check the ratio conditions in case of the doublets:

				#if line['name'] in self.doublets_info:

					#ratio_info = self.doublets_info[line['name']]
					#amplitude_ratio = ratio_info['amplitude_ratio']
					#sigma_ratio = ratio_info['sigma_ratio']

					#if abs(a0 / sd - amplitude_ratio / sigma_ratio) > self.noise_scale:
						##print('This the invalid parameter combnination ratio check')
						##print(abs(a0 / sd - amplitude_ratio / sigma_ratio))
						#return False  # Invalid parameter combination
		return True


	def main(self, theta0, logprob, args):
		sampler = emcee.EnsembleSampler(self.nwalkers, self.ndim, logprob, args=args, moves=[
			(emcee.moves.DEMove(), 0.8),
			(emcee.moves.DESnookerMove(), 0.2),])

		#theta0 = theta0.reshape(self.nwalkers, self.ndim)
		#std = 0.1 * np.ones_like(theta0)
		#initial_positions = emcee.utils.sample_ball(theta0, std=std, size=self.nwalkers)

		##print("Running burn-in...")
		p0, prob, state = sampler.run_mcmc(theta0, 100)
		sampler.reset()

		##print("Running production...")
		pos, prob, state = sampler.run_mcmc(p0, self.niter)

		return sampler, pos, prob, state

	def mcmc_chains(self):
		args = [*self.data]
		# #print('Finding the parameters for the best model ')
		# Running the chains:
		sampler, pos, prob, state = self.main(self.theta0, self.log_probability, args)
		samples = sampler.flatchain
		# Best model:
		theta_max  = samples[np.argmax(sampler.flatlnprobability)]
		# Errors in the model:
		theta_errors = np.std(samples, axis=0)
		# Saving the ENTIRE samples in a file:
		self.CHAINS = pd.DataFrame(samples)
		self.theta_max = np.array(theta_max)
		self.continuum = self.theta_max[self.ndim-3:]
		self.theta_errors = np.array(theta_errors)
		self.samples = samples

		return self.theta_max, self.theta_errors

	def goodness(self):
		x2 = self.chi2(self.theta_max, *self.data)  # Use the correct log-likelihood method
		x2r = x2 / (self.nlambda - self.ndim)
		BIC = self.bic(x2, self.ndim, self.nlambda)
		self.goodness = {'Chi-Squared': x2,
			'Reduced Chi-Squared': x2r,
			'BIC': BIC,
			}
		return self.goodness

	def calculate_line_errors(self):
		print('~ Calculating fit errors ...')
		# Create a dictionary to store the results
		results_dict = {
			'Line Name': [],
			'Component': [],
			'err_Centroid': [],
			'err_Amplitude': [],
			'err_Sigma': [],
			'err_Sigma (km/s)': []
		}

		param_start = 0
		for line in self.linelist:
			Ncomp = line['Ncomp']
			params_per_line = Ncomp*3
			LINE = self.CHAINS.iloc[:, param_start:param_start+params_per_line]
			for j in range(Ncomp):
				results_dict['Line Name'].append(line['name'])
				results_dict['Component'].append(j+1)
				# All chains parameters for this component
				centroids = LINE.iloc[:, j*3]
				amplitudes = LINE.iloc[:, j*3 + 1]
				sigmas = LINE.iloc[:, j*3 + 2]
				# Calculations of st. deviations:
				err_centroid = np.std(centroids)
				err_amplitude = np.std(amplitudes)
				err_sigma = np.std(sigmas)
				sigmas_kms = spl.vel_correct_sigma(sigmas, centroids)
				err_sigma_kms = np.std(sigmas_kms)
				# Append the values:
				results_dict['err_Centroid'].append(err_centroid)
				results_dict['err_Amplitude'].append(err_amplitude)
				results_dict['err_Sigma'].append(err_sigma)
				results_dict['err_Sigma (km/s)'].append(err_sigma_kms)
			param_start += params_per_line

		self.lines_errors = pd.DataFrame(results_dict)
		# print('This errors dataframe')
		# print(len(self.lines_errors))

		return self.lines_errors


	def dataframe_results(self):
		'''
		Final dataframe of results
		'''
		# Create a dictionary to store the results
		results_dict = {
			'Line Name': [],
			'Component': [],
			'Centroid': [],
			'err_Centroid': [],
			'Amplitude': [],
			'err_Amplitude': [],
			'Sigma': [],
			'err_Sigma': [],
			'Sigma (km/s)': [],
			'err_Sigma (km/s)': [],
		}

		# Populate the dictionary with the results
		param_start = 0
		for i, line in enumerate(self.linelist):
			Ncomp = line['Ncomp']
			params_per_line = Ncomp*3
			for j in range(Ncomp):
				results_dict['Line Name'].append(line['name'])
				results_dict['Component'].append(j+1)
				# Extracting parameters:
				centroid = self.theta_max[param_start + j*3]
				amplitude = self.theta_max[param_start + j*3 + 1]
				sigma = self.theta_max[param_start + j*3 + 2]
				sigma_kms = spl.vel_correct_sigma(sigma, centroid)
				# Extracting errors:
				err_centroid = self.theta_errors[param_start + j*3]
				err_amplitude = self.theta_errors[param_start + j*3 + 1]
				err_sigma = self.theta_errors[param_start + j*3 + 2]
				err_sigma_kms = spl.vel_correct_sigma(err_sigma, centroid)
				# Appending results:
				results_dict['Centroid'].append(centroid)
				results_dict['Amplitude'].append(amplitude)
				results_dict['Sigma'].append(sigma)
				results_dict['Sigma (km/s)'].append(sigma_kms)
				results_dict['err_Centroid'].append(err_centroid)
				results_dict['err_Amplitude'].append(err_amplitude)
				results_dict['err_Sigma'].append(err_sigma)
				results_dict['err_Sigma (km/s)'].append(err_sigma_kms)
			param_start += params_per_line

		self.lines_parameters = pd.DataFrame(results_dict)

		# Append the errors:

		goodness_df = pd.DataFrame({
			'Line Name': ['Chi-squared', 'Reduced Chi-squared', 'BIC'],
			'Component': [self.goodness['Chi-Squared'], self.goodness['Reduced Chi-Squared'],
				self.goodness['BIC']],
			'Centroid': ['', '', ''],
			'err_Centroid': ['a (continuum)', 'l_0 (continuum)', 'b (continuum)'],
			'Amplitude': [self.continuum[0], self.continuum[1], self.continuum[2]],
			'err_Amplitude': ['', '', ''],
			'Sigma': ['', '', ''],
			'err_Sigma': ['', '', ''],
			'Sigma (km/s)': ['', '', ''],
			'err_Sigma (km/s)': ['', '', ''],
			})
#
		#Concatenate the header DataFrame and the results DataFrame
		self.results_df = pd.concat([self.lines_parameters, goodness_df], ignore_index=True)

		return self.results_df

	def calculate_line_properties(self):
		results = []
		param_start = 0
		x_mock = np.linspace(min(self.data[0]), max(self.data[0]), 10000)
		for i, line in enumerate(self.linelist):
			line_name = line['name']
			Ncomp = line['Ncomp']
			params_per_line = Ncomp*3
			line_params = self.theta_max[param_start:param_start+params_per_line]
			flux = spl.intflux_gauss(line_params, Ncomp)
			# Calculate the line profile using the composed_gaussians function
			l = line_params.flatten()
			#print('Flatten', l)
			line_profile = spl.composed_gaussians(Ncomp, x_mock, line_params.flatten())

			# Calculate integrated flux using intflux_gauss
			# flux0 = spi.simps(line_profile, x_mock)
			# Calculate FWHM based on the line profile
			max_amplitude = np.max(line_profile)
			half_max = max_amplitude / 2.0
			indexes = np.where(line_profile >= half_max)[0]
			fwhm = x_mock[indexes[-1]] - x_mock[indexes[0]]

			# Calculate Equivalent Width (EW)
			# EW = Flux / (peak amplitude)
			peak_amplitude = max_amplitude
			ew = flux / peak_amplitude

			# Calculate the mean centroid:
			loc_c = spl.mean_loc(line_params, Ncomp)

			# Calculate radial velocity, W80, W90, A, and K:
			v_r = spl.V_flux(Ncomp, line_params, 0.5, loc_c)[0] # radial velocity
			v_90 = spl.V_flux(Ncomp, line_params, 0.9, loc_c)[0]
			v_10 = spl.V_flux(Ncomp, line_params, 0.1, loc_c)[0]
			v_95 = spl.V_flux(Ncomp, line_params, 0.95, loc_c)[0]
			v_05 = spl.V_flux(Ncomp, line_params, 0.05, loc_c)[0]
			W80, W90 = v_90 - v_10, v_95 - v_05
			A = ((v_90 - v_r) - (v_r - v_10)) / W80
			K = W90 / (1.397 * fwhm)

			results.append({
				'Line Name': line_name,
				'Flux': flux,
				'FWHM (Angstrom)': fwhm,
				'Equivalent Width (Angstrom)': ew,
				'V_med (km/s)': v_r,
				'W80 (km/s)': W80,
				'A': A,
				'K': K
			})

			param_start += params_per_line

		# Convert the list of dictionaries into a Pandas DataFrame
		self.properties_df = pd.DataFrame(results)

		return self.properties_df

# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# ·········································  SCIPY CLASS   ·········································
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

class SpELFiG_Scipy:

	def __init__(self, spectrum, linelist, prev_dataframe=None, continuum0=None,
		line_first_estimations=True):
		self.spec = spectrum
		self.data = (self.spec[:,0], self.spec[:,1], self.spec[:,2])
		self.minspec = min(self.data[0])
		self.maxspec = max(self.data[0])
		self.maxflux = max(self.data[1])

		self.linelist = linelist
		if continuum0 is not None:
			self.continuum = continuum0
		if line_first_estimations is True:
			self.dfparams, self.continuum = self.estimate_parameters()
		elif line_first_estimations is False:
			self.dfparams = prev_dataframe
		self.fit_spectral_model()
		self.results_df, self.lines_parameters = self.dataframe_results()
		self.lines_properties_df = self.calculate_line_properties()

	def estimate_parameters(self):

		def continuum_params(x_continuum, y_continuum):
			a_guess = abs(np.mean(y_continuum))
			loc0_guess = np.mean(x_continuum)
			g_guess = -1.0
			p0 = [a_guess, loc0_guess, g_guess]
						# Parameter bounds
			bounds = ([0.0, self.minspec, -2.0], [self.maxflux, self.maxspec, 2.0])
			params, _ = curve_fit(spl.continuum_function, x_continuum, y_continuum, p0=p0,
				bounds=bounds)
			return params

		def narrow__component_est(x, y, dy, loc0, max_width, threshold=0.4):
			# Find local maxima
			x = np.squeeze(x)  # Ensure x is 1-D
			y = np.squeeze(y)  # Ensure y is 1-D
			peaks, _ = find_peaks(y, height=threshold * np.max(y))
			# Emergency values in case the above does not work:
			sigma_max = max_width / (2.0 * np.sqrt(2.0 * np.log(2.0)))
			if not peaks.size:
				# Handle the case where no valid maxima are found
				return [loc0, np.max(y), random.uniform(0.5, 1.)*sigma_max]

			# Find the index of the maximum closest to loc0
			max_idx = peaks[np.argmin(np.abs(x[peaks] - loc0))]
			if max_idx >= len(x) or max_idx < 0:
				# Handle the case where the index is out of bounds
				return [loc0, np.max(y), random.uniform(0.5, 1.)*sigma_max]

			a_max = y[max_idx]
			loc = x[max_idx]

			# Define the range for FWHM or sigma calculation
			range_left = max(0, max_idx - 1)
			range_right = min(len(x) - 1, max_idx + 1)
			# Iterate outward based on the initial range and zero leve
			while (range_left > 0 and y[range_left] > a_max * 0.5 and
				x[range_right] - x[range_left] <= max_width):
				range_left -= 1
			while (range_right < len(x) - 1 and y[range_right] > a_max * 0.5 and
				x[range_right] - x[range_left] <= max_width):
				range_right += 1

			# Determine whether to use FWHM or sigma based on symmetry and maximum width
			fwhm = x[range_right] - x[range_left]
			if fwhm > max_width:
				fwhm = max_width
			sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))

			return [loc, a_max, sigma]

		initial_params = np.zeros((len(self.linelist), 1, 3))
		continuum_data = np.copy(self.spec)
		Ncomps = []
		line_names = []
		min_locs = []
		max_locs = []
		minfluxes = []
		maxfluxes = []
		min_sds = []
		max_sds = []

		for i, line in enumerate(self.linelist):
			# Define the restframe wavelength.
			center_wavelength = line['wavelength']
			line_name = line['name']
			line_names.append(line_name)
			Ncomps.append(line['Ncomp'])
			min_locs.append(line['min_loc'])
			max_locs.append(line['max_loc'])
			minfluxes.append(0.0)
			maxfluxes.append(line['max_flux'])
			min_sds.append(line['min_sd'])
			max_sds.append(line['max_sd'])

			wline = line['max_sd']*(2.0 * np.sqrt(2.0 * np.log(2.0)))
			selec_line = np.where((self.spec[:, 0] >= center_wavelength - wline*0.5) &
				(self.spec[:, 0] <= center_wavelength + wline*0.5))

			x = self.spec[selec_line, 0].flatten()
			y = self.spec[selec_line, 1].flatten()
			dy = self.spec[selec_line, 2]
			line['max_flux'] = max(y)
			continuum_data[selec_line, 1] -= y
			# Estimate Gaussian parameters for non-continuum lines
			lineparams = narrow__component_est(x, y, dy, center_wavelength, wline)
			initial_params[i] = lineparams

		def array_to_dataframe(data_array):
			# Flatten the array
			flat_array = data_array.reshape((len(data_array), -1))
			# Define column names without the 'err' columns
			column_names = ['Centroid', 'Amplitude', 'Sigma']
			# Create a DataFrame
			df = pd.DataFrame(flat_array, columns=column_names)
			# Display the resulting DataFrame
			return df

		priors_dataframe = array_to_dataframe(initial_params)
		priors_dataframe['Component'] = Ncomps
		priors_dataframe['Line Name'] = line_names
		priors_dataframe['min_loc'] = min_locs
		priors_dataframe['max_loc'] = max_locs
		priors_dataframe['minflux'] = minfluxes
		priors_dataframe['maxflux'] = maxfluxes
		priors_dataframe['min_sd'] = min_sds
		priors_dataframe['max_sd'] = max_sds

		continuum0 = continuum_params(continuum_data[:, 0], continuum_data[:, 1])

		return priors_dataframe, continuum0


	def spectral_model(self, x, *params):
		# Initialize the model with an empty array
		model = np.zeros_like(x)
		# Start index for unpacking theta
		param_start = 0
		# Case where the params object does not have the desired input form:
		if not len(params) > 1:
			params = np.concatenate(params)
			params = params.flatten()
		for i, line in enumerate(self.linelist):
			Ncomp = line['Ncomp']

			# So, the counter of parameters for this line:
			params_per_line = Ncomp*3

			# Extract parameters for the current component
			component_params = params[param_start:param_start + params_per_line]

			# Add the Gaussian component to the model
			model += spl.composed_gaussians(Ncomp, x, component_params)

			# Move the start index to the next component
			param_start += params_per_line

		# Add the continuum to the model
		model += spl.continuum_function(x, *self.continuum)
		self.model = model

		return self.model

	def fit_spectral_model(self):

		# self.ncomp_plus = int(self.ncomp + self.n_extra)
		x = self.data[0]
		y = self.data[1]

		# Get initial guesses from the parameter prior dataframe:
		initial_guess = self.dfparams[['Centroid', 'Amplitude', 'Sigma']].values.flatten()
		bounds_min = self.dfparams[['min_loc', 'minflux', 'min_sd']].values.flatten()
		bounds_max = self.dfparams[['max_loc', 'maxflux', 'max_sd']].values.flatten()
		bounds = (bounds_min, bounds_max)


		# VERIFY PARAMETER AND BOUNDS ONE BY ONE:
		for i, param in enumerate(initial_guess):
			if not (bounds_min[i] <= param <= bounds_max[i]):
				print('&&&&& H E R E &&&&&&&')
				print(f'{bounds_min_flat[i]},  {param},  {bounds_max_flat[i]}')

		# PERFORMING THE FIT: ººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººº
		s_ftol = 0.01  # 1% error tolerance for the cost function (chi-squared)
		s_xtol = 0.01  # 1% error tolerance for the parameters

		# Calculating the spectral curve:
		popt, _ = curve_fit(self.spectral_model, x, y, p0=initial_guess, bounds=bounds,
			method='trf', ftol=s_ftol, xtol=s_xtol) # maxfev=6000 try later

		# Calculating residuals and goodness parameters:
		residuals = y - self.spectral_model(x, *popt)
		ssr = np.sum(residuals ** 2)
		n = len(x)
		p = len(popt)
		mse = ssr / (n - p)
		goodness = {
			'Reduced Chi-Squared': mse,
			'BIC': mse + p * np.log(n),
		}

		# Organize results into a DataFrame
		# self.linelist = results['Line Name'].unique()
		self.goodness = goodness
		self.best_fit = popt
		# ººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººººº

	def dataframe_results(self):
		# Create a dictionary to store the results
		results_dict = {
			'Line Name': [],
			'Component': [],
			'Centroid': [],
			'Amplitude': [],
			'Sigma': [],
			}

		# Populate the dictionary with the results from popt
		p_index = 0  # Initialize the index for popt
		for line in self.linelist:
			Ncomp = line['Ncomp']
			# Now the storage loop:
			for j in range(Ncomp):
				# Extract values from popt
				centroid = self.best_fit[p_index]
				amplitude = self.best_fit[p_index + 1]
				sigma = self.best_fit[p_index + 2]

				# Append results to the dictionary
				results_dict['Line Name'].append(line['name'])
				results_dict['Component'].append(j + 1)
				results_dict['Centroid'].append(centroid)
				results_dict['Amplitude'].append(amplitude)
				results_dict['Sigma'].append(sigma)

				# Increment the popt index
				p_index += 3

		# Create a pandas DataFrame from the dictionary
		results_df = pd.DataFrame(results_dict)
		params_ranges_df = results_df.copy()

		# Add the ranges which were used for the fit:
		params_ranges_df['min_loc'] = self.dfparams['min_loc']
		params_ranges_df['max_loc'] = self.dfparams['max_loc']
		params_ranges_df['minflux'] = self.dfparams['minflux']
		params_ranges_df['maxflux'] = self.dfparams['maxflux']
		params_ranges_df['min_sd'] = self.dfparams['min_sd']
		params_ranges_df['max_sd'] = self.dfparams['max_sd']

		return results_df, params_ranges_df

	def calculate_line_properties(self):
		results = []
		param_start = 0
		x_mock = np.linspace(min(self.data[0]), max(self.data[0]), 10000)
		for i, line in enumerate(self.linelist):
			line_name = line['name']
			Ncomp = line['Ncomp']
			params_per_line = Ncomp*3
			line_params = self.best_fit[param_start:param_start+params_per_line]
			flux = spl.intflux_gauss(line_params, Ncomp)
			# Calculate the line profile using the composed_gaussians function
			l = line_params.flatten()
			#print('Flatten', l)
			line_profile = spl.composed_gaussians(Ncomp, x_mock, line_params.flatten())

			# Calculate integrated flux using intflux_gauss
			# flux0 = spi.simps(line_profile, x_mock)
			# Calculate FWHM based on the line profile
			max_amplitude = np.max(line_profile)
			half_max = max_amplitude / 2.0
			indexes = np.where(line_profile >= half_max)[0]
			fwhm = x_mock[indexes[-1]] - x_mock[indexes[0]]

			# Calculate Equivalent Width (EW)
			# EW = Flux / (peak amplitude)
			peak_amplitude = max_amplitude
			ew = flux / peak_amplitude

			# Calculate the mean centroid:
			loc_c = spl.mean_loc(line_params, Ncomp)

			# Calculate radial velocity, W80, W90, A, and K:
			v_r = spl.V_flux(Ncomp, line_params, 0.5, loc_c)[0] # radial velocity
			v_90 = spl.V_flux(Ncomp, line_params, 0.9, loc_c)[0]
			v_10 = spl.V_flux(Ncomp, line_params, 0.1, loc_c)[0]
			v_95 = spl.V_flux(Ncomp, line_params, 0.95, loc_c)[0]
			v_05 = spl.V_flux(Ncomp, line_params, 0.05, loc_c)[0]
			W80, W90 = v_90 - v_10, v_95 - v_05
			A = ((v_90 - v_r) - (v_r - v_10)) / W80
			K = W90 / (1.397 * fwhm)

			results.append({
				'Line Name': line_name,
				'Flux': flux,
				'FWHM (Angstrom)': fwhm,
				'Equivalent Width (Angstrom)': ew,
				'V_med (km/s)': v_r,
				'W80 (km/s)': W80,
				'A': A,
				'K': K
				})

			param_start += params_per_line

		# Convert the list of dictionaries into a Pandas DataFrame
		self.properties_df = pd.DataFrame(results)

		return self.properties_df
