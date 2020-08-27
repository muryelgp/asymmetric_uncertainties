import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate
from scipy.optimize import curve_fit
import random
import sys

class aufloat:
	"""
	A class for dealing with measurments with asymmetric uncertanties

	e.g. Value = 10 (-1/+2)
	"""
	def __init__(self, nv=None, err_n=None, err_p=None, confidence= None, delta_chi2 = None, data=[], mod=None):
		"""
		Instantiates the class aufloat.
	   
		Parameters
		----------
		nv (nomianal value) : float
			The measured value, i.e. the mode of the distribution.

		err_n : float
			The negative error, the minus error (insert only the value, without the (-) sign).

		err_p : float
			The positive error, the plus error (insert only the value, without the (+) sign).

		confidence : int
			The statistical significance of the errors given by the percentile (%) confience intervals.

			The implemented options are 68, 90 and 95:

				Confidence(%)	Δχ2
				68%	1.00
				90%	2.71
				95%	4.00

		delta_chi2: float
			The user can either define the error as a pre-defined confidence interval (confidence parameter) or to define it using a Δχ2 value
			Values tested range from 0.5 < delta_chi2 < 5.0.

		data: list or np.array
			Another way to Instantiate the class is to input a data list with a sampled PDF. 
			Used to propagate uncertanties. See documetation for details.

		Returns
		-------
		Nothing.
		"""		
		import warnings


		warnings.filterwarnings("ignore", message="divide by zero encountered")
		warnings.filterwarnings("ignore", message="invalid value encountered")

		if not any(data):
			self.data = np.asarray([])
			self._creation_method = 'by_parameters'

			if ((err_n < 0) or (err_n is None)):
				raise RuntimeError("the 'err_n' value needs to be positive")
		   
			if ((err_n < 0) or (err_n is None)):
				raise RuntimeError("the 'err_p' value needs to be positive")

		elif((any(data)) & (nv is None) & (err_n is None) & (err_p is None)):
			self.data = np.asarray(data)
			self._creation_method = 'by_data'
			if(len(data) < 100):
				raise RuntimeError("'data' list needs to have more than 100 elements for a minimuum relliability of the Monte Carlo simulation. At least 1000 elements is more advisible though")
		else:
			raise RuntimeError('Something went wrong! You have to either insert the parameters (nm, err_n, err_p) or a list of values (data). You can not insert both!')



		try:
			self.delta_chi2 = self._get_delta_chi2(confidence)
			self._confidence = confidence
		except:
			try: 
				if delta_chi2 is not None:
					self.delta_chi2 = delta_chi2
					self._confidence = confidence
				else:
					print('Confidence interval value needs to be chosen among 68, 90, 95, or you have to insert a Δχ2 value. Assuming 68% (Δχ2 = 1.0)!.')
					confidence = 68
					self.delta_chi2 = self._get_delta_chi2(confidence)
					self._confidence = confidence
			except:
				pass


		self.nv= nv
		self.err_n = err_n
		self.err_p = err_p
		

		
		if(str(self._creation_method) == 'by_parameters'):
			self._N = int(1e4)
			delta = 4*np.max([self.err_n,self.err_p])/np.sqrt(self.delta_chi2)
			self.x_lim =  [self.nv - delta, self.nv + delta]
			self.x_values = np.linspace(self.x_lim[0], self.x_lim[1], self._N)
			

			self.norm = 1.0
			self.norm = self._calculate_norm()

			
			self.pdf_values = np.asarray(self.pdf(self.x_values))
			self.cdf_values = np.asarray(self._calculate_cdf_values())
			self.delta_chi2_curve = np.asarray(self.delta_chi2_func(self.x_values))
			
		
		elif str(self._creation_method) == 'by_data':
			self._N = self.data.size
			

			self._fit(guess = mod)
			#self.sigma_n, self.sigma_p = self.estimate()

			delta = 4*np.max([self.err_n,self.err_p])/np.sqrt(self.delta_chi2)
			self.x_lim =  [self.nv - delta, self.nv + delta]
			self.x_values = np.linspace(self.x_lim[0], self.x_lim[1], self._N)
			

			self.norm = 1.0
			self.norm = self._calculate_norm()

			
			self.pdf_values = np.asarray(self.pdf(self.x_values))
			self.cdf_values = np.asarray(self._calculate_cdf_values())
			self.delta_chi2_curve = np.asarray(self.delta_chi2_func(self.x_values))

			self.data = np.asarray([])
			

	def __str__(self):
		if self._confidence is None:
			output = "{:.2e} (-{:.2e}/+{:.2e}) [Δχ2 = {:.2f}]"
			return output.format(self.nv, self.err_n, self.err_p, self.delta_chi2)
		else:
			output = "{:.2e} (-{:.2e}/+{:.2e})\n{:.0f}% confidence level (Δχ2={:.2f})"
			return output.format(self.nv, self.err_n, self.err_p, self._confidence, self._get_delta_chi2(self._confidence))


	def __repr__(self):
		if self._confidence is None:
			output = "{:.2e} (-{:.2e}/+{:.2e}) [Δχ2 = {:.2e}]"
			return output.format(self.nv, self.err_n, self.err_p, self.delta_chi2)
		else:
			output = "{:.2e} (-{:.2e}/+{:.2e})\n{:.0f}% confidence level (Δχ2={:.2f})"
			return output.format(self.nv, self.err_n, self.err_p, self._confidence, self._get_delta_chi2(self._confidence))


	
	def _get_delta_chi2(self,conf_inter ):
		c_dic = {
				'68': 1.00,
				'90': 2.71,
				'95': 4.00}
		const = c_dic[str(conf_inter)]
		return const

	def _integrate(self):
		delta_x = self.x_lim[1] - self.x_lim[0]
		c = delta_x / (self._N - 1)
		# x_values = np.linspace(self.x_limits[0], self.x_limits[1], self.N, dtype=float)
		area = np.sum(c * self.pdf(self.x_values))
		return area

	def _calculate_norm(self):
		area = self._integrate()
		norm = 1/area
		return norm

	def pdf(self, x):
		"""
		Measures the Probability Density Function (CDF) for a given (x) value.
		(assumming the Variable Width Gaussian PDF, from R. Barlow's 2004 paper "Asymmetric Statistical Errors") 

		Returns
		--------
		pdf_x : float
			Value of PDF at x.
		"""
		c = np.sqrt(self.delta_chi2)
		
		par_1 = (2.0 * self.err_p * self.err_n) / ((self.err_p + self.err_n)*c)
		par_2 = (self.err_p  - self.err_n) / ((self.err_p  + self.err_n)*c)
		par_3= (-1.0/2.0) * ((self.nv - x)/(par_1 + par_2*(x - self.nv)))**2.0
		
		#par_4 = self.norm / (2.0 * np.pi)**0.5
		par_4 = self.norm
		value = par_4 * np.exp(par_3)
		#print(v1, v2)
		return value
	
	def delta_chi2_func(self, x):
		"""
		Measures the Δχ2 for a given (x) value.
		(assumming the Variable Width Gaussian PDF, from R. Barlow's 2004 paper "Asymmetric Statistical Errors") 
		
		Returns
		--------
		delta_chi2 : float
			Value of Δχ2 curve at x.
		"""

		c = np.sqrt(self.delta_chi2)
		par_1 = (2.0 * self.err_p * self.err_n) / ((self.err_p + self.err_n)*c)
		par_2 = (self.err_p  - self.err_n) / ((self.err_p  + self.err_n)*c)
		value = ((self.nv - x)/(par_1 + par_2*(x - self.nv)))**2.0
		return value

	def _calculate_cdf_values(self):
		
		cdf_values = np.asarray([])
		for i in range(self._N):
			area = np.trapz(self.pdf_values[:i], x = self.x_values[:i])
			cdf_values = np.append(cdf_values, area)
		return cdf_values

	def cdf(self, x):
		"""
		Measures the Cumulative Density Function (CDF) for a given (x) value.
		
		Returns
		--------
		cdf_x : float
			Value of CDF at x.
		"""
		cdf = interpolate.interp1d(self.x_values, self.cdf_values, kind='nearest')
		return cdf(x)

	def _calculate_inverse_cdf(self, x):
		inverse_cdf = interpolate.interp1d(self.cdf_values, self.x_values, kind='nearest')
		return inverse_cdf(x)


	def gen_random(self):
		"""
		Generates a random number using the probability density function of the object (see .plot_pdf() function).
		
		Returns
		--------
		n : float
			Random number.
		"""
		rnd_prob = np.random.uniform(self.cdf_values[0], self.cdf_values[-1])
		n_rand = self._calculate_inverse_cdf(rnd_prob)
		return n_rand

	def gen_data(self, N=int(1e4)):
		rnd_prob = np.random.uniform(self.cdf_values[0],self.cdf_values[-1], N)
		return (self._calculate_inverse_cdf(rnd_prob))

	def rescale_uncertainties(self, confidence = None, delta_chi2 = None):
		"""
		Rescale the confidence level of the errors. 

		Implemented options are 68, 90 and 95% (e.g confidence = 90)

				Confidence(%)	Δχ2
				68%	1.00
				90%	2.71
				95%	4.00

		or any other value for Δχ2 (e.g. delta_chi2 = 4.7)


		Returns
		--------
		errors : list
			2-size list with the error values at the new confidence interval.
		"""


		if (confidence is None) & (delta_chi2 is None):
			raise RuntimeError("You have to select a 'confidence' value or a Δχ2 (delta_chi2) value!")


		if confidence is not None:
			try:
				self.delta_chi2 = self._get_delta_chi2(confidence)
				self._confidence = confidence
				pass
			except:
				raise RuntimeError('Confidence interval value needs to be chosen among 68, 90 and 95')	
		else:
			self.delta_chi2 = delta_chi2


		flag_up = self.x_values > self.nv
		dif =  np.abs( self.delta_chi2_curve[flag_up] - self.delta_chi2)
		self.err_p = (self.x_values[flag_up])[np.argmin(dif)] - self.nv

		flag_low = self.x_values < self.nv
		dif =  np.abs( self.delta_chi2_curve[flag_low] - self.delta_chi2)
		self.err_n = self.nv - (self.x_values[flag_low])[np.argmin(dif)]


		print(self)

		return [self.err_n, self.err_p]

	@staticmethod
	def fit_func(x, norm, me, err_n, err_p, c):
		
		par_1 = (2.0 * err_p * err_n) / ((err_p + err_n)*c)
		par_2 = (err_p - err_n) / ((err_p + err_n)*c)
		par_3 = (-1.0 / 2.0) * ((me - x) / (par_1 + par_2 * (x - me))) ** 2.0
		#par_4 = norm / (2.0 * np.pi) ** 0.5
		par_4 = norm 
		value = par_4 * np.exp(par_3)
		
		return value

	def _fit(self, expected=None, guess = None):
		plt.ioff()
		y, x, _ = plt.hist(self.data, bins=50, density =True, visible = False)
		plt.ion()

		x = (x[1:] + x[:-1]) / 2  # for len(x)==len(y)
		
		if guess is not None:
			mod = guess
		else:
			mod =  np.mean(x[np.where(y==max(y))])


		norm = (max(y))
	

		
		min_data = min(self.data)
		max_data = max(self.data)
		
		c = np.sqrt(self.delta_chi2)
		if not expected:
			expected = np.array([1, mod, ((mod - min_data) * 0.4), ((max_data - mod) * 0.4), c])

		#expected = [expected_values[0], expected_values[1], expected_values[2], expected_values[3], expected_values[4])
		
		x_scale= [1, 10**int(np.log10(abs(mod))), 10**int(np.log10(abs(((mod - min_data) * 0.4)))), 10**int(np.log10(abs(((max_data - mod) * 0.4)))),  1e-10]
	
		bounds = ([0, -np.inf, 0, 0, c], [np.inf, np.inf, np.inf, np.inf, c+1e-10])
		params, cov = curve_fit(self.fit_func, x, y/norm, expected,  bounds = bounds, x_scale=x_scale, method='trf')
		#print(norm)
		self.norm = params[0]*norm 
		self.nv = params[1]
		
		if params[2] > 0.0:
			self.err_n = (params[2])
			self.err_p = (params[3])
		else:
			self.err_n = (params[3])
			self.err_p = (params[2])

	def plot_pdf(self, show=True, save=False):

		plt.clf()
		if self._confidence is None:
			output = "{:.2e} (-{:.2e}/+{:.2e}) [Δχ2 = {:.2f}]"
			plt.plot(self.x_values, self.pdf_values, color="blue", label = output.format(self.nv, self.err_n, self.err_p, self.delta_chi2))
		else:
			output = "{:.2e} (-{:.2e}/+{:.2e})\n{:.0f}% confidence level (Δχ2={:.2f})"
			plt.plot(self.x_values, self.pdf_values, color="blue", label = output.format(self.nv, self.err_n, self.err_p, self._confidence, self._get_delta_chi2(self._confidence)))


		
		ymin,ymax = plt.ylim()
		plt.ylim(ymin, 1.4*ymax)
		plt.vlines(self.nv + self.err_p, ymin =0,ymax=self.pdf(self.nv + self.err_p), color = 'black', linestyle = '--', linewidth = 0.8)
		plt.vlines(self.nv - self.err_n, ymin =0,ymax=self.pdf(self.nv - self.err_n), color = 'black',  linestyle = '--', linewidth =  0.8)
		plt.vlines(self.nv, ymin =0,ymax=self.pdf(self.nv), color = 'black',  linestyle = '--', linewidth = 1.5)
		plt.axhline(y=0, color="black", linewidth =  0.8)
		plt.legend(fontsize = 12)
		plt.xlabel("x", fontsize=12)
		plt.ylabel("Normalized Probability", fontsize=12)

		if save:
			plt.savefig("plot_pdf.png", dpi=300)

		if show:
			plt.show()

	def plot_delta_chi2_curve(self, show=True, save=False):
		plt.clf()
		plt.plot(self.x_values, self.delta_chi2_curve)
		plt.ylim([-0.5, 6.5])
		

		plt.xlabel("x", fontsize=14)
		plt.ylabel(r'$\Delta\chi^{2}$', fontsize=14)

		plt.axhline(y=self.delta_chi2, color="black", ls="--", label = r'$\Delta\chi^{2} = $'+str(self.delta_chi2))
		plt.axhline(y=0, color="black", label = r'$\Delta\chi^{2} = 0$')
		plt.legend()
		if save:
			plt.savefig("plot_delta_chi2_curve.png", dpi=300)

		if show:
			plt.show()

	def plot_cdf(self, show=True, save=False):
		plt.plot(self.x_values, self.cdf(self.x_values))

		if save:
			plt.savefig("plot_cdf.png", dpi=300)

		if show:
			plt.show()


	def __add__(self, other):
		if isinstance(other, self.__class__):
			add = self.gen_data(N=10000) + other.gen_data(N=10000)
			mod = self.nv + other.nv
			#print(len(add))
		elif isinstance(other, (int, float)):
			add = self.gen_data() + float(other)
			mod = self.nv + float(others)
		else:
			print("Unindentified input type! ({}, {})".format(other, type(other)))
			sys.exit()
		temp_obj = aufloat(data=add, confidence =self._confidence, delta_chi2 = self.delta_chi2, mod = mod)
		return temp_obj

	def __radd__(self, other):
		if isinstance(other, self.__class__):

			add = other.gen_data(N=10000) + self.gen_data(N=10000)
			mod = other.nv + self.nv
		elif isinstance(other, (int, float)):
			add = float(other) + self.gen_data(N=10000)
			mod = float(othr) + self.nv
		else:
			print("Unindentified input type! ({}, {})".format(other, type(other)))
			sys.exit()
		temp_obj = aufloat(data=add, confidence =self._confidence, delta_chi2 = self.delta_chi2, mod = mod)
		return temp_obj

	def __sub__(self, other):
		if isinstance(other, self.__class__):
	
			add = self.gen_data(N=10000) - other.gen_data(N=10000)
			mod - self.nv - other.nv
		elif isinstance(other, (int, float)):
			add = self.gen_data(N=10000) - float(other)
			mod = self.nv - float(other)
		else:
			print("Unindentified input type! ({}, {})".format(other, type(other)))
			sys.exit()
		temp_obj = aufloat(data=add, confidence =self._confidence, delta_chi2 = self.delta_chi2, mod = mod)
		return temp_obj

	def __rsub__(self, other):
		if isinstance(other, self.__class__):
			
			add = other.gen_data(N=10000) - self.gen_data(N=10000)
			mod = other.gen_data(N=10000) - self.nv
		elif isinstance(other, (int, float)):
			add = float(other) - self.gen_data(N=10000)
			mod = float(other) - self.nv
		else:
			print("Unindentified input type! ({}, {})".format(other, type(other)))
			sys.exit()
		temp_obj = aufloat(data=add, confidence =self._confidence, delta_chi2 = self.delta_chi2, mod = mod)
		return temp_obj

	def __mul__(self, other):
		if isinstance(other, self.__class__):
			
			add = self.gen_data(N=10000) * other.gen_data(N=10000)
			mod = self.nv * other.nv
		elif isinstance(other, (int, float)):
			add = self.gen_data(N=10000) * float(other)
			mod = self.nv * float(other)
		else:
			print("Unindentified input type! ({}, {})".format(other, type(other)))
			sys.exit()
		temp_obj = aufloat(data=add, confidence =self._confidence, delta_chi2 = self.delta_chi2, mod = mod)
		return temp_obj

	def __rmul__(self, other):
		if isinstance(other, self.__class__):
			
			add = other.gen_data(N=10000) * self.gen_data(N=10000)
			mod = other.nv * self.nv
		elif isinstance(other, (int, float)):
			add = float(other) * self.gen_data(N=10000)
			mod = float(other) * self.nv
		else:
			print("Unindentified input type! ({}, {})".format(other, type(other)))
			sys.exit()
		temp_obj = aufloat(data=add, confidence =self._confidence, delta_chi2 = self.delta_chi2, mod = mod)
		return temp_obj

	def __truediv__(self, other):
		if isinstance(other, self.__class__):
			
			add = self.gen_data(N=10000) / other.gen_data(N=10000)
			mod = self.nv / other.nv
		elif isinstance(other, (int, float)):
			add = self.gen_data(N=10000) / float(other)
			mod = self.nv / float(other)
		else:
			print("Unindentified input type! ({}, {})".format(other, type(other)))
			sys.exit()

		
		temp_obj = aufloat(data=add, confidence =self._confidence, delta_chi2 = self.delta_chi2, mod = mod)
		return temp_obj

	def __rtruediv__(self, other):
		if isinstance(other, self.__class__):
			
			add = other.gen_data(N=10000) / self.gen_data(N=10000)
			mod = other.nv / self.nv
		elif isinstance(other, (int, float)):
			add = float(other) / self.gen_data(N=10000)
			mod = float(other) / self.nv
		else:
			print("Unindentified input type! ({}, {})".format(other, type(other)))
			sys.exit()
		
		temp_obj = aufloat(data=add, confidence =self._confidence, delta_chi2 = self.delta_chi2, mod = mod)
		return temp_obj

	def __pow__(self, other):
		if isinstance(other, self.__class__):
			
			add = self.gen_data(N=10000) ** other.gen_data(N=10000)
			mod = self.nv ** other.nv
		elif isinstance(other, (int, float)):
			add = self.gen_data() ** float(other)
			mod = self.nv ** float(other)
		else:
			print("Unindentified input type! ({}, {})".format(other, type(other)))
			sys.exit()

		temp_obj = aufloat(data=add, confidence =self._confidence, delta_chi2 = self.delta_chi2, mod = mod)
		return temp_obj

	def __rpow__(self, other):
		if isinstance(other, self.__class__):
			
			add = other.gen_data(N=10000) ** self.gen_data(N=10000)
			mod = other.nv ** self.nv
		elif isinstance(other, (int, float)):
			add = float(other) ** self.gen_data(N=10000)
			mod = float(other) ** self.nv
		else:
			print("Unindentified input type! ({}, {})".format(other, type(other)))
			sys.exit()

		
		temp_obj = aufloat(data=add, confidence =self._confidence, delta_chi2 = self.delta_chi2, mod = mod)
		return temp_obj	



if __name__ == "__main__":
	

	examples:

	a = aufloat(10,.5,1, confidence = 68)

	b = aufloat(10,.5,1, delta_chi2 = 2.7)

	c = a + b

	d = aufloat(data = a.gen_data, confidence = 90)


