import numpy as np
import hazel
import matplotlib.pyplot as pl
import pdb
import pymultinest

class bayesHazel(object):
	
	def __init__(self, fileObs, noise):
		
		stokes = np.loadtxt(fileObs, skiprows=1)
		
		left = 100
		right = 195
		self.lambdaAxisInput = np.asarray([stokes[left,0],stokes[right,0]])
		self.nLambdaInput = stokes[left:right,0].shape[0]
		self.wavelength = stokes[left:right,0]
		self.stokes = stokes[left:right,1:].T
		self.sigmaNoise = noise
		
		self.synModeInput = 0
		self.nSlabsInput = 1
		self.hInput = 3.e0
		self.boundaryInput  = np.asarray([0.0,0.0,0.0,0.0])
		self.transInput = 1
		self.atomicPolInput = 1
		self.anglesInput = np.asarray([90.0,0.0,90.0])
		
		self.nbarInput = np.asarray([0.0,0.0,0.0,0.0])
		self.omegaInput = np.asarray([0.0,0.0,0.0,0.0])
		
		hazel.init()

		self.nParams = 7
		self.lower = [4.0,-3.0,0.0,0.0,0.0,1e-5,1e-5]
		self.upper = [10.0,3.0,200.0,180.0,180.0,1e-1,1e-2]
		
	def synthesize(self, x):
		
		#x = [doppler,v,b,theta,chi,sigmaI,sigmaPol]		
		self.B1Input = np.asarray(x[2:5])
		self.B2Input = np.asarray([0.0,0.0,0.0])
		self.tau1Input = 0.e0
		self.tau2Input = 0.e0
		self.dopplerWidthInput = x[0]
		self.dopplerWidth2Input = 0.e0
		self.dampingInput = 0.0
		self.dopplerVelocityInput = x[1]
		self.dopplerVelocity2Input = 0.e0
		self.ffInput = 0.e0
		
		
		[self.l, stokes, etaOutput, epsOutput] = hazel.synth(self.synModeInput, self.nSlabsInput, self.B1Input, self.B2Input, self.hInput, 
			self.tau1Input, self.tau2Input, self.boundaryInput, self.transInput, self.atomicPolInput, self.anglesInput, 
			self.lambdaAxisInput, self.nLambdaInput, self.dopplerWidthInput, self.dopplerWidth2Input, self.dampingInput, 
			self.dopplerVelocityInput, self.dopplerVelocity2Input, self.ffInput, self.nbarInput, self.omegaInput)
		
		return stokes
		
	def prior(self, cube, ndim, nparams):		
		for i in range(ndim):
			cube[i] = cube[i] * (self.upper[i]-self.lower[i]) + self.lower[i]
			
	def logLike(self, cube, ndim, nparams):
		for i in range(self.nParams):
			self.pars[i] = cube[i]		
		
		stokesSynth = self.synthesize(self.pars[0:6])
		
# Log-likelihood
		logL = 0.0
		logLI = -self.nLambdaInput * np.log(self.pars[-2]) - 0.5 * np.sum((self.stokes[0,:]-stokesSynth[0,:])**2 / self.pars[-2]**2)
		logLPol = -3.0 * self.nLambdaInput * np.log(self.pars[-1]) - 0.5 * np.sum((self.stokes[1:,:]-stokesSynth[1:,:])**2 / self.pars[-1]**2)
			
# Log-prior
		logLPrior = -np.log(self.pars[-2]) - np.log(self.pars[-1])
		
		logL = logLI + logLPol + logLPrior
		
		print self.pars, logL
		return logL
		
	def sample(self):
				
		self.pars = np.zeros(self.nParams)

# run MultiNest
		pymultinest.run(self.logLike, self.prior, self.nParams, importance_nested_sampling = False, resume = False, verbose = True, 
			sampling_efficiency = 'parameter', n_live_points = 100)
		
	def analyze(self):
		
		which = [0,2,3,4,5,6]		
		pars = [r'$v_\mathrm{Doppler}$',r'$B$',r'$\theta_B$',r'$\chi_B$',r'$\sigma_I$ [x10$^{-2}$]', r'$\sigma_{QUV}$ [x10$^{-2}$]']
		self.chains = pymultinest.Analyzer(n_params = self.nParams)
		samples = self.chains.get_equal_weighted_posterior()
		samples[:,-3:] *= 1e2
		fig, ax = pl.subplots(nrows=2,ncols=3, figsize=(16,8))
		ax = ax.flatten()

		for i in range(self.nParams-1):
			n, bins, patches = ax[i].hist(samples[:,which[i]], bins=20)
			pl.setp(patches, 'facecolor', 'g', 'alpha', 0.75)		
			ax[i].set_xlabel(pars[i])
		pl.savefig("bayesHazel.png")
		
h = bayesHazel('prominence_ApJ_642_554.prof', 0.0003)
h.sample()
h.analyze()