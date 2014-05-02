import numpy as np
import hazel
import matplotlib.pyplot as p

class bayesHazel(object):
	
	def __init__(self, wavelength, stokes):
		
		self
		
		self.synModeInput = 5
		self.nSlabsInput = 1
		self.hInput = 3.e0
		self.boundaryInput  = [0.0,0.0,0.0,0.0]
		self.transInput = 1
		self.atomicPolInput = 1
		self.anglesInput = [90.0,0.0,90.0]
		
		lambdaAxisInput = [-1.5e0,2.5e0]
		nLambdaInput = 150
		self.nbar = [0.0,0.0,0.0,0.0]
		self.omega = [0.0,0.0,0.0,0.0]
		
		
	def synth(self, x):
				
		self.B1Input = [3.0,80.0,41.0]
		self.tau1Input = 1.e0
		self.dopplerWidthInput = 6.e0
		self.dampingInput = 0.e0
		self.dopplerVelocityInput = 0.e0
		
