import math
from scipy.integrate import quad
import numpy as np 
from sympy import symbols, diff,integrate, Function, Derivative, summation, oo
import variables
import constants

np.prod(variables.xa([variables.i], variables.tau, variables.sigma)) == diff(variables.L)/diff(diff(variables.xa(variables.i), variables.tau)) == 1/(2*math.pi*variables.SlopeRegge)*variables.GWorldSheetMetric(variables.sigma, variables.sigma) * diff(variables.xa(variables.i), variables.i) == variables.pa(variables.pos)/variables.l*variables.xa(variables.i, variables.tau) #for eta in particle example, see page 19 & associated hamiltonian
