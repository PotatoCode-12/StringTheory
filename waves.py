import math
from scipy.integrate import quad
import numpy as np 
from sympy import symbols, diff,integrate, Function, Derivative, summation, oo
import variables
import constants

n= symbols('n') #for summation
diff(diff(variables.xa(variables.i), variables.tau), variables.tau) == constants.c**2*diff(diff(variables.xa(variables.i), variables.sigma), variables.sigma) # reference page 19 for conditions
variables.xa(variables.i, variables.tau, variables.sigma) == variables.xas(variables.i) + variables.pa(variables.i)/(variables.pa(variables.pos))*variables.tau +variables.i(variables.SlopeRegge)**(1/2)*sum((1/n*variables.TensionString(variables.i, n)*math.e**(-math.pi*variables.tau.imag/variables.l) * math.cos((math.pi*n*variables.sigma)/variables.l)) for n in range(-math.inf, math.inf) if n != 0) #solution to earlier wave equation with bounds, check page 20
variables.TensionString(variables.i, -n) == np.conjugate(variables.TensionString(variables.i, n).T) #conjugate of the tension string, needed for above solution
