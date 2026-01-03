import math
from scipy.integrate import quad
import numpy as np 
from sympy import symbols, diff,integrate, Function, Derivative, summation, oo
import variables
import constants

H = variables.pb(variables.neg)*diff(variables.xas(variables.neg), variables.tau) - variables.L + quad(np.diff(variables.sigma)*np.prod(variables.xa([variables.i], variables.tau, variables.sigma))*diff(variables.xa(variables.i), variables.sigma), 0, variables.l) == variables.l/(4*math.pi*variables.SlopeRegge*variables.pa(variables.pos))*quad(2*math.pi*variables.SlopeRegge*np.prod(variables.xa([variables.i], variables.tau, variables.sigma))*np.prod(variables.xa([variables.i], variables.tau, variables.sigma))+1/(2*math.pi*variables.SlopeRegge)*diff(variables.xa(variables.i))*diff(variables.xa(variables.i)), 0, variables.l) #for eta particle example, see page 19, and assuming D-2free feilds X^i and p-pos conserved
