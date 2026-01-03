import math
from scipy.integrate import quad
import numpy as np 
from sympy import symbols, diff,integrate, Function, Derivative, summation, oo
import variables
import constants

variables.ActionWorldLine == 1/2(integrate(np.diff(variables.tau)(variables.eta**(-1))*variables.xadiff(variables.mu)*variables.xbdiff(variables.mu)-variables.eta*variables.m**2))
