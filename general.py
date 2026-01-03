import math
from scipy.integrate import quad
import numpy as np 
from sympy import symbols, diff,integrate, Function, Derivative, summation, oo
import variables
import actions

variables.chi == 1/(4*math.pi)*quad(np.diff(variables.tau)*np.diff(variables.sigma)*(-variables.GWorldSheetMetric)**(-1/2)*variables.ScalarRicci, variables.WorldSheet, math.inf)#additional term beyond original Sp scale (p15)
