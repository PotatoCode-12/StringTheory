import math
from scipy.integrate import quad
import numpy as np 
from sympy import symbols, diff,integrate, Function, Derivative, summation, oo
import variables
import actions

(-variables.NWorldSheetMetric)**(1/2)*variables.ScalarRicciP == (-variables.GWorldSheetMetric)**(1/2)*(variables.ScalarRicci - 2* variables.rescaled_laplacian(variables.omega)) #rescaled laplacian for ricci scalar definition
