import math
from scipy.integrate import quad
import numpy as np 
from sympy import symbols, diff,integrate, Function, Derivative, summation, oo
import variables
import tensions

variables.h(variables.a, variables.b) == 1/2(variables.WorldSheetMetric(variables.a*variables.b)*variables.GWorldSheetMetric**(variables.c*variables.d)*variables.h(variables.c, variables.d)) #for motion equation
variables.WorldSheetMetric(variables.a*variables.b)*(variables.NambuGotoVar/(variables.NambuGotoVar*variables.WorldSheetMetric(variables.a*variables.b))) == 0
(tensions.T(variables.a))**variables.a == 0
