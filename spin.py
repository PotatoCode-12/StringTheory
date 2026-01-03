import math
from scipy.integrate import quad
import numpy as np 
from sympy import symbols, diff,integrate, Function, Derivative, summation, oo
import variables
import constants

n = symbols('n') # integer for summation
variables.S(variables.i, variables.j) == -variables.i*summation(1/n*(variables.TensionStringab(variables.i, -n)*variables.TensionStringab(variables.j, n) - variables.TensionStringab(variables.j, -n)*variables.TensionString(variables.i, n)), (n, 1, math.inf)) # transverse direction spin check conditions pg 24
variables.S(23) <= 1 +variables.SlopeRegge*variables.m**2
