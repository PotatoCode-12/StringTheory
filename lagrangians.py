import math
from scipy.integrate import quad
import numpy as np 
from sympy import symbols, diff,integrate, Function, Derivative, summation, oo
import variables
import constants

#125.22
def zeroPointConstant():
    (variables.Lab(variables.m, 0) + variables.Lab(variables.g, 0)) * variables.dirac(variables.psi, variables.down) == 0
    (variables.Ha(variables.m) + variables.Ha(variables.g)) * variables.dirac(variables.psi, variables.down) == 0
