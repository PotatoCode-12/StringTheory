import math
from scipy.integrate import quad
import numpy as np 
from sympy import symbols, diff,integrate, Function, Derivative, summation, oo
import variables
import constants

variables.E == -(variables.c + variables.ctilde)/24 # radial 73.17
variables.E == -(math.pi*(variables.c + variables.ctilde))/(12* variables.l) # general
