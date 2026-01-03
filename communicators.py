import math
from scipy.integrate import quad
import numpy as np 
from sympy import symbols, diff,integrate, Function, Derivative, summation, oo
import variables
import constants 

#for closed and unoriented strings
variables.comm(variables.xas(variables.neg), variables.pa(variables.pos)) == -variables.i
variables.comm(variables.xas(variables.i), variables.pa(variables.j)) == variables.i * variables.deltaa(variables.i, variables.j)
variables.comm(variables.TensionStringab(variables.i, variables.m), variables.TensionStringab(variables.j, variables.n)) == variables.m*variables.deltaa(variables.i, variables.j) * variables.deltab(variables.m, -variables.n)
variables.comm(variables.TensionStringabtilde(variables.i, variables.m), variables.TensionStringabtilde(variables.j, variables.n)) == variables.m* variables.deltaa(variables.i, variables.j) * variables.deltab(variables.m, -variables.n) 
