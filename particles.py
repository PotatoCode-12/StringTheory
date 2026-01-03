import math
from scipy.integrate import quad
import numpy as np 
from sympy import symbols, diff,integrate, Function, Derivative, summation, oo
import variables
import constants

n = symbols('n', integer=True)  # for summation limits
variables.m**2 == 2*variables.pa(variables.pos)* variables.H - variables.pa(variables.i) * variables.pa(variables.i) == 1/(variables.SlopeRegge) * ((variables.N + 2-variables.d)/24) # point particle pg 22
variables.N == summation(summation(n*variables.N(variables.i, n) (n, 1, math.inf)), (variables.i, 2, variables.D-1))
