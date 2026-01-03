import math
from scipy.integrate import quad
import numpy as np 
from sympy import symbols, diff,integrate, Function, Derivative, summation, oo
import variables
import constants

#for vacuum expansion of hilbert space reference page 21
variables.D == 26
variables.A == -1
n = symbols('n', integer=True) #for summation limits
variables.hilbert == variables.hilbert1 + variables.hilbert2
variables.H == variables.pa(variables.i)*variables.pa(variables.i)/(2*variables.pa(variables.pos)) + 1/(2*variables.pa(variables.pos)+variables.SlopeRegge)*(summation(variables.TensionStringab(variables.i, -n)*variables.TensionStringab(variables.i, n)+variables.A, (n, 1, math.inf)))
variables.A == (variables.D-2)/2*summation(n, (n, 1, math.inf)) #for hamiltonian
summation(n,(n, 1, math.inf)) == (2*variables.l*variables.pa(variables.pos)*variables.SlopeRegge)/(variables.epsilon**2*math.pi)-1/12+variables.O(variables.epsilon)
