import math
from scipy.integrate import quad
import numpy as np 
from sympy import symbols, diff,integrate, Function, Derivative, summation, oo, ln
import variables
import constants

variables.xa(variables.mu) == variables.weight(0, 0)
diff(variables.xa(variables.mu), variables.z) == variables.weight(1, 0) 
diff(variables.xa(variables.mu), variables.zl) == variables.weight(0, 1)
diff(diff(variables.xa(variables.mu), variables.z), variables.z) == variables.weight(2,0)
variables.point(variables.ea(variables.i, variables.k*variables.X)) == variables.weight((variables.SlopeRegge*variables.ka(2))/4, (variables.SlopeRegge*variables.ka(2))/4)
variables.point([diff(variables.xa(variables.mub(variables.i)), variables.mb(variables.i)) for variables.i in variables.xa(variables.mub(variables.i))], [diff(variables.xa(variables.vb(variables.j)), variables.nb(variables.j)) for variables.i in variables.nb(variables.j)] * variables.ea(variables.i,variables.k*variables.X)) == variables.weight((variables.SlopeRegge*variables.ka(2))/4 + summation(variables.mb(variables.i), (variables.i, 1, math.inf)), (variables.SlopeRegge*variables.ka(2))/4 + summation(variables.nb(variables.j), (variables.j, 1, math.inf))) # conditions page 46-7

variables.chi == variables.chitilde + 1/4 * variables.nb(variables.c) # if
math.exp(-variables.lamda * variables.chitilde) == math.exp(-variables.lamda*variables.chi + variables.lamda*variables.nb(variables.c)/4)
