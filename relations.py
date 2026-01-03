import math
from scipy.integrate import quad
import numpy as np 
from sympy import symbols, diff,integrate, Function, Derivative, summation, oo
import variables
import constants

x = symbols('x') #placeholder for interval notation
variables.EquMotion == variables.WorldSheetMetric(variables.a*variables.b) * (variables.WorldSheetMetric(variables.a*variables.b))**(variables.a*variables.b)*variables.NambuGotoVar*variables.GWorldSheetMetric == -variables.NambuGotoVar*variables.GWorldSheetMetric*variables.WorldSheetMetric(variables.a*variables.b)*(variables.WorldSheetMetric(variables.a*variables.b))**(variables.a*variables.b)


#for heisenberg operators in wave function page 20
if variables.xas(variables.neg) <= x <= variables.pa(variables.pos):
    variables.i*variables.etaa(variables.neg, variables.pos) == -variables.i
if variables.xa(variables.i, variables.sigma) <= x <= variables.ConjMomentum(variables.j, variables.sigmap):
    (variables.i*variables.NambuGotoVar(variables.i, variables.j) * variables.NambuGotoVar(variables.sigma - variables.sigmap))
#fourier compoents
if variables.xas(variables.i) <= x <= variables.pa(variables.j):
    variables.i*variables.NambuGotoVar(variables.i, variables.j)
if variables.TensionStringab(variables.i, variables.m) <= x <= variables.TensionStringab(variables.j, variables.n):
    variables.m*variables.NambuGotoVar(variables.i, variables.j)*variables.EquMotion(variables.m, -variables.n)
#reference normalization bounds on page 20
variables.m >> 0
variables.pa(variables.pos)*variables.dirac(0, variables.k) == variables.k(variables.pos)*variables.dirac(0, variables.k)
variables.pa(variables.i)*variables.dirac(0, variables.k) == variables.k(variables.i)*variables.dirac(0, variables.k)
variables.TensionStringab(variables.i, variables.m)*variables.dirac(0, variables.k) == 0
variables.dirac(variables.N, variables.k) == variables.dirac(0, variables.k) #check parameters page 21
