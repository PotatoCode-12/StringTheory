import math
from scipy.integrate import quad
import numpy as np 
from sympy import symbols, diff,integrate, Function, Derivative, summation, oo
import sympy as sp
import variables
import constants

quad(variables.d2z(variables.d2va(variables.v, variables.z) + variables.d2va(variables.v, variables.zl)), variables.R) == variables.i*quad(variables.va(variables.z)*diff(variables.zl) - variables.va(variables.zl)*diff(variables.z), (variables.tau, 0, 2*math.pi)) # with complex cordinates pg 34

#theorems pg 42
quad(diff(variables.A)*variables.nb(variables.a)*variables.ja(variables.a)*variables.AO(variables.sigmab(variables.d)), (diff(variables.R))) == 2*math.pi/(variables.i*variables.epsilon)*variables.delta*variables.AO(variables.sigmab(variables.D0))
#when in two flat dimensions
quad((variables.j*variables.d*variables.z - variables.jtilde*variables.d*variables.zl)*variables.AO(variables.zb(variables.D0), variables.zlb(variables.D0)), (variables.tau, 0, 2*math.pi)) == 2*math.pi/variables.epsilon*variables.delta*variables.AO(variables.zb(variables.D0), variables.zlb(variables.D0))

#theta-based divergence
summation((variables.n - variables.theta), (variables.n, 1, math.inf)) == 1/24 -1/8*(2*variables.theta -1)**2
