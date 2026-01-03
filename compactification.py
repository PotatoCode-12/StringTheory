import math
from scipy.integrate import quad
import numpy as np 
from sympy import symbols, diff,integrate, Function, Derivative, summation, oo, ln, Abs
import sympy as sp
import variables
import constants

#narin compactification, check ommited vertex OPE operator approximation 251.16
np.dot(variables.lb(variables.L), variables.lbp(variables.L) ) - np.dot(variables.lb(variables.R), variables.lbp(variables.R)) == variables.l(variables.lp(variables.ZB))
2 * (variables.l(variables.lp)) == ((variables.l + variables.lp)(variables.l + variables.lp)) - variables.l(variables.l) - variables.lp(variables.lp) # to 2 * variables.ZB
variables.Zbf(variables.Rho, variables.tau) == Abs(variables.etaf(variables.tau))**(-2 * variables.k) * summation(math.exp(math.pi * variables.i * variables.tau * variables.lab(2, variables.L) - math.pi * variables.i * variables.taul * variables.lab(2, variables.R)), (variables.l, variables.Rho))
summation(variables.deltaf(variables.l - variables.lp), (variables.lp, variables.Rho)) == variables.Vab(-1, variables.Rho) * summation(math.exp(2 * math.pi * variables.i * variables.lpp(variables.l)), (variables.lpp, variables.star(variables.Rho)))
variables.Zbf(variables.Rho, variables.tau) == variables.Vab(-1, variables.Rho) * Abs(variables.etaf(variables.tau))**(-2 * variables.k) * summation(integrate(variables.da(2 * variables.k) * variables.l * math.exp(2 * math.pi * variables.i * variables.lpp(variables.l) + math.pi * variables.i * variables.tau * variables.lab(2, variables.L) - math.pi * variables.i * variables.taul * variables.lab(2, variables.R))), (variables.lpp, variables.star(variables.Rho))) == variables.Vab(-1, variables.Rho)(variables.tau * variables.taul)**(-variables.k/2) * summation(math.exp(-math.pi * variables.i * (variables.lppab(2, variables.L))/variables.tau + math.pi * variables.i * (variables.lppab(2, variables.R))/variables.taul)) == variables.Vab(-1, variables.Rho) * variables.Zbf(variables.star(variables.Rho, -1/variables.tau))
variables.Rho == variables.star(variables.Rho)
#check condition 252.25
variables.Rhop == variables.Lamda * variables.Rho
variables.lb((variables.L, variables.R)) == (variables.n)/(variables.r) + variables.posneg * (variables.m * variables.r)/2
variables.lbp(variables.L) == variables.lb(variables.L) * math.cos(variables.h * variables.lamda) + variables.lb(variables.R) * math.sin(variables.h * variables.lamda)
variables.lbp(variables.R) == variables.lb(variables.L) * math.sin(variables.h * variables.lamda) + variables.lb(variables.R) * math.cos(variables.h * variables.lamda)
(2 * variables.k * (2 * variables.k -1))/2 - variables.k *(variables.k -1) == variables.ka(2)
variables.Lamda * variables.Rhob(0) == variables.Lamdap * variables.Lamda * variables.Lamdapp * variables.Rhob(0)
variables.Lamdap == np.cross(variables.Of(variables.k, variables.RB), variables.Of(variables.k, variables.RB))
variables.Lamdapp == variables.Of(variables.k, variables.k, variables.ZB)
# check space condition 253.31
#periodicty
variables.xaps(variables.m) == variables.Lab(variables.m, variables.n) * variables.xas(variables.n)
variables.bb(variables.m * variables.n) == variables.bb(variables.m * variables.n) + variables.Nb(variables.m * variables.n)
# skipped example, 255.36+

