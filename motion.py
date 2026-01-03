import math
from scipy.integrate import quad
import numpy as np 
from sympy import symbols, diff,integrate, Function, Derivative, summation, oo, ln
import variables
import actions
import constants 

n = symbols('n') #summation variable
diff_resultb = variables.xb(variables.mu)
for _ in range(int(variables.l)):
    diff_resultb = diff(diff_resultb, variables.mu)
b = diff_resultb
def nambuGotoAction():
    variables.EquMotion * actions.AppNambuGoto == -1/(4*math.pi*variables.SlopeRegge) * integrate(np.diff(variables.tau)*np.diff(variables.sigma)*(-variables.GWorldSheetMetric)**(1/2)*variables.EquMotion*(variables.h(variables.a,variables.b)-1/2*variables.GWorldSheetMetric**(variables.a*variables.b)*(variables.h(variables.a, variables.b)-((1/2)*variables.WorldSheetMetric(variables.a*variables.b)*variables.GWorldSheetMetric**(variables.c*variables.d)*variables.h(variables.c, variables.d)), variables.WorldSheet, math.inf)), variables.M)

variables.eta**2 == -variables.xap(variables.mu)*variables.xbp(variables.mu)/(variables.m**2)
variables.NambuGotoVar*actions.AppNambuGoto == 1/(2*math.pi*variables.SlopeRegge)*quad(np.diff(variables.tau) * quad(np.diff(variables.sigma)*(-variables.GWorldSheetMetric)**(-1/2)*variables.NambuGotoVar*variables.xa(variables.mu)*variables.FreeStringPropagation(variables.tau)*variables.xb(variables.mu) -1/(2*math.pi*variables.SlopeRegge)*quad(diff(variables.tau)*(-variables.GWorldSheetMetric)**(1/2)*variables.NambuGotoVar*variables.xa(variables.mu)*diff_resultb,-math.inf, math.inf), 0, variables.l), -math.inf, math.inf)
# follow page 19 and associated hamiltonian
diff(variables.xas(variables.neg), variables.tau) == diff(variables.H)/diff(variables.pb(variables.neg)) == variables.H/(variables.pa(variables.pos))
diff(variables.xa(variables.i), variables.tau) == variables.NambuGotoVar*variables.H/(variables.NambuGotoVar*np.prod(variables.xa([variables.i], variables.tau, variables.sigma))) == 2*math.pi*variables.SlopeRegge*constants.c*np.prod(variables.xa([variables.i], variables.tau, variables.sigma))
diff(variables.pa(variables.pos), variables.tau) == diff(variables.H)/diff(variables.xas(variables.neg)) == 0
diff(np.prod(variables.xa([variables.i], variables.tau, variables.sigma)), variables.tau) == -diff(variables.H)/diff(variables.xa(variables.i)) == variables.c/(2*math.pi*variables.SlopeRegge)*diff(diff(variables.xa(variables.i), variables.sigma), variables.sigma)

#centermass variable limits follow page 20
variables.xas(variables.i, variables.tau) == 1/(variables.l)*quad(np.diff(variables.sigma)*variables.xa(variables.i, variables.tau, variables.sigma), 0, variables.l)
variables.pa(variables.i, variables.tau) == quad(np.diff(variables.sigma)* np.prod(variables.xa([variables.i], variables.tau, variables.sigma)), 0, variables.l) == variables.pa(variables.pos)/variables.l*quad(np.diff(variables.sigma)*diff(variables.xa(variables.i, variables.tau, variables.sigma), variables.tau), 0, variables.l) 

variables.xa(variables.i, variables.tau, variables.sigma) == variables.xas(variables.i) + variables.pa(variables.i)/variables.pa(variables.pos)*variables.tau + variables.i*(variables.SlopeRegge/2)**(1/2) * (summation((variables.TensionStringab(variables.i, n))/n*math.exp(-(2*math.pi*variables.i*n(variables.sigma+constants.c*variables.tau))/variables.l) + (variables.TensionStringabtilde(variables.i, n))/n* math.exp((2*math.pi*variables.i*n(variables.sigma-constants.c*variables.tau))/variables.l), (n, -math.inf, -1)) + (summation((variables.TensionStringab(variables.i, n))/n*math.exp(-(2*math.pi*variables.i*n(variables.sigma+constants.c*variables.tau))/variables.l) + (variables.TensionStringabtilde(variables.i, n))/n* math.exp((2*math.pi*variables.i*n(variables.sigma-constants.c*variables.tau))/variables.l), (n, 1, math.inf)))) #pg 26

diff(diff(variables.xa(variables.mu, variables.z, variables.zl), variables.zl), variables.z) == 0 # classical equation of motion pg 34
diff(diff(variables.xa(variables.mu), variables.zl), variables.z) == diff(diff(variables.xa(variables.mu), variables.z), variables.zl) == 0

variables.point(variables.PI) == math.exp(variables.SlopeRegge/4*integrate(variables.d2z1*variables.d2z2*ln(variables.z12)**(variables.D2)*variables.delta/(variables.delta*variables.xa(variables.mu, variables.z1, variables.zl1))*(variables.delta/(variables.delta*variables.xb(variables.mu, variables.z2, variables.zl2)))))*variables.PI

diff(variables.cf(variables.z), variables.zl) == diff(variables.bf(variables.z), variables.zl) == 0
diff(variables.bf(variables.z), variables.zl) * variables.cf(variables.D0) == 2*math.pi*variables.deltaa(variables.D2, variables.z, variables.zl) # operator equations pg 50
variables.point(variables.bf(variables.z1)*variables.cf(variables.z2)) == variables.bf(variables.z1) * variables.cf(variables.z2) - 1/(variables.z12)
diff(1/variables.z, variables.zl) == diff(1/variables.zl, variables.z) == 2*math.pi*variables.deltaa(variables.D2, variables.z, variables.zl)
variables.bf(variables.z1)*variables.bf(variables.z2) == variables.O(variables.z12) # products
variables.cf(variables.z1) * variables.bf(variables.z2) == variables.O(variables.z12)

#momentum for spacetime current 
variables.pa(variables.mu) == 1/(2*math.pi*variables.i)*integrate(diff(variables.z) * variables.ja(variables.mu) - diff(variables.zl) * variables.jatilde(variables.mu), 0, 2*math.pi) == (2/variables.SlopeRegge)**(1/2) * variables.alphaab(variables.mu, 0) * (2/variables.SlopeRegge)**(1/2) * variables.alphaabtilde(variables.mu, 0) # contour
variables.xa(variables.mu, variables.z, variables.zl) == variables.xas(variables.mu) - variables.i*variables.SlopeRegge/2*variables.pa(variables.mu) * ln(variables.z)**(variables.D2) + variables.i*(variables.SlopeRegge/2)**(1/2) * summation(1/variables.m*((variables.alphaab(variables.mu, variables.m))/variables.za(variables.m) + (variables.alphaabtilde(variables.mu, variables.m))/variables.zla(variables.m)), (variables.m, 1, math.inf)) + summation(1/variables.m*((variables.alphaab(variables.mu, variables.m))/variables.za(variables.m) + (variables.alphaabtilde(variables.mu, variables.m))/variables.zla(variables.m)), (variables.m, -math.inf, -1)) # post-integration
variables.brackets(variables.alphaab(variables.mu, variables.m), variables.alphaab(variables.v, variables.n)) == variables.brackets(variables.alphabtilde(variables.mu, variables.m), variables.alphaabtilde(variables.v, variables.n)) == variables.m*variables.deltab(variables.m, -variables.n) * variables.etaa(variables.mu, variables.v)
variables.brackets(variables.xas(variables.mu), variables.pa(variables.v)) == variables.i*variables.etaa(variables.mu, variables.v)
variables.Lb(variables.D0) == (variables.SlopeRegge*variables.pa(variables.D2))/4 + summation(variables.alphaab(variables.mu, -variables.n) * variables.alphab(variables.mu, variables.n), (variables.n, 0, math.inf)) + variables.aa(variables.X) # normal ordering constant
2*variables.Lb(variables.D0)*variables.dirac(0, 0) == (variables.Lb(variables.D1)*variables.Lb(-variables.D1) - variables.Lb(-variables.D1) * variables.Lb(variables.D1)) * variables.dirac(0, 0) == 0 # virasaro
variables.aa(variables.X) == 0 #then
variables.Lb(variables.m) == 1/2(summation(variables.ANO(variables.alphaab(variables.mu, variables.m-variables.n) * variables.alphab(variables.mu, variables.n))), (variables.n, -math.inf, math.inf))
variables.xa(variables.mu, variables.z, variables.zl) * variables.xa(variables.v, variables.zp, variables.zpl) == variables.ANO(variables.xa(variables.mu, variables.z, variables.zl) * variables.xa(variables.v, variables.zp, variables.zpl)) + variables.SlopeRegge/2 * variables.etaa(variables.mu, variables.v) * variables.brackets(-ln(variables.z)** variables.D2 + summation(1/variables.m * ((variables.zpa(variables.m))/(variables.za(variables.m))) + (variables.zpla(variables.m))/variables.zla(variables.m), (variables.m, 1, math.inf))) == variables.ANO(variables.xa(variables.mu, variables.z, variables.zl) * variables.xa(variables.v, variables.zp, variables.zpl)) - variables.SlopeRegge/2*variables.etaa(variables.mu, variables.v) * ln(variables.z - variables.zp)**variables.D2 # w/ communicators 60.11
variables.ANO(variables.xa(variables.mu, variables.z, variables.zl) * variables.xa(variables.v, variables.zp, variables.zpl)) == variables.point(variables.xa(variables.mu, variables.z, variables.zl) * variables.xa(variables.v, variables.zp, variables.zpl))
variables.brackets1(variables.xa(variables.mu, variables.z, variables.zl) * variables.xa(variables.v, variables.zp, variables.zpl)) == variables.brackets2(variables.xa(variables.mu, variables.z, variables.zl) * variables.xa(variables.v, variables.zp, variables.zpl)) + variables.etaa(variables.mu, variables.v) * variables.deltaf(variables.z, variables.zl, variables.zp, variables.zpl)

# vertex operator momenta 250.14
variables.kb(variables.r, variables.L) == variables.eba(variables.r, variables.m) * (variables.vb(variables.m * variables.L))/variables.SlopeRegge 
variables.kb(variables.r * variables.R) == variables.eba(variables.r, variables.m) * (variables.vb(variables.m * variables.R))/variables.SlopeRegge 
variables.ma(2) == 1/2 * (variables.kb(variables.r * variables.L) * variables.kb(variables.r * variables.L) + variables.kb(variables.r * variables.R) * variables.kb(variables.r * variables.R)) + 2/variables.SlopeRegge * (variables.N + variables.Ntilde - 2)
0 == variables.SlopeRegge * (variables.kb(variables.r * variables.L) * variables.kb(variables.r * variables.L) - variables.kb(variables.r * variables.R) * variables.kb(variables.r *variables.R)) + 4 * (variables.N - variables.Ntilde)
