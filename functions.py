import math
from scipy.integrate import quad
import numpy as np 
from sympy import symbols, diff,integrate, Function, Derivative, summation, oo, ln
import variables
import constants

diff_resultb = variables.xa(variables.mu, variables.z2, variables.zl2)
for _ in range(int(variables.k)):
    diff_resultb = diff(diff_resultb, variables.zl)
b = diff_resultb
diff_resulta = variables.xa(variables.mu, variables.z2, variables.zl2)
for _ in range(int(variables.k)):
    diff_resulta = diff(diff_resulta, variables.z)
a = diff_resulta
variables.xa(variables.mu, variables.z1, variables.zl1) * variables.xa(variables.v, variables.z2, variables.zl2) == -variables.SlopeRegge/2*variables.etaa(variables.mu*variables.v)*ln(variables.z12)**variables.D2 + variables.point(variables.xa(variables.v)*variables.xa(variables.mu, variables.z2, variables.zl2)) + summation(1/math.factorial(variables.k)*(variables.z12**variables.k * variables.point(variables.xa(variables.v)*a) + variables.zl12**variables.k*variables.point(variables.xa(variables.v)*b)), (variables.k, 1, math.inf)) # check free taylor expansion pg 38

#theta functions 214.31
variables.thetavar(variables.v, variables.tau) == summation(math.exp(math.pi * variables.i * variables.n**2 * variables.tau + 2 * math.pi * variables.i * variables.n * variables.v), (variables.n, -math.inf, math.inf))
#periodicity
variables.thetavar(variables.v + 1, variables.tau) == variables.thetavar(variables.v, variables.tau)
variables.thetavar(variables.v + variables.tau, variables.tau) == math.exp(-math.pi * variables.i * variables.tau - 2 * math.pi * variables.i * variables.v) * variables.thetavar(variables.v, variables.t)
variables.thetavar(variables.v, variables.tau + 1) == variables.thetavar(variables.v + 1/2, variables.tau)
variables.thetavar(variables.v/variables.tau, -1/variables.tau) == (-variables.i * variables.tau)**(1/2) * math.exp(-math.pi * variables.i * variables.va(2)/(variables.tau)) * variables.thetavar(variables.v, variables.tau)
variables.thetavar(variables.v, variables.tau) == np.prod((1-variables.pa(variables.m))*(1 + variables.z * variables.qa(variables.m - 1/2)) * (1 + variables.za(-1) * variables.pa(variables.m-1/2)), (variables.m, 1, math.inf)) # so long as
variables.q == math.exp(2 * math.pi * variables.i * variables.tau)
variables.z == math.exp(2 * math.pi * variables.i * variables.v)
variables.jacobithetavar(variables.a, variables.b, variables.v, variables.tau) == math.exp(math.pi * variables.i * variables.aa(2) * variables.tau + 2 * math.pi * variables.i * variables.a(variables.v + variables.b)) * variables.thetavar(variables.v + variables.a * variables.tau + variables.b, variables.tau) == summation(math.exp(math.pi * variables.i * (variables.n + variables.a)**2 * variables.tau + 2 * math.pi * variables.i * (variables.n + variables.a) * (variables.v + variables.b)), (variables.n, -math.inf, math.inf))
variables.thetavarb((0, 0), variables.v, variables.tau) == variables.thetavarb(3, variables.v, variables.tau) == variables.jacobithetavar(0, 0, variables.v, variables.tau) == summation(variables.qa((variables.na(2))/2) * variables.za(variables.n), (variables.n, -math.inf, math.inf))
variables.thetavarb((0, 1), variables.v, variables.tau) == variables.thetavarb(3, variables.v, variables.tau) == variables.jacobithetavar(0, 1/2, variables.v, variables.tau) == summation((-1)**variables.n * variables.qa(variables.na(2)/2) * variables.za(variables.n), (variables.n, -math.inf, math.inf))
variables.thetavarb((1, 0), variables.v, variables.tau) == variables.thetavarb(2, variables.v, variables.tau) == variables.jacobithetavar(1/2, 0, variables.v, variables.tau) == summation(variables.qa((variables.n - 1/2)**2/2) * variables.za(variables.n - 1/2), (variables.n, -math.inf, math.inf))
variables.thetavarb((1, 1), variables.v, variables.tau) == -variables.thetavarb(1, variables.v, variables.tau) == variables.jacobithetavar(1/2, 1/2, variables.v, variables.tau) == -variables.i * summation((-1)**variables.n * variables.qa((variables.n - 1/2)**2/2) * variables.za(variables.n - 1/2), (variables.n, -math.inf, math.inf))
#product representation
variables.thetavarb((0, 1), variables.v, variables.tau) == np.prod((1-variables.qa(variables.m))*(1 + variables.z * variables.qa(variables.m-1/2)) * (1 + variables.za(-1) * variables.qa(variables.m-1/2)), (variables.m, 1, math.inf))
variables.thetavarb((0, 1), variables.v, variables.tau) == np.prod((1 - variables.qa(variables.m)) * (1 - variables.z * variables.qa(variables.m-1/2)) * (1 - variables.za(-1) * variables.qa(variables.m-1/2)), (variables.m, 1, math.inf))
variables.thetavarb((1, 0), variables.v, variables.tau) == 2 * math.exp(math.pi * variables.i * variables.tau/4) * math.cos(math.pi * variables.v) * np.prod((1-variables.qa(variables.m)) * (1 + variables.z * variables.qa(variables.m)) * (1 + variables.za(-1) * variables.qa(variables.m)), (variables.m, 1, math.inf))
variables.thetavarb((1, 1), variables.v, variables.tau) == -2 * math.exp(math.pi * variables.i * variables.tau/4) * math.sin(math.pi * variables.v) * np.prod((1 - variables.qa(variables.m)) * (1 - variables.z * variables.qa(variables.m)) * (1 - variables.za(-1) * variables.qa(variables.m)), (variables.m, 1, math.inf))
# modular transformation
variables.thetavarb((0, 0), variables.v, variables.tau + 1) == variables.thetavarb((0, 1), variables.v, variables.tau)
variables.thetavarb((0, 1), variables.v, variables.tau + 1) == variables.thetavarb((0, 0), variables.v, variables.tau)
variables.thetavarb((1, 0), variables.v, variables.tau + 1) == math.exp(math.pi * variables.i/4) * variables.thetavarb((1, 0), variables.v, variables.tau)
variables.thetavarb((1, 1), variables.v, variables.tau + 1) =math.exp(variables.math.pi * variables.i/4) * variables.thetavarb((1, 1), variables.v, variables.tau)
variables.thetavarb((0, 0), variables.v/variables.tau, -1/variables.tau) == (-variables.i * variables.tau)**(1/2) * math.exp(math.pi * variables.i * variables.va(2)/variables.tau) * variables.thetavarb((0, 0), variables.v, variables.tau)
variables.thetavarb((0, 1), variables.v/variables.tau, -1/variables.tau) == (- variables.i * variables.tau)**(1/2) * math.exp(math.pi * variables.i * variables.va(2)/variables.tau) * variables.thetavarb((1, 0), variables.v, variables.tau)
variables.thetavarb((1, 0), variables.v/variables.tau, -1/variables.tau) == (- variables.i * variables.tau)**(1/2) * math.exp(math.pi * variables.i * variables.va(2)/variables.tau) * variables.thetavarb((0, 1), variables.v, variables.tau)
variables.thetavarb((1, 1), variables.v/variables.tau, -1/variables.tau) == -variables.i * (- variables.i * variables.tau)**(1/2) * math.exp(math.pi * variables.i * variables.va(2)/variables.tau) * variables.thetavarb((1, 1), variables.v, variables.tau)
variables.jacobiriemann(4, (0, 0), 0, variables.tau) - variables.jacobiriemann(4, (0, 1), 0, variables.tau) - variables.jacobiriemann(4, (1,0), 0, variables.tau) == 0 # quartic identity of a Riemann
variables.thetavarb((1, 1), 0, variables.tau) == 0
variables.etaf(variables.tau) == variables.qa(1/24) * np.prod(1 - variables.qa(variables.m), (variables.m, 1, math.inf)) == variables.brackets((diff(variables.thetavarb((1, 1), 0, variables.tau), variables.v))/(-2 * math.pi))**(1/3)
#modular transofrmations
variables.etaf(variables.tau + 1) == math.exp(variables.i * math.pi /12) * variables.etaf(variables.tau)
variables.etaf(-1/variables.tau) == (-variables.i * variables.tau)**(1/2) * variables.etaf(variables.tau)
