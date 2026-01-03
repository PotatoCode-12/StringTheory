import math
from scipy.integrate import quad
import numpy as np 
from sympy import symbols, diff,integrate, Function, Derivative, summation, oo, Abs
import variables
import constants

variables.MOb(variables.r) == variables.gOb(variables.r)/(variables.holderb(variables.diffxweyl, variables.r)) # 150.1
variables.MOb(variables.r, variables.n) == (variables.gOb(variables.r) * variables.MOa(variables.n))/(variables.holderb(variables.diffxweyl, variables.r))
variables.delta * variables.gb(variables.a * variables.b) == -2 * variables.holderb(variables.P1 * variables.delta * variables.sigma, variables.a * variables.b) + (2 * variables.delta * variables.omegas - variables.Delta * variables.delta * variables.sigma) * variables.gb(variables.a * variables.b)
0 == integrate(variables.d2sigma * variables.g**(1/2) * variables.deltap * variables.gb(variables.a * variables.b)) * variables.brackets(-2 * variables.holdera(variables.P1 * variables.delta * variables.sigma, variables.a * variables.b) + (2 * variables.delta * variables.omegas * variables.Delta * variables.delta * variables.sigma) * variables.gb(variables.a * variables.b)) == integrate(variables.d2sigma * variables.g**(1/2)) * variables.brackets(-2 * variables.holderb(variables.P1a(variables.T) * variables.deltap * variables.g, variables.a) * variables.delta * variables.sigmaa(variables.a) + variables.deltap * variables.gb(variables.a * variables.b) * variables.ga (variables.a * variables.b) * (2 * variables.delta * variables.omegas - variables.Delta * variables.delta * variables.sigma))
variables.ga(variables.a * variables.b) * variables.deltap * variables.gb(variables.a * variables.b) == 0
variables.holderb(variables.P1a(variables.T) * variables.deltap * variables.g, variables.a) == 0
variables.holderb(variables.P1 * variables.delta * variables.sigma, variables.a * variables.b)# confromal killing equation
diff(variables.deltap * variables.gb(variables.z *( variables.z)), variables.zl) == diff(variables.deltap * variables.gb(variables.zl * variables.zl), variables.z) ==  0
diff(variables.delta * variables.z, variables.zl) == diff(variables.delta * variables.zl, variables.z) == 0
variables.mu - variables.keta == - 3 * variables.chi
integrate(variables.d2sigma * variables.g**(1/2) * variables.holderb(variables.P1 * variables.delta * variables.sigma, variables.a * variables.b) * variables.holdera(variables.P1 * variables.delta * variables.sigma, variables.a * variables.b)) == integrate(variables.d2sigma * variables.g **(1/2) * variables.delta * variables.sigmab(variables.a) * variables.holdera(variables.P1a(variables.T) * variables.P1 * variables.delta * variables.sigma, variables.a)) == integrate(variables.d2sigma * variables.g**(1/2) * (1/2 * variables.Deltab(variables.a) * variables.delta * variables.sigmab(variables.b) * variables.Deltaa(variables.a) * variables.delta * variables.sigmaa(variables.b) - variables.R/4 * variables.delta * variables.sigmab(variables.a) * variables.delta * variables.sigmaa(variables.a)))
# rules
variables.chi > 0
variables.keta == 3* variables.chi # then
variables.mu == 0
variables.chi < 0
variables.keta == 0 # then
variables.mu == -3 * variables.chi
#153.12
variables.sigmaab(variables.a, variables.m) == variables.fabf(variables.a, variables.m * variables.n, variables.sigmab(variables.n))
variables.zab(variables.a, variables.m) == variables.fabf(variables.a, variables.m *variables.n, variables.zb(variables.n))
# conditions
variables.w == variables.w + 2 *math.pi == variables.w + 2 * math.pi * variables.tau 

# 155.2
variables.brackets(diff(variables.g)) * diff(variables.sigma, 1, variables.n * 2) == variables.brackets(diff(variables.zeta)) * diff(variables.t, 1, variables.mu) * diff(variables.sigma, 1, 2*variables.n - variables.keta)
1 == variables.Deltabf(variables.FP, (variables.g, variables.sigma)) * integrate(diff(variables.mu) * variables.t * integrate(variables.brackets(diff(variables.zeta)) * variables.delta  *(variables.g - variables.hat(variables.g) * variables.t**variables.zeta) * np.prod(variables.delta * (variables.sigmaab(variables.a, variables.i) - variables.hat(variables.sigmaab(variables.zeta * variables.a, variables.i))), ((variables.a, variables.i), variables.f)), variables.diffxweyl), variables.F)
variables.delta * variables.gb(variables.a * variables.b) == summation(variables.delta * variables.ta(variables.k) * diff(variables.hat(variables.gb(variables.a * variables.b)), variables.ta (variables.k)) - 2 * variables.holderb(variables.hat(variables.P1) * variables.delta * variables.sigma, variables.a * variables.b) + (2 * variables.delta * variables.omegas - variables.hat(variables.Delta) * variables.delta * variables.sigma) * variables.hat(variables.gb(variables.a * variables.b)))
variables.Deltabf(variables.FP, (variables.hat(variables.g), variables.hat(variables.sigma))**(-1)) == variables.nb(variables.R) * integrate(variables.da(variables.mu) * variables.delta * variables.t)* variables.brackets(variables.d * (variables.delta * variables.omegas) * variables.d * (variables.delta * variables.sigma)) * variables.delta * (variables.delta * variables.gb(variables.a * variables.b)) * np.prod(variables.delta * (variables.delta * variables.sigmaa(variables.a) * variables.hat(variables.sigmab(variables.i))), ((variables.a, variables.i), variables.f)) == variables.nb(variables.R) * integrate(variables.da(variables.mu) * variables.delta * variables.t * variables.da(variables.keta) * variables.x) * variables.brackets(variables.d * variables.Betap * variables.d * variables.delta * variables.sigma) * math.exp(2 * math.pi * variables.i * (variables.Betap, 2 * variables.hat(variables.P1) * variables.delta * variables.sigma - variables.delta * variables.ta(variables.k) * diff(variables.hat(variables.g), variables.k)) + variables.i * 2 * math.pi * summation(variables.xbs(variables.a * variables.i) * variables.delta * variables.sigmaa(variables.a) * (variables.hat(variables.sigmab(variables.i))), ((variables.a * variables.i), variables.f)))
variables.delta * variables.sigmaa(variables.a) == variables.ca(variables.a)
variables.Betabp(variables.a * variables.b) == variables.bb(variables.a * variables.b)
variables.xbs(variables.a * variables.i) == variables.etab(variables.a * variables.i)
variables.delta * variables.ta(variables.k) == variables.zetaa(variables.k)
variables.Deltabf(variables.FP, (variables.hat(variables.g), variables.hat(variables.sigma))) == 1/(variables.nb(variables.R)) * integrate(variables.brackets(diff(variables.b) * diff(variables.c)) * variables.daf(variables.mu, variables.zeta * variables.daf(variables.eta, variables.keta))) * math.exp(-1/(4 * math.pi) * (variables.b, 2* variables.P1 * variables.c - variables.zetaa(variables.k) * diff(variables.hat(variables.g), variables.k)) + summation(variables.etab(variables.a * variables.i) * variables.ca(variables.a) * variables.hat(variables.sigmab(variables.i)), ((variables.a, variables.i), variables.f))) == 1/(variables.nb(variables.R)) * integrate(variables.brackets(diff(variables.b) * diff(variables.c)) * math.exp(-variables.Sb(variables.g)) * np.prod(1/(4*math.pi) * (variables.b, diff(variables.hat(variables.g), variables.k)) * np.prod(variables.ca(variables.a) * variables.hat(variables.sigmab(variables.i)), ((variables.a, variables.i), variables.f)), (variables.k, 1, variables.mu)))

# surfaces 166.1, sphere
variables.u == 1/variables.z
diff(variables.s)**2 == math.exp(2 * variables.omegas (variables.z, variables.zl)) * diff(variables.z * diff(variables.zl))
diff(variables.s)**2 == (4 * variables.r**2 * diff(variables.z) * diff(variables.zl))/(1+variables.z * variables.zl)**2 == (4 * variables.r**2 * diff(variables.u) * diff(variables.ul))/(1+variables.u * variables.ul)**2
variables.delta * variables.u == (diff(variables.u, variables.z))/(diff(variables.z, variables.z))**(-2) * variables.delta * variables.gb(variables.z * variables.z) == variables.z**4 * variables.delta * variables.gb(variables.z * variables.z)
#general CKV
variables.delta * variables.z == variables.ab(0) + variables.ab(1) * variables.z + variables.ab(2) ** variables.za(2)
variables.delta * variables.zl == variables.abp(0) + variables.abp(1) * variables.zl + variables.abp(2) * variables.zla(2)
variables.zp == (variables.alpha * variables.z + variables.Beta)/(variables.gamma * variables.z + variables.delta)
#for disk
variables.zp == 1/(variables.zl)
variables.zp == variables.zl
#projectile plane
variables.zp == -1/variables.zl 
#reference prior terms x2

#Riemann surfaces, torus, 206.1
variables.w == variables.w + 2 * math.pi == variables.w + 2  * math.pi * variables.tau 
(variables.sigmaa(1), variables.sigmaa(2)) == (variables.sigmaa(1) + 2 * math.pi, variables.sigmaa(2)) == (variables.sigmaa(1) + 2 * math.pi * variables.taub(1), variables.sigmaa(2) + 2 * math.pi * variables.taub(2))
variables.z == variables.z * math.exp(-2 * math.pi * variables.i * variables.tau)
1 <= Abs(variables.z) <= math.exp(2 * math.pi * variables.taub(2))
# cylinder 207.5
0 <= variables.Re * variables.w <= math.pi
variables.w == variables.w + 2 * math.pi * variables.i * variables.t 
variables.wp == -variables.wl
(variables.sigmaa(1), 0) == (variables.sigmaa(1), 2 * math.pi * variables.t)
#klein bottle 207.8
variables.w == variables.w + 2 * math.pi == - variables.wl + 2 * math.pi * variables.i * variables.t 
(variables.sigmaa(1), variables.sigmaa(2)) == (variables.sigmaa(1) + 2 * math.pi, variables.sigmaa(2)) == (-variables.sigmaa(1), variables.sigmaa(2) + 2 * math.pi * variables.t)
variables.wp == - variables.wl + 2 * math.pi * variables.i * variables.t 
# mobius strip
0 <= variables.Re * variables.w <= math.pi
variables.w == - variables.wl + math.pi + 2 * math.pi * variables.i * variables.t 
variables.wp == -variables.wl
variables.wp == variables.w + math.pi * (2 * variables.i * variables.t + 1)
