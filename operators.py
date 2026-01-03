import math
from scipy.integrate import quad
import numpy as np 
from sympy import symbols, diff,integrate, Function, Derivative, summation, oo, ln, Abs
import sympy as sp
import variables
import constants

NO1 = variables.order(variables.xa(variables.mu, variables.z, variables.zl)) == variables.xa(variables.mu, variables.z, variables.zl)
NO2 = variables.order(variables.xa(variables.mu, variables.z1, variables.zl1)*variables.xa(variables.v, variables.z2, variables.zl2)) == variables.xa(variables.mu, variables.z1, variables.zl1) * variables.xa(variables.v, variables.z2, variables.zl2) + variables.SlopeRegge/2*variables.etaa(variables.mu, variables.v)*ln(variables.z12)*(variables.D2)
variables.zb(variables.i, variables.j) == variables.zb(variables.i) - variables.z(variables.j) # normal ordering and parameters check pg 36
diff(diff(variables.point(variables.xa(variables.mu, variables.z1, variables.zl1) * variables.xa(variables.v, variables.z2, variables.zl2)), variables.zl1), variables.z1)
diff(diff(ln(variables.z)(variables.D2), variables.zl), variables.z) == 2*math.pi* variables.delta(variables.D2, variables.z, variables.zl)

#local operators pg 37
variables.AO(variables.i, variables.sigmab(variables.D1)) * variables.AO(variables.j, variables.sigmab(variables.D2)) == summation(variables.cab(variables.k, variables.i*variables.j)*(variables.sigmab(variables.D1) - variables.sigmab(variables.D2))*variables.AO(variables.k, variables.sigmab(variables.D2)), (variables.k))

variables.point(variables.NO) * variables.point(variables.GO) == math.exp(-variables.SlopeRegge/2*integrate(variables.d2z1*variables.d2z2*ln(variables.z12)**variables.D2*variables.delta/(variables.delta*variables.xab(variables.mu, variables.F, variables.z1, variables.zl1)) * variables.delta/(variables.delta*variables.xb(variables.G, variables.mu, variables.z2, variables.zl2)))) * variables.point(variables.NO*variables.GO) #39 expansion of contradicting pairs OPE

diff_resulta = variables.sigma
for _ in range(int(variables.d)):
    diff_resulta = diff(diff_resulta, variables.d)
a = diff_resulta
variables.cab(variables.k, variables.i, variables.j)*(variables.sigmab(variables.D1) - variables.sigmab(variables.D2))*variables.bbb(variables.s, variables.y, variables.m) == +variables.cab(variables.k, variables.j, variables.i)(variables.sigmab(variables.D2) - variables.sigmab(variables.D1))*variables.bbb(variables.s, variables.y, variables.m) == -variables.cab(variables.k, variables.j, variables.i)(variables.sigmab(variables.D2) - variables.sigmab(variables.D1))*variables.bbb(variables.s, variables.y, variables.m) 
variables.delta*variables.AO(variables.sigmab(variables.D0)) + variables.epilson/(2*math.pi*variables.i) * quad(a*variables.g**(1/2)*variables.Deltab(variables.a)*variables.ja(variables.a)*variables.AO(variables.sigmab(variables.D0)), variables.R) == 0 #42
#residues, pg 42
variables.res(variables.z, variables.zb(variables.D0))*variables.jf(variables.z) * variables.AO(variables.zb(variables.D0), variables.zlb(variables.D0)) + variables.resl(variables.zl, variables.zlb(variables.D0)) * variables.jtilde(variables.z)*variables.AO(variables.zb(variables.D0), variables.zlb(variables.D0)) == 1/(variables.i*variables.epsilon)*variables.delta*variables.AO(variables.zb(variables.D0), variables.zlb(variables.D0))
#after spacetime translation and example case scenario
diff_resultb = variables.xa(variables.mu)
for _ in range(int(variables.a)):
    diff_resultb = diff(diff_resultb)
b = diff_resultb
variables.delta*variables.S == (variables.epsilon*variables.ab(variables.mu))/(2*math.pi*variables.SlopeRegge)*integrate(variables.d2sigma*b*diff(variables.rho, variables.a))
variables.jab(variables.mu, variables.a) == variables.i/(variables.SlopeRegge)*variables.xa(variables.mu)

#check conditions pg 46
variables.Tf(variables.z)*variables.OO(variables.D0, variables.D0) == variables.h/(variables.z**variables.D2)*variables.OO(variables.D0, variables.D0) + 1/variables.z*diff(variables.OO(variables.D0, variables.D0), variables.z) # increase terms for precision

variables.AOb(variables.i, variables.zb(variables.D1), variables.zlb(variables.D1)) * variables.AOb(variables.j, variables.zb(variables.D2), variables.zlb(variables.D2)) == summation(variables.zab(variables.hb(variables.k) - variables.hb(variables.i) - variables.hb(variables.j), variables.D12)*variables.zlab(variables.hbtilde(variables.k)-variables.hbtilde(variables.i)-variables.hbtilde(variables.j), variables.D12)* variables.cab(variables.k, variables.i, variables.j) * variables.AOb(variables.k, variables.z2, variables.zl2), variables.k) #47

#OPE of energy momentum tensor
variables.Tf(variables.z)*variables.Tf(variables.D0) == (variables.etaab(variables.mu, variables.mu))/(2*variables.za(variables.D4)) - 2/(variables.SlopeRegge*variables.za(variables.D2)) * variables.point(diff(variables.xa(variables.mu, variables.z), variables.z) * diff(variables.xb(variables.mu, variables.D0), variables.z)) + variables.point(variables.Tf(variables.z)*variables.Tf(variables.D0)) # check excluded approximate equal pg 47

#general 60.7
variables.brackets1(variables.GO) == math.exp(1/2*integrate(variables.d2z*variables.d2zp*variables.deltaf(variables.z, variables.zl, variables.zp, variables.zpl) * variables.delta/(variables.delta*variables.xa(variables.mu, variables.z, variables.zl)) * variables.delta/(variables.delta*variables.xb(variables.mu, variables.zp, variables.zpl)))) * variables.brackets2(variables.GO)

#corresponding operators 65.7
variables.alphaab(variables.mu, -variables.m) == variables.i* (2/variables.SlopeRegge)**(1/2) * 1/(math.factorial(variables.m-1)) * diff(variables.xa(variables.mu, variables.D0), variables.m, variables.z)
variables.alphaabtilde(variables.mu, -variables.m) == variables.i* (2/variables.SlopeRegge)**(1/2) * 1/(math.factorial(variables.m-1)) * diff(variables.xa(variables.mu, variables.D0), variables.m, variables.zl)
variables.xabs == variables.xa(variables.mu, 0, 0)
variables.dirac(0, variables.k) == variables.point(math.exp(variables.i*variables.k*variables.X(0, 0)))
variables.bb(variables.m) * variables.dirac(1) == 0
variables.m >= -1
variables.cb(variables.m) * variables.dirac(1) == 0
variables.m >= 2
variables.dirac(1) == variables.bb(-1) * variables.dirac(variables.down)
#translation of rising operators
variables.bb(-variables.m) == 1/(math.factorial(variables.m - 2)) * diff(variables.bf(0), (variables.m-2), variables.z)
variables.m >=2
variables.cb(-variables.m) == 1/(math.factorial(variables.m+1)) * diff(variables.cf(0), (variables.m + 1), variables.z)
diff(variables.w, variables.z)*variables.jfb(variables.w, variables.w) == variables.jfb(variables.z, variables.z) + variables.qb(variables.D0) * (diff(variables.w, 2, variables.z))/(diff(variables.w, variables.z)) == variables.jfb(variables.z, variables.z) - (variables.qb(variables.D0))/variables.z #ghost transfer current
#frame charge
variables.Qa(variables.g) == 1/(2*math.pi*variables.i) * integrate(diff(variables.z) * variables.jb(variables.z)) == variables.Na(variables.g) + variables.qb(variables.D0) # conic integration

#69.2
variables.AOb(variables.i*variables.j, variables.z, variables.zl) == summation(variables.za(variables.hb(variables.k) - variables.hb(variables.i) - variables.hb(variables.j)) * variables.zla(variables.hbtilde(variables.k)-variables.hbtilde(variables.i)-variables.hbtilde(variables.j)) * variables.cab(variables.k, variables.i*variables.j) * variables.AOb(variables.k), variables.k) # eigenstates
#for 
variables.Lb(variables.m) * variables.dirac(variables.AO) == variables.Lb(variables.m) * variables.AO(variables.D0, variables.D0)
variables.Tf(variables.z) * variables.AO(variables.D0, variables.D0) == summation(variables.za(-variables.m-2) * variables.Lb(variables.m) * variables.AO(variables.D0, variables.D0), (variables.m, -math.inf, math.inf))
variables.delta * variables.AO(variables.z, variables.zl) == -variables.epsilon * summation(1/(math.factorial(variables.n)) * variables.brackets(diff(variables.vf(variables.z), variables.z) * variables.Lb(variables.n-1) + (diff(variables.vf(variables.z), variables.n, variables.z) * variables.trans * variables.Lbtilde(variables.n-1))) * variables.AO(variables.z, variables.zl), (variables.n, 0, math.inf))
# result
variables.Lb(-variables.D1) * variables.AO == diff(variables.AO, variables.z)
variables.Lbtilde(-variables.D1) * variables.AO == diff(variables.AO, variables.zl)
variables.Lb(variables.D0) * variables.AO == variables.h*variables.AO
variables.Lbtilde(variables.D0) * variables.AO == variables.htilde * variables.AO
variables.Lb(variables.D0) * variables.dirac(variables.OO) == variables.h * variables.dirac(variables.OO)
variables.Lbtilde(variables.D0) * variables.dirac(variables.OO) == variables.htilde * variables.dirac(variables.OO)
variables.Lb(variables.m) * variables.dirac(variables.OO) == variables.Lbtilde(variables.m) * variables.dirac(variables.OO) == 0
variables.m > 0
variables.Lb(variables.m) * variables.dirac(variables.D1) == variables.Lbtilde(variables.m) * variables.dirac(variables.D1) == 0
variables.m >= -1 # check full algebra 71.9 7 10

#closed string tachyon state 101.1
variables.Vb(variables.D0) == 2*variables.gb(variables.c) * integrate(variables.d2sigma * variables.g**(1/2) * math.e**(variables.i*variables.k*variables.X)) == variables.gb(variables.c) * integrate(variables.d2z * variables.point(math.e**(variables.i*variables.k*variables.X)))
variables.m**2 == -variables.k**2 == -4/variables.SlopeRegge #tensor condition
2*variables.gb(variables.c)/(variables.SlopeRegge) * integrate(variables.d2z * variables.point(diff(variables.xa(variables.mu), variables.z) * diff(variables.xa(variables.v), variables.zl) * variables.math.e**(variables.i*variables.k*variables.X))) # vertex operator 102.3
variables.h == variables.htilde == 1 + variables.SlopeRegge*variables.k**2/4
#renormalized 102.5
variables.RO == math.exp(1/2 * integrate(variables.d2sigma * variables.d2sigmap * variables.Deltaf(variables.sigma, variables.sigmap) * variables.delta/(variables.delta * variables.xa(variables.mu, variables.sigma)) * variables.delta/(variables.delta * variables.xb(variables.mu, variables.sigmap)))) * variables.FO
variables.Deltaf(variables.sigma, variables.sigmap) == variables.SlopeRegge/2 * ln(variables.d2f(variables.sigma, variables.sigmap))
variables.deltab(variables.W) * variables.RO == variables.bracketsb(variables.deltab(variables.W) * variables.FO, variables.r) + 1/2 * integrate(variables.d2sigma * variables.d2sigmap * variables.deltab(variables.W) * variables.Deltaf(variables.sigma, variables.sigmap) * variables.delta/(variables.delta * variables.xa(variables.mu, variables.sigma)) *  variables.delta/(variables.delta * variables.xb(variables.mu, variables.sigmap)) * variables.RO)
variables.deltab(variables.W) * variables.Vb(variables.D0) == 2 * variables.gb(variables.c) * integrate (variables.d2sigma * variables.g**(1/2) * (2 * variables.delta * variables.omegasf(variables.sigma) - variables.k**2/2 * variables.deltab(variables.W) * variables.Deltaf(variables.sigma, variables.sigma)) * variables.bracketsb(variables.math.e**(variables.i*variables.k*variables.Xf(variables.sigma)), variables.r))
variables.d2f(variables.sigma, variables.sigmap) == (variables.sigma - variables.sigmap)**2 * math.exp(2 * variables.omegasf(variables.sigma)) # for short distance
variables.Deltaf(variables.sigma, variables.sigmap) == variables.SlopeRegge * variables.omegasf(variables.sigma) + variables.SlopeRegge/2 * ln(variables.sigma - variables.sigmap) **2
variables.deltab(variables.W) * variables.Deltaf(variables.sigma, variables.sigma) == variables.SlopeRegge* variables.delta * variables.omegasf(variables.sigma) # weyl variation 103.11

#massless closed string vertex 104.14
variables.Vb(variables.D1) == variables.gb(variables.c)/(variables.SlopeRegge) * integrate(variables.d2sigma * variables.g**(1/2) * variables.bracketc((variables.ga(variables.a*variables.b) * variables.sb(variables.mu*variables.v) + variables.i*variables.epsilona(variables.a*variables.b) * variables.ab(variables.mu *variables.v)) * variables.bracketsb(diff(variables.xa(variables.mu), variables.a) * diff(variables.xa(variables.v), variables.b), variables.r) + variables.SlopeRegge * variables.phi * variables.R * variables.bracketsb(math.e**variables.i*variables.k*variables.X, variables.r)))
variables.deltab(variables.W) * variables.Vb(variables.D1) == variables.gb(variables.c)/2 * integrate(variables.d2sigma * variables.g**(1/2) * variables.delta * variables.omegas * variables.bracketc((variables.ga(variables.a*variables.b) * variables.Sb(variables.mu* variables.v) + variables.i * variables.epsilona(variables.a*variables.b) * variables.Ab(variables.mu*variables.v)) * variables.bracketsb(diff(variables.xa(variables.mu), variables.a) * diff(variables.xa(variables.v), variables.b) * math.e**(variables.i*variables.k*variables.X), variables.r) + variables.SlopeRegge*variables.F * variables.R * variables.bracketsb(math.e**(variables.i*variables.k*variables.X), variables.r))) # 105.16
variables.Sb(variables.mu * variables.v) == -variables.k**2 * variables.sb(variables.mu * variables.v) + variables.kb(variables.v) * variables.ka(variables.omegas) * variables.sb(variables.mu * variables.omegas) + variables.kb(variables.mu) * variables.ka(variables.omegas) * variables.sb(variables.v * variables.omegas) - (1 + variables.gamma) * variables.kb(variables.mu) * variables.kb(variables.v) * variables.sab(variables.omegas, variables.omegas) + 4 * variables.kb(variables.mu) * variables.kb(variables.v) * variables.phi
variables.Ab(variables.mu * variables.v) == -variables.k**2 * variables.ab(variables.mu * variables.v) + variables.kb(variables.v) * variables.ka(variables.omegas) * variables.ab(variables.mu * variables.omegas) - variables.kb(variables.mu) * variables.ka(variables.omegas) * variables.ab(variables.v * variables.omegas)
variables.F == (variables.gamma - 1) * variables.k**2 * variables.phi + 1/2 * variables.gamma * variables.ka(variables.mu) * variables.ka(variables.v) * variables.sb(variables.mu * variables.v) - 1/4 * (variables.gamma) * (1 + variables.gamma) * variables.k**2 * variables.sab(variables.v, variables.v)
variables.bracketsb(variables.Deltaa(variables.D2) * variables.xa(variables.mu) * math.e**(variables.i*variables.k*variables.X), variables.r) == variables.i * (variables.SlopeRegge * variables.gamma)/4 * variables.ka(variables.mu) * variables.R * variables.bracketsb(math.e**(variables.i*variables.k*variables.X), variables.r)
variables.Sb(variables.mu * variables.v) == variables.Ab(variables.mu * variables.v) == variables.F == 0 # weyl invariance 106.19
variables.na(variables.mu) * variables.sb(variables.mu * variables.v) == variables.na(variables.mu) * variables.ab(variables.mu * variables.v) == 0 # check restriction states 106.21
variables.k**2 == 0
variables.ka(variables.v) * variables.sb(variables.mu * variables.v) == variables.ka(variables.mu) * variables.ab(variables.mu * variables.v) == 0
variables.phi == (1+variables.gamma)/4 * variables.sab(variables.mu, variables.mu) # excluded later polarization framework
variables.bracketsb(math.e**(variables.i*variables.k*variables.X), variables.DR) == variables.bracketsb(math.e**(variables.i*variables.k*variables.X), variables.r)
variables.bracketsb(diff(variables.xa(variables.mu), variables.a) * math.e**(variables.i*variables.k*variables.X), variables.DR) == variables.bracketsb(diff(variables.xa(variables.mu), variables.a) * math.e**(variables.i*variables.k*variables.X), variables.r)
variables.bracketsb(diff(variables.xa(variables.mu), variables.a) * diff(variables.xa(variables.v), variables.b) *math.e**(variables.i*variables.k*variables.X), variables.DR) == variables.bracketsb(diff(variables.xa(variables.mu), variables.a) * diff(variables.xa(variables.v), variables.b) * math.e**(variables.i*variables.k*variables.X), variables.r) - variables.SlopeRegge/12 * variables.gb(variables.a*variables.b) * variables.etaa(variables.mu * variables.v) * variables.R * variables.bracketsb(math.e**(variables.i*variables.k*variables.X), variables.r)
variables.bracketsb(variables.Deltab(variables.a) * diff(variables.xa(variables.mu), variables.b)*math.e**(variables.i*variables.k*variables.X), variables.DR) == variables.bracketsb(variables.Deltab(variables.a) * diff(variables.xa(variables.mu), variables.b) * math.e**(variables.i*variables.k*variables.X), variables.r) + variables.i* variables.SlopeRegge/12 * variables.gb(variables.a*variables.b) * variables.ka(variables.mu) * variables.R * variables.bracketsb(math.e**(variables.i*variables.k*variables.X), variables.r)
#open string vertex operators 107.25
variables.gb(variables.o) * integrate(diff(variables.s) * variables.bracketsb(math.e**(variables.i*variables.k*variables.X), variables.r), 0, 2*math.pi) # contour
-variables.i*(variables.gb(variables.o))/(2*variables.SlopeRegge)**(1/2) * variables.eb(variables.mu) * integrate(diff(variables.s) * variables.bracketsb(diff(variables.xa(variables.mu), variables.tau) * math.e**(variables.i*variables.k*variables.X), variables.r), 0, 2*math.pi)

# 169.1
variables.Zf(variables.brackets(variables.J)) == variables.bracket(math.exp(variables.i * integrate(variables.d2sigma * variables.Jf(variables.sigma) * variables.Xf(variables.sigma))))
variables.xa(variables.mu, variables.sigma) == summation(variables.xabs(variables.mu, variables.I) * variables.Xbfn(variables.I, variables.sigma), variables.I)
variables.Deltaa(2) * variables.Xbn(variables.I) == -variables.omegaabs(2, variables.I) * variables.Xbn(variables.I)
integrate(variables.d2sigma * variables.g**(1/2) * variables.Xbn(variables.I) * variables.Xbn(variables.Ip), variables.M) == variables.deltab(variables.I * variables.Ip)
variables.Zf(variables.brackets(variables.J)) == np.prod(integrate(diff(variables.xabs(variables.mu, variables.I) * math.exp(-(variables.omegaabs(2, variables.I) * variables.xbas(variables.I, variables.mu) * variables.xbs((variables.I, variables.mu))/(4 * math.pi * variables.SlopeRegge) + variables.i * variables.xabs(variables.mu, variables.I) * variables.Jb((variables.I, variables.mu)))))), variables.I, variables.mu)
variables.Jab(variables.mu, variables.I) == integrate(variables.d2sigma * variables.Jaf(variables.mu, variables.sigma)* variables.Xbfn(variables.I, variables.sigma))
variables.Xbn(0) == (integrate(variables.sigma* variables.g**(1/2)))**(2)
variables.Zf(variables.brackets(variables.J)) == variables.i * variables.holdera(2 * math.pi, variables.a) * variables.deltaaf(variables.d, variables.Jb(0)) * (np.prod(variables.holdera((4 * math.pi**2 * variables.SlopeRegge)/(variables.omegaabs(2, variables.I)), (variables.d/2)) * math.exp(-(math.pi * variables.SlopeRegge * variables.Jb(variables.I * variables.Jb(variables.I))/(variables.omegaabs(2, variables.I)))), (variables.I, -math.inf, -1)) + np.prod(variables.holdera((4 * math.pi**2 * variables.SlopeRegge)/(variables.omegaabs(2, variables.I)), (variables.d/2)) * math.exp(-(math.pi * variables.SlopeRegge * variables.Jb(variables.I * variables.Jb(variables.I))/(variables.omegaabs(2, variables.I)))), (variables.I, 1, math.inf))) == variables.i * variables.holdera(2 * math.pi, variables.d) * variables.deltaaf(variables.d, variables.Jb(0)) * variables.holdera(variables.detp((-variables.Deltaa(2))/(4 * math.pi * variables.SlopeRegge)), (-variables.d/2)) * math.exp(-1/2 * integrate(variables.d2sigma * variables.d2sigmap * variables.Jf(variables.sigma) * variables.Jf(variables.sigmap) * variables.Gfp(variables.sigma, variables.sigmap)))
variables.Gfp(variables.sigmab(1), variables.sigmab(2)) == summation((2 * math.pi * variables.SlopeRegge)/(variables.omegaabs(2, variables.I) * variables.Xbfn(variables.I, variables.sigmab(1)) * variables.Xbfn(variables.I, variables.sigmab(2))), (variables.I, -math.inf, -1)) +  summation((2 * math.pi * variables.SlopeRegge)/(variables.omegaabs(2, variables.I) * variables.Xbfn(variables.I, variables.sigmab(1)) * variables.Xbfn(variables.I, variables.sigmab(2))), (variables.I, 1, math.inf))
-1/( 2 * math.pi * variables.SlopeRegge) * variables.Deltaa(2) * variables.Gfp(variables.sigmab(1), variables.sigmab(2)) == summation(variables.Xbfn(variables.sigmab(1)) * variables.Xbfn(variables.I, variables.sigmab(2)), (variables.I, -math.inf, -1)) + summation(variables.Xbfn(variables.sigmab(1)) * variables.Xbfn(variables.I, variables.sigmab(2)), (variables.I, -1, math.inf)) == variables.g**(-1/2) * variables.deltaaf(2, variables.sigmab(1) - variables.sigmab(2)) - variables.Xabn(2, 0)
# ZE SPHERE 170.9
variables.Gfp(variables.sigmab(1), variables.sigmab(2))  == -variables.SlopeRegge/2 * ln(variables.z12)**2 + variables.f(variables.z1 + variables.zl1) + variables.f(variables.z2, variables.zl2)
variables.f(*variables.z, variables.zl) == (variables.SlopeRegge * variables.Xabn(2, 0))/4  * integrate(variables.d2zp * math.exp(variables.brackets(2 * variables.omegasf(variables.zp, variables.zpl))) * ln(variables.z-variables.zp)**2 +variables.k)
#check tachyon vertex operators 171.11
variables.Jf(variables.sigma)  == summation(variables.kb(variables.i) * variables.deltaaf(2, variables.sigma - variables.sigmab(variables.i)), (variables.i, 1, variables.n))
if variables.i < variables.j:
    variables.Aabf(variables.n, variables.Sb(2), (variables.k, variables.sigma)) == variables.i * variables.Cab(variables.X, variables.Sb(2)) * variables.holdera(2 * math.pi, variables.d) * variables.deltaa(variables.d) * (summation(variables.kb(variables.i)), variables.i) * math.exp(-summation(variables.kb(variables.i) * variables.kb(variables.j) * variables.Gfp(variables.sigmab(variables.i), variables.sigmab(variables.j)) - 1/2 * summation(variables.kab(2, variables.i) * variables.Gbfp(variables.r, (variables.sigmab(variables.i), variables.sigmab(variables.i))), (variables.i, 1, variables.n))), ((variables.i, variables.j), 1, variables.n))
variables.Cab(variables.X, variables.Sb(2)) == variables.Xabn(-variables.d, 0) * variables.holderab(variables.detp((-variables.Deltaa(2))/(4 * math.pi**2 * variables.SlopeRegge)), (-variables.d/2), variables.Sb(2))
variables.Gbfp(variables.r, (variables.sigma, variables.sigmap)) == variables.Gfp(variables.sigma, variables.sigmap) + variables.SlopeRegge/2 * ln(variables.d2(variables.sigma, variables.sigmap))
variables.Gbfp(variables.sigma, variables.sigma)  == 2 * variables.f(variables.z, variables.zl) + variables.SlopeRegge * variables.omegasf(variables.z, variables.zl)
#path integral
variables.Aabf(variables.n, variables.Sb(2), (variables.k, variables.sigma)) == variables.i * variables.Cab(variables.X, variables.Sb(2)) * variables.holdera(2 * math.pi, variables.d) * variables.deltaa(variables.d) * (summation(variables.kb(variables.i), variables.i)) * math.exp(-variables.SlopeRegge/2 * summation(variables.kab(2, variables.i) * variables.omegaf(variables.sigmab(variables.i)), variables.i)) * np.prod(Abs(variables.zb(variables.i * variables.j)**variables.SlopeRegge * variables.kb(variables.i) * variables.kb(variables.j)), ((variables.i, variables.j), 1, variables.n, variables.i <variables.j)) # double check bounds on product expansion
#reference ignored additions included in final calculations 172.18-19
variables.vaf(variables.mu, variables.z) == - variables.i  * variables.SlopeRegge/2 * summation((variables.kab(variables.mu, variables.i))/(variables.z - variables.zb(variables.i)), (variables.i, 1, variables.n))
variables.vabtilde(variables.mu, variables.zl) == -variables.i * variables.SlopeRegge/2 * summation((variables.kab(variables.mu, variables.i))/(variables.zl - variables.zlb(variables.i)), (variables.i, 1, variables.n))
#172.21 - 173.25 excluded for example reasons/additional terms
variables.Aabf(variables.n, variables.Sb(2), (variables.k, variables.sigma)) * summation(variables.kab(variables.mu, variables.i), (variables.i, 1, variables.n)) == 0
variables.pa(variables.mu) == 1/(2 * math.pi * variables.i) * integrate(diff(variables.z) * variables.jab(variables.mu, variables.z) - diff(variables.zl) * variables.jab(variables.mu, variables.zl), variables.C) # contour integral
# final result OPE 174.28
- variables.i * variables.SlopeRegge/2 * variables.Aabf(variables.n, variables.Sb(2), (variables.k, variables.sigma)) * (variables.kab(variables.mu, 1)/(variables.z - variables.zb(1)) + summation((variables.kab(variables.mu, variables.i))/(variables.zb(1) - variables.zb(variables.i)) + variables.Of(variables.z - variables.zb(1)), (variables.i, 2, variables.n)))
variables.i * variables.kb(1) * diff(variables.Xf(variables.z)) * variables.point(math.e**(variables.i * variables.kb(1) * variables.Xf(variables.zb(1), variables.zlb(1)))) == (variables.SlopeRegge * variables.kab(2, 1))/(2 * (variables.z - variables.zb(1))) * variables.point(math.e**(variables.i * variables.kb(1) * variables.Xf(variables.zb(1), variables.zlb(1)))) + diff(variables.point(math.e**(variables.i * variables.kb(1) * variables.Xf(variables.zb(1), variables.zlb(1)))), variables.zb(1)) + variables.Of(variables.z - variables.zb(1))
diff(variables.Aabf(variables.n, variables.Sb(2), (variables.k, variables.sigma)), variables.zb(1)) == variables.SlopeRegge/2 * variables.Aabf(variables.n, variables.Sb(2), (variables.k, variables.sigma)) * summation((variables.kb(1) * variables.kb((variables.i)))/(variables.zb(1 * variables.i)), (variables.i, 2, variables.n))
# voided normalization 174.31
# for disk
variables.Gfp(variables.sigmab(1), variables.sigmab(2)) == -variables.SlopeRegge/2 * ln(Abs(variables.z1 - variables.z2)**2) - variables.SlopeRegge/2 * ln (Abs(variables.z1 - variables.zl2)**2)
if variables.i < variables.j:
    variables.bracketb(np.prod(variables.point(math.e**(variables.i * variables.kb(variables.i) * variables.Xf(variables.zb(variables.i), variables.zlb(variables.i))))), variables.Db(2)) == variables.i * variables.Cab(variables.X, variables.Db(2)) * variables.holdera(2 * math.pi, variables.d) * variables.deltaa(variables.d) * (summation(variables.kb(variables.i)), variables.i) * np.prod(Abs(variables.zb(variables.i) - variables.zlb(variables.i))**(variables.SlopeRegge * (variables.kab(2, variables.i))/2), (variables.i, 1, variables.n)) * np.prod(Abs(variables.zb(variables.i) - variables.zb(variables.j))**(variables.SlopeRegge * variables.kb(variables.i) * variables.kb(variables.j)) * Abs(variables.zb(variables.i) - variables.zlb(variables.j))**(variables.SlopeRegge * variables.kb(variables.i) * variables.kb(variables.j)), ((variables.i, variables.j), 1, variables.n))
variables.point(variables.xa(variables.mu, (variables.yb(1))) * variables.xa(variables.v, variables.yb(2))) == variables.xa(variables.mu, variables.yb(1)) * variables.xa(variables.v, variables.yb(2)) + 2 * variables.SlopeRegge * variables.etaa(variables.mu * variables.v) * ln(Abs(variables.yb(1) - variables.yb(2)))
if variables.i< variables.J : 
    variables.bracketb(np.prod(variables.point(math.e**(variables.i * variables.kb(variables.i) * variables.Xf(variables.yb(variables.i))))), (variables.i, 1, variables.n), variables.Db(2)) == variables.i * variables.Cab(variables.X, variables.Db(2)) * variables.holdera(2 * math.pi, variables.d) * variables.deltaa(variables.d) * summation(variables.kb(variables.i), variables.i) * np.prod(Abs(variables.yb(variables.i) - variables.yb(variables.j)**(2 * variables.SlopeRegge * variables.kb(variables.i) * variables.kb(variables.j))), ((variables.i, variables.j), 1, variables.n))
variables.bracketb(np.prod(variables.point(math.e**(variables.i * variables.kb(variables.i) * variables.Xf(variables.yb(variables.i)))), (variables.i, 1, variables.n)) * np.prod(diff(variables.xa(variables.mu * variables.j, variables.ybp(variables.j)), variables.y), (variables.j, 1, variables.p)), variables.Db(2)) == variables.i * variables.Cab(variables.X, variables.Db(2)) * variables.holdera(2 * math.pi, variables.d) * variables.deltaa(variables.d) * summation(variables.kb(variables.i), variables.i) * np.prod(Abs(variables.yb(variables.i * variables.j))**(2 * variables.SlopeRegge * variables.kb(variables.i) * variables.kb(variables.j)) * variables.bracketb(np.prod(variables.brackets(variables.vaf(variables.mu * variables.j, variables.ybp(variables.j)) + variables.qaf(variables.mu * variables.j, variables.ybp(variables.j)), (variables.j, 1, variables.p)), variables.Db(2))), ((variables.i, variables.j), 1, variables.n, variables.i<variables.j)) # check for prod conditions
variables.vaf(variables.mu, variables.y) == -2 * variables.i * variables.SlopeRegge * summation((variables.kab(variables.mu, variables.i))/(variables.y - variables.yb(variables.i)), (variables.i, 1, variables.n))
variables.Gfp(variables.sigmab(1), variables.sigmab(2)) == - variables.SlopeRegge/2 * ln(Abs(variables.z1 - variables.z2)**2) - variables.SlopeRegge/2 * ln(Abs(1 + variables.z1 * variables.zl2)**2)
variables.bracketb(np.prod(variables.point(math.e**(variables.i * variables.kb(variables.i) * variables.Xf(variables.zb(variables.i), variables.zlb(variables.i)))), (variables.i, 1,variables.n)), variables.R * variables.Pb(2)) == variables.i * variables.Cab(variables.X, variables.R * variables.Pb(2)) * variables.holdera(2 * math.pi, variables.d) * variables.deltaa(variables.d) * summation(variables.kb(variables.i), variables.i) * np.prod(Abs(1 + variables.zb(variables.i) * variables.zlb(variables.i)**(variables.SlopeRegge * (variables.kab(2, variables.i))/2)), (variables.i, 1, variables.n)) * np.prod(Abs(variables.zb(variables.i) - variables.zb(variables.j))**(variables.SlopeRegge * variables.kb(variables.i) * variables.kb(variables.j)) * Abs(1 + variables.zb(variables.i) * variables.zlb(variables.j))**(variables.SlopeRegge ** variables.kb(variables.i) * variables.kb(variables.j)), ((variables.i, variables.j), 1, variables.n, variables.i<variables.j))

#bc CFT, check ordering 176.3
variables.Can(variables.z) == 1, variables.z, variables.za(2)
variables.Can(variables.zl) == 0
variables.Can(variables.z) == 0
variables.Can(variables.zl) == 1, variables.zl, variables.zla(2)
variables.Cab(variables.g, variables.Sb(2)) * sp.Determinant(sp.Matrix([[1, 1, 1],
                                                                        [variables.zb(1), variables.zb(2), variables.zb(3)],
                                                                        [variables.zab(2, 1), variables.zab(2, 2), variables.zab(2, 3)]])) * sp.Determinant([[1, 1, 1], 
                                                                                                                                                             [variables.zlb(4), variables.zlb(5), variables.zlb(6)], 
                                                                                                                                                             [variables.zlab(2, 4), variables.zlab(2, 5), variables.zlab(2, 6)]]) == variables.Cab(variables.g, variables.Sb(2)) * variables.zb(12) * variables.zb(13) * variables.zb(23) * variables.zlb(45) * variables.zlb(46) * variables.zlb(56)
# ignored path integral bc permutations 177.5
integrate((diff(variables.z))/(2 * math.pi * variables.i) * variables.jb(variables.z), variables.C) == - integrate((diff(variables.u))/(2 * math.pi * variables.i) * variables.jb(variables.u) + 3, variables.C) == 3 # contour integrals
variables.bftilde(variables.zl) == variables.bf(variables.zp)
variables.cftilde(variables.zl) == variables.cf(variables.zp)
variables.zp == variables.zl 
variables.i * variables.z > 0 # for im permutations
variables.bracketb(variables.cf(variables.zb(1)) * variables.cf(variables.zb(2)) * variables.cf(variables.zb(3)), variables.Db(2)) == variables.Cab(variables.g, variables.Db(2)) * variables.zb(12) * variables.zb(13) * variables.zb(23)
variables.bracketb(variables.cf(variables.zb(1)) * variables.cf(variables.zb(2)) * variables.cbtilde(variables.zlb(3)), variables.Db(2)) == variables.bracketb(variables.cf(variables.zb(1)) * variables.cf(variables.zb(2)) * variables.cf(variables.zbp(3)), variables.Db(2)) == variables.Cab(variables.g, variables.Db(2)) * variables.zb(12) * (variables.zb(1) - variables.zlb(3)) * (variables.zb(2) - variables.zlb(3))
variables.bftilde(variables.zl) == ((diff(variables.zp, variables.z))/(diff(variables.zl, variables.z)))**2 * variables.bf(variables.zp) == variables.zpa(4) * variables.bf(variables.zp) # with next term 178.12
variables.cftilde(variables.zl) == ((diff(variables.zp, variables.z))/(diff(variables.zl, variables.z)))**(-1) * variables.cf(variables.zp) == variables.zpa(-2) * variables.cf(variables.zp)
variables.bracketb(variables.cf(variables.zb(1)) * variables.cf(variables.zb(2)) * variables.cf(variables.zb(3)), variables.R * variables.Pb(2)) == variables.Cab(variables.g, variables.R * variables.Pb(2)) * variables.zb(12) * variables.zb(13) * variables.zb(23)
variables.bracketb(variables.cf(variables.zb(1)) * variables.cf(variables.zb(2)) * variables.cftilde(variables.zlb(3)), variables.R * variables.Pb(2)) == variables.zabp(-2, 3) * variables.bracketb(variables.cf(variables.zb(1)) * variables.cf(variables.zb(2)) * variables.cf(variables.zbp(3)), variables.R * variables.Pb(2)) == variables.Cab(variables.g, variables.R * variables.Pb(2)) * variables.zb(12) * (1 + variables.zb(1) * variables.zlb(3)) * (1 + variables.zb(2) * variables.zlb(3))

#201.19, operator calculations, 
variables.bracketb(variables.pointp(variables.ea(variables.i * variables.kb(4) * variables.Xf(math.inf, math.inf))) * variables.point(variables.ea(variables.i * variables.kb(1) * variables.Xf(variables.zb(1), variables.zlb(1)))) * variables.point(variables.ea(variables.i * variables.kb(2) * variables.Xf(variables.zb(2), variables.zlb(2)))) * variables.point(variables.ea(variables.i * variables.kb(3) * variables.Xf(0, 0))), variables.Sb(2)) == variables.diracbs(variables.dirac(0, variables.kb(4)), variables.Tf(variables.point(variables.ea(variables.i * variables.kb(1) * variables.xb(1))) * variables.point(variables.ea(variables.i * variables.kb(2) * variables.xb(2)))), variables.dirac(0, variables.kb(3)))
variables.point(variables.ea(variables.i * variables.k * variables.X)) == variables.ea(variables.i * variables.k * variables.xb(variables.C)) * variables.ea(variables.i * variables.k * variables.xb(variables.A))
variables.xab(variables.mu, variables.C, (variables.z, variables.zl)) == variables.xas(variables.mu) - variables.i * (variables.SlopeRegge/2)**(1/2) * summation(1/variables.m * (variables.alphaab(variables.mu, -variables.m) * variables.za(variables.m) + variables.alphaabtilde(variables.mu, -variables.m) * variables.zla(variables.m)), (variables.m, 1, math.inf))
variables.xab(variables.mu, variables.A, (variables.z, variables.zl)) == variables.i * variables.SlopeRegge/2 * variables.pa(variables.mu) * ln(Abs(variables.z)**2) + variables.i * (variables.SlopeRegge/2)**(1/2) * summation(1/variables.m * ((variables.alphaab(variables.mu, variables.m))/(variables.za(variables.m)) + (variables.alphaabtilde(variables.mu, variables.m))/(variables.zla(variables.m))), (variables.m, 1, math.inf))
variables.diracbs(variables.dirac(0, variables.kb(4)), variables.ea(variables.i * variables.kb(1) * variables.xb((1, variables.C))) * variables.ea(variables.i * variables.kb(1) * variables.xb((1, variables.A))) * variables.ea(variables.i * variables.kb(2) * variables.xb((2, variables.C))) * variables.ea(variables.i * variables.kb(2) * variables.xb((2 * variables.A))), variables.dirac(0, variables.kb(3))) # matric element
variables.ea(np.dot(variables.i * variables.kb(1), variables.xb((1, variables.A)))) * variables.ea(np.dot(variables.i * variables.kb(2), variables.xb((2, variables.C)))) == variables.ea(np.dot(variables.i * variables.kb(2), variables.xb((2, variables.C)))) * variables.ea(np.dot(variables.i * variables.kb(1), variables.xb((1, variables.A)))) * variables.ea(-variables.brackets(np.dot(variables.kb(1), variables.xb((1, variables.A))), np.dot(variables.kb(2), variables.xb((2 * variables.C))))) == variables.ea(np.dot(variables.i * variables.kb(2), variables.xb((2, variables.C)))) * variables.ea(np.dot(variables.i * variables.kb(1), variables.xb((1, variables.A)))) * Abs(variables.zb(12))**(variables.SlopeRegge * variables.kb(1) * variables.kb(2))
Abs(variables.zb(12))**(np.dot(variables.SlopeRegge * variables.kb(1), variables.kb(2))) * variables.diracbs(variables.dirac(0, variables.kb(4)), variables.ea(np.dot(variables.i * variables.kb(1), variables.xb((1, variables.C))) + np.dot(variables.i * variables.kb(2), variables.xb((2, variables.C)))) * variables.ea(np.dot(variables.i * variables.kb(1), variables.xb((1, variables.A))) + np.dot(variables.i * variables.kb(2), variables.xb((2, variables.A)))), variables.dirac(0, variables.kb(3))) == Abs(variables.zb(12))**(np.dot(variables.SlopeRegge * variables.kb(1), variables.kb(2))) * variables.diracbs(variables.dirac(0, variables.kb(4)), variables.ea(np.dot(variables.i * (variables.kb(1) + variables.kb(2)), variables.x)) * variables.ea(variables.SlopeRegge * (variables.kb(1) * ln(variables.zb(1)) + variables.kb(2) * ln(variables.zb(2)))), variables.dirac(0, variables.kb(3))) == Abs(variables.zb(12))**(np.dot(variables.SlopeRegge * variables.kb(1), variables.kb(2))) * Abs(variables.zb(1))**(np.dot(variables.SlopeRegge * variables.kb(1), variables.kb(3))) * Abs(variables.zb(2))**(np.dot(variables.SlopeRegge * variables.kb(2), variables.kb(3))) * variables.diracbs(variables.dirac(0, variables.kb(1) + variables.kb(2) + variables.kb(4)), variables.dirac(0, variables.kb(3))) == variables.i * variables.Cab(variables.X, variables.Sb(2)) * (2 * math.pi)**(variables.d) * variables.deltaa(variables.d) * summation(variables.kb(variables.i), variables.i) * Abs(variables.zb(12))**(np.dot(variables.SlopeRegge * variables.kb(1), variables.kb(2))) * Abs(variables.zb(1))**(np.dot(variables.SlopeRegge * variables.kb(1), variables.kb(2))) * Abs(variables.zb(2))**(np.dot(variables.SlopeRegge * variables.kb(2), variables.kb(3)))
#relationship between inner products
variables.diracbs(variables.dirac(0, variables.k), variables.dirac(0, variables.l)) == variables.i * variables.Cab(variables.X, variables.Sb(2)) * variables.holdera(2 * math.pi, variables.d) * variables.deltaa(variables.d) * variables.k + variables.l 
variables.bracket(variables.dirac(0, variables.k), variables.dirac(0, variables.l)) == variables.holdera(2 * math.pi, variables.d) * variables.deltaa(variables.d) * (variables.l - variables.k)
variables.diracbs(0, variables.k) == variables.i * variables.Cab(variables.X, variables.Sb(2)) * variables.diracb(0, -variables.k)
variables.AOl(variables.p) == variables.dagger(variables.AOp(variables.pp))
variables.OO(variables.z) == variables.ia(variables.h) * summation((variables.OOb(variables.n))/(variables.za(variables.n + variables.h)), (variables.n, -math.inf, math.inf))
variables.dagger(variables.OO(variables.z)) == variables.ia(-variables.h) * summation((variables.dagger(variables.OOb(variables.n)))/(variables.zla(variables.n + variables.h)), (variables.n, -math.inf, math.inf))
variables.OOl(variables.z) == variables.ia(-variables.h) * (-variables.za(-2))**(variables.h) * summation((variables.dagger(variables.OOb(variables.n)))/(variables.za(-variables.n-variables.h)), (variables.n, -math.inf, math.inf)) == variables.ia(variables.h) * summation((variables.dagger(variables.OOb(-variables.n)))/(variables.za(variables.n + variables.h)), (variables.n, -math.inf, math.inf))
variables.diracbs(variables.AObl(variables.i)) == variables.K * variables.diracb((variables.AOb(variables.i)))
variables.diracbs(variables.il, variables.j) == variables.diracbs(1, variables.AObl(variables.i, (math.inf, math.inf)), variables.j) == variables.K * variables.bracket(1, variables.AObl(variables.i, (math.inf, math.inf)), variables.j) == variables.K * variables.star(variables.bracket(1, variables.dagger(variables.AOb(variables.i, (0, 0))), 1)) == variables.K * variables.star(variables.bracket(variables.j, variables.i)) == variables.K * variables.bracket(variables.i, variables.j)
#cft on torus
2/variables.SlopeRegge * diff(diff(variables.Gpf(variables.semigroup(variables.w, variables.wl), variables.semigroup(variables.wp, variables.wpl)), variables.z), variables.zl) == -2 * math.pi * variables.deltaa(2) * (variables.w - variables.wp) + 1/(4 * math.pi * variables.taub(2))
variables.Gpf(variables.semigroup(variables.w, variables.wl), variables.semigroup(variables.zp, variables.wpl)) == -variables.SlopeRegge/2 * ln(variables.thetavarb(1, (variables.w - variables.wp)/(2 * math.pi), variables.tau))**2 + variables.SlopeRegge * variables.brackets(variables.Im * (variables.w - variables.wp))**2/(4 * math.pi * variables.taub(2)) + variables.kf(variables.tau, variables.taul)
variables.bracketb(np.prod(variables.point(variables.ea(np.dot(variables.i * variables.kb(variables.i), variables.Xf(variables.zb(variables.i), variables.zlb(variables.i))))), (variables.i, 1, variables.n)), variables.Ta(2)) == np.cross(variables.i * variables.Cab(variables.X, variables.Ta(2)) * (variables.tau) * (2 * math.pi)**variables.d * variables.deltaa(variables.d) * summation(variables.kb(variables.i), variables.i), np.prod(Abs((2 * math.pi)/(diff(variables.thetavarb(1, 0, variables.tau), variables.v)) * variables.thetavarb(1, (variables.wb(variables.i * variables.j))/(2 * math.pi), variables.tau) * math.exp(variables.brackets(-(variables.Im * variables.wb(variables.i * variables.j))**2/(4 * math.pi * variables.taub(2)))))**(np.dot(variables.SlopeRegge * variables.kb(variables.i), variables.kb(variables.j))), (variables.i < variables.j)))
variables.Zf(variables.tau) == variables.Trf(variables.brackets(math.exp(2 * math.pi * variables.i * variables.taub(1) * variables.P - 2 * math.pi * variables.taub(2) * variables.H))) == (variables.q * variables.ql)**(-variables.d/24) * variables.Trf(variables.qa(variables.Lb(0)) * variables.qa(variables.Lbtilde(0)))
summation(variables.qa(variables.n * variables.N), (variables.N, 0, math.inf)) == (1- variables.qa(variables.n)) **(-1)
variables.Zf(variables.tau) == variables.i * variables.Vb(variables.d) * variables.Zbf(variables.X, variables.tau**variables.d)
variables.Zbf(variables.X, variables.tau) == (4 * math.pi**2 * variables.SlopeRegge * variables.taub(2))**(-1/2) * Abs(variables.etaf(variables.tau))**(-2)
variables.etaf(variables.tau) == variables.qa(1/24) * np.prod(1- variables.qa(variables.n), (variables.n, 1, math.inf))
variables.ds2 == diff(variables.w) * diff(variables.wl) + variables.star(variables.epsilon) * diff(variables.wa(2)) + variables.epsilon * diff(variables.wla(2)) == (1 + variables.star(variables.epsilon) + variables.epsilon) * diff(variables.w + variables.epsilon * (variables.wl - variables.w)) * diff(variables.wl + variables.star(variables.epsilon) * (variables.w-variables.wl)) + variables.Of(variables.epsilona(2))
variables.wp == variables.wp + 2 * math.pi == variables.wp + 2 * math.pi * (variables.tau - 2 * variables.i * variables.taub(2) * variables.epsilon)
variables.delta * variables.tau == -2 * variables.i * variables.taub(2) * variables.epsilon
variables.delta * variables.Zf(variables.tau) == -1/(2 * math.pi) * integrate(variables.d2w) * variables.brackets(variables.delta * variables.gb(variables.wl * variables.wl) * variables.bracket(variables.Tbf(variables.w * variables.w, variables.w)) + variables.delta * variables.gb(variables.w * variables.w) * variables.bracket(variables.Tbf(variables.wl * variables.wl, variables.wl))) == -2 * math.pi * variables.i * variables.brackets(variables.delta * variables.tau * variables.bracket(variables.Tbf(variables.w * variables.w, 0)) - variables.delta * variables.taul * variables.bracket(variables.Tbf(variables.wl * variables.wl, 0)))
diff(variables.xa(variables.mu, variables.w), variables.w) * diff(variables.xb(variables.mu, 0), variables.w) == - (variables.SlopeRegge * variables.d)/(2 * variables.wa(2)) - variables.SlopeRegge * variables.Tbf(variables.w * variables.w, 0) + variables.Of(variables.w)
variables.Zf(variables.tau)**(-1) * variables.bracket(diff(variables.xa(variables.mu, variables.w), variables.w) * diff(variables.xb(variables.mu, 0), variables.w)) == variables.d * diff(diff(variables.Gpf(variables.semigroup(variables.w, variables.wl), variables.semigroup(variables.wp, variables.wpl)), variables.wp), variables.w) == (variables.SlopeRegge * variables.d)/2 * (variables.thetavarb(1) * diff(diff(variables.thetavarb(1), variables.w), variables.w) - diff(variables.thetavarb(1), variables.w) ( variables.thetavarb(1), variables.w))/(variables.thetavarab(2, 1)) + (variables.SlopeRegge * variables.d)/(8 * math.pi * variables.taub(2)) # evaluate from wp=0 to inf
variables.bracket(variables.Tbf(variables.w * variables.w, 0)) == (-variables.d/6 * (diff(variables.thetavarb(1, 0, variables.tau), variables.w, 3))/(diff(variables.thetavarb(1, 0, variables.tau), variables.w)) - variables.d/(8 * math.pi * variables.taub(2))) * variables.Zf(variables.tau)
diff(ln(variables.Zf(variables.tau)), variables.tau) == (math.pi * variables.i * variables.d)/3 * (diff(variables.thetavarb(1, 0, variables.tau), variables.w, 3))/(diff(variables.thetavarb(1, 0, variables.tau), variables.w)) + (variables.i * variables.d)/(4 * variables.taub(2))
diff(variables.thetavarb(1, (variables.w)/(2 * math.pi), variables.tau), variables.w, 2) == variables.i/(math.pi) * diff(variables.thetavarb(variables.w/(2 * math.pi), variables.tau), variables.tau)
diff(ln(variables.Zf(variables.tau)), variables.tau) == -variables.d/3 * diff(ln(diff(variables.thetavarb(1, 0, variables.tau), variables.w)), variables.tau) + (variables.i * variables.d)/(4 * variables.taub(2))
variables.Zf(variables.tau) == Abs(diff(variables.thetavarb(1, 0, variables.tau), variables.tau))**(-2 * variables.d/3) * variables.tauab((-variables.d/2), 2)
# bc cft 212.23
variables.Trf(math.exp(2 * math.pi * variables.i * variables.taub(1) * variables.P - 2 * math.pi * variables.taub(2) * variables.H)) == (variables.q * variables.ql)**(13/12) * variables.Trf(variables.qa(variables.Lb(0)) * variables.qla(variables.Lbtilde(0))) == 4 * (variables.q * variables.ql)**(1/12) * np.prod(Abs(1 + variables.qa(variables.n))**4, (variables.n, 1, math.inf))
variables.Zf(variables.tau) == variables.Trf(variables.brackets((-1)**(variables.F) * math.exp(2 * math.pi * variables.i * variables.taub(1) * variables.P - 2 * math.pi * variables.taub(2) * variables.H))) == 0
(variables.q * variables.ql)**(1/12) * np.prod(Abs(1 - variables.qa(variables.n)**4), (variables.n, 1, math.inf)) == Abs(variables.etaf(variables.tau))**4
variables.Zf(variables.tau) == summation(variables.qa(variables.hb(variables.i) - variables.c/24) * variables.qla(variables.hbtilde(variables.i) - variables.ctilde/24) * (-1)**(variables.Fb(variables.i)), variables.i)

#240.26 check approx bounds on DDF operators
variables.Vaf(variables.i, (variables.n * variables.kb(0), variables.z)) == diff(variables.Xaf(variables.i, variables.z), variables.z)* variables.ea(variables.i * variables.n * variables.kb(0) * variables.Xaf(variables.pos, variables.z)) * (2/variables.SlopeRegge)**(1/2)
#approx .28
variables.Aab(variables.i,variables.n) == integrate((diff(variables.z))/(2 * math.pi) * variables.Vaf(variables.i, (variables.n * variables.kb(0), variables.z)), 2 * math.pi) # contour integral
variables.brackets(variables.Aab(variables.i, variables.m), variables.Aab(variables.j, variables.n)) == variables.m * variables.deltaa(variables.i * variables.j) * variables.deltab(variables.m, -variables.n) * (variables.SlopeRegge * variables.kb(0) * variables.pa(variables.pos))/2

# 2.10.2
variables.dirac(1) == variables.dirac(0)
variables.psiab(variables.mu, -variables.r) == 1/(math.factorial(variables.r - 1/2)) * diff(variables.psiaf(variables.mu, 0), variables.z, variables.r-1/2)
variables.deltab(variables.eta) * variables.AO(variables.z, variables.zl) == - variables.epsilon * summation(1/(math.factorial(variables.n)) * np.cross(variables.brackets(diff(variables.etaf(variables.z), variables.z, variables.n) * variables.Gb(variables.n-1/2) + variables.star(diff(variables.etaf(variables.z), variables.z, variables.n)) * variables.Gbtilde(variables.n-1/2)), variables.AO(variables.z, variables.zl)), (variables.n, 0 , math.inf))
# check 2.11.5&6a
variables.ea(variables.i * variables.Hf(variables.z)) * variables.ea(variables.i * variables.Hf(0)) == variables.Of(variables.z)
variables.ea(-variables.i * variables.Hf(variables.z)) * variables.ea(-variables.i * variables.Hf(0)) == variables.Of(variables.z)
variables.bracketb(np.prod(variables.ea(variables.i * variables.epsilonb(variables.i) * variables.Hf(variables.zb(variables.i))), (variables.i)), variables.Sb(2)) == np.prod(variables.zab(variables.epsilonb(variables.i) * variables.epsilonb(variables.j), variables.i * variables.j), (variables.i<variables.j)) 
summation(variables.epsilonb(variables.i), variables.i) == 0
variables.psi == 2**(-1/2) * (variables.psia(1) + variables.i * variables.psia(2))
variables.psil == 2**(-1/2)  * (variables.psia(1) - variables.i * variables.psia(2))
# no 2.11.9a
variables.psif(variables.z) * variables.pisf(0) == variables.Of(variables.z)
variables.psilf(variables.z) * variables.psilf(0) == variables.Of(variables.z)
variables.psif(variables.z) == variables.ea(variables.i * variables.Hf(variables.z))
variables.psilf(variables.z) == variables.ea(-variables.i * variables.Hf(variables.z))
variables.psiftilde(variables.zl) == variables.ea(variables.i * variables.Hftilde(variables.z))
variables.psilftilde(variables.zl) == variables.ea(-variables.i * variables.Hftilde(variables.z))
variables.ea(variables.i * variables.Hf(variables.z)) * variables.ea(-variables.i * variables.Hf(-variables.z)) == 1/(2 * variables.z) + variables.i * diff(variables.Hf(0), variables.z) + 2 * variables.z * variables.Tabf(variables.H, variables.B, 0) + variables.Of(variables.za(2))
variables.psif(variables.z) * variables.psilf(-variables.z) == 1/(2 * variables.z) + variables.psi * variables.psilf(0) + 2 * variables.z * variables.Tabf(variables.psi, variables.B, 0) + variables.Of(variables.za(2))
variables.psi * variables.psil == variables.i * diff(variables.H, variables.z)
variables.Tab(variables.psi, variables.B) == variables.Tab(variables.H, variables.B)
variables.psif(variables.z) == variables.point(variables.ea(variables.i * variables.Hf(variables.z)))
variables.point(variables.ea(variables.i * variables.Hf(variables.z))) * variables.point(variables.ea(variables.i * variables.Hf(variables.zp))) == math.exp(variables.bracketc(-variables.brackets(variables.Hf(variables.z), variables.Hf(variables.zp)))) * variables.point(variables.ea(variables.i * variables.Hf(variables.zp))) * variables.point(variables.ea(variables.i * variables.Hf(variables.z))) == - variables.point(variables.ea(variables.i * variables.Hf(variables.zp))) * variables.point(variables.ea(variables.i * variables.Hf(variables.z)))
variables.Hf(variables.z) * variables.point(variables.ea(variables.i * variables.Hf(variables.zp))) == variables.point(variables.ea(variables.i * variables.Hf(variables.zp))) * (variables.Hf(variables.z) + variables.i * variables.brackets(variables.Hf(variables.z), variables.Hf(variables.zp))) == variables.point(variables.ea(variables.i * variables.Hf(variables.zp))) * (variables.Hf(variables.z) - math.pi * variables.sign(variables.sigmab(1), variables.sigmabp(1)))
variables.Trf(variables.qa(variables.Lb(0))) == (summation(variables.qa((variables.kab(2, variables.L))/2), (variables.kb(variables.L), variables.ZB))) * np.prod((1 - variables.qa(variables.n))**(-1), (variables.n, 1, math.inf))
variables.Trf(variables.qa(variables.Lb(0))) == np.prod((1 + variables.qa(variables.n-1/2))**(2), (variables.n, 1, math.inf))
variables.psif(variables.w + 2 * math.pi) == math.exp(2 * math.pi * variables.i * variables.v) * variables.psif(variables.w)
variables.psif(variables.z) == summation((variables.psib(variables.r))/(variables.za(variables.r + 1/2)), (variables.r, variables.ZB + variables.v)) 
variables.psilf(variables.z) == summation((variables.psibl(variables.s))/(variables.za(variables.s + 1/2)), (variables.s, variables.ZB - variables.v))
variables.bracketc(variables.psib(variables.r, variables.psibl(variables.s))) == variables.deltab(variables.r, -variables.s)
# check reference state condition 2.14.23
variables.psif(variables.z) * variables.AObf(variables.v, 0) == variables.Of(variables.za(-variables.v + 1/2))
variables.psilf(variables.z) * variables.AObf(variables.v, 0) == variables.Of(variables.za(variables.v-1/2))
math.exp(variables.brackets(variables.i * (-variables.v + 1/2) * variables.H)) == variables.AOb(variables.v)
variables.bracketb(0, variables.v+1) == variables.psibl(-variables.v) * variables.diracb(0, variables.v)
variables.dirac(variables.s) == variables.ea(variables.i * variables.s * variables.H) 
variables.s == variables.posneg * 1/2
2**(-1/2) * (variables.posneg * variables.psia(0) + variables.psia(1)) == variables.ea(variables.posneg * variables.i * variables.Ha(0))
2**(-1/2) * (variables.psia(2 * variables.a) + variables.posneg * variables.i * variables.psia(2 * variables.a + 1)) == variables.ea(variables.posneg * variables.i * variables.Ha(variables.a)) # check domain for a 2.14.28b
variables.Thetab(variables.sB) == math.exp(variables.brackets(variables.i * summation(variables.sb(variables.a) * variables.Ha(variables.a), variables.a)))
variables.Tab((variables.lamda), variables.B) == variables.Tab((1/2), variables.B) - (variables.lamda - 1/2) * diff(variables.b * variables.c, variables.z)
variables.Tab((variables.lamda), variables.B) == variables.Tab(variables.H, variables.B) - variables.i * (variables.lamda - 1/2) * diff(variables.H, variables.z, 2)
variables.b == variables.ea(variables.i * variables.H) 
variables.c == variables.ea(-variables.i * variables.H)
variables.c == 1-3* (2 * variables.lamda -1)**2 == 1 + 12 * variables.Va(2)
variables.H == variables.i * variables.rho 
variables.c == variables.ea(variables.rho) 
variables.b == variables.ea(-variables.rho)
