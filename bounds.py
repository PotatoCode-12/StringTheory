import math
from scipy.integrate import quad
import numpy as np 
from sympy import symbols, diff,integrate, Function, Derivative, summation, oo, ln, Abs
from sympy.diffgeom import WedgeProduct
import sympy as sp
import variables
import constants

diff_resulta = variables.xa(variables.mu, variables.tau, variables.l)
for _ in range(int(variables.sigma)):
    diff_resulta = diff(diff_resulta, variables.mu)
a = diff_resulta
diff_resultb = variables.xa(variables.mu, variables.tau, 0)
for _ in range(int(variables.sigma)):
    diff_resultb = diff(diff_resultb, variables.mu)
b = diff_resultb
-math.inf < variables.tau < math.inf
0 < variables.sigma < variables.l
variables.bnds*diff(variables.xb, variables.a) == 0 #def for diff(M), pg 14
variables.xa(variables.mu, variables.tau, variables.l) == variables.xa(variables.mu, variables.tau, 0) #cont
a == b #cont
variables.WorldLineMetric(variables.tau, variables.l) == variables.WorldLineMetric(variables.tau, 0) #periodic field restriction for string propagation w/o sources
variables.cords(variables.pb(variables.i), variables.xa(variables.j)) == -variables.i*variables.NambuGotoVar(variables.i, variables.j)
variables.cords(variables.pb(variables.neg), variables.xa(variables.neg)) == -variables.i

#boundry normalization 62.25
variables.alphaab(variables.mu, variables.D0) == (2*variables.SlopeRegge)**(1/2) * variables.pa(variables.mu)
variables.xa(variables.mu, variables.z, variables.zl) == variables.xas(variables.mu) - variables.i*variables.SlopeRegge*variables.pa(variables.mu) * ln(variables.z) **variables.D2 + variables.i(variables.SlopeRegge/2)**(1/2) * (summation((variables.alphaab(variables.mu, variables.m))/variables.m * (variables.za(-variables.m)) + variables.zla(-variables.m), (variables.m, -math.inf, -1)) + summation((variables.alphaab(variables.mu, variables.m))/variables.m * (variables.za(-variables.m)) + variables.zla(-variables.m), (variables.m, 1, math.inf)))
variables.Lb(variables.D0) == variables.SlopeRegge*variables.pa(variables.D2) + summation(variables.alphaab(variables.mu, -variables.n) * variables.alphab(variables.mu, variables.n), (variables.n, 1, math.inf))
variables.brackets(variables.alphaab(variables.mu, variables.m), variables.alphaab(variables.v, variables.n)) == variables.m*variables.deltab(variables.m, -variables.n) * variables.etaa(variables.mu, variables.v) #communicators
variables.brackets(variables.xas(variables.mu), variables.pa(variables.v)) == variables.i*variables.etaa(variables.mu, variables.v)
variables.cf(variables.z) == variables.cftilde(variables.zl)
variables.bf(variables.z) == variables.bftilde(variables.zl)
variables.i*variables.z == 0 # all imaginary permutations
variables.cf(variables.z) == variables.cftilde(variables.zpl)
variables.bf(variables.z) == variables.bftilde(variables.zpl)
variables.i*variables.z <= 0 # double check 63.30, sure wrong
variables.zp == variables.zl

variables.xa(variables.mu, 'upper', variables.sigma) == variables.xa(variables.mu, 'lower', math.pi-variables.sigma)

#95.24
variables.Tab(variables.a, variables.a) == variables.ab(variables.D1) * variables.R + variables.ab(variables.D2)
variables.Sb(variables.c*variables.t) == variables.b * integrate(variables.d2sigma * variables.g**(1/2))
variables.Tab(variables.a, variables.a) == variables.ab(variables.D1) * variables.R + (variables.ab(variables.D2) + 4 * math.pi * variables.b)
variables.Sb(variables.c*variables.t) * variables.integrate (variables.d2sigma * variables.g**(1/2) * variables.bb(variables.D1), variables.M) + integrate (diff(variables.s) * (variables.bb(variables.D2) + variables.bb(variables.D3) * variables.k), 0, 2*math.pi) # contour, 96. 28
variables.deltab(variables.W) * variables.Sb(variables.c * variables.t) == 2 * integrate(variables.d2sigma * variables.g**(1/2) * variables.bb(variables.D1) * variables.delta * variables.omegas, variables.M) + variables.integrate (diff(variables.s) * (variables.bb(variables.D2) + variables.bb(variables.D3) * variables.na (variables.a) * diff(1, variables.a)) * variables.delta * variables.omegas, 0, math.pi * 2) # contour
variables.Tabf(variables.a, variables.a, variables.sigma) == -variables.GO(variables.sigma)/12 * variables.Rf(variables.sigma) # lorentz invariance

#291.1 mapping
variables.gamma
(variables.zp - variables.U)/(variables.zp - variables.V) == variables.K * (variables.z - variables.U)/(variables.z - variables.V)
Abs(variables.K)**(-1/2) < Abs((variables.z - variables.U)/(variables.z - variables.V)) < Abs(variables.K)**(1/2) 
diff(variables.s**2) == (diff(variables.z) * diff(variables.zl))/(variables.Im * variables.z)**2
# check defined path 292.4 
# period matrix 293.5
diff(variables.omegabs(variables.z), variables.zl) == 0
variables.omegabs(variables.zl) == 0
variables.dim * variables.ker * variables.Pab(variables.T, 0) - variables.dim * variables.ker * variables.Pb(0) == -variables.chi == 2 * variables.g - 2
integrate(diff(variables.z) * variables.omegabs(variables.z * variables.j), variables.Ab(variables.i)) == variables.deltab(variables.i * variables.j) # contour
variables.taub(variables.i * variables.j) == integrate(diff(variables.z) * variables.omegabs(variables.z * variables.j), variables.Bb(variables.i))
variables.ya(2) == np.prod(variables.z - variables.Zb(variables.i), (variables.i, 1, 2 * variables.k))
#9.3 cont
variables.zb(1) * variables.zb(2) == variables.q 
(1 - variables.epsilon)**(-1) * Abs(variables.q)**(1/2) > Abs(variables.zb(1)) > (1 - variables.epsilon)* Abs(variables.q)**(1/2)
variables.M == variables.Mbf(1, math.inf) * variables.Mbf(2, (variables.zb(1), variables.zb(2), variables.q))
ln(Abs(variables.q)) < variables.Im * variables.w < 0 
variables.w == variables.w + 2 * math.pi 
variables.za((0)) == (variables.a * variables.z)/(variables.z - 2)
variables.za((1)) == (variables.a * (variables.z - 1))/(variables.z + 1)
variables.za((math.inf)) == variables.a/(2 * variables.z - 1) 
variables.zb(4) == ((variables.q - variables.aa(2))**2)/(4 * variables.aa(2) * variables.q)
variables.r == 3 * variables.g + variables.n - 2 
#9.4 cont, check .1 &2 for unincluded states and terms
variables.Gan(variables.i * variables.j) * variables.gOb(variables.j * variables.k) == variables.deltaab(variables.i, variables.k)
variables.Gan(variables.i * variables.j) == variables.gOa(variables.i * variables.j) 
variables.AOab((variables.zbp(2)), variables.j) == variables.qa(variables.hb(variables.j)) * variables.qla(variables.hbtilde(variables.j)) * variables.AOab((variables.zb(2)), variables.j)
# removed 302.6 & 7
variables.zbp(1) == variables.zb(1) + summation(variables.epsilonb(variables.n) * variables.zab(variables.n + 1, 1), (variables.n, -math.inf, math.inf))
variables.zbp(2) == variables.zb(2) - summation(variables.epsilonb(-variables.n) * variables.qa(-variables.n) * variables.zab(variables.n + 1, 2), (variables.n, -math.inf, math.inf))
# removed 303.9
2 * math.pi * (variables.m, variables.n)
(variables.sigmaa(1), variables.sigmaa(2)) == (variables.sigmaap(1), variables.sigmaap(2)) * np.array([[variables.m, variables.n], 
                                                                                                       variables.q, variables.p])
(variables.sigmaap(1), variables.sigmaap(2)) == 2 * math.pi * (1, 0)

#2.164.1
1/2 * variables.bracketc(np.array([[variables.Qb(variables.alpha)],
                                   [variables.Qbtilde(variables.alpha)]]), np.array([[variables.dagger(variables.Qb(variables.Beta)), variables.dagger(variables.Qbtilde(variables.Beta))]])) == variables.M * variables.deltab(variables.alpha * variables.Beta) * np.array([[1, 0],
                                                                                                                                                                                                                                                                             [0, 1]]) + (variables.Lb(1))/(2 * math.pi * variables.SlopeRegge) * variables.holderb(variables.Rhoa(0) * variables.Rhoa(1), variables.alpha * variables.Beta) * np.array([[variables.p, variables.q/variables.g],
                                                                                                                                                                                                                                                                                                                                                                                                                                                        [variables.q/variables.g, -variables.p]])
variables.M + variables.posneg * variables.Lb(1) * ((variables.pa(2) + variables.qa(2)/variables.ga(2))**(1/2))/(2 * math.pi * variables.SlopeRegge) 
variables.M/(variables.Lb(1)) >= variables.Lb(1) * ((variables.pa(2) + variables.qa(2)/variables.ga(2))**(1/2))/(2 * math.pi * variables.SlopeRegge)
variables.taub((0, 1)) + variables.taub((1, 0)) == (variables.ga(-1) + 1)/(2 * math.pi * variables.SlopeRegge)
variables.taub((1, 1)) == ((variables.qa(-2) + 1)**(1/2))/(2 * math.pi * variables.SlopeRegge)
variables.Rhoa(0) * variables.Rhoa(1) * variables.Q == variables.Q # left moving
variables.Rhoa(0) * variables.Rhoa(1) * variables.Qtilde == -variables.Qtilde # right moving
variables.SBb(1) == -variables.Tb(1) * integrate(variables.d2(variables.zeta) * variables.ea(-variables.Phi) * variables.brackets(-sp.Determinant(variables.Gb(variables.a * variables.b) + variables.Bb(variables.a * variables.b) + 2 * math.pi * variables.SlopeRegge * variables.Fb(variables.a * variables.b)))**(1/2))
variables.taub((1, 0)) + variables.taub((0, 1)) - variables.taub((1, 1)) == (1 - variables.Of(variables.g))/(2 * math.pi * variables.SlopeRegge)
integrate(diff(variables.s) * variables.delta * variables.xa(variables.mu) * (1/(2 * math.pi * variables.SlopeRegge) * diff(variables.xb(variables.mu), variables.n) + variables.i * variables.Fb(variables.mu * variables.v) * diff(variables.xa(variables.v), variables.t)), diff(variables.M, variables.z)) # contour integral
diff(variables.xb(variables.mu), variables.n) + 2 * math.pi * variables.SlopeRegge * variables.i * variables.Fb(variables.mu * variables.v) * diff(variables.xa(variables.v), variables.t) == 0
variables.taub((1, 1)) + variables.taub((0, 1)) == ((variables.ga(-2) + 1)**(1/2) + variables.ga(-1))/(2 * math.pi * variables.SlopeRegge) == (2 * variables.ga(-1) + variables.g/2 + variables.Of(variables.ga(3)))/(2 * math.pi * variables.SlopeRegge)
variables.taub((1, 2)) == ((4 * variables.ga(-2) + 1)**(1/2))/(2 * math.pi * variables.SlopeRegge) == (2 * variables.ga(-1) + variables.g/4 + variables.Of(variables.ga(3)))/(2 * math.pi * variables.SlopeRegge)
# check proprtional matric 2.157.13
variables.Wf(variables.Phi) == variables.Trf(variables.Phib(1) * variables.brackets(variables.Phib(2), variables.Phib(3))) + variables.m * variables.Trf(variables.Phib(variables.i) * variables.Phib(variables.i))
1/2 * variables.bracketc(variables.brackets(np.array([[variables.Qb(variables.alpha)],
                                                      [variables.Qbtilde(variables.alpha)]])), variables.brackets(np.array([[variables.dagger(variables.Qb(variables.Beta)), variables.dagger(variables.Qbtilde(variables.Beta))]]))) == variables.M * np.array([[1, 0],
                                                                                                                                                                                                                                                                [0, 1]]) * variables.deltab(variables.alpha * variables.Beta) + np.array([[0, variables.Zb(variables.alpha * variables.gamma)],
                                                                                                                                                                                                                                                                                                                                          [-variables.dagger(variables.Zb(variables.alpha * variables.gamma)), 0]]) * variables.Rhoab(0, variables.gamma * variables.Beta)
variables.Z == variables.taub(0) + variables.taub(variables.p) * variables.Vb(*variables.p) * variables.Beta
variables.Ma(2) * np.array([[1, 0],
                            [0, 1]]) >= np.array([[0, variables.Z],
                                                  [-variables.dagger(variables.Z)]]) * variables.Rhoa(0) * np.array([[0, variables.Z],
                                                                                                                     [-variables.dagger(variables.Z), 0]]) * variables.Rhoa(0) == np.array([[variables.Z * variables.dagger(variables.Z), 0],
                                                                                                                                                                                           [0, variables.dagger(variables.Z) * variables.Z]])
variables.Z * variables.dagger(variables.Z) == variables.tauab(2, 0) + variables.taub(0) * variables.taub(variables.p) * variables.Vb(variables.p) * (variables.Beta + variables.Betap) + variables.tauab(2, variables.p) * variables.Vab(2, variables.p) * variables.Beta * variables.dagger(variables.Beta)
variables.M >= variables.taub(0) + variables.taub(variables.p) * variables.Vb(variables.p)
variables.M >= (variables.tauab(2, 0) + variables.tauab(2, variables.p) * variables.Vab(2, variables.p))**(1/2)
variables.taub(0) + (variables.taub(0) + (variables.pab(2, 9))/(2 * variables.taub(0)))
2 * variables.taub(0) + (variables.pab(2, 9))/(4 * variables.taub(0))
integrate(variables.Fb(2), variables.D2) == 2 * math.pi 
variables.i * variables.mub(2) * integrate((variables.Cb(3) + WedgeProduct(2 * math.pi * variables.SlopeRegge * variables.Fb(2), variables.Cb(1))))
summation(variables.qa(variables.n) * variables.Db(variables.n), (variables.n, 0, math.inf)) == 2**8 * np.prod(((1 + variables.qa(variables.k))/(1 - variables.qa(variables.k)))**8, (variables.k, 1, math.inf))
(variables.gab(2, variables.D0))/4 * summation((variables.dagger(variables.chib(variables.i)) * variables.sigmaab(variables.A, variables.i * variables.j) * variables.chib(variables.j))**2, (variables.A, 1, 3)) + summation(((variables.xb(variables.i) - variables.xbp(variables.i))**2)/((2 * math.pi * variables.SlopeRegge)**2) * variables.dagger(variables.chi) * variables.chi, (variables.i, 1, 5))
variables.Xb(variables.i) - variables.xbp(variables.i) == 0
variables.chi != 0
variables.Xb(variables.i) - variables.xbp(variables.i) != 0 
variables.chi == 0
variables.Da(variables.A) == variables.dagger(variables.chib(variables.i)) * variables.sigmaab(variables.A, variables.i * variables.j) * variables.chib(variables.j) == 0
variables.dagger(variables.chib(variables.i * variables.a)) * variables.sigmaab(variables.A, variables.i * variables.j) * variables.chib(variables.j * variables.a) == 0
variables.chib(variables.i * variables.a) == variables.v * variables.deltab(variables.i * variables.a)
variables.chib(variables.i * variables.a) == variables.v * variables.Ub(variables.i * variables.a)
variables.starf(variables.Fb(2)) == variables.posneg * variables.Fb(2)
variables.delta * variables.lamda == variables.Fb(variables.M * variables.N) * variables.Rhoa(variables.M * variables.N) * variables.zeta 
(variables.B4, variables.B2, variables.B1) + (variables.B4p, variables.B1, variables.B2)
1/2 * (2 * math.pi * variables.SlopeRegge)**2 * variables.mub(4) * integrate(WedgeProduct(variables.Cb(1), variables.Trf(WedgeProduct(variables.Fb(2), variables.Fb(2)))))
integrate(variables.Trf(WedgeProduct(variables.Fb(2), variables.Fb(2))), variables.D4) == 8 * math.pi**2 
