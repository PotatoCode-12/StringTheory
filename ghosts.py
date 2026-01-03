import math
from scipy.integrate import quad
import numpy as np 
from sympy import symbols, diff,integrate, Function, Derivative, summation, oo
import sympy as sp
import variables
import constants
import states

def noGhostTheorem():
    variables.lb(variables.m) == variables.Lab(variables.X, variables.m) + variables.Lab(variables.K, variables.m) + variables.Lab(variables.g, variables.m) #137.1
    variables.dirac((variables.N, variables.I), variables.k), variables.dirac((variables.N, variables.Ntilde, variables.I), variables.k) #general state
    -summation(variables.kb(variables.mu) * variables.ka(variables.mu), (variables.mu, 0, variables.d-1)) == variables.m**2
    variables.SlopeRegge * variables.m**2 == summation(variables.n * (variables.Nb(variables.b*variables.n) + variables.Nb(variables.c * variables.n) + summation(variables.Nb(variables.mu * variables.n), (variables.mu, 0, variables.d-1))) + variables.Lab(variables.K, 0) - 1, (variables.n, 1, math.inf))
    -summation(variables.kb(variables.mu) * variables.ka(variables.mu), (variables.mu, 0, variables.d-1)) == variables.m**2 == variables.mtilde**2
    variables.SlopeRegge/4 * variables.m**2 == summation(variables.n * (variables.Nb(variables.b*variables.n) + variables.Nb(variables.c * variables.n) + summation(variables.Nb(variables.mu * variables.n), (variables.mu, 0, variables.d-1))) + variables.Lab(variables.K, 0) - 1, (variables.n, 1, math.inf))
    variables.SlopeRegge/4 * variables.mtilde**2 == summation(variables.n * (variables.Nbtilde(variables.b*variables.n) + variables.Nbtilde(variables.c * variables.n) + summation(variables.Nbtilde(variables.mu * variables.n), (variables.mu, 0, variables.d-1))) + variables.Labtilde(variables.K, 0) - 1, (variables.n, 1, math.inf))
    variables.bracket(variables.dirac((0, variables.I), variables.kb), variables.dirac((0, variables.Ip), variables.kbp)) == variables.bracket(variables.dirac((0, 0, variables.I), variables.kb), variables.dirac((0, 0, variables.Ip), variables.kbp)) == 2*variables.ka(0) * (2*math.pi)**(variables.d-1) * variables.deltaa(variables.d-1) * (variables.kb- variables.kbp) * variables.deltab(variables.I*variables.Jp)
    variables.alphaab(variables.posneg, variables.m) == 2**(-1/2) * (variables.alphaab(0, variables.m)  + variables.posneg*variables.alphaab(1, variables.m))
    variables.brackets(variables.alphaab(variables.pos, variables.m), variables.alphaab(variables.neg, variables.n)) == -variables.m * variables.deltab(variables.m, -variables.n)
    variables.brackets(variables.alphaab(variables.pos, variables.m), variables.alphaab(variables.pos, variables.n)) == variables.brackets(variables.alphaab(variables.neg, variables.m), variables.alphaab(variables.neg, variables.n)) == 0
    variables.Na(variables.l*variables.c) == summation(1/variables.m * variables.cords(variables.alphaab(variables.pos, -variables.m) * variables.alphaab(variables.neg, variables.m)))
    variables.Qb(variables.B) == variables.Qb(1) + variables.Qb(0) + variables.Qb(-1)
    variables.brackets(variables.Na(variables.l * variables.c), variables.Qb(variables.j)) == variables.j * variables.Qb(variables.j)
    variables.brackets(variables.Na(variables.g), variables.Qb(variables.j)) == variables.Qb(variables.j)
    (variables.Qab(2, 1)) + (variables.bracketc(variables.Qb(1), variables.Qb(0))) + (variables.bracketc(variables.Qb(1), variables.Qb(-1)) + variables.Qab(2, 0)) + (variables.bracketc(variables.Qb(0), variables.Qb(-1))) + (variables.Qab(2, -1)) == 0
    variables.Qb(1) == -(2*variables.SlopeRegge)**(1/2) * variables.ka(variables.pos) * (summation(variables.alphaab(variables.neg, -variables.m) * variables.cb(variables.m), (variables.m, -math.inf, -1)) + summation(variables.alphaab(variables.neg, -variables.m) * variables.cb(variables.m), (variables.m, 1, math.inf)))
    variables.R == 1/((2*variables.SlopeRegge)**(1/2) * variables.ka(variables.pos)) * (summation(variables.alphaab(variables.pos, -variables.m) * variables.bb(variables.m), (variables.m, -math.inf, -1)) + summation(variables.alphaab(variables.pos, -variables.m) * variables.bb(variables.m), (variables.m, 1, math.inf)))
    variables.S == variables.bracketc(variables.Qb(1), variables.R) == summation(variables.m*variables.bb(-variables.m) * variables.cb(variables.m) + variables.m * variables.cb(-variables.m) * variables.bb(variables.m) - variables.alphaab(variables.pos, -variables.m) * variables.alphaab(variables.neg, variables.m) - variables.alphaab(variables.neg, -variables.m) * variables.alphaab(variables.pos, variables.m)) == summation(variables.m *(variables.Nb(variables.b*variables.m) + variables.Nb(variables.c * variables.m) + variables.Nab(variables.pos, variables.m) + variables.Nab(variables.neg, variables.m)), (variables.m, 1, math.inf))
    variables.dirac(variables.psi) == 1/variables.s * variables.bracketc(variables.Qb(1), variables.R) * variables.dirac(variables.psi) == 1/variables.s * variables.Qb(1) * variables.R * variables.dirac(variables.psi)
    0 == variables.Qb(1) * variables.S * variables.dirac(variables.psi) == variables.S * variables.Qb(1) * variables.dirac(variables.psi)
    variables.S + variables.U == variables.bracketc(variables.Qb(variables.B), variables.R)
    # check 140.19
    variables.bracket(variables.psi, variables.psip) == variables.bracketc(variables.psib(0), variables.psibp(0))
    variables.hilbertb(variables.O * variables.C * variables.Q) == variables.hilbertb(variables.B * variables.R * variables.S * variables.T) == variables.hilbertb(variables.lightcone)
    variables.Qb(variables.B) * variables.dirac(variables.psi, variables.down) == summation(variables.cb(-variables.n) * (variables.Lab(variables.m, variables.n) - variables.deltab(variables.n, 0)) * variables.dirac(variables.psi, variables.down)) == 0
    variables.dirac(variables.psi, variables.down)  - variables.dirac(variables.psip, variables.down, (variables.n, 0, math.inf)) == variables.Qb(variables.B) * variables.dirac(variables.chi)
    variables.dirac(variables.psi, variables.down) - variables.dirac(variables.psip, variables.down) == summation(variables.cb(variables.m) * variables.Lab(variables.m, -variables.n) * variables.bb(-variables.n) * variables.dirac(variables.chib(variables.n), variables.down), ((variables.m, variables.n), 1, math.inf)) == summation(variables.Lab(variables.m, -variables.n) * variables.dirac(variables.chib(variables.n), variables.down), (variables.n, 1, math.inf))

# fields 157.10
variables.caf(variables.a, variables.sigma) == summation(variables.cb(variables.J) * variables.Cabnf(variables.a, variables.J, variables.sigma), variables.J)
variables.bbf(variables.a * variables.b, variables.sigma) == summation(variables.bb(variables.K) * variables.Bbnf(variables.K * variables.a * variables.b, variables.sigma), variables.K)
variables.Sb(variables.g) == 1/(2 * math.pi) * (variables.b, variables.P1 * variables.c) == 1/(2 * math.pi) * (variables.P1a(variables.T) * variables.b, variables.c)
variables.P1a(variables.T) * variables.P1 * variables.Cabn(variables.a, variables.J) == variables.vabp(2, variables.J) * variables.Cabn 
variables.P1 * variables.P1a(variables.T) * variables.Bbn(variables.K * variables.a * variables.b) == variables.vab(2, variables.K) * variables.Bbn(variables.K * variables.b *variables.a)
(variables.Cbn(variables.J), variables.Cbn(variables.Jp)) == integrate(variables.d2sigma * variables.g**(1/2) * variables.Cabn(variables.a, variables.J) * variables.Cbn(variables.Jp * variables.a)) == variables.deltab(variables.J * variables.Jp)
(variables.Bbn(variables.K), variables.Bbn(variables.Kp)) == integrate(variables.d2sigma * variables.g**(1/2) * variables.Bbn(variables.K * variables.b *variables.a) * variables.Babn(variables.a * variables.b, variables.Kp)) == variables.deltab(variables.K * variables.Kp)
(variables.P1 * variables.P1a(variables.T)) * variables.P1 * variables.Cbn(variables.J) == variables.P1 *(variables.P1a(variables.T) * variables.P1) * variables.Cbn(variables.J) == variables.vabp(2, variables.J) * variables.P1 * variables.Cbn(variables.J)
variables.Bbn(variables.J * variables.a * variables.b) == 1/(variables.vb(variables.J)) *  variables.holderb(variables.P1 * variables.Cbn(variables.J), variables.a *variables.b)
variables.vb(variables.J) == variables.vbp(variables.J) != 0
# check Delta FP 158.16
variables.Deltab(variables.FP)  == integrate(np.prod(diff(variables.bb(0 * variables.k)), (variables.k, 1, variables.mu))) * np.prod(summation(variables.bb((0*variables.kpp)/(4 * math.pi)) * (variables.Bbn(0 * variables.kpp), diff(variables.hat(variables.g), variables.kp)), (variables.kpp, 1, variables.mu)), (variables.kp, 1, variables.mu)) * integrate(np.prod(variables.d * variables.cb(0, variables.j), (variables.j, 1, variables.keta))) * np.prod(summation(variables.cb(0, variables.jp) * variables.Cabnf(variables.a, (0, variables.jp), variables.sigmab(variables.i)), (variables.jp, 1, variables.keta)), ((variables.a, variables.i), variables.f)) * integrate(np.prod(variables.d * variables.bb(variables.J) * variables.d * variables.cb(variables.J) * math.exp((-variables.vb(variables.J) * variables.bb(variables.J) * variables.cb(variables.J))/(2 * math.pi)), variables.J))
variables.Deltab(variables.FP) == sp.Determinant((variables.Bbn(0, variables.k), diff(variables.hat(variables.g), variables.kp))/(4 * math.pi)) * sp.Determinant(variables.Cabnf(variables.a, (0, variables.j), variables.sigmab(variables.i))) * variables.detp(((variables.P1a(variables.T) * variables.P1)/(4 * math.pi**2))**(1/2))
#158.19
variables.Deltab(variables.a) * variables.ja(variables.a) == (1 - 2* variables.lamda)/4 * variables.R
(variables.delta * (variables.brackets(diff(variables.phi))) * math.exp(-variables.S))/(variables.brackets(diff(variables.phi)) * math.exp(-variables.S)) == (variables.i * variables.epsilon)/(2*math.pi) * integrate(variables.d2sigma * variables.g**(1/2) * variables.Deltab(variables.a) * variables.ja(variables.a)) # check arrow for transformation 159.20
1/2 * (variables.keta - variables.mu) == 1/2 * (variables.dim * variables.ker * variables.P1 - variables.dim * variables.ker * variables.P1a(variables.T))
variables.dim * variables.ker * variables.Pb(variables.n) - variables.dim * variables.ker * variables.Pab(variables.T, variables.n) == (2 * variables.n + 1) * variables.chi

#2.16.1
variables.Betaf(variables.z) == summation((variables.Betab(variables.r))/(variables.za(variables.r + 3/2)), (variables.r, variables.ZB + variables.v))
variables.gammaf(variables.z) == summation((variables.gammab(variables.r))/(variables.za(variables.r - 1/2)), (variables.r, variables.ZB + variables.v))
variables.bf(variables.z) == summation((variables.bb(variables.m))/(variables.za(variables.m + 2)) (variables.m, -math.inf, math.inf))
variables.cf(variables.z) == summation((variables.cb(variables.m))/(variables.za(variables.m-1)), (variables.m, -math.inf, math.inf))
variables.brackets(variables.gammab(variables.r), variables.Betab(variables.s)) == variables.deltab(variables.r, -variables.s)
variables.bracketc(variables.bb(variables.m), variables.cb(variables.n)) == variables.deltab(variables.n, -variables.m)
variables.Betab(variables.r) * variables.diracb(0, variables.NS) == 0 
variables.r >= 1/2 
variables.gammab(variables.r) * variables.diracb(0, variables.NS) == 0

variables.Betab(variables.r) * variables.diracb(0, variables.R) == 0 
variables.r >= 0 
variables.gammab(variables.r) * variables.diracb(0, variables.R) == 0
variables.r >= 1 
variables.bb(variables.m) * variables.diracb(0, (variables.NS, variables.R)) == 0
variables.m >= 0
variables.cb(variables.m) * variables.diracb(0, (variables.NS, variables.R)) == 0
variables.m >= 1
variables.Lab(variables.g, variables.m) == summation((variables.m + variables.n) * variables.point(variables.bb(variables.m-variables.n) * variables.cb(variables.n)), (variables.n, variables.ZB)) + summation(1/2 * (variables.m + 2 * variables.r) * variables.point(variables.Betab(variables.m-variables.r) * variables.gammab(variables.r)) + variables.aa(variables.g) * variables.deltab(variables.m, 0), (variables.r, variables.ZB + variables.v))
variables.Gab(variables.g, variables.r) == - summation(variables.brackets(1/2 * (2 * variables.r + variables.n) * variables.Betab(variables.r-variables.n) * variables.cb(variables.n) + 2 * variables.bb(variables.n) * variables.gammab(variables.r-variables.n)), (variables.n, variables.ZB))
# R:
variables.aa(variables.g) == -5/8
# NS:
variables.aa(variables.g) == -1/2 
variables.Betab(variables.r) * variables.dirac(1) == 0 
variables.r >= -1/2 
variables.gammab(variables.r) * variables.dirac(1) == 0
variables.r >= 3/2 
variables.gammaf(variables.z) * variables.deltaf(variables.gammaf(0)) == variables.Of(variables.z)
variables.Betaf(variables.z) * variables.deltaf(variables.gammaf(0)) == variables.Of(variables.za(-1))
# voided 2.17.8-11
variables.Betaf(variables.z) * variables.betaf(0) == variables.Of(variables.za(0))
variables.Betaf(variables.z) * variables.gammaf(0) == variables.Of(variables.za(-1))
variables.gammaf(variables.z) * variables.gammaf(0) == variables.Of(variables.za(0))
variables.Betaf(variables.z) == variables.ea(-variables.phif(variables.z)) * diff(variables.zetaf(variables.z), variables.z)
variables.gamma == variables.ea(variables.phif(variables.z)) * variables.etaf(variables.z)
variables.etaf(variables.z) * variables.etaf(0) == variables.Of(variables.z) 
diff(variables.zetaf(variables.z), variables.z) * diff(variables.zetaf(0), variables.z) == variables.Of(variables.z)
# check OPE 2.18.15
variables.Tab(variables.phi, variables.B) == -1/2 * diff(variables.phi, variables.z) * diff(variables.phi, variables.z) + 1/2 * (1 - 2 * variables.lamdap) * diff(variables.phi, variables.z, 2)
variables.Tab(variables.eta * variables.zeta, variables.B) == - variables.eta * diff(variables.zeta, variables.z)
variables.Tab(variables.Beta * variables.gamma, variables.B) == variables.Tab(variables.phi, variables.B) + variables.Tab(variables.eta * variables.zeta, variables.B)
variables.eta == variables.ea(-variables.chi) 
variables.zeta == variables.ea(variables.chi)
variables.Beta == variables.ea(-variables.phi + variables.chi) * diff(variables.chi, variables.z)
variables.gamma == variables.ea(variables.phi - variables.chi)
variables.Tb(variables.B) == -1/2 * diff(variables.phi, variables.z) * diff(variables.phi, variables.z) +1/2 * diff(variables.chi, variables.chi) * diff(variables.chi, variables.z) + 1/2 * (1 - variables.lamdap) * diff(variables.phi, variables.z, 2) + 1/2 * diff(variables.chi, variables.z, 2)
variables.delta * (variables.gamma) == variables.ea(-variables.phi)
variables.h == 1/2 
variables.ea(-variables.phi) 
variables.ea(-variables.phi) * variables.ea(variables.posneg * variables.i * variables.Ha(variables.a))
variables.Betaf(variables.z) * variables.Sigmaf(0) == variables.Of(variables.za(-1/2))
variables.gammaf(variables.z) * variables.Sigma(0) == variables.Of(variables.za(1/2)) 
variables.Sigma == variables.ea(-variables.phi/2)
variables.h == 3/8 
variables.VOb(variables.sB) * variables.ea(-variables.phi/2) * variables.Thetab(variables.sB)
math.exp(np.cross(variables.i * variables.kb(variables.L), variables.Hb(variables.L)) + np.cross(variables.i * variables.kb(variables.R), variables.Hb(variables.R)))
variables.Cbf(variables.k, variables.alphab(0)) * variables.point(math.exp(np.cross(variables.i * variables.kb(variables.L), variables.Hb(variables.L)) + np.cross(variables.i * variables.kb(variables.R), variables.Hb(variables.R))))
variables.Cbf(variables.k, variables.alphab(0)) == math.exp(math.pi * variables.i * summation(variables.nb(variables.alpha) * variables.alphab(0 * variables.Beta) * variables.kbf(variables.alpha, variables.kb(variables.Beta)), (variables.alpha > variables.Beta))) 
