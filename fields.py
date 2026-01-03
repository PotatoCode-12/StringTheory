import math
from scipy.integrate import quad
import numpy as np 
from sympy import symbols, diff,integrate, Function, Derivative, summation, oo, ln, Abs
import variables
import constants

1/(math.pi*variables.SlopeRegge)*diff(diff(variables.xa(variables.mu, variables.z, variables.zl) * variables.xa(variables.v, variables.zp, variables.zpl), variables.zl), variables.z) == -variables.etaa(variables.mu, variables.v)* variables.delta(variables.D2, (variables.z- variables.zp), (variables.zl - variables.zpl)) # limiting field conditions, reference 36
variables.point(variables.xa(variables.mub(variables.D1), variables.z1, variables.zl1)), variables.xa(variables.mub(variables.n), variables.zb(variables.n), variables.zlb(variables.n)) == variables.xa(variables.mub(variables.D1), variables.z1, variables.zl1), variables.xa(variables.mub(variables.n), variables.zb(variables.n), variables.zlb(variables.n)) + summation(variables.subs)

variables.phipb(variables.TensionString, variables.sigma) == variables.phib(variables.TensionString, variables.sigma) + variables.delta*variables.phib(variables.TensionString, variables.sigma) # field symmetry

#unitary CFT
variables.dagger(variables.Lb(variables.m)) == variables.Lb(-variables.m)
variables.dagger(variables.Lbtilde(variables.m)) == variables.Lbtilde(-variables.m)
variables.dirac([variables.D0, variables.k], [variables.D0, variables.kp]) == 2*math.pi*variables.delta(variables.k - variables.kp) #spacelike product mapping
variables.dagger(variables.alphab(variables.m)) == variables.alphab(-variables.m)
variables.dagger(variables.alphabtilde(variables.m)) == variables.alphabtilde(-variables.m)

#129.19
def pointParticleExample():
    variables.deltab(variables.taub(1)) * variables.xa(variables.mu, variables.tau) == -variables.delta * (variables.tau - variables.taub(1)) * diff(variables.xa(variables.mu, variables.tau), variables.tau)
    variables.deltab(variables.taub(1)) * variables.ef(variables.tau) == -diff(variables.brackets(variables.delta * (variables.tau - variables.taub(1) * variables.ef(variables.tau))), variables.tau)
    variables.brackets(variables.deltab(variables.taub(1)), variables.deltab(variables.taub(2))) * variables.xa(variables.mu, variables.tau) == - variables.brackets(variables.delta * (variables.tau - variables.taub(1)) * diff(variables.delta * (variables.tau - variables.taub(2)), variables.tau) * diff(variables.delta * (variables.tau - variables.taub(1)), variables.tau)) * diff(variables.xa(variables.mu, variables.tau), variables.tau) == integrate(diff(variables.taub(3)) * variables.fab(variables.taub(3), variables.taub(1) * variables.taub(2)) * variables.deltab(variables.taub(3)) * variables.xa(variables.mu, variables.tau))
    variables.fab(variables.taub(3), variables.taub(1) * variables.taub(2)) == variables.delta * (variables.taub(3) - variables.taub(1)) * diff(variables.delta * (variables.taub(3) - variables.taub(2)), variables.taub(3)) - variables.delta * (variables.taub(3) - variables.taub(2)) * diff(variables.delta * (variables.taub(3) - variables.taub(1)), variables.taub(3))
    variables.deltab(variables.B) * variables.xa(variables.mu) == variables.i * variables.epsilon * variables.c * diff(variables.xa(variables.mu), variables.tau)
    variables.deltab(variables.B) * variables.e == variables.i * variables.epsilon * diff(variables.c * variables.e, variables.tau)
    variables.deltab(variables.B) * variables.B == 0
    variables.deltab(variables.B) * variables.b == variables.epsilon * variables.B
    variables.deltab(variables.B) * variables.c == variables.i* variables. epsilon * variables.c * diff(variables.c, variables.tau)
    # gauge fixed action
    variables.S == integrate(diff(variables.tau) * (1/2 * variables.e**(-1) * diff(variables.xa(variables.mu), variables.tau) * diff(variables.xb(variables.mu), variables.tau) + 1/2 * variables.e * variables.m**2 + variables.i * variables.B * (variables.e - 1) - variables.e * diff(variables.b, variables.tau)))
    variables.S == integrate(diff(variables.tau) * (1/2 * diff(variables.xa(variables.mu), variables.tau) * diff(variables.xb(variables.mu), variables.tau) + 1/2 * variables.m**2 - diff(variables.b, variables.tau) * variables.c))
    variables.deltab(variables.B) * variables.xa(variables.mu) == variables.i * variables.epsilon * variables.c * diff(variables.xa(variables.mu), variables.tau)
    variables.deltab(variables.B) * variables.b == variables.i* variables.epsilon(-1/2 * diff(variables.xa(variables.mu), variables.tau) * diff(variables.xb(variables.mu), variables.tau) + 1/2 * variables.m**2 - diff(variables.b, variables.tau) * variables.c)
    variables.deltab(variables.B) * variables.c == variables.i * variables.epsilon * variables.c * diff(variables.c)
    variables.brackets(variables.pa(variables.mu), variables.xa(variables.v)) == - variables.i * variables.etaa(variables.mu * variables.v)
    variables.bracketc(variables.b, variables.c) == 1
    variables.Qb(variables.B) == variables.c * variables.Hn
    variables.pa(variables.mu) * variables.dirac(variables.k, variables.down) == variables.ka(variables.mu) * variables.dirac(variables.k, variables.down)
    variables.pa(variables.mu) * variables.dirac(variables.k, variables.up) == variables.ka(variables.mu) * variables.dirac(variables.k, variables.up)
    variables.b * variables.dirac(variables.k, variables.down) == 0
    variables.b * variables.dirac(variables.k, variables.up) == variables.dirac(variables.k, variables.down)
    variables.c * variables.dirac(variables.k, variables.down) == variables.dirac(variables.k, variables.up)
    variables.c * variables.dirac(variables.k, variables.up) == 0
    variables.Qb(variables.B) * variables.dirac(variables.k, variables.down) == 1/2 * (variables.k**2 + variables.m **2) * variables.dirac(variables.k, variables.up)
    variables.Qb(variables.B) * variables.dirac(variables.k, variables.up) == 0
    # reference closed states
    variables.b * variables.dirac(variables.psi) == 0

# open strings 263.1
variables.Abf(25, variables.xas(variables.M)) == -(variables.theta)/(2 * math.pi * variables.R) == - variables.i * variables.Lamdaa(-1) * (diff(variables.Lamda, variables.z))/(diff(variables.xas(25), variables.z))
variables.Lamdaf(variables.xas(25)) == math.exp(-(variables.i * variables.theta * variables.xas(25))/(2 * math.pi * variables.R))
variables.Wb(variables.q) == math.exp(variables.i * variables.q * integrate(diff(variables.xas(25)) * variables.Ab(25), 2 * math.pi)) == math.exp(-variables.i * variables.q * variables.theta) # contour integral
variables.S == integrate(diff(variables.tau) * (1/2 * diff(variables.xa(variables.M), variables.tau) * diff(variables.xb(variables.M), variables.tau) + (variables.ma(2))/2 - variables.i * variables.q * variables.Ab(variables.M) * diff(variables.xa(variables.M), variables.tau)))
variables.pb(25) == - (diff(variables.L), variables.z)/(diff(variables.va(25), variables.z)) == variables.va(25) - (variables.q * variables.theta)/(2 * math.pi * variables.R)
variables.vb(25) == (2 * math.pi * variables.l + variables.q * variables.theta)/(2 * math.pi * variables.R)
variables.H == 1/2 * (variables.pb(variables.mu) * variables.pa(variables.mu) + variables.vab(2, 25) + variables.ma(2))
variables.vb(25) == variables.pb(25) == (2 * math.pi * variables.l + variables.q * variables.theta)/(2 * math.pi * variables.R)
# check gauge transforamtion, 264.8
variables.vb(25) == (2 * math.pi * variables.l - variables.thetab(variables.j) + variables.thetab(variables.i))/(2 * math.pi * variables.R)
variables.ma(2) == ((2 * math.pi * variables.l - variables.thetab(variables.j) + variables.thetab(variables.i))**2)/(4 * math.pi**2 * variables.R**2) + 1/variables.SlopeRegge * (variables.N - 1)
variables.ma(2) == ((variables.thetab(variables.j) - variables.thetab(variables.i))**(2))/(4 * math.pi**2 * variables.R**2)
# ommitted gauge symetry condition 265.12
# T - duality
variables.xap(25, (variables.z, variables.zl)) == variables.xab(25, variables.L, variables.z) - variables.xab(25, variables.R, variables.zl)
diff(variables.xa(25), variables.n) == - variables.i * diff(variables.xap(25), variables.t)
variables.xap(25, math.pi) - variables.xap(25, 0) == quad(diff(variables.sigmaa(1)) * diff(variables.xap(25), 1), (0, math.pi)) == -variables.i * quad(diff(variables.sigmaa(1)) * diff(variables.xa(25), 2), 0, math.pi) == -2 * math.pi * variables.SlopeRegge * variables.va(25) == -(2 * math.pi * variables.SlopeRegge * variables.l)/(variables.R) == -2 * math.pi * variables.l * variables.R 
variables.Delta * variables.xap(25) == variables.xap(25, math.pi) - variables.xap(25, 0) == -(2 * math.pi * variables.l  - variables.thetab(variables.j) + variables.thetab(variables.i))*variables.Rp 
variables.xap(25) == variables.thetab(variables.i) * variables.Rp == -2 * math.pi * variables.SlopeRegge * variables.Ab((25, variables.i *variables.i))
variables.xap(25, (variables.z, variables.zl)) == variables.thetab(variables.i) * variables.Rp - (variables.i * variables.Rp)/(2 * math.pi) * (2 * math.pi * variables.l - variables.thetab(variables.j) + variables.thetab(variables.i)) * ln(variables.z/variables.zl) + variables.i * (variables.SlopeRegge/2)**(1/2) * (summation((variables.alphaab(25, variables.m))/variables.m * (variables.za(-variables.m) - variables.zla(-variables.m)), (variables.m, -math.inf, -1)) + summation((variables.alphaab(25, variables.m))/variables.m * (variables.za(-variables.m) - variables.zla(-variables.m)), (variables.m, 1, math.inf))) == variables.thetab(variables.i) * variables.Rp + (variables.sigmaa(1))/math.pi * variables.Delta * variables.xap(25) - (2 *variables.SlopeRegge)**(1/2) * (summation((variables.alphaab(25, variables.m))/variables.m * math.exp(-variables.m * variables.sigmaa(2)) * math.sin(variables.m * variables.sigmaa(1)), (variables.m, -math.inf, -1)) +summation((variables.alphaab(25, variables.m))/variables.m * math.exp(-variables.m * variables.sigmaa(2)) * math.sin(variables.m * variables.sigmaa(1)), (variables.m, 1, math.inf)))
variables.ma(2) == ((variables.Delta * variables.xap(25))/(2 * math.pi * variables.SlopeRegge))**2 + 1/variables.SlopeRegge * (variables.N - 1)

#2.103.1
variables.thetaa(2) == variables.thetala(2) == variables.bracketc(variables.theta, variables.thetal) == 0
diff(1, variables.z) == (diff(variables.zp, variables.z))/(diff(variables.z, variables.z)) * diff(1, variables.zp) + (diff(variables.zpl, variables.z))/(diff(variables.z, variables.z)) * diff(1, variables.zpl)
variables.Db(variables.theta) == diff(1, variables.theta) == variables.theta * diff(1, variables.z)
variables.Db(variables.thetal) == diff(1, variables.thetal) + variables.thetal * diff(1, variables.zl)
variables.Dab(2, variables.theta) == diff(1, variables.z)
variables.Dab(2, variables.thetal) == diff(1, variables.zl)
variables.bracketc(variables.Db(variables.theta), variables.Db(variables.thetal)) == 0
variables.Db(variables.theta) == variables.Db(variables.theta) * variables.thetap * diff(1, variables.thetap) + variables.Db(variables.theta) * variables.zp * diff(1, variables.zp) + variables.Db(variables.theta) * variables.thetapl * diff(1, variables.thetapl) + variables.Db(variables.theta) * variables.zpl * diff(1, variables.zpl)
variables.Db(variables.theta) * variables.thetapl == variables.Db(variables.theta) * variables.zpl == 0
variables.Db(variables.theta) * variables.zp == variables.thetap * variables.Db(variables.theta) * variables.thetap 
variables.Db(variables.theta) == (variables.Db(variables.theta) * variables.thetap) * variables.Db(variables.thetap)
diff(variables.zp, variables.zl) == diff(variables.zp, variables.thetal) == diff(variables.thetap, variables.zl) == diff(variables.thetap, variables.thetal) == 0
variables.zpf(variables.z, variables.theta) == variables.ff(variables.z) + variables.theta * variables.gf(variables.z) * variables.hf(variables.z)
variables.thetapf(variables.z, variables.theta) == variables.gf(variables.z) + variables.theta * variables.hf(variables.z)
variables.hf(variables.z) == variables.ponseg * variables.brackets(diff(variables.ff(variables.z), variables.z) + variables.gf(variables.z) * diff(variables.gf(variables.z), variables.z))**(1/2)
variables.delta * variables.z == variables.epsilon * variables.brackets(variables.vf(variables.z) - variables.i * variables.theta * variables.etaf(variables.z)) 
variables.delta * variables.theta == variables.epsilon * variables.brackets(- variables.i * variables.etaf(variables.z) + 1/2 * variables.theta * diff(variables.vf(variables.z), variables.z))
(variables.Db(variables.theta) * variables.thetap)**(2 * variables.h) * (variables.Db(variables.thetal) * variables.thetapl)**(2 * variables.htilde) * variables.phifpB(variables.zBp, variables.zbpl) == variables.phifB(variables.zB, variables.zBl)
variables.delta * variables.phifB(variables.zB, variables.zBl) == - variables.epsilon * variables.brackets(2 * variables.h * variables.theta * diff(variables.etaf(variables.z), variables.z) + variables.etaf(variables.z) * variables.Qb(variables.theta) + 2 * variables.htilde * variables.thetal * diff(variables.etafl(variables.zl), variables.zl) + variables.etafl(variables.zl) * variables.Qb(variables.thetal)) * variables.phifB(variables.zB, variables.zBl)
variables.phifB(variables.z) == variables.OOf(variables.z) + variables.theta * variables.Psif(variables.z)
variables.delta * variables.OO == -variables.epsilon * variables.eta * variables.Psi 
variables.delta * variables.Psi == -variables.epsilon * variables.brackets(2 * variables.h * diff(variables.eta * variables.OO, variables.z) + variables.eta * diff(variables.OO, variables.z))
np.dot(variables.Gb(-1/2), variables.OO) == variables.Psi 
np.dot(variables.Gb(variables.r), variables.OO) == 0
variables.R >= 1/2
np.dot(variables.Gb(-1/2), variables.Psi) == diff(variables.OO, variables.z) 
np.dot(variables.Gb(1/2), variables.Psi) == 2 * variables.h * variables.OO 
np.dot(variables.Gb(variables.r), variables.Psi) == 0
variables.r >= 3/2
diff(variables.zp) * diff(variables.thetap) == diff(variables.z) * diff(variables.theta) * variables.Db(variables.theta) * variables.thetap 
variables.S == 1/(4 * math.pi) * integrate(variables.d2z * variables.d2(variables.theta) * variables.Db(variables.thetal) * variables.XaB(variables.mu) * variables.Db(variables.theta) * variables.XbB(variables.mu))
variables.XafB(variables.mu, (variables.zB, variables.zBl)) == variables.xa(variables.mu) + variables.i * variables.theta * variables.psia(variables.mu) + variables.i * variables.thetal * variables.psiatilde(variables.mu) + variables.theta * variables.thetal * variables.Fa(variables.mu)
variables.S == 1/(4 * math.pi) * integrate(variables.d2z * (diff(variables.xa(variables.mu), variables.zl)) * diff(variables.xb(variables.mu), variables.z) + variables.psia(variables.mu) * diff(variables.psib(variables.mu), variables.zl) + variables.psiatilde(variables.mu) * diff(variables.psibtilde(variables.mu), variables.z) + variables.Fa(variables.mu) * variables.Fb(variables.mu))
variables.Db(variables.theta) * variables.Db(variables.thetal) * variables.XafB(variables.mu, (variables.zB, variables.zBl)) == 0
variables.Sb(variables.B * variables.C) == 1/(2 * math.pi) * integrate(variables.d2z * variables.d2(variables.theta) * variables.B * variables.Db(variables.thetal) * variables.C)
variables.Db(variables.thetal) * variables.B == variables.Db(variables.thetal) * variables.C == 0
variables.Bf(variables.zBb(1)) == variables.Betaf(variables.z) + variables.theta * variables.bf(variables.z)
variables.Cf(variables.zB) == variables.cf(variables.z) + variables.theta * variables.gammaf(variables.z)
variables.S == 1/(4 * math.pi) * integrate(variables.d2z * variables.d2(variables.theta) * variables.brackets(variables.Gbf(variables.mu * variables.v, variables.XB) + variables.Bbf(variables.mu * variables.v, variables.XB)) * variables.Db(variables.thetal) * variables.XaB(variables.v) * variables.Db(variables.theta) * variables.XaB(variables.mu)) == 1/(4 * math.pi) * integrate(variables.d2z * variables.bracketc(variables.brackets(variables.Gbf(variables.mu * variables.v, variables.X) + variables.Bbf(variables.mu * variables.v, variables.X)) * diff(variables.xa(variables.mu), variables.z) * diff(variables.xa(variables.v), variables.zl) + variables.Gbf(variables.mu * variables.v, variables.X) * (variables.psia(variables.mu) * variables.DOb(variables.zl) * variables.psia(variables.v) + variables.psiatilde(variables.mu) * variables.DOb(variables.z) * variables.psiatilde(variables.v)) + 1/2 * variables.Rbf((variables.mu * variables.v, variables.rho * variables.sigma), variables.X) * variables.psia(variables.mu) * variables.psia(variables.v) * variables.psiatilde(variables.rho) * variables.psiatilde(variables.sigma)))
variables.DOb(variables.zl) * variables.psia(variables.v) == diff(variables.psia(variables.v), variables.zl) + variables.brackets(variables.Rhoabf(variables.v, variables.rho * variables.sigma, variables.X) + 1/2 * variables.Habf(variables.v, variables.rho * variables.sigma, variables.X)) * diff(variables.xa(variables.rho), variables.zl) * variables.psia(variables.sigma)
variables.DOb(variables.z) * variables.psiatilde(variables.v) == diff(variables.psiatilde(variables.v), variables.z) + variables.brackets(variables.Rhoabf(variables.v, variables.rho * variables.sigma, variables.X) - 1/2 * variables.Habf(variables.v, variables.rho * variables.sigma, variables.X)) * diff(variables.xa(variables.rho), variables.z) * variables.psiatilde(variables.sigma)
variables.XaB(variables.mu) == variables.xa(variables.mu) + variables.i * variables.thetal * variables.psiatilde(variables.mu)
variables.lamdaaB(variables.A) == variables.lamdaa(variables.A) + variables.thetal * variables.Ga(variables.A)
variables.S == 1/(4 * math.pi) * integrate(variables.d2z * diff(variables.thetal) * variables.bracketc(variables.brackets(variables.Gbf(variables.mu * variables.v, variables.XB) + variables.Bbf(variables.mu * variables.v, variables.XB)) * diff(variables.XaB(variables.mu), variables.z) * variables.Db(variables.thetal) * variables.XaB(variables.v) - variables.lamdaaB(variables.A) * variables.DOb(variables.thetal) * variables.lamdaaB(variables.A))) == 1/(4 * math.pi) * integrate(variables.d2z * variables.bracketc(variables.brackets(variables.Gbf(variables.mu * variables.v, variables.X) + variables.Bbf(variables.mu * variables.v, variables.X)) * diff(variables.xa(variables.mu), variables.z) * diff(variables.xa(variables.v), variables.zl) + variables.Gbf(variables.mu * variables.v, variables.X) * variables.psiatilde(variables.mu) * variables.DOb(variables.z) * variables.psiatilde(variables.v) + variables.lamdaa(variables.A) * variables.DOb(variables.zl) * variables.lamdaa(variables.A) + variables.i/2 * variables.Fabf(variables.A * variables.B, variables.rho * variables.sigma, variables.X) * variables.lamdaa(variables.A) * variables.lamdaa(variables.B) * variables.psiatilde(variables.rho) * variables.psiatilde(variables.sigma)))
variables.DOb(variables.thetal) * variables.lamdaaB(variables.A) == variables.Db(variables.thetal) * variables.lamdaaB(variables.A) - variables.i * variables.Aabf(variables.A * variables.B, variables.mu, variables.XB) * variables.Db(variables.thetal) * variables.XaB(variables.mu) * variables.lamdaaB(variables.B)
variables.DOb(variables.zl) * variables.lamdaa(variables.A) == diff(variables.lamdaa(variables.A), variables.zl) - variables.i * variables.Aabf(variables.A * variables.B, variables.mu, variables.X) * diff(variables.xa(variables.mu), variables.zl) * variables.lamdaa(variables.B)
variables.delta * variables.Aab(variables.A * variables.B, variables.mu) == variables.Db(variables.mu) * variables.chia(variables.A * variables.B) 
variables.delta * variables.lamdaa(variables.A) == variables.i * variables.chia(variables.A * variables.B) * variables.lamdaa(variables.B)
variables.Aabf(variables.A * variables.B, variables.zl, (variables.z, variables.zl)) == 1/(2 * math.pi) * variables.Aabf(variables.A * variables.B, variables.mu, variables.X) * diff(variables.xa(variables.mu), variables.zl)
variables.delta * variables.Za(variables.brackets(variables.A)) == 1/(8 * math.pi) * integrate(variables.d2z * variables.Trbf(variables.v, variables.brackets(variables.chif(variables.X) * variables.Fbf(variables.mu * variables.v, variables.X))) * diff(variables.xa(variables.mu), variables.z) * diff(variables.xa(variables.v), variables.zl))
variables.delta * variables.Bb(variables.mu * variables.v) == 1/2 * variables.Trbf(variables.v, variables.chi * variables.Fb(variables.mu * variables.v))
(variables.kappaab(2, 10))/(variables.gab(2, 10)) == 1/2 == variables.SlopeRegge/4
(variables.kappaa(2))/(variables.gab(2, variables.Y * variables.M)) == (variables.ea(2 * variables.Phi) * variables.kapaab(2, 10))/(variables.ea(1 * variables.Phi) * variables.gab(2, 10)) == variables.SlopeRegge/4
variables.deltaf(variables.gamma) * variables.deltaf(variables.gammatilde) == variables.ea(-variables.phi - variables.phitilde)
variables.VOa(0, 0) == np.dot(variables.Gb(-1/2) * variables.Gbtilde(-1/2), variables.OO)
variables.VOa(-1, -1) == variables.ea(-variables.phi - variables.phitilde) * variables.OO 
variables.psiab(variables.mu, -1/2) * variables.psiabtilde(variables.v, -1/2) * variables.diracb((0, variables.k), variables.NS)
variables.VOa(-1, -1) == variables.gb(variables.c) * variables.ea(-variables.phi - variables.phitilde) * variables.psia(variables.mu) * variables.psiatilde(variables.v) * variables.ea(np.dot(variables.i * variables.k, variables.X))
variables.Gb(-1/2) * variables.Gbtilde(-1/2) * variables.psiab(variables.mu, -1/2) * variables.psiabtilde(variables.v, -1/2) * variables.diracb((0, variables.k), variables.NS) == - (variables.alphaab(variables.mu, -1) + np.dot(variables.alphab(0), variables.psib(-1/2) * variables.psiab(variables.mu, -1/2))) * (variables.alphaabtilde(variables.v, -1) + np.dot(variables.alphab(0), variables.psibtilde(-1/2) * variables.psiabtilde(variables.v, -1/2))) * variables.diracb((0, variables.k), variables.NS)
variables.VOa(0, 0) == -(2 * variables.gb(variables.c))/variables.SlopeRegge * (variables.i * diff(variables.xa(variables.mu), variables.z) + np.dot(1/2 * variables.SlopeRegge * variables.k, variables.psi * variables.psia(variables.mu))) * (variables.i * diff(variables.xa(variables.v), variables.zl) + np.dot(1/2 * variables.SlopeRegge * variables.k, variables.psitilde * variables.psiatilde(variables.v))) * variables.ea(np.dot(variables.i * variables.k, variables.X))
variables.VOa(-1) == variables.gb(variables.o) * variables.ea(-variables.psi) * variables.ta(variables.a) * variables.psia(variables.mu) * variables.ea(np.dot(variables.i * variables.k, variables.X))
variables.VOa(0) == variables.gb(variables.o) * (2 * variables.SlopeRegge)**(-1/2) * variables.ta(variables.a) * (variables.i * diff(variables.xa(variables.mu), variables.tau) + np.dot(2 * variables.SlopeRegge * variables.k, variables.psi * variables.psia(variables.mu))) * variables.ea(np.dot(variables.i * variables.k, variables.X))
# type I
variables.gb(variables.o) == variables.gb(variables.Y * variables.M) * (2 * variables.SlopeRegge)**(1/2) 
variables.gb(variables.Y * variables.M) == variables.gb(10) * variables.ea(variables.Phi/2)
# heterotic
variables.gb(variables.c) == (variables.kappa)/(2 * math.pi) == (variables.SlopeRegge**(1/2) * variables.gb(variables.Y * variables.M))/(4 * math.pi) 
variables.kappa == variables.kappab(10) * variables.ea(variables.Phi) 
variables.gb(variables.Y * variables.M) == variables.gb(10) * variables.ea(variables.Phi)
# type I/II
variables.gb(variables.c) == variables.kappa/(2 * math.pi) 
variables.kappa == variables.kappab(10) * variables.ea(variables.Phi)
 
# 2.386.1
summation(variables.lamdaa(variables.pos * variables.K) * variables.lamdaa(variables.neg * variables.K), (variables.K, 4, 8))
variables.lamdaa(variables.A) * variables.lamdaa(variables.B) 
variables.ThetabB(variables.B16) * math.exp(3**(1/2) * variables.i * variables.H/2) 
variables.ThetabB(variables.B16l) * math.exp(-3**(1/2) * variables.i * variables.H /2) 
variables.i * diff(variables.H, variables.z)
variables.B45 + variables.B16 + variables.B16l + variables.B1 
variables.VO == variables.lamdaa(variables.A) * variables.PhiaB(variables.pos, variables.pos)
variables.ThetabB(16) * variables.PhiafB((variables.pos, variables.pos), 1 == -1/2) # reference 2.387.5
variables.h == variables.Q/2 + (variables.Qpa(2) - variables.Qa(2))/6
variables.htilde == (variables.Qtilde)/2 
variables.PhiafB((variables.pos, variables.pos), 1 == -2)
np.dot(variables.Gab(variables.neg, -1/2), variables.PhiaB(variables.pos, variables.pos)) 
np.dot(variables.Gab(variables.neg, -1/2) * variables.Gab(variables.neg, -1/2), variables.PhiaB(variables.pos, variables.posneg)) == 0
np.dot(variables.Gab(variables.neg, -1/2) * variables.Gab(variables.pos, -1/2), variables.PhiaB(variables.neg, variables.posneg)) == np.dot((2 * variables.Lb(-1) - variables.Gab(variables.pos, -1/2) * variables.Gab(variables.neg, -1/2)), variables.PhiaB(variables.neg, variables.posneg)) == 2 * diff(variables.PhiaB(variables.neg, variables.posneg), variables.z)
(variables.Gab(variables.neg, -1/2))**2 == 0
np.dot(variables.Gab(variables.neg, -1/2), variables.PhiaB(variables.neg, variables.posneg)) == 0
variables.Gbp(variables.A * variables.Bl) == math.exp(variables.brackets(variables.kappaa(2) * (variables.Kb(2) - variables.Kb(1))/3)) * variables.Gb(variables.A * variables.Bl)
variables.Wf(variables.Phi) == variables.Phiab(variables.A, variables.xl) * variables.Phiab(variables.B, variables.yl) * variables.Phiab(variables.C, variables.zl) * variables.da(variables.xl * variables.yl * variables.zl) * diff(diff(diff(variables.Fbf(variables.T), variables.C), variables.B), variables.A) 
variables.Gbp(variables.a * variables.bl) == math.exp(variables.brackets(variables.kappaa(2) * (variables.Kb(1) - variables.Kb(2))/3)) * variables.Gb(variables.a * variables.bl)
variables.Wf(variables.chi) == variables.chiab(variables.a, variables.x) * variables.chiab(variables.b, variables.y) * variables.chiab(variables.c, variables.z) * variables.da(variables.x * variables.y * variables.z) * diff(diff(diff(variables.Fbf(2, variables.Z), variables.c), variables.b), variables.a)
# cont. 19.4, 2.391.1
variables.c == 3 - 6/(variables.k + 2) == (3 * variables.k)/(variables.k + 2)
# NS:
variables.h == (variables.l * (variables.l + 2) - variables.qa(2))/(4 * (variables.k + 2))
variables.Q == variables.q/(variables.k + 2)
# R:
variables.h == (variables.l * (variables.l + 2) - (variables.q + variables.posneg * 1)**2)/(4 * (variables.k + 2)) + 1/8
variables.Q == (variables.q + variables.posneg * 1)/(variables.k + 2) + variables.negpos * 1/2 
# -
variables.ja(variables.pos) == variables.psib(1) * math.exp(variables.brackets(variables.i * (2/variables.k)**(1/2) * variables.H))
variables.ja(variables.neg) == variables.dagger(variables.psib(1)) * math.exp(variables.brackets(-variables.i * (2/variables.k)**(1/2) * variables.H))
variables.Tab(variables.pos, variables.F) == variables.psib(1) * math.exp(variables.brackets(variables.i * ((variables.k + 2)/(variables.k))**1/2 * variables.H)) 
variables.Tab(variables.neg, variables.F) == variables.dagger(variables.psib(1)) * math.exp(variables.brackets(-variables.i * ((variables.k + 2)/(variables.k))**(1/2) * variables.H))
1 - 1/2 * (2/variables.k) + 1/2 * (variables.k + 2)/(variables.k) == 3/2
variables.OOab(variables.j, variables.m) == variables.psiab(variables.j, variables.m) * math.exp(variables.brackets(variables.i * variables.m * (2/variables.k)**(1/2) * variables.H))
variables.OOpab(variables.j, variables.m) == variables.psiab(variables.j, variables.m) * math.exp(variables.brackets(variables.i * (2 * variables.m)/(variables.ka(1/2) * (variables.k + 2)**(1/2)) * variables.H))
variables.h == (variables.j * (variables.j + 1))/(variables.k + 2) - (variables.ma(2))/(variables.k) + (2 * variables.ma(2))/(variables.k * (variables.k + 2)) == (variables.j * (variables.j + 1) - variables.ma(2))/(variables.k + 2)
variables.j == variables.i * variables.brackets(variables.k/(variables.k + 2))**(1/2) * diff(variables.H, variables.z)
variables.psiab(variables.j, variables.m) * math.exp(variables.brackets(variables.i * (2 * variables.m + variables.posneg * variables.k/2)/(variables.ka(1/2) * (variables.k + 2)**(1/2)) * variables.H))
variables.Wf(variables.PhiB) == variables.PhiaB(variables.k + 2)
variables.sigma == variables.lamda * variables.sigma 
variables.Phi == variables.lamdaa(variables.omegas) * variables.Phi 
variables.psi == variables.lamdaa(variables.omegas - 1/2) * variables.psi 
variables.F == variables.lamdaa(variables.omegas - 1) * variables.F 
variables.lamdaa(2 - 1 + (variables.k + 2 * variables.omegas))
variables.hb(variables.Phi) == variables.hbtilde(variables.Phi) == 1/(2 * (variables.k + 2))
variables.Qb(variables.Phi) == variables.Qbtilde(variables.Phi) == 1/(variables.k + 2)
diff(variables.Wf(variables.Phi), variables.Phi) == (variables.k + 2) * variables.Phia(variables.k + 1) == 0
variables.Q == variables.c/3 == variables.k/(variables.k + 2) 
variables.PhiB == math.exp((2 * math.pi * variables.i)/(variables.k + 2)) * variables.PhiB 
math.exp(2 * math.pi * variables.i * variables.Q) 
# gepner models, 2.395.1 
summation((variables.kb(variables.i))/(variables.kb(variables.i) + 2), variables.i) == 3
variables.l + variables.sb(0) + variables.sb(1) + variables.Q == 2 * variables.ZB 
variables.gb(variables.q) == math.exp(math.pi * variables.i * variables.s + 2 * math.pi * variables.i * variables.Q) == math.exp(math.pi * variables.i * variables.s)  * np.prod(math.exp(2 * math.pi * variables.i * variables.Qb(variables.i)), variables.i)
(variables.N * variables.k)/(variables.k + 2) == 3 
variables.ka(variables.N) == 1**(9), 2**(6), 3**(5), 6**(4)
np.prod(variables.psiab(variables.jb(variables.i), variables.mb(variables.i)) * variables.psiabtilde(variables.jb(variables.i), variables.mb(variables.i)) * math.exp(variables.brackets(variables.i * (2 * variables.m * (variables.Hb(variables.i) + variables.Hbtilde(variables.i)))/(variables.ka(1/2) * (variables.k + 2)**(1/2)))), (variables.i, 1, variables.N))
math.exp(variables.brackets(variables.i * variables.l * ((variables.k + 2)/(variables.k))**(1/2) * variables.Hb(variables.i)))
math.exp(variables.brackets(variables.i * variables.n * (variables.k/(variables.k + 2))**(1/2) * variables.Hb(variables.i)))
np.prod(variables.psiab(variables.jb(variables.i), variables.mb(variables.i)) * variables.psiabtilde(variables.jb(variables.i), variables.mb(variables.i)) * math.exp(variables.brackets(variables.i * ((2 * variables.mb(variables.i) + variables.n * variables.k) * variables.Hb(variables.i) + 2 * variables.mb(variables.i) * variables.Hbtilde(variables.i))/(variables.ka(1/2) * (variables.k + 2)**(1/2)))), (variables.i, 1, variables.N))
variables.Q == 1/(variables.k + 2) * summation((2 * variables.mb(variables.i) + variables.n * variables.k), (variables.i, 1, variables.N)) == 3 * variables.n  + 2/(variables.k + 2) * summation(variables.mb(variables.i), (variables.i, 1, variables.N))
variables.Lb(0) - variables.Lbtilde(0) == 1/(2 * variables.k * (variables.k + 2)) * summation(variables.bracketc((2 * variables.mb(variables.i) + variables.n * variables.k)**2 - (2 * variables.mb(variables.i))**2), (variables.i, 1, variables.N)) == (3 * variables.na(2))/2 + (2 * variables.n)/(variables.k + 2) * summation(variables.mb(variables.i), (variables.i, 1, variables.N))
variables.htilde - variables.Qtilde/2 == summation((variables.jb(variables.i) * (variables.jb(variables.i) + 1) - variables.mb(variables.i) * (variables.mb(variables.i) + 1))/(variables.k + 2), (variables.i,1 , variables.N))
summation(variables.jb(variables.i), variables.i) =(variables.k + 2)/(2) 
Abs(variables.jb(variables.i)) <= variables.k/2 
variables.mb(variables.i) == -variables.mbtilde(variables.i) == -variables.jb(variables.i) 
# check the numbers for models, 2.397.15
variables.SemiDirect(variables.Sb(5), variables.ZabB(4, 5))
math.exp(2 * math.pi * variables.i * variables.Qb(variables.i)) 
variables.i == 2, 3, 4, 5 
variables.gf(variables.z) == variables.zab(5, 1) + variables.zab(5, 2) + variables.zab(5, 3) + variables.zab(5, 4) + variables.zab(5, 5)
variables.zb(variables.i) == math.exp(2 * math.pi * variables.i * (variables.nb(variables.i))/5) * variables.zb(variables.i) 
variables.nb(1) == 0
summation(variables.PhiabB(5, variables.i), (variables.i, 1, 5))
variables.W == variables.P * variables.Gf(variables.PhiB)
variables.qb(variables.PhiB) == 1
variables.qb(variables.P) == - 5
variables.U == Abs(variables.Gf(variables.Phi))**2 == Abs(variables.p)**2 * summation(Abs((diff(variables.G, variables.z))/(diff(variables.Phib(variables.i))))**2, (variables.i, 1, 5)) + (variables.ea(2))/2 * (variables.r + 5 * abs(variables.p)**2 - summation(Abs(variables.Phib(variables.i))**2, (variables.i, 1, 5)))**2 + (variables.Aab(2, 2) + variables.Aab(2, 3)) * (25 * Abs(variables.p)**2 + summation(Abs(variables.Phib(variables.i))**2, (variables.i, 1, 5)))
(diff(variables.G, variables.z))/diff(variables.Phib(variables.i), variables.z) == 0
variables.p == 0 
summation(Abs(variables.Phib(variables.i))**2, (variables.i, 1, 5)) == variables.r 
variables.Gf(variables.Phi) == 0
variables.Rab(2, variables.c) == variables.r 
Abs(variables.p)**2 == variables.r/5 
variables.Phib(variables.i) == 0
variables.Ab(2) == variables.Ab(3) == 0
variables.W == variables.bracket(variables.p) * variables.Gf(variables.PhiB)
variables.p == variables.p 
variables.Phib(variables.i) == math.exp(2 * math.pi * variables.i/5) * variables.Phib(variables.i)
variables.i * (variables.theta)/(2 * math.pi) * integrate(variables.Fb(2))
variables.Fb(12) == variables.theta /(2 * math.pi)
