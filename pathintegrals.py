import math
from scipy.integrate import quad
import numpy as np 
from sympy import symbols, diff,integrate, Function, Derivative, summation, oo, Abs, ln
from sympy.diffgeom import WedgeProduct
import variables
import constants

variables.PI(variables.X) == integrate(diff(variables.X)*math.exp(-variables.S)*variables.PI(variables.X)) # 34
integrate(diff(variables.X)*variables.NambuGotoVar/(variables.NambuGotoVar(variables.xb(variables.mu, variables.z, variables.zl)))*math.exp(-variables.S)) == -integrate(diff(variables.X)*math.exp(-variables.S)*(variables.NambuGotoVar*variables.S)/(variables.NambuGotoVar*variables.xb(variables.z, variables.zl))) == -variables.bracket((variables.NambuGotoVar*variables.S)/(variables.NambuGotoVar*variables.xb(variables.mu, variables.z, variables.zl))) == 1/(math.pi * variables.SlopeRegge) * variables.bracket(diff(diff(variables.xa(variables.mu, variables.z, variables.zl), variables.zl), variables.z)) == 0 # more formal application with Grassmann variables
variables.bracket(diff(diff(variables.xa(variables.mu, variables.z, variables.zl), variables.zl), variables.z)) == 0 
diff(diff(variables.xahat(variables.mu, variables.z, variables.zl), variables.zl), variables.z) == 0 # follows Hilbert space formallism
integrate(diff(variables.X)*variables.NambuGotoVar/(variables.NambuGotoVar*variables.xb(variables.mu, variables.z, variables.zl))*(math.exp(-variables.S)*variables.xa(variables.v, variables.zp, variables.zpl))) == integrate(diff(variables.X)* math.exp(-variables.S)*(variables.etaa(variables.mu, variables.v)*variables.delta(variables.D2, (variables.z - variables.zp), (variables.zl - variables.zpl)) + 1/(math.pi*variables.SlopeRegge)(diff(diff(variables.xa(variables.mu, variables.z, variables.zl) * variables.xa(variables.v, variables.zp, variables.zpl), variables.zl), variables.z)))) == variables.etaa(variables.mu, variables.v) * variables.bracket(variables.delta(variables.D2, (variables.z - variables.zp), (variables.zl - variables.zpl))) + 1/(math.pi*variables.SlopeRegge)*diff(diff(variables.bracket(variables.xa(variables.mu, variables.z, variables.zl) * variables.xa(variables.v, variables.zp, variables.zpl)), variables.zl), variables.z)

#66.17
variables.psib(variables.AO(variables.phib(variables.b))) == integrate(variables.bracketsb(diff(variables.phib(variables.i))), variables.phib(variables.b) * math.exp(-variables.S*variables.brackets(variables.phib(variables.i))) * variables.AO(variables.D0))
variables.Psib(variables.D1) * variables.brackets(variables.xb(variables.b)) == integrate(variables.bracketsb(diff(variables.xb(variables.i)), variables.xb(variables.b)) * math.exp(-1/(2*math.pi * variables.SlopeRegge) * integrate(variables.d2z * diff(variables.X, variables.z) * diff(variables.X, variables.zl))))
variables.xb(variables.i) == variables.xb(variables.c*variables.l) + variables.xbp(variables.i) # evalutation with gausian method
variables.xb(variables.c*variables.l, variables.z, variables.zl) == variables.xb(variables.D0) + summation(variables.za(variables.n) * variables.xb(variables.n) + variables.zla(variables.n) * variables.xb(-variables.n), (variables.n, 1, math.inf))
variables.Psib(variables.D1) * variables.brackets(variables.xb(variables.b)) == math.exp(-variables.Sb(variables.c*variables.l)) * integrate(variables.bracketsb(diff(variables.xbp(variables.i)), variables.xb(variables.b)) * math.exp(-1/(2*math.pi*variables.SlopeRegge) * integrate(variables.d2z * diff(variables.Xp, variables.z) * diff(variables.Xp, variables.zl))))
variables.xb(variables.b) == 0
variables.Sb(variables.c*variables.l) == 1/(2*math.pi*variables.SlopeRegge) * summation(variables.m*variables.n*variables.xb(variables.m) * variables.xb(-variables.n) * quad(variables.d2z * variables.za(variables.m-1) * variables.zla(variables.n-1), Abs(variables.z) < 1), ((variables.m, variables.n), 1, math.inf)) == 1/(variables.SlopeRegge) * summation(variables.m*variables.xb(variables.m)*variables.xb(-variables.m), (variables.m, 1, math.inf))
variables.Psib(variables.D1) * variables.brackets(variables.xb(variables.b)) == math.exp(-1/(variables.SlopeRegge) * summation(variables.m*variables.xb(variables.m)*variables.xb(-variables.m), (variables.m, 1, math.inf)))
variables.alphab(variables.n) == (variables.i*variables.n)/(2*variables.SlopeRegge)**(1/2) * variables.xb(-variables.n) - variables.i*(variables.SlopeRegge/2)**(1/2) * diff(1/(diff(variables.xb(variables.n), variables.z)), variables.z)
variables.alphabtilde(variables.n) == -(variables.i*variables.n)/(2*variables.SlopeRegge)**(1/2) * variables.xb(-variables.n) - variables.i*(variables.SlopeRegge/2)**(1/2) * (diff(1/(diff(variables.xb(-variables.n), variables.z))), variables.z)
variables.alphab(variables.n) * variables.Psib(variables.D1) * variables.brackets(variables.xb(variables.b)) == variables.alphabtilde(variables.n) * variables.Psib(variables.D1) * variables.brackets(variables.xb(variables.b)) == 0
variables.n >= 0
variables.dirac(variables.D1) == variables.dirac(variables.D0, variables.D0)
variables.dirac(diff(variables.X, variables.k, variables.z)) == math.factorial(variables.k) * variables.xb(variables.k) * variables.Psib(variables.D1) == -variables.i * (variables.SlopeRegge/2)**(1/2) * math.factorial(variables.k-1) * variables.alphab(-variables.k) * variables.dirac(variables.D0, variables.D0)

#82.4
integrate(variables.brackets(diff(variables.eta) * diff(variables.X))*  math.exp(variables.i/2 * integrate(diff(variables.tau) * (variables.eta**(-1) * diff(variables.xa(variables.mu)) * diff(variables.xb(variables.mu)) - variables.eta * variables.m**2)))) # check conditions for redefinition of integral

# 121.1
def pathintegralchange():
    variables.delta * variables.bracket(variables.f, variables.i) == -1/(4*math.pi) * integrate(variables.d2sigma * variables.gf(variables.sigma)**(1/2) * variables.delta*variables.gbf(variables.a*variables.b, variables.sigma) * variables.bracket(variables.f, variables.Taf(variables.a*variables.b, variables.sigma), variables.i))
    variables.bracket(variables.psi, variables.Taf(variables.a*variables.b, variables.sigma), variables.psip)
    variables.Tb(variables.a*variables.b) == variables.Tab(variables.X, variables.a*variables.b) + variables.Tab(variables.g, variables.a*variables.b)
    variables.Tb(variables.a*variables.b) == variables.Tab(variables.m,variables.a*variables.b) + variables.Tab(variables.g, variables.a*variables.b)
    (variables.Lab(variables.m, variables.n) + variables.A * variables.deltab(variables.n, variables.D0)) * variables.dirac(variables.psi) == 0 # for
    variables.n >=0
    variables.n < 0 #then
    variables.bracket(variables.psi, variables.Lab(variables.m, variables.n), variables.psip) == variables.bracket(variables.Lab(variables.m, -variables.n) * variables.psi, variables.psip) == 0
    variables.dagger(variables.Lab(variables.m, variables.n)) == variables.Lab(variables.m, -variables.n) #needed
    variables.dirac(variables.chi) == summation(variables.Lab(variables.m, -variables.n) * variables.dirac(variables.chib(variables.n)), (variables.n, 1, math.inf))
    variables.dirac(variables.psi) == variables.dirac(variables.psi) + variables.dirac(variables.chi)
    variables.hilbertb(variables.O, variables.C, variables.Q) == variables.hilbertb(variables.phys)/(variables.hilbertb(variables.null))

#145.1
def euclidianPathIntegralEvaluation():
    diff(variables.taup)/(diff(variables.tau)) == variables.ef(variables.tau)
    variables.taupf(variables.tau) == quad(diff(variables.taupp) * variables.ef(variables.taupp), 0, variables.tau)
    variables.taupf(1) == quad(diff(variables.tau) * variables.ef(variables.tau), 0, variables.l) == variables.l 
    variables.ep == variables.l 
    0 <= variables.tau <= 1
    # or
    variables.ep == 1
    0 <= variables.tau <= variables.l 
    # for region
    0 <= variables.sigmaa(1) <= 2*math.pi 
    0 <= variables.sigmaa(2) <= 2*math.pi 
    variables.point(variables.sigmaa(1), variables.sigmaa(2)) == variables.point(variables.sigmaa(1), variables.sigmaa(2) + 2*math.pi * variables.point(variables.m, variables.n))
    diff(variables.sa(2)) == variables.linea(diff(variables.sigmaa(1)) + variables.tau * diff(variables.sigmaa(2)), 2)
    variables.sigmaatilde(variables.a) == variables.sigmaatilde(variables.a) + 2 * math.pi * (variables.m * variables.ua(variables.a) + variables.n * variables.va(variables.a))
    variables.w == variables.w + 2*math.pi * (variables.m + variables.n * variables.tau)
    variables.taup == variables.tau + 1 # T transofrmation 148.12
    variables.taup == -1/variables.tau # S
    variables.taup == (variables.a * variables.tau + variables.b)/(variables.c * variables.tau + variables.d)
    variables.a*variables.d - variables.b * variables.c == 1
    np.array([variables.sigmaa(1)], [variables.sigmaa(2)]) == np.array([variables.d, variables.b], [variables.c, variables.a]) * np.array([variables.sigmaap(1)], [variables.sigmaap(2)])
    # check boundry conditions at 148.15
    variables.sigmaa(variables.a) = variables.sigmaa(variables.a) + variables.va(variables.a)

#200.10
variables.diracbs(variables.i, variables.j) == variables.bracketb(variables.AObp(variables.i, (math.inf, math.inf)) * variables.AOb(variables.j, (0, 0)), variables.Sb(2))
variables.diracbs(variables.i, variables.j) == variables.posneg * variables.diracbs(variables.j, variables.i)
variables.bracketb(variables.AObp(variables.i, (math.inf, math.inf)) * variables.AOb(variables.k, (1, 1)) * variables.AOb(variables.j, (0, 0)), variables.Sb(2)) == variables.diracbs(variables.i, variables.hat(variables.AOb(variables.k, (1, 1))), variables.j)
summation(variables.cab(variables.l, variables.k * variables.j) * variables.bracketb(variables.AObp(variables.i, (math.inf, math.inf)) * variables.AOb(variables.l, (0, 0)), variables.Sb(2)), variables.l) == variables.cb(variables.i * variables.k * variables.j)
variables.bracketb(variables.AObp(variables.i, (math.inf, math.inf)) * variables.AOb(variables.k, (variables.z(variables.zb(1), variables.zlb(1)))) * variables.AOb(variables.j, (0, 0)), variables.Sb(2)) == variables.zab(variables.hb(variables.i) - variables.hb(variables.k) - variables.hb(variables.j), 1) * variables.zlb(variables.hbtilde(variables.i) - variables.hbtilde(variables.k) - variables.hbtilde(variables.j)) * variables.cb(variables.i * variables.k * variables.j)
variables.bracketb(variables.AObp(variables.i, (math.inf, math.inf)) * variables.AOb(variables.k, (variables.zb(1), variables.zlb(1))) * variables.AOb(variables.l, (variables.zb(2), variables.zlb(2))) * variables.AOb(variables.j, (0, 0)), variables.Sb(2)) == variables.diracbs(variables.i, variables.Tf(variables.hat(variables.AOb(variables.k, (variables.zb(1), variables.zlb(1)))) * variables.hat(variables.AOb(variables.l, (variables.zb(2), variables.zlb(2))))), variables.j)
1 == variables.dirac(variables.m) * variables.gOa(variables.m * variables.n) * variables.diracbs(variables.n)
summation(variables.zab(variables.hb(variables.i) - variables.hb(variables.k) - variables.hb(variables.m), 1) * variables.zlab(variables.hbtilde(variables.i) - variables.hbtilde(variables.k) - variables.hbtilde(variables.m), 1) * variables.zab(variables.hb(variables.m) - variables.hb(variables.l) - variables.hb(variables.j)) * variables.zlab(variables.hbtilde(variables.m) - variables.hbtilde(variables.l) - variables.hbtilde(variables.j), 2) * variables.cb(variables.i * variables.k *variables.m) * variables.cab(variables.m, variables.l * variables.j), variables.m) # four point amplitude

# 311.1
integrate(variables.brackets(diff(variables.Psi)) * math.exp(variables.i * variables.S * variables.brackets(variables.Psi)))
variables.delta * variables.Psi == variables.Qb(variables.B) * variables.Lamda 
variables.Sb(0) == 1/2 * variables.bracket(variables.dirac(variables.Psi, variables.Qb(variables.B), variables.Psi))
variables.Qb(variables.B) * variables.Psi == 0
variables.Psi * variables.brackets(variables.X, variables.c, variables.ctilde) == summation(variables.Phibf(variables.i, variables.x) * variables.Psib(variables.i) * variables.brackets(variables.Xp, variables.c, variables.ctilde), variables.i)
integrate(variables.brackets(diff(variables.Psi))) == np.prod(integrate(variables.brackets(diff(variables.Phib(variables.i)))), variables.i)
#ommitted 312.6
variables.delta * variables.Psi == variables.Qb(variables.B) * variables.Lamda + variables.starf(variables.g * variables.Psi, variables.Lamda) - variables.starf(variables.g * variables.Lamda, variables.Psi)
variables.Qb(variables.B) * (variables.starf(variables.psib(1), variables.Psib(2))) == variables.starf(variables.Qb(variables.B) * variables.Psib(1), variables.Psib(2)) + variables.starf(variables.Psib(1), variables.Qb(variables.B) * variables.Psib(2))
integrate(variables.starf(variables.Psib(1), variables.Psib(2))) == integrate(variables.starf(variables.Psib(2), variables.Psib(1)))
variables.S == 1/2 * integrate(variables.starf(variables.Psi, variables.Qb(variables.B) * variables.Psi) + (2 * variables.g)/3 * integrate(variables.starf(variables.Psi, variables.Psi, variables.Psi)))
variables.S == 1/2 * variables.bracketb(np.dot(variables.VOb(variables.Psi) * variables.Qb(variables.B), variables.VOb(variables.Psi)), variables.Db(2)) + (2 * variables.g)/3 * variables.bracketb(variables.VOb(variables.Psi) * variables.VOb(variables.Psi) * variables.VOb(variables.Psi), variables.Db(2))
quad(diff(variables.t) * math.exp(-variables.t * variables.Lb(0)), 0, math.inf) == variables.Lab(-1, 0)
variables.Sbp(0) == 1/2 * variables.bracket(variables.dirac(variables.Psi, (variables.cb(0) - variables.cbtilde(0)) * variables.Qb(variables.B), variables.Psi))
# 9.7, 315.1
quad(diff(variables.y) * math.exp(variables.brackets(- (variables.ya(2))/(math.factorial(2)) - variables.lamda * ((variables.ya(3))/(math.factorial(3))))), -math.inf, math.inf) == variables.lamdaa(-1) * quad(diff(variables.z) * math.exp(variables.brackets(-(1)/(variables.lamdaa(2) * ((variables.za(2))/(math.factorial(2)) + (variables.za(3))/(math.factorial(3)))))), -math.inf, math.inf)
summation(((-variables.lamda)**variables.n)/(6**(variables.n) * math.factorial(variables.n)) * quad(diff(variables.y) * variables.ya(3 * variables.n) * math.exp(-(variables.ya(2))/2), -math.inf, math.inf), (variables.n, 0, math.inf)) == (2 * math.pi)**(1/2) * summation(variables.lamdaa(2 * variables.k) * variables.Cb(2 * variables.k), (variables.k, 0, math.inf))
variables.Cb(2 * variables.k) == (2**variables.k * variables.Rho(3 * variables.k + 1/2))/(math.pi**(1/2) * 3**(2 * variables.k) * math.factorial(2 * variables.k))
variables.Cb(variables.k * 2) == variables.k **variables.k == math.factorial(variables.k)
summation(variables.lamdaa(2 * variables.k) * math.factorial(variables.k) * variables.fb(2 *variables.k), (variables.k, 0 , math.inf))
math.exp(-variables.Of(1/(variables.lamdaa(2))))
quad(diff(variables.t) * math.exp(-variables.t) * summation((variables.t * variables.lamdaa(2))**variables.k * variables.fb(2 * variables.k), (variables.k, 0, math.inf)), 0, math.inf)
summation(variables.gab(2 * variables.k, variables.o) * variables.Of(math.factorial(variables.k)), (variables.k, 0, math.inf)) == summation(variables.gab(variables.k, variables.c) * variables.Of(math.factorial(variables.k)), (variables.k, 0, math.inf))
math.exp(-variables.Of(1/(variables.gb(variables.c))))

#hard scattering, 9.8, 317.1
(8 * math.pi * variables.i * variables.gab(2, variables.c))/(variables.SlopeRegge) * (2 * math.pi)**26 * variables.deltaaf(26, summation(variables.kb(variables.i), variables.i)) * integrate(variables.d2(variables.zb(4)) * Abs(variables.zb(4))**(-variables.SlopeRegge * variables.u / 2 - 4) * Abs(1 - variables.zb(4))**(-variables.SlopeRegge * variables.t/2 - 4))
(diff(1, variables.z))/(diff(variables.zb(4), variables.z)) * (variables.u * ln(Abs(variables.zb(4))**2) + variables.t * ln(Abs(1 - variables.zb(4))**2)) == (variables.u)/(variables.zb(4)) + variables.t/(variables.zb(4) - 1) == 0
variables.S == math.exp(variables.brackets(-variables.SlopeRegge * (variables.s * ln(variables.s) + variables.t * ln(variables.t) + variables.u * ln(variables.u))/2))
math.exp(variables.brackets(-summation(np.dot(variables.kb(variables.i), variables.kb(variables.j) * variables.Gpf(variables.sigmab(variables.i), variables.sigmab(variables.j))), (variables.i < variables.j))))
variables.zab(variables.N, 2) == ((variables.zb(1) - variables.ab(1)) * (variables.zb(1 - variables.ab(2))))/((variables.zb(1) - variables.ab(3)) * (variables.zb(1) - variables.ab(4)))
variables.Sb(variables.g) == math.exp(variables.brackets(- variables.SlopeRegge * (variables.s * ln(variables.s) + variables.t * ln(variables.t) + variables.u * ln(variables.u))/(2 * (variables.g + 1)))) # double check 319.7
math.exp(variables.brackets(variables.SlopeRegge * variables.t * ln(variables.s)/(2 * (variables.g + 1))))
variables.bracket(0, variables.bracketsa(variables.xa(1, variables.sigma) - variables.xas(1), 2), 0) == summation((variables.SlopeRegge)/(2 * variables.na(2)) * variables.bracket(0, variables.alphab(variables.n) * variables.alphab(-variables.n) + variables.alphabtilde(variables.n) * variables.alphabtilde(-variables.n), 0), (variables.n, 1, math.inf)) == variables.SlopeRegge * summation(1/variables.n, (variables.n, 1, math.inf)) 
# ommitted 230.10
quad(diff(variables.m) * variables.nf(variables.m) * math.exp(-variables.SlopeRegge * math.pi * variables.ma(2) * variables.ls), 0, math.inf) == math.exp(4 * math.pi /variables.ls)
quad(diff(variables.m) * math.exp(4 * math.pi * variables.m * variables.SlopeRegge**(1/2) * math.exp(-variables.m/variables.t)), 0, math.inf)
variables.Tb(variables.H) == 1/(4 * math.pi * variables.SlopeRegge**(1/2))
variables.Ff(variables.T, variables.ma(2)) == variables.T * integrate(variables.da(variables.d-1) * variables.kB/((2 * math.pi**(variables.d-1))) * ln(Abs(1 - math.exp((-variables.omegabs(variables.k))/variables.T)))) == - quad((diff(variables.t))/variables.t * (2 * math.pi* variables.t)**(-variables.d/2) * summation(math.exp(-(variables.ma(2) * variables.t)/2 - (variables.ra(2))/(2 * variables.Ta(2) * variables.t))), 0, math.inf)
variables.Ff(variables.T) == - integrate((diff(variables.tau) * diff(variables.taul))/(2 * variables.taub(2)) * (4 * math.pi**2 * variables.SlopeRegge * variables.taub(2))**(-13) * Abs(variables.etaf(variables.tau))**(-48) * summation(math.exp(-(variables.ra(2))/(4 * math.pi * variables.Ta(2) * variables.SlopeRegge * variables.taub(2))), (variables.r, 1, math.inf)), variables.R)
variables.R 
variables.taub(1) <= 1/2 
Abs(variables.taub(2)) > 0
variables.etaf(variables.i * variables.taub(2)) == variables.etaf(variables.i/(variables.taub(2))) * variables.tauab(-1/2, 2) == math.exp(-math.pi/(12 * variables.taub(2))) # as taub(2) -> 0
variables.Fftilde(variables.T) == variables.Ff(variables.T) + variables.rhob(0)
1/variables.T * variables.Fftilde(variables.T) == variables.T/(4 * variables.Tab(2, variables.H)) * variables.Fftilde(4 * variables.Tab(2, variables.H)/variables.T)
variables.Fftilde(variables.T) == (variables.Ta(2))/(4 * variables.Tab(2, variables.H)) * variables.rhob(0)
#323.1, 9.9
variables.Gbf(variables.mu * variables.v, variables.x) == variables.etab(variables.mu * variables.v)
variables.Bbf(variables.mu * variables.v, variables.x) == 0
variables.Phif(variables.x) == variables.Vb(variables.mu) * variables.xas(variables.mu)
variables.Vb(variables.mu) == variables.deltaba(variables.mu, 1) * ((26 - variables.D)/(6 * variables.SlopeRegge))**(1/2)
-diff(diff(variables.Tf(variables.x), variables.mu), 1, variables.mu) + 2 * variables.Va(variables.mu) * diff(variables.Tf(variables.x), variables.mu) - 4/(variables.SlopeRegge) * variables.Tf(variables.x) == 0
variables.Tf(variables.x) == math.exp(np.cross(variables.q, variables.x))
(variables.q - variables.V)**2 == (2 - variables.D)/(6 * variables.SlopeRegge) 
variables.qb(1) == variables.alphab(variables.posneg) == ((26 - variables.D)/(6 * variables.SlopeRegge))**(1/2) + variables.posneg * ((2 - variables.D)/(6 * variables.SlopeRegge))**(1/2) 
variables.Sb(variables.sigma) == 1/(4 * math.pi * variables.SlopeRegge) * integrate(variables.d2sigma * variables.g**(1/2) * variables.brackets(variables.ga(variables.a * variables.b) * variables.etab(variables.mu * variables.v) * diff(variables.xa(variables.mu), variables.a) * diff(variables.xa(variables.v), variables.b) + variables.SlopeRegge * variables.R * variables.Vb(1) * variables.xa(1) + variables.Tb(0) * math.exp(variables.alpha_xa1)), variables.M)
variables.S == 1/(4 * math.pi * variables.SlopeRegge) * integrate(variables.d2sigma * variables.ga(1/2) * (variables.ga(variables.a * variables.b) * diff(variables.xa(variables.mu), variables.a) * diff(variables.xb(variables.mu), variables.b) + variables.mu), variables.M)
variables.gbf(variables.a * variables.b, variables.sigma) == math.exp(2 * variables.phif(variables.sigma)) * variables.hat(variables.gbf(variables.a * variables.b, variables.sigma))
variables.ga(1/2) * variables.R == variables.hat(variables.ga(1/2)) *(variables.hat(variables.R) - 2 * variables.hat(variables.Deltaa(2)) * variables.phi)
variables.S == 1/(4 * math.pi * variables.SlopeRegge) * integrate(variables.d2sigma * variables.hat(variables.ga(1/2)) * variables.brackets(variables.hat(variables.ga(variables.a * variables.b)) * diff(variables.xa(variables.mu), variables.a) * diff(variables.xb(variables.mu), variables.b) + variables.mu * math.exp(2 * variables.phi) + (13 * variables.SlopeRegge)/3 * (variables.hat(variables.ga(variables.a * variables.b)) * diff(variables.phi, variables.a) * diff(variables.phi, variables.b) + variables.hat(variables.R) * variables.phi)), variables.M)
variables.hat(variables.gbf(variables.a * variables.b, variables.sigma)) == math.exp(2 * variables.omegasf(variables.sigma)) * variables.hat(variables.gbf(variables.a * variables.b, variables.sigma))
variables.phif(variables.sigma) == variables.phif(variables.sigma) - variables.omegasf(variables.sigma)

# 2.136.1
variables.Xabf(9, variables.R, variables.zl) == - variables.Xabf(9, variables.R, variables.zl)
variables.psiafptilde(9, variables.zl) == - variables.psiaftilde(9, variables.zl)
variables.VObfp(variables.alpha, variables.z) == variables.VObf(variables.alpha, variables.z)
variables.VObfptilde(variables.alpha, variables.zl) == variables.Betaab(9, variables.alpha * variables.Beta) * variables.VObftilde(variables.Beta, variables.zl)
variables.Cb(9) == variables.C 
(variables.Cb(variables.mu), variables.Cb(variables.mu * variables.v, 9)) == (variables.Cb(variables.mu, 9), variables.Cb(variables.mu * variables.v))
variables.Cb(variables.mu * variables.v * variables.lamda) == variables.Cb(variables.mu * variables.v, variables.lamda * 9)
np.prod(variables.Betaa(variables.m), variables.m) 
variables.Betaa(variables.m) * variables.Betaa(variables.n) == math.exp(math.pi * variables.i * variables.FBtilde) * variables.Betaa(variables.n) * variables.Betaa(variables.m)
variables.psiab(variables.mu, -1/2) * variables.diracb(variables.k, variables.NS) 
variables.psiab(9, -1/2) * variables.diracb(variables.k, variables.NS) 
variables.diracb((variables.alpha, variables.k), variables.R)
variables.Qbp(variables.alpha) + variables.holderb(variables.Betaa(9) * variables.Qptilde, variables.alpha)
1/(2 * math.pi * variables.SlopeRegge) * integrate(diff(variables.s) * diff(variables.Xpa(9), variables.n), diff(variables.M, variables.z))
integrate(diff(variables.s) * variables.VObp(variables.alpha), diff(variables.M, variables.z)) == - integrate(diff(variables.s) * variables.holderb(variables.Betaa(9) * variables.VOptilde, variables.alpha), diff(variables.M, variables.z))
variables.Qbp(variables.alpha) + variables.holderb(variables.Betaa(variables.perp) * variables.Qptilde, variables.alpha)
variables.Betaa(variables.perp) == np.prod(variables.Betaa(variables.m), variables.m)
integrate(variables.Cb(variables.p + 1))
variables.Fb(2) == variables.d * (variables.Cb(1)) 
WedgeProduct(variables.d, variables.starf(variables.d * variables.Cb(1))) == 1
variables.starf(variables.Fb(2)) == variables.holderb(variables.starf(variables.F), 9) == variables.d * variables.Cb(7) 
WedgeProduct(variables.d, variables.starf(variables.d * variables.Cb(7))) == 0
variables.bracketc(variables.Qb(variables.alpha), variables.Qlb(variables.Beta)) == - 2 * variables.brackets(variables.Pb(variables.M) + (2 * math.pi * variables.SlopeRegge)**(-1) * variables.Qab(variables.NS, variables.M)) * variables.Rhoab(variables.M, variables.alpha * variables.Beta) 
variables.bracketc(variables.Qbtilde(variables.alpha), variables.Qlbtilde(variables.Beta)) == - 2 * variables.brackets(variables.Pb(variables.M) - (2 * math.pi * variables.SlopeRegge)**(-1) * variables.Qab(variables.NS, variables.M)) * variables.Rhoab(variables.M, variables.alpha * variables.Beta)
