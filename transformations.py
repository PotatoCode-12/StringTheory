import math
from scipy.integrate import quad
import numpy as np 
from sympy import symbols, diff,integrate, Function, Derivative, summation, oo
import variables
import constants
import states

t = variables.t
def NambuGotoSymetricalTransformations():
    variables.tetrad(t)*diff(t, t)= diff(variables.tetrad(t), diff(t,t)) #with worldine reparameterization invariance
    variables.n**2 == -diff(variables.xa(variables.mu)*variables.xb(variables.mu), t)/variables.m #for transformation to motion equation

variables.xap,(variables.mup, variables.taup, variables.sigmap) == variables.xa(variables.mu, variables.t, variables.sigma) #two dimesnional coordinate invariance translation
if variables.sigmap == variables.l - variables.sigma:
    if variables.taup == variables.tau:
        if states.string() == "open" :
            variables.bomega * variables.TensionStringab(variables.i, variables.n) *variables.bomega**(-1) == (-1)**variables.n * variables.TensionString(variables.i, variables.n)
        elif states.string() == "closed":
            variables.bomega * variables.TensionStringab(variables.i, variables.n) *variables.bomega**(-1) == variables.TensionStringabtilde(variables.i, variables.n)
            variables.bomega*variables.TensionStringabtilde(variables.i, variables.n) *variables.bomega**(-1) == variables.TensionStringab(variables.i, variables.n)
        variables.bomega*variables.dirac(variables.N, variables.k) == (-1)**(variables.N) * variables.dirac(variables.N, variables.k)
        variables.bomega*variables.dirac(variables.N, variables.Ntilde, variables.k) == variables.dirac(variables.Ntilde, variables.N, variables.k)

diff_resulta = variables.sigma
for _ in range(int(variables.d)):
    diff_resulta = diff(diff_resulta, variables.d)
a = diff_resulta
variables.PIC(diff(variables.phip)) * math.exp(-variables.S*variables.PIC(variables.phip)) == variables.PIC(diff(variables.phi)) * math.exp(-variables.S*variables.PIC(variables.phi)) * variables.PIC((1 + variables.i*variables.epsilon)/(2*math.pi)* integrate(a*variables.g**(1/2)*variables.ja(variables.a, variables.sigma)* diff(variables.rho(variables.sigma), variables.a) + variables.O(variables.epsilon**2))) # 41

diff_resultc = variables.xb(variables.mu)
for _ in range(int(variables.c)):
    diff_resultc = diff(diff_resultc)
c = diff_resultc
#for world sheet translation example page 43
variables.jb(variables.a) == variables.i*variables.va(variables.b)*variables.Tb(variables.a, variables.b)
variables.Tb(variables.a, variables.b) == -1/(variables.SlopeRegge) *variables.point(diff(variables.xa(variables.mu), variables.a)*diff(variables.xb(variables.mu), variables.b) - 1/2(variables.EquMotion(variables.a, variables.b)*diff(variables.xa(variables.mu), variables.mu)*c))
variables.delta*variables.xa(variables.mu) == -variables.epsilon* variables.vf(variables.z) * diff(variables.xa(variables.mu), variables.z) - variables.epsilon*variables.vf(variables.z)*variables.trans(diff(variables.xa(variables.mu), variables.zl)) #44
variables.xap(variables.mu, variables.zp, variables.zpl) == variables.zp == variables.f(variables.z) #45


#conformal invariance 46
n = symbols('n') # summation stepping
diff_resultb = variables.vf(variables.z)
for _ in range(int(n)):
    diff_resultb = diff(diff_resultb, variables.z)
b = diff_resultb
diff_resultd = variables.vf(variables.z)
for _ in range(int(n)):
    diff_resultd = diff(diff_resultd, variables.zl)
d = diff_resultd
variables.delta*variables.AO(variables.z, variables.zl) == -variables.epsilon * summation(1/(math.factorial(n)) * variables.brackets(b*variables.AOa(n, variables.z, variables.zl)+d*variables.trans(variables.AOatilde(n, variables.z, variables.zl))))
variables.AOp(variables.zp, variables.zpl) == variables.zetaa(-variables.h)*variables.zetaal(-variables.htilde) * variables.AO(variables.z, variables.zl) #eigenstates
variables.OOp(variables.zp, variables.zpl) == diff(variables.zp, variables.z)**(-variables.h)*diff(variables.zpl, variables.zl)**(-variables.htilde)*variables.OO(variables.z, variables.zl) #general conformal transformation

variables.epsilon**(-1)*variables.delta*variables.Tf(variables.z) == -variables.D/12*diff(variables.vf(variables.z), variables.z, 3) - 2*diff(variables.vf(variables.z), variables.z)*variables.Tf(variables.z) - variables.vf(variables.z)*diff(variables.Tf(variables.z), variables.z) #law established pg 48
variables.epsilon**(-1)*variables.delta*variables.Tf(variables.z) == -constants.c/12 * diff(variables.vf(variables.z), variables.z, 3) - 2*diff(variables.vf(variables.z), variables.z) * variables.Tf(variables.z) - variables.vf(variables.z) * diff(variables.Tf(variables.z), variables.z) # transformation to general CFT
(diff(variables.zp, variables.z))**2*variables.Tfp(variables.zp) == variables.Tf(variables.z) - constants.c/12*variables.bracketc(variables.zp, variables.z)
variables.bracketc(variables.f, variables.z) == (2*diff(variables.f, variables.z, 3)*diff(variables.f, variables.z) - 3*diff(variables.f, variables.z, 2)*diff(variables.f, variables.z, 2))/(2*diff(variables.f, variables.z)*diff(variables.f, variables.z))

variables.delta*variables.xa(variables.mu) == variables.epsilon*variables.v*diff(variables.xa(variables.mu), variables.z) - variables.epsilon*variables.v*variables.trans*diff(variables.xa(variables.mu), variables.zl) - variables.epsilon/2 * variables.SlopeRegge*variables.Va(variables.mu) * variables.brackets(diff(variables.v, variables.z) + variables.trans(diff(variables.v, variables.z))) # for Ward identity and central charge defined pg 49

variables.epsilon**(-1)*variables.delta*variables.j == -variables.v*diff(variables.j) - variables.j*diff(variables.v) + (2*variables.lamda - 1)/2 * variables.diff(variables.v, 2) # transformation law page 51, check rules
diff(variables.zp, variables.z) * variables.jb(variables.zp, variables.zp) == variables.jb(variables.z, variables.z) + (2*variables.lamda - 1)/2 * diff(variables.zp, variables.z, 2)/diff(variables.zp, variables.z)
variables.psi == 2**(-1/2)* (variables.psib(variables.D1) + variables.i*variables.psib(variables.D2)) #for earlier conditions
variables.psil == 2**(-1/2)* (variables.psib(variables.D1) - variables.i*variables.psib(variables.D2))
variables.S == 1/(4*math.pi) * integrate(variables.d2z * (variables.psib(variables.D1)) * diff(variables.psib(variables.D1), variables.zl) + variables.psib(variables.D2) * diff(variables.psib(variables.D2), variables.zl))
variables.TensionString == -1/2*variables.psib(variables.D1)*diff(variables.psib(variables.D1), variables.z) - 1/2*variables.psib(variables.D2) * diff(variables.psi(variables.D2), variables.z)

variables.gab(variables.zeta, variables.a*variables.b, variables.sigmap) == math.exp(variables.brackets(2*variables.omegasf(variables.sigma))) * diff(variables.sigmaa(variables.c))/(diff(variables.sigmaap(variables.a))) * diff(variables.sigmaa(variables.d))/(diff(variables.sigmaap(variables.b))) * variables.gbf(variables.c*variables.d, variables.sigma)
1 == variables.Deltabf(variables.FP, variables.g) * integrate((diff(variables.zeta)) * variables.delta * (variables.g - variables.hat(variables.ga(variables.zeta))))
variables.Zf(variables.bracket(variables.hat(variables.g))) == integrate((variables.brackets(diff(variables.zeta) * diff(variables.X) * diff(variables.g)))/(variables.Vb(variables.dw)) * variables.Deltabf(variables.FP, variables.g) * variables.delta(variables.g-variables.hat(variables.ga(variables.zeta))) * math.exp(-variables.S * variables.brackets(variables.X, variables.g)))
variables.Zf(variables.hat(variables.g)) == integrate(variables.brackets(diff(variables.zeta) * diff(variables.xa(variables.zeta)))/variables.Vb(variables.dw) * variables.Deltabf(variables.FP, variables.hat(variables.ga(variables.zeta))) * math.exp(-variables.S * variables.brackets(variables.xa(variables.zeta), variables.hat(variables.ga(variables.zeta)))))
variables.Zf(variables.brackets(variables.hat(variables.g))) == integrate(variables.brackets(diff(variables.zeta) * diff(variables.X))/(variables.Vb(variables.dw)) * variables.Deltabf(variables.FP, variables.hat(variables.g)) * math.exp(-variables.S * variables.brackets(variables.X, variables.hat(variables.g))))
variables.Zf(variables.hat(variables.g)) == integrate(variables.brackets(diff(variables.X)) * variables.Deltabf(variables.FP, variables.hat(variables.g)) * math.exp(-variables.S * variables.brackets(variables.X, variables.g)))
variables.delta*variables.gb(variables.a*variables.b) == 2*variables.delta*variables.omegas*variables.gb(variables.a*variables.b) - variables.Deltab(variables.a) * variables.delta * variables.sigmab(variables.b) - variables.Deltab(variables.b) * variables.delta * variables.sigmab(variables.a) == (2*variables.delta*variables.omegas - variables.Deltab(variables.c) * variables.delta * variables.sigmaa(variables.c)) * variables.gb(variables.a*variables.b) - 2* variables.holderb(variables.Pb(1) * variables.delta*variables.sigma, variables.a*variables.b)
variables.holderb(variables.Pb(1) * variables.delta*variables.sigma, variables.a*variables.b) == 1/2 * (variables.Deltab(variables.a) * variables.delta*variables.sigmab(variables.b) + variables.Deltab(variables.b) * variables.delta * variables.sigmab(variables.a) - variables.gb(variables.a*variables.b) * variables.Deltab(variables.c) * variables.delta * variables.sigmaa(variables.c))
variables.Deltabf(variables.FP, variables.hat(variables.g) **(-1)) == integrate(variables.brackets(diff(variables.delta *variables.omegas) * diff(variables.delta*variables.sigma)) * variables.delta(variables.brackets(-(2*variables.delta*variables.omegas - variables.hat(variables.Delta) * variables.delta* variables.sigma) * variables.hat(variables.g) + 2*variables.hat(variables.P1) * variables.delta * variables.sigma))) == integrate(variables.brackets(diff(variables.delta*variables.omegas) * diff(variables.Beta) * diff(variables.delta * variables.sigma)) * math.exp(variables.bracketc(2*math.pi*variables.i * integrate(variables.d2sigma * variables.hat(variables.g)**(1/2) * variables.Betaa(variables.a*variables.b) * variables.bracketsb(-(2*variables.delta*variables.omegas - variables.hat(variables.Delta) * variables.delta * variables.sigma) * variables.hat(variables.g) + 2*variables.hat(variables.P1) * variables.delta * variables.sigma), variables.a*variables.b)))) == integrate(variables.brackets(diff(variables.Betap) * diff(variables.delta * variables.sigma)) * math.exp(4*math.pi * variables.i * integrate(variables.d2sigma * (variables.hat(variables.g))**(1/2) * variables.Betapa(variables.a*variables.b) * variables.holderb(variables.hat(variables.P1) * variables.delta * variables.sigma, variables.a*variables.b)))) # 88.18
#for grassman ghost field 88.19
variables.Deltabf(variables.FP, variables.hat(variables.g)) == integrate(variables.brackets(diff(variables.b) * diff(variables.c)) * math.exp(-variables.Sb(variables.g)))
variables.Sb(variables.g) == 1/(2*math.pi) * integrate(variables.d2sigma * variables.hat(variables.g) ** (1/2) * variables.bb(variables.a*variables.b) * variables.hat(variables.Deltaa(variables.a)) * variables.ca(variables.b)) == 1/(2*math.pi) * integrate(variables.d2sigma * variables.hat(variables.g) ** (1/2) * variables.bb(variables.a*variables.b) * variables.holdera(variables.hat(variables.P1) * variables.c, variables.a*variables.b)) #normalization
variables.Zf(variables.brackets(variables.hat(variables.g))) == integrate(variables.brackets(diff(variables.X) * diff(variables.b) * diff(variables.c)) * math.exp(-variables.Sb(variables.x) - variables.Sb(variables.g)))
variables.Zf(variables.brackets(variables.hat(variables.g))) == ((variables.hat(variables.Delta**2)).det())**(-variables.D/2) * (variables.hat(variables.P1)).det() 
# gauge
variables.Sb(variables.g) == 1/(2*math.pi) * integrate(variables.d2z * (variables.bb(variables.z*variables.z) * variables.Deltab(variables.zl) * variables.ca(variables.z) + variables.bb(variables.zl*variables.zl*variables.Deltab(variables.z) * variables.ca(variables.zl)))) == 1/(2*math.pi * integrate(variables.d2z * (variables.bb(variables.z * variables.z) * diff(variables.ca(variables.z), variables.zl) + variables.bb(variables.zl* variables.zl) * diff(variables.ca(variables.zl), variables.z))))
#vanishing coomponent
variables.nb(variables.a) * variables.delta * variables.sigmaa(variables.a) == 0
variables.nb(variables.a) * variables.ca(variables.a) == 0
0 == integrate(diff(variables.s) * variables.na(variables.a) * variables.bb(variables.a*variables.b) *variables.delta * variables.ca(variables.b), 0, 2*math.pi)# contour integral
variables.nb(variables.a) * variables.tb(variables.b) * variables.ba(variables.a* variables.b) == 0

variables.epsilon**(-1) * variables.delta * variables.Tbf(variables.z* variables.z , variables.z) == -variables.c/12 * diff(variables.vaf(variables.z, variables.z), 3, variables.z) - 2*(diff(variables.vaf(variables.z, variables.z), variables.z) * variables.Tbf(variables.z*variables.z, variables.z)) - variables.vaf(variables.z, variables.z) * diff(variables.Tbf(variables.z * variables.z, variables.z), variables.z) # Weyl transformation 93.13
variables.deltab(variables.W) * variables.Tb(variables.z, variables.z) == -variables.c/6 * diff(variables.delta*variables.omegas, 2, variables.z) #for leading around a flat space 93.14
variables.c == -12 * variables.ab(variables.D1)
variables.Tab(variables.a, variables.a) == -variables.c/12 * variables.R

#126.1
def gaugeTransformationSatisfactions():
    variables.brackets(variables.deltab(variables.alpha), variables.deltab(variables.Beta)) == variables.fab(variables.gamma, variables.alpha*variables.Beta) * variables.deltab(variables.gamma)
    variables.Faf(variables.A, variables.phi) == 0 # gauge condition fixing
    integrate(variables.brackets(diff(variables.phib(variables.i)))/variables.Vb(variables.gauge) * math.exp(-variables.Sb(-variables.D1))) == integrate(variables.brackets(diff(variables.phib(variables.i)) * diff(variables.Bb(variables.A)) * diff(variables.bb(variables.A)) * diff(variables.ca(variables.alpha)))) * math.exp(-variables.Sb(variables.D1) - variables.Sb(variables.D2) - variables.Sb(variables.D3))
    variables.Sb(variables.D2) == -variables.i*variables.Bb(variables.A) * variables.Faf(variables.A, variables.phi) # gauge fixed action
    variables.Sb(variables.D3) == variables.bb(variables.A) * variables.ca(variables.alpha) * variables.deltab(variables.alpha) * variables.Faf(variables.A, variables.phi)
    variables.deltab(variables.B) * variables.phib(variables.i) == -variables.i*variables.epsilon*variables.ca(variables.alpha) * variables.deltab(variables.alpha) * variables.deltab(variables.alpha) * variables.phib(variables.i)
    variables.deltab(variables.B) * variables.Bb(variables.A) == 0
    variables.deltab(variables.B) * variables.bb(variables.A) == variables.epsilon * variables.Bb(variables.A)
    variables.deltab(variables.B) * variables.ca(variables.alpha) == variables.i/2 * variables.epsilon * variables.fpb(variables.alpha, variables.Beta * variables.gamma) * variables.ca(variables.Beta) * variables.ca(variables.gamma)
    # 127.7
    variables.deltab(variables.B) * (variables.bb(variables.A) * variables.Fa(variables.A))  == variables.i*variables.epsilon * (variables.Sb(variables.D2) + variables.Sb(variables.D3))
    variables.epsilon * variables.delta * variables.bracket( variables.f, variables.i) == variables.i * variables.bracket(variables.f, variables.deltab(variables.B) * (variables.bb(variables.A) * variables.delta * variables.Fa(variables.A)), variables.i) == -variables.epsilon * variables.bracket(variables.f, variables.bracketc(variables.Qb(variables.B), variables.bb(variables.A), variables.delta * (variables.Fa(variables.A))), variables.i) # ghost action product
    variables.bracket(variables.psi, variables.bracketc(variables.Qb(variables.B), variables.bb(variables.A) * variables.delta * variables.Fa(variables.A)), variables.psip) == 0
    variables.Qb(variables.B) * variables.dirac(variables.psi) == variables.Qb(variables.B) * variables.dirac(variables.psip) == 0
    variables.epsilona(-1) * variables.deltab(variables.B) * (variables.bb(variables.A) * variables.Bb(variables.B) * variables.Ma(variables.A*variables.B)) == -variables.Bb(variables.A) * variables.Bb(variables.B) * variables.Ma(variables.A8variables.B)
    #hamiltonian
    0 == variables.brackets(variables.Qb(variables.B), variables.bracketc(variables.Qb(variables.B), variables.bb(variables.A) * variables.delta*variables.Fa(variables.A))) == variables.Qab(2, variables.B) * variables.bb(variables.A) * variables.delta * variables.Fa(variables.A) - variables.Qb(variables.B) * variables.bb(variables.A) * variables.delta* variables.Fa(variables.A) * variables.Qb(variables.B) + variables.Qb(variables.B) * variables.bb(variables.A) * variables.delta * variables.Fa(variables.A) * variables.Qb(variables.B) - variables.bb(variables.A) * variables.delta * variables.Fa(variables.A) * variables.Qab(2, variables.B) == variables.brackets(variables.Qab(2, variables.B), variables.bb(variables.A) * variables.delta * variables.Fa(variables.A))
    variables.Qab(2, variables.B) == 0
    variables.deltab(variables.B) * (variables.deltabp(variables.B) * variables.ca(variables.alpha)) == -1/2 * variables.epsilon * variables.epsilonp * variables.fab(variables.alpha, variables.Beta * variables.gamma) * variables.fab(variables.gamma, variables.delta * variables.epsilon) * variables.ca(variables.Beta) * variables.ca(variables.delta) * variables.ca(variables.epsilon) == 0
    variables.bracket(variables.psi, (variables.Qb(variables.B) * variables.dirac(variables.chi))) == variables.bracket(variables.psi, variables.Qb(variables.B), variables.chi) == 0
    variables.dirac(variables.psip) == variables.dirac(variables.psi) + variables.Qb(variables.B) * variables.dirac(variables.chi)
    variables.hilbertb(variables.BRST) == variables.Hb(variables.closed)/(variables.Hb(variables.exact))

#131.1
def BRSTQuantisation():
    variables.deltab(variables.B) * variables.xa(variables.mu) == variables.i * variables.epsilon *(variables.c * diff(1, variables.z) + variables.ctilde * diff(1, variables.zl)) * variables.xa(variables.mu)
    variables.deltab(variables.B) * variables.b == variables.i * variables.epsilon *(variables.Ta(variables.X) + variables.Ta(variables.g))
    variables.deltab(variables.B)  * variables.btilde == variables.i* variables.epsilon * (variables.Tatilde(variables.X) + variables.Tatilde(variables.g))
    variables.deltab(variables.B) * variables.c == variables.i * variables.epsilon * variables.c * diff(variables.c, variables.z)
    variables.deltab(variables.B) * variables.ctilde == variables.i * variables.epsilon * variables.ctilde * diff(variables.ctilde, variables.zl)
    variables.i/(4*math.pi) * integrate(variables.d2sigma * variables.g**(1/2) * variables.Ba(variables.a*variables.b) (variables.deltab(variables.a*variables.b) - variables.gb(variables.a*variables.b)))
    variables.jb(variables.B) == variables.c * variables.Ta(variables.m) + 1/2 * variables.point(variables.c * variables.ta(variables.g)) + 3/2 * diff(variables.c, variables.z, 2) == variables.c * variables.T + variables.point(variables.b * variables.c * diff(variables.c, variables.z)) + 3/2 * diff(variables.c, variables.z, 2)
    #check ghost fields
    variables.Qb(variables.B) == 1/(2 * math.pi * variables.i) * integrate(variables.diff(variables.z) * variables.jb(variables.B) - diff(variables.zl) * variables.jbtilde(variables.B), 0, 2*math.pi) # contour
    variables.bracketc(variables.Qb(variables.b), variables.bb(variables.m)) == variables.Lab(variables.m, variables.m) + variables.Lab(variables.g, variables.m)
    variables.Qb(variables.B) == summation(variables.cb(variables.n) * variables.Lab(variables.m, -variables.n) + variables.cbtilde(variables.n) * variables.Labtilde(variables.m, -variables.n), (variables.n, -math.inf, math.inf)) + summation((variables.m - variables.n)/2 * variables.point(variables.cb(variables.m) * variables.cb(variables.n) * variables.bb(-variables.m-variables.n) + variables.cbtilde(variables.m) * variables.cbtilde(variables.n) * variables.bbtilde(-variables.m-variables.n)) + variables.aa(variables.B) * (variables.cb(0) + variables.cbtilde(0)), ((variables.m, variables.n), -math.inf, math.inf)) # double check 'point'
    variables.bracketc(variables.Qb(variables.B), variables.bb(0)) == variables.Lab(variables.m, 0) + variables.Lab(variables.g, 0)
    if variables.ca(variables.m) == 26:
        variables.bracketc(variables.Qb(variables.B), variables.Qb (variables.B)) == 0
    variables.brackets(variables.Gb (variables.I), variables.Gb(variables.J)) == variables.i * variables.Gab(variables.K, variables.I * variables.J) * variables.Gb(variables.K)
    variables.bracketc(variables.ca(variables.I), variables.bb(variables.J)) == variables.deltab(variables.I, variables.J)
    variables.bracketc(variables.ca(variables.I), variables.ca(variables.J)) == variables.bracketc(variables.bb(variables.I), variables.bb(variables.J)) == 0
    variables.Qb(variables.B) == variables.ca(variables.I) * variables.Gab(variables.m, variables.I) - variables.i/2 * variables.gab(variables.K, variables.I* variables.J) * variables.ca(variables.I) * variables.ca(variables.J) * variables.bb(variables.K) == variables.ca(variables.I) *(variables.Gab(variables.m, variables.I) + 1/2 * variables.Gab(variables.g, variables.I))
    variables.Gab(variables.g, variables.I) == -variables.i * variables.gab(variables.K, variables.I * variables.J) * variables.ca(variables.J) * variables.bb(variables.K)
    variables.Qab(2, variables.B) == 1/2 * variables.bracketc(variables.Qb(variables.B), variables.Qb(variables.B)) == -1/2 * variables.gab(variables.K, variables.I * variables.J) * variables.gab(variables.M, variables.K, variables.Ls) * variables.ca(variables.I) * variables.ca(variables.J) * variables.ca(variables.L) * variables.bb(variables.M) == 0
    variables.dagger(variables.alphaab(variables.mu, variables.m)) == variables.alphaab(variables.mu, -variables.m)
    variables.dagger(variables.alphaabtilde(variables.mu, variables.m)) == variables.alphaabtilde(variables.mu, -variables.m)
    variables.dagger(variables.bb(variables.m)) == variables.bb(-variables.m)
    variables.dagger(variables.bbtilde(variables.m)) == variables.bbtilde(-variables.m)
    variables.dagger(variables.cb(variables.m)) == variables.cb(-variables.m)
    variables.dagger(variables.cbtilde(variables.m)) == variables.cbtilde(-variables.m)
    # for open string
    variables.bracket(variables.dirac(0, variables.k), variables.cb(0), variables.dirac(0, variables.kp)) == (2*math.pi)**26 * variables.delta**26 * (variables.k-variables.kp)
    # closed string
    variables.bracket(variables.dirac(0, variables.k), variables.cbtilde(0) * variables.cb(0), variables.dirac(0, variables.kp)) == variables.i *(2*math.pi)**26 * variables.delta**26 * (variables.k-variables.kp)
    variables.bb(0) * variables.dirac(variables.psi) == 0
    variables.Lb(0) * variables.dirac(variables.psi) == variables.bracketc(variables.Qb(variables.B), variables.bb(0)) * variables.dirac(variables.psi) == 0
    variables.Lb(0) == variables.alphap * (variables.pa(variables.mu) * variables.pb(variables.mu) + variables.m**2)
    variables.alphap * variables.m**2 == summation(variables.n * (variables.Nb(variables.b *variables. n) + variables.Nb(variables.c * variables.n) + summation(variables.Nb(variables.mu, variables.n), (variables.mu, 0, 25))) - 1, (variables.n, 1, math.inf))
    -variables.k**2 == -1/(variables.SlopeRegge) # 135.23
    variables.Qb(variables.B) * variables.dirac(0, variables.kb) == 0
    variables.dirac(variables.psib(1)) == (variables.e * variables.alphab(-1) + variables.Beta * variables.bb(-1) + variables.gamma * variables.cb(-1)) * variables.dirac(0, variables.kb)
    -variables.k**2 == 0
    # check norm state 135.26
    0 == variables.Qb(variables.B) * variables.dirac(variables.psib(1)) == (2 * variables.SlopeRegge)**(1/2) * (variables.cb(-1) * variables.k * variables.alphab(1) + variables.cb(1) * variables.k * variables.alphab(-1)) * variables.dirac(variables.psib(1)) == (2 * variables.SlopeRegge)**(1/2) * (variables.k * variables.e * variables.cb(-1) + variables.Beta* variables.k * variables.alphab(-1)) * variables.dirac(0, variables.kb)
    variables.Qb(variables.B) * variables.dirac(variables.chi) == (2 * variables.SlopeRegge)**(1/2) * (variables.k * variables.ep * variables.cb(-1) + variables.Betap * variables.k * variables.alphab(-1)) * variables.dirac(0, variables.kb)
    variables.bb(0) * variables.dirac(variables.psi) == variables.bbtilde(0) * variables.dirac(variables.psi) == 0
    variables.Lb(0) * variables.dirac(variables.psi) == variables.Lbtilde(0) * variables.dirac(variables.psi) == 0
    variables.Lb(0) == variables.SlopeRegge/4 * (variables.p**2 + variables.m**2)
    variables.Lbtilde(0) == variables.SlopeRegge/4 * (variables.p**2 + variables.mtilde**2)
    variables.SlopeRegge/4 * variables.m**2 == summation(variables.n * (variables.Nb(variables.b * variables.n) + variables.Nb(variables.c*variables.n) + summation(variables.Nb(variables.mu * variables.n), (variables.mu, 0, 25))) - 1, (variables.n, 1, math.inf))
    variables.SlopeRegge/4 * variables.mtilde**2 == summation(variables.n * (variables.Nbtilde(variables.b * variables.n) + variables.Nbtilde(variables.c*variables.n) + summation(variables.Nbtilde(variables.mu * variables.n), (variables.mu, 0, 25))) - 1, (variables.n, 1, math.inf))

#2.259.1
variables.deltab(variables.s) * variables.z == variables.epsilon * variables.z 
variables.deltab(variables.s) * variables.gb(variables.a * variables.b) == 2 * variables.epsilon * variables.gb(variables.a * variables.b)
-variables.epsilon/(2 * math.pi) * integrate(variables.d2sigma * variables.Tabf(variables.a, variables.a , variables.sigma))
variables.Tab(variables.a, variables.a) == diff(variables.KOa(variables.a), variables.a)
variables.epsilona(-1) * variables.deltab(variables.s) * variables.bracket(np.prod(variables.AObf(variables.ib(variables.m), variables.sigmab(variables.m)), variables.m)) == -1/(2 * math.pi) * integrate(variables.d2sigma * variables.bracket(variables.Tabf(variables.a, variables.a, variables.sigma) * np.prod(variables.AObf(variables.ib(variables.m), variables.sigmab(variables.m)), variables.m))) - summation(variables.Deltaba(variables.ib(variables.n), variables.j) * variables.bracket(variables.AObf(variables.ib(variables.n), variables.j) * np.prod(variables.AObf(variables.ib(variables.m), variables.sigmab(variables.m)), variables.m != variables.n)), variables.n)
variables.epsilona(-1) * variables.deltab(variables.s) * variables.AObf(variables.i, variables.sigma) == - variables.Deltaba(variables.i, variables.j) * variables.AObf(variables.j, variables.sigma)
integrate(variables.da(variables.d) * variables.sigma * variables.Tab(variables.a, variables.a)) == - 2 * math.pi * summation(variables.ga(variables.i) * integrate(variables.da(variables.d) * variables.sigma * variables.AObf(variables.i, variables.sigma)), (variables.i, variables.i != 0))
variables.epsilona(-1) * variables.deltab(variables.s) * variables.bracket(np.prod(variables.AObf(variables.ib(variables.m), variables.sigmab(variables.m)), variables.m)) == - summation((variables.Betaaf(variables.i, variables.g) * (diff(1, variables.z))/(diff(variables.ga(variables.i), variables.z)) * variables.bracket(np.prod(variables.AObf(variables.ib(variables.m), variables.sigmab(variables.m)), variables.m))), (variables.i, variables.i != 0)) - summation(variables.Deltaba(variables.ib(variables.n), variables.j) * variables.bracket(variables.AObf(variables.j, variables.sigmab(variables.n)) * np.prod(variables.AObf(variables.ib(variables.m), variables.sigmab(variables.m)), (variables.m != variables.n))), variables.n)
variables.Ff(variables.ra(2)) == variables.za(4) * variables.bracket(variables.Tbf(variables.z * variables.z, (variables.z, variables.zl)) * variables.Tbf(variables.z * variables.z, (0, 0)))
variables.Gf(variables.ra(2)) == 4 * variables.za(3) * variables.zl * variables.bracket(variables.Tbf(variables.z * variables.z, (variables.z, variables.zl)) * variables.Tbf(variables.z * variables.zl, (0, 0)))
variables.Hf(variables.ra(2)) == 16 * variables.za(2) * variables.zla(2) * variables.bracket(variables.Tbf(variables.z * variables.zl, (variables.z, variables.zl)) * variables.Tbf(variables.z * variables.zl, (0, 0)))
4 * diff(variables.F, variables.tau) + diff(variables.G, variables.tau) - 3 * variables.G == 0 
4 * diff(variables.G, variables.tau) + 4 * variables.G + diff(variables.H, variables.tau) - 2 * variables.H == 0
variables.C == 2 * variables.F - variables.G - 3/8 * variables.H
diff(variables.C, variables.tau) == -3/4 
variables.S == variables.Sb(0) + variables.lamdaa(variables.i) * integrate(variables.d2z * variables.OOb(variables.i))
-variables.Tbf(variables.z * variables.z, (variables.z, variables.zl)) * variables.lamdaa(variables.i) * integrate(variables.d2(variables.w) * variables.OObf(variables.i, (variables.w, variables.wl)))
diff(variables.Tbf(variables.z * variables.z, variables.z), variables.zl) * variables.OObf(variables.i, (variables.w, variables.wl)) == diff(variables.brackets((variables.z - variables.w)**(-2) * variables.hb(variables.i) + (variables.z - variables.w)**(-1) * diff(1, variables.w)), variables.zl) * variables.OObf(variables.i, (variables.w, variables.wl)) == -2 * math.pi * variables.hb(variables.i) * diff(variables.deltaa(2), variables.z) * (variables.z - variables.w) * variables.OObf(variables.i, (variables.w, variables.wl)) + 2 * math.pi * variables.deltaa(2) * (variables.z - variables.w) * diff(variables.OObf(variables.i, (variables.w, variables.wl)), variables.w)
diff(variables.Tbf(variables.z * variables.z, (variables.z, variables.zl), variables.z)) == 2 * math.pi * variables.lamdab(variables.i) * (variables.h - 1) * diff(variables.OObf(variables.i, (variables.z, variables.zl)), variables.z) 
diff(variables.Tb(variables.z * variables.z), variables.zl) + diff(variables.Tb(variables.zl * variables.z), variables.z) == 0
variables.Tb(variables.zl, variables.z) == 2 * math.pi * variables.lamdaa(variables.i) * (1 - variables.hb(variables.i)) * variables.OObf(variables.i, (variables.z, variables.zl)) 
variables.Betaa(variables.i) == 2 * (variables.hb(variables.i) - 1) * variables.lamdaa(variables.i)
1/2 * integrate(variables.d2z * variables.OObf(variables.i, (variables.z, variables.zl))) * integrate(variables.d2(variables.w) * variables.OObf(variables.j, (variables.w, variables.wl)))
2 * math.pi * integrate((diff(variables.r))/variables.r * variables.cab(variables.k, variables.i * variables.j)) * integrate(variables.d2(variables.w) * variables.OObf(variables.k, (variables.w, variables.wl))) 
variables.delta * variables.lamdaa(variables.k) == - 2 * math.pi * variables.epsilon * variables.cab(variables.k, variables.i * variables.j) * variables.lamdaa(variables.i) * variables.lamdaa(variables.j)
variables.Betaa(variables.k) == 2 * math.pi * variables.cab(variables.k, variables.i * variables.j) * variables.lamdaa(variables.i) * variables.lamdaa(variables.j)
variables.cab(variables.k, variables.i * variables.j) * variables.lamdaa(variables.i) * variables.lamdaa(variables.j) == 0
variables.Betaa(variables.i) == 2 * (variables.hb(variables.i) - 1) * variables.lamdaa(variables.i) + 2 * math.pi * variables.cab(variables.k, variables.i * variables.j) * variables.lamdaa(variables.i) * variables.lamdaa(variables.j)
diff(variables.C, variables.tau) == -12 * math.pi**2 * variables.Betaa(variables.i) * variables.Betaa(variables.j) * variables.Gb(variables.i * variables.j)
variables.Gb(variables.i * variables.j) == variables.za(2) * variables.zla(2) * variables.bracket(variables.OObf(variables.i, (variables.z, variables.zl)) * variables.OObf(variables.j, (0, 0)))
variables.Betaa(variables.i) == (diff(1, variables.z))/(diff(variables.lamdab(variables.i), variables.z)) * variables.Uf(variables.lamda)
variables.Uf(variables.lamda) == (variables.hb(variables.i) - 1) * variables.lamdaa(variables.i) * variables.lamdab(variables.i) + (2 * math.pi)/3 * variables.cb(variables.i, variables.j * variables.k) * variables.lamdaa(variables.i) * variables.lamdaa(variables.j) * variables.lamdaa(variables.k)
diff(variables.C, variables.tau) == 24 * math.pi**2 * variables.Betab(variables.j) * diff(variables.lamdaa(variables.j), variables.tau) == 24 * math.pi **2 * diff(variables.U, variables.tau)
variables.C == variables.c + 24 * math.pi**2 * variables.U 
diff(variables.lamda, variables.tau) == (1 - variables.h)*variables.lamda - math.pi * variables.cb(1, 1, 1) * variables.lamdaa(2)
variables.lamdap == (1 - variables.h)/(math.pi * variables.cb(1, 1, 1))
variables.cp == variables.c - 8 * ((1 - variables.h)**(3))/(variables.cab(2, (1, 1, 1)))
# statistical mechanics 2.266.1
variables.Z == integrate(variables.brackets(diff(variables.q)) * math.exp(-variables.Beta * variables.H))
variables.Z == integrate(variables.brackets(diff(variables.phi) * math.exp(-variables.S/variables.hl)))
variables.H == -summation(variables.sigmab(variables.i) * variables.sigmab(variables.ip), variables.links)
# check defined states 2.267.7
variables.h == 0 
variables.h == 1/16 
variables.h == 1/2 
# -
variables.sigmab(variables.i) == variables.sigmaf(variables.z, variables.zl) == variables.OObf((1, 2), variables.z) * variables.OObftilde((1, 2), variables.zl)
variables.bracket(variables.sigmaf(variables.z, variables.zl) * variables.sigmaf(0, 0)) == (variables.z * variables.zl)**(-2 * variables.h) == (variables.z * variables.zl)**(-1/8)
variables.OObf((1, 2), variables.z) * variables.OObftilde((1, 3), variables.zl) 
variables.brackets(variables.OOb(1, 1) * variables.OObtilde(1, 1)) + variables.brackets(variables.OOb(1,2) * variables.OObtilde(1, 2)) + variables.brackets(variables.OOb(1, 3) * variables.OObtilde(1, 3))
variables.H == - summation(variables.Ref(variables.sigmab(variables.i) * variables.star(variables.sigmab(variables.ip))), variables.links)
variables.LO == diff(variables.phi, variables.z) * diff(variables.phi, variables.zl) + variables.lamdab(1) * variables.phia(2) + variables.lamdab(2) * variables.phia(4)
variables.OOb(2, 2) * variables.OOb(2, 2) == variables.brackets(variables.OOb(1, 1)) + variables.brackets(variables.OOb(3, 1)) + variables.brackets(variables.OOb(3, 3)) + variables.brackets(variables.OOb(1, 3))
variables.phia(variables.n) == variables.OOb(variables.n + 1, variables.n + 1)
0 <= variables.n <= variables.m - 2 
variables.phia(variables.m - 1 + variables.n) == variables.OOb(variables.n + 1, variables.n + 2) 
0 <= variables.n <= variables.m - 3
variables.m * variables.lamdab(variables.m) * variables.phia(2 * variables.m - 3) == diff(diff(variables.phi, variables.zl), variables.z) == np.dot(variables.Lb(-1) * variables.Lbtilde(-1), variables.phi)
variables.OOb(variables.m -1, variables.m - 2) == variables.OOb(1, 3) 
variables.h == 1 - 2/(variables.m + 1)
variables.OOb(1, 3) * variables.OOb(1, 3) == variables.brackets(variables.OOb(1, 1)) + variables.brackets(variables.OOb(1, 3)) + variables.brackets(variables.OOb(1, 5))
variables.cp == variables.c - 12/(variables.ma(3))
