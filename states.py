import math
from scipy.integrate import quad
import numpy as np 
import sympy as sp
from sympy import symbols, diff,integrate, Function, Derivative, summation, oo, Abs, ln
from sympy.diffgeom import WedgeProduct
import variables
import constants

n = symbols('n') #placeholder for steps
def String():
    if variables.point == variables.dirac(0, variables.k): # check point variable
        variables.m**2 == (2-variables.D)/(24*variables.SlopeRegge)
    elif variables.point == variables.TensionStringab(variables.i, -1)*variables.dirac(1, variables.k):
        variables.m**2 == (26-variables.D)/(24*variables.SlopeRegge)
    elif variables.point == variables.dirac(variables.N, variables.Ntilde, variables.k):
        variables.point == np.prod([np.prod([(variables.TensionStringab(variables.i, -n)**variables.Nb(variables.i, n)*variables.TensionStringabtilde(variables.i, -n)**(variables.Nbtilde(variables.i, n)))/((n**(variables.Nb(variables.i, n))*math.factorial(variables.Nb(variables.i, n))*n**(variables.Nbtilde(variables.i, n))*math.factorial(variables.Nbtilde(variables.i, n)))**(1/2)) for n in np.linspace(1, math.inf)]) for variables.i in np.linspace(2, variables.D-1)]) * variables.dirac(0, 0, variables.k)
        variables.m == 2*variables.pa(variables.pos)*variables.H - variables.pa(variables.i)* variables.pa(variables.i) == 2/variables.SlopeRegge * summation(variables.TensionStringab(variables.i, -n)*variables.TensionStringab(variables.i, n) + variables.TensionStringabtilde(variables.i, -n)*variables.TensionStringabtilde(variables.i, n) + variables.A + variables.Atilde, (n, 1, math.inf)) == 2/variables.SlopeRegge(variables.N + variables.Ntilde + variables.A + variables.Atilde)
        variables.A == variables.Atilde == (2-variables.D)/24
        variables.P == -quad(np.diff(variables.sigma)*np.prod(diff(variables.xa[variables.i]), variables.sigma), 0, variables.l) == -2*math.pi/variables.l(summation(variables.TensionStringab(variables.i, -n)*variables.TensionStringab(variables.i, n) - variables.TensionStringabtilde(variables.i, -n)*variables.TensionStringabtilde(variables.i, n)) +variables.A - variables.Atilde, n, 1, math.inf) == -2*math.pi/variables.l(variables.N - variables.Ntilde)
        variables.N == variables.Ntilde
    elif variables.point == variables.dirac(0, 0, variables.k):
        variables.m**2 == (26-variables.D)/(6*variables.SlopeRegge) #check excited states pg 27
        variables.A =variables.Atilde == -1
        variables.D == 26
        
#oscilator ground state valueation
variables.bb0(variables.dirac(variables.down)) == 0
variables.bb0(variables.dirac(variables.up)) == variables.dirac(variables.down)
variables.cb0(variables.dirac(variables.down)) == variables.dirac(variables.up)
variables.cb0(variables.dirac(variables.up)) == 0
variables.bbn(variables.dirac(variables.down)) == variables.bbn(variables.dirac(variables.up)) == variables.cbn(variables.dirac(variables.down)) == variables.cbn(variables.dirac(variables.up)) == 0 # check bounds on 'n' 61.18

# groundstate moving forawrd
variables.dirac(1) == variables.dirac(0, 0)
variables.m >= 1
variables.alphaab(variables.mu, -variables.m) == (2/(variables.SlopeRegge))**(1/2) * integrate(diff(variables.z)/(2*math.pi) * variables.za(-variables.m) * diff(variables.xa(variables.mu, variables.z), variables.m, variables.z)) == (2/variables.SlopeRegge)**(1/2) * variables.i/(math.factorial(variables.m-1)) * diff(variables.xa(variables.mu, variables.D0), variables.m) =variables.alphaab(variables.mu, -variables.m) * variables.dirac(variables.D1) #conic integration
variables.alphaabtilde(variables.mu, -variables.m) * variables.dirac(1) == (2/(variables.SlopeRegge))**(1/2) * variables.i/(math.factorial(variables.m-1)) * diff(variables.xa(variables.mu, variables.D0), variables.m, variables.zl) # double check partial diff

def nextLevelStates():
    # 123.12
    variables.dirac(variables.e, variables.k) == variables.e * variables.alphab(-variables.D1) * variables.dirac(0, variables.k)
    variables.bracket(variables.dirac(variables.e, variables.k), variables.dirac(variables.e, variables.kp)) == variables.bracket(variables.dirac(0, variables.k), (variables.star(variables.e)  * variables.alphab(1) * variables.e * variables.alphab(-1)), variables.dirac(0, variables.kp)) 
    # D states
    variables.dirac(variables.e, variables.k) == variables.e * variables.alpha(-1) * variables.dirac(0, variables.k)
    variables.bracket(variables.dirac(variables.e, variables.k), variables.dirac(variables.e, variables.kp)) == variables.brackets(variables.dirac(0, variables.k), variables.e*variables.alphab(1) * variables.e * variables.alphab(-1), variables.dirac(0, variables.kp)) == variables.brackets(variables.dirac(0, variables.k), (variables.star(variables.e) * variables.e + variables.star(variables.e) * variables.alphab(-1) * variables.e * variables.alphab(1)), variables.dirac(0, variables.kp)) == variables.star(variables.ea(variables.mu)) * variables.eb(variables.mu) * (2*math.pi)**(variables.D) * variables.deltaa(variables.D) * (variables.k-variables.kp)
    variables.dagger(variables.alphaab(variables.mu, variables.n)) == variables.alphaab(variables.mu, -variables.n)
    variables.bracket(variables.dirac(0, variables.k), variables.dirac(0, variables.kp)) == (2*math.pi)**(variables.D) * variables.deltaa(variables.D) * (variables.k - variables.kp)
    variables.m**2 == (1+variables.A)/variables.SlopeRegge
    variables.e * variables.k * variables.dirac(0, variables.k) == 0
    # check case scenarios for following condition check 124.18
    variables.Lab(variables.m, -1) * variables.dirac(0, variables.k) == (2*variables.SlopeRegge)**(1/2) * variables.k * variables.alphab(-1) * variables.dirac(0, variables.k)
    #at all mass levels
    variables.hilbertb(variables.O*variables.C*variables.Q) == variables.hilbertb(variables.lightcone)
    # check first two levels 125.21

#241.1, mass shell conditions
variables.ma(2) == -variables.ka(variables.mu) * variables.kb(variables.mu) == (variables.kab(25, variables.L))**2 + 4/variables.SlopeRegge * (variables.N - 1) == (variables.kab(25, variables.R))**2 +4/variables.SlopeRegge * (variables.Ntilde - 1)
variables.ma(2) == (variables.na(2))/(variables.Ra(2)) + (variables.wa(2) * variables.Ra(2))/(variables.SlopeRegge**2) +2/variables.SlopeRegge * (variables.N + variables.Ntilde -2) 
0= variables.n * variables.w + variables.N - variables.Ntilde 
diff(variables.xa(variables.mu), variables.z) * diff(variables.xa(25), variables.zl) - diff(variables.xa(25), variables.z) * diff(variables.xa(variables.mu), variables.zl) == diff((variables.xa(25) * diff(variables.xa(variables.mu), variables.z)), variables.zl) - diff(variables.xa(25) * diff(variables.xa(variables.mu), variables.zl), variables.z)
# check all freefloating expressions 241-242 
(variables.n + variables.w)**2 + 4 * variables.nN == (variables.n - variables.w)**2 + 4 * variables.Ntilde == 4 
variables.n == variables.w == variables.posneg * 1
variables.N == 0
variables.Ntilde == 1
variables.n == -variables.w == variables.posneg * 1
variables.N == 1
variables.Ntilde == 0
variables.n == variables.posneg * 2
variables.w == variables.N == variables.Ntilde == 0
variables.w == variables.posneg * 2
variables.n == variables.N == variables.Ntilde == 0
variables.jaf(1, variables.z) == variables.point(math.cos(variables.brackets(2 * variables.SlopeRegge**(-1/2) * variables.Xabf(25, variables.L, variables.z))))
variables.jaf(2, variables.z) == variables.point(math.sin(variables.brackets(2 * variables.SlopeRegge**(-1/2) * variables.Xabf(25, variables.L, variables.z))))
variables.jaf(3, variables.z) =variables.i * diff(variables.Xabf(25, variables.L, variables.z), variables.z)/(variables.SlopeRegge**(1/2))
variables.jaf(variables.i, variables.z) == summation((variables.jab(variables.i, variables.m))/(variables.za(variables.m + 1)), (variables.m, -math.inf, math.inf))
variables.brackets(variables.jab(variables.i, variables.m), variables.jab(variables.j, variables.n)) == (variables.m)/2 * variables.deltab(variables.m, -variables.n) * variables.deltaa(variables.i * variables.j) + variables.i * variables.epsilona(variables.i * variables.j * variables.k) * variables.jab(variables.k, variables.m + variables.n)
variables.gab(2, 25) == (2 * variables.ketaab(2, 25))/variables.SlopeRegge
variables.gab(2, 4) == (2 * variables.ketaab(2, 4))/variables.SlopeRegge
variables.gabf(2, (variables.G, 4), variables.E) == variables.ketaab(2, 4) * variables.Ea(2)
variables.gab(2, 5) == 2 * math.pi * variables.rhob(5) * variables.gab(2, 4)
variables.hat(variables.gab(2, 5)) == variables.gab(2, 5) * variables.E == 2 * math.pi * variables.rhob(5) * variables.E * variables.gab(2, 4)
variables.m == (Abs(variables.Ra(2) - variables.SlopeRegge))/(variables.R * variables.SlopeRegge) == 2/variables.SlopeRegge * Abs(variables.R - variables.SlopeRegge**(1/2))
# check 246.22
variables.Uf(variables.m) == (diff(variables.Uf(variables.M), variables.z))/(diff(variables.Mb(variables.i * variables.j), variables.z)) == 0
variables.Mb(11) * variables.Mb(22) * variables.Mb(33) == variables.Mb(11) * variables.Mb(22) == variables.Mb(11) * variables.Mb(33) == variables.Mb(22) * variables.Mb(33) =0
variables.ma(2) == (variables.na(2))/(variables.Ra(2)) + (variables.wa(2) * variables.Ra(2))/(variables.SlopeRegge**2) + 2/variables.SlopeRegge * (variables.N + variables.Ntilde + 2)
variables.R == variables.Rp == variables.SlopeRegge/variables.R 
variables.n == variables.w 
variables.pab(25, variables.R) == - variables.pab(25, variables.R)
variables.Xap(25, (variables.z, variables.zl)) == variables.Xabf(25, variables.L, (variables.z)) - variables.Xabf(25, variables.R, variables.zl)
variables.Rb(variables.selfdual) == variables.Rb(np.cross(variables.S * variables.Uf(2), variables.S * variables.Uf(2))) == variables.SlopeRegge**(1/2)
variables.rhop == variables.SlopeRegge/variables.rho 
variables.kappap == (variables.SlopeRegge**(1/2))/variables.rho * variables.kappa 
variables.ea(variables.Phip) == (variables.SlopeRegge**(1/2))/variables.rho * variables.ea(variables.Phi)

# 249.1, compactification of several dimensions, periodic conditions:
variables.xa(variables.m) == variables.xa(variables.m) + 2 * math.pi * variables.R 
26 - variables.k <= variables.m <= 25
variables.SB == ((2 * math.pi * variables.R)**variables.k)/(2 * variables.kappaab(2, 0)) * integrate(variables.da(variables.d) * variables.x * (-variables.Gb(variables.d))**(1/2) * variables.ea(-2 * variables.Phib(variables.d)) * variables.brackets(variables.RBb(variables.d) + 4 * diff(variables.Phib(variables.d), variables.mu) * diff(variables.Phib(variables.d), 1, variables.mu) - 1/4 * variables.Ga(variables.m * variables.n) * variables.Ga(variables.p * variables.q) * (diff(variables.Gb(variables.m * variables.p), variables.mu) * diff(variables.Gb(variables.n * variables.q), 1, variables.mu) + diff(variables.Bb(variables.m * variables.p), variables.mu) * diff(variables.Bb(variables.n * variables.q), 1, variables.mu)) - 1/4 * variables.Gb(variables.m * variables.n) * variables.Fab(variables.m, variables.mu * variables.v) * variables.Fa(variables.n * variables.mu * variables.v) - 1/4 * variables.Ga(variables.m * variables.n) * variables.Hb(variables.m * variables.mu * variables.v) * variables.Hab(variables.mu * variables.v, variables.n) - 1/12 * variables.Hb(variables.mu * variables.v * variables.lamda) * variables.Ha(variables.mu * variables.v * variables.lamda)))
variables.Phib(variables.d) == variables.Phi - 1/4 * ln(sp.Determinant(variables.Gb(variables.m * variables.n)))
# string spectrum
variables.Bb(variables.m * variables.n) * diff(variables.ga(1/2) * variables.epsilona(variables.a * variables.b) * variables.xa(variables.m) * diff(variables.xa(variables.n), variables.b), variables.a)
variables.Xaf(variables.m, (variables.sigmaa(1), variables.sigmaa(2))) == variables.xafs(variables.m, (variables.sigmaa(2))) + variables.wa(variables.m) * variables.R * variables.sigmaa(1)
variables.L == 1/(2 * variables.SlopeRegge) * variables.Gb(variables.m * variables.n) * (diff(variables.xas(variables.m), variables.tau) * diff(variables.xas(variables.n), variables.tau) + variables.wa(variables.m) * variables.wa(variables.n) * variables.Ra(2)) - variables.i/variables.SlopeRegge * variables.Bb(variables.m * variables.n) * diff(variables.xas(variables.m), variables.tau) * variables.wa(variables.n) * variables.R
variables.pb(variables.m) == -(diff(variables.L, variables.z))/(diff(variables.va(variables.m), variables.z)) == 1/variables.SlopeRegge * (variables.Gb(variables.m * variables.n) * variables.va(variables.n) + variables.Bb(variables.m * variables.n) * variables.wa(variables.n) * variables.R)
variables.vb(variables.m) == variables.SlopeRegge * (variables.nb(variables.m))/(variables.R) - variables.Bb(variables.m * variables.n) * variables.wa(variables.n) * variables.R 
1/(2 * variables.SlopeRegge) * variables.Gb(variables.m * variables.n) * (variables.va(variables.m) * variables.va(variables.n) + variables.wa(variables.m) * variables.wa(variables.n) * variables.Ra(2))
variables.ma(2) == 1/(2 * variables.SlopeRegge**2) * variables.Gb(variables.m * variables.n) * (variables.vab(variables.m, variables.L) * variables.vab(variables.n, variables.L) + variables.vab(variables.m, variables.R) * variables.vab(variables.n, variables.R)) + 2/variables.SlopeRegge * (variables.N + variables.Ntilde -2)
variables.vab(variables.m, (variables.L, variables.R)) == variables.va(variables.m) + variables.posneg * variables.wa(variables.m) * variables.R 
0 == variables.Gb(variables.m * variables.n) * (variables.vab(variables.m, variables.L) * variables.vab(variables.n, variables.L) - variables.vab(variables.m, variables.R) * variables.vab(variables.n, variables.R)) + 4 * variables.SlopeRegge * (variables.N - variables.Ntilde) == 4 * variables.SlopeRegge * (variables.nb(variables.m) * variables.wa(variables.m) + variables.N - variables.Ntilde)
variables.xa(variables.m) == (variables.wab(variables.m, 1) * variables.sigmaa(1) + variables.wab(variables.m, 2) * variables.sigmaa(2)) * variables.R
2 * math.pi * variables.i * variables.bb(variables.m * variables.n) * variables.wab(variables.m, 1) * variables.wab(variables.n, 2) 
variables.Gb(variables.m * variables.n) == variables.eba(variables.m, variables.r) * variables.eba(variables.n, variables.r) 

# 268.1, D-banes
variables.alphaab(variables.mu, -1) * variables.dirac(variables.k, variables.i * variables.i)
variables.VO == variables.i * diff(variables.xa(variables.mu), variables.t)
variables.alphaab(25, -1) * variables.dirac(variables.k, variables.i * variables.i)
variables.VO == variables.i * diff(variables.xa(25), variables.t) == diff(variables.xap(25), variables.n)
#action, 270.2
variables.SBb(variables.p) == - variables.Tb(variables.p) * integrate(variables.da(variables.p + 1) * variables.zeta * variables.ea(-variables.Phi) * variables.brackets(-sp.Determinant(variables.Gb(variables.a * variables.b) + variables.Bb(variables.a * variables.b) + 2 * math.pi * variables.SlopeRegge * variables.Fb(variables.a * variables.b)))**(1/2))
variables.Gbf(variables.a * variables.b, variables.zeta) == (diff(variables.xa(variables.mu), variables.z))/(diff(variables.zetaa(variables.a), variables.z)) * (diff(variables.xa(variables.v), variables.z))/(diff(variables.zetaa(variables.b), variables.z)) * variables.Gbf(variables.mu * variables.v, variables.Xf(variables.zeta))
variables.Bbf(variables.a * variables.b, variables.zeta) == (diff(variables.xa(variables.mu), variables.z))/(diff(variables.zetaa(variables.a), variables.z)) * (diff(variables.xa(variables.v), variables.z))/(diff(variables.zetaa(variables.b), variables.z)) * variables.Bbf(variables.mu * variables.v, variables.Xf(variables.zeta))
variables.xap(2) == -2 * math.pi * variables.SlopeRegge * variables.xa(1) * variables.Fb(12) 
integrate(diff(variables.xa(1)) * variables.brackets(1 + (diff(variables.xap(2), 1))**2)**(1/2)) == integrate(diff(variables.xa(1)) * variables.brackets(1 + (2 * math.pi * variables.SlopeRegge * variables.Fb(12))**2)**(1/2))
#check world-sheet action definition, 271.6
variables.delta * variables.Ab(variables.mu) == diff(variables.lamda, variables.mu)
variables.delta * variables.Bb(variables.mu * variables.v) == diff(variables.zetab(variables.v), variables.mu) - diff(variables.zetab(variables.mu), variables.v) 
variables.delta * variables.Ab(variables.mu) == (-variables.zetab(variables.mu))/(2 * math.pi * variables.SlopeRegge)
variables.Bb(variables.mu * variables.v) + 2 * math.pi * variables.Fb(variables.mu * variables.v) == 2 * math.pi * variables.SlopeRegge * variables.FOb(variables.mu * variables.v)
# check field strength leading term approximation 272.11
variables.SBb(variables.p) == - variables.Tb(variables.p) * integrate(variables.da(variables.p + 1) * variables.zeta * variables.Trf(variables.bracketc(variables.ea(-variables.Phi) * variables.brackets(-sp.Determinant(variables.Gb(variables.a * variables.b) + variables.Bb(variables.a * variables.b) + 2 * math.pi * variables.SlopeRegge * variables.Fb(variables.a * variables.b)))**(1/2) + variables.Of(variables.brackets(variables.xa(variables.m), variables.xa(variables.n))**2))))
variables.Tb(variables.p) * variables.ea(-variables.Phi) * np.prod(2 * math.pi * variables.Rb(variables.i), (variables.i, 1, variables.p))
# check like terms, 273.13-15
variables.Tb(variables.p) == (variables.Tb(variables.p -1))/(2 * math.pi * variables.SlopeRegge)**(1/2)
variables.AO == variables.i * variables.Vb(variables.p + 1) * quad((diff(variables.t))/variables.t * (8 * math.pi **2 * variables.SlopeRegge * variables.t)**(-(variables.p + 1)/2) * math.exp(-variables.t * variables.y**2/(2 * math.pi * variables.SlopeRegge)) * variables.etaf(variables.i * variables.t)**(-24), 0, math.inf) 
(variables.i * variables.Vb(variables.p +1))/((8 * math.pi **2 * variables.SlopeRegge)**((variables.p + 1)/2)) * quad(diff(variables.t) * variables.t**((21-variables.p)/2) * math.exp(-variables.t * variables.ya(2)/(2 * math.pi * variables.SlopeRegge))) # second part of 275.17, unifished expasnion of cross product
variables.AO == variables.i * variables.Vb(variables.p + 1) * (24/(2**12)) * (4 * math.pi**2 * variables.SlopeRegge)**(11-variables.p) * math.pi**((variables.p-23)/2) * variables.Rho * ((23 - variables.p)/2) * Abs(variables.y)**(variables.p - 23) == variables.i * variables.Vb(variables.p + 1) * (24 * math.pi)/(2**10) * (4 * math.pi**2 * variables.SlopeRegge)**(11-variables.p) * variables.Gbf(25 - variables.p, variables.y)
variables.SB == 1/(2 *variables.kappa**2) * integrate(variables.da(26) * variables.X * (-variables.Gtilde)**(1/2) * (variables.RBtilde - 1/6 * variables.Deltab(variables.mu) * variables.Phitilde * variables.Deltaatilde(variables.mu) * variables.Phitilde)) 
variables.SBb(variables.p) == -variables.taub(variables.p) * integrate(variables.da(variables.p + 1) * variables.zeta * math.exp((variables.p - 11)/12 * variables.Phitilde) * (-sp.Determinant(variables.Gbtilde(variables.a * variables.b)))**(1/2))
variables.Fb(variables.v) == diff(variables.hb(variables.mu * variables.v), 1, variables.hat(variables.mu)) - 1/2 * diff(variables.hab(variables.hat(variables.mu), variables.mu), variables.v) == 0
variables.SB == -1/(8 * variables.kappaa(2)) * integrate(variables.da(26) * variables.X * (diff(variables.hb(variables.v * variables.lamda), variables.mu) * diff(variables.ha(variables.hat(variables.v) * variables.hat(variables.lamda)), 1, variables.hat(variables.mu)) - 1/2 * (diff(variables.hab(variables.hat(variables.v), variables.v), variables.mu) * diff(variables.hab(variables.hat(variables.lamda), variables.lamda), 1, variables.hat(variables.mu))) + 2/3 * diff(variables.Phitilde, variables.mu) * diff(variables.Phitilde, 1, variables.hat(variables.mu))))
variables.bracket(variables.Phitilde * variables.Phitilde) == - ((variables.D - 2) * variables.i * variables.kappaa(2))/(4 * variables.ka(2))
variables.bracket(variables.hb(variables.mu * variables.v) * variables.hb(variables.sigma * variables.rho)) == -((2 * variables.i * variables.kappaa(2)))/(variables.ka(2)) * (variables.etab(variables.mu * variables.sigma) * variables.etab(variables.v * variables.rho) + variables.etab(variables.mu * variables.rho) * variables.etab(variables.v * variables.sigma) - 2/(variables.D -2) * variables.etab(variables.mu * variables.v) * variables.etab(variables.sigma * variables.rho))
variables.SBb(variables.p) == -variables.taub(variables.p) * integrate(variables.da(variables.p + 1) * variables.zeta * ((variables.p -11)/12 * variables.Phitilde - 1/2 * variables.hb(variables.a * variables.a)))
variables.AO == (variables.i * variables.kappaa(2) * variables.tauab(2, variables.p))/(variables.kab(2, variables.perp)) * variables.Vb(variables.p + 1) * variables.bracketc(6 * variables.brackets((variables.p-11)/12)**2 + 1/2 * variables.brackets(2 * (variables.p + 1) - 1/12 * (variables.p + 1)**2)) == (6 * variables.i * variables.kappaa(2) * variables.tauab(2, variables.p))/(variables.kab(2, variables.perp)) * variables.Vb(variables.p +1)
variables.tauab(2, variables.p) == math.pi/(256 * variables.kappaa(2)) * (4 * math.pi**2 * variables.SlopeRegge)**(11-variables.p)
# check expanded action 276.27
(variables.gab(2, variables.o))/(variables.gb(variables.c)) == (4 * math.pi * variables.SlopeRegge * variables.gabp(2, variables.o))/variables.kappa == 2 **18 * math.pi**25/2 * variables.SlopeRegge**6

#T- duality of unoriented theories 277.1
variables.omega 
variables.xab(variables.M, variables.L, variables.z) == variables.xab(variables.M, variables.R, variables.z)
variables.omega
variables.xap(variables.m, (variables.z, variables.zl)) == - variables.xap(variables.m, (variables.zl, variables.z))
variables.xa(variables.mu, (variables.z, variables.zl)) == variables.xa(variables.mu, (variables.zl, variables.z))
#tensor conditions 277.3
variables.Gbf(variables.mu * variables.v, variables.xp) == variables.Gbf(variables.mu * variables.v, variables.x)
variables.Bbf(variables.mu * variables.v, variables.xp) == -variables.Bbf(variables.mu * variables.v, variables.x)
variables.Gbf(variables.mu * variables.n, variables.xp) == - variables.Gbf(variables.mu * variables.n, variables.x)
variables.Bbf(variables.mu * variables.n, variables.xp) == variables.Bbf(variables.mu * variables.n, variables.x)
variables.Gbf(variables.m * variables.n, variables.xp) == variables.Gbf(variables.m * variables.n, variables.x)
variables.Bbf(variables.m * variables.n, variables.xp) == -variables.Bbf(variables.m * variables.n, variables.x)
# check open string eigenvalue pairing 279.4
2**(12 - variables.k) * variables.Tb(variables.p) * integrate(variables.da(variables.p + 1) * variables.zeta * variables.ea(-variables.Phi) * (-sp.Determinant(variables.Gb(variables.a * variables.b)))**(1/2))

# OCQ,2.20.1
variables.Lab(variables.m, variables.n) * variables.dirac(variables.psi) == 0
variables.n > 0
variables.Gab(variables.m, variables.r) * variables.dirac(variables.psi) == 0
variables.r >= 0
variables.Lab(variables.m, variables.n) * variables.dirac(variables.chi) == 0
variables.n < 0
variables.Gab(variables.m, variables.r) * variables.dirac(variables.chi) == 0
variables.r < 0 
variables.Lb(0) * variables.dirac(variables.psi) == variables.Hf(variables.dirac(variables.psi)) == 0
variables.H == variables.SlopeRegge * variables.pa(2) + variables.N - 1/2 # NS
variables.H == variables.SlopeRegge * variables.pa(2) + variables.N # R
8 * (-1/24 - 1/48) == -1/2 # transverse mode NS
8 * (-1/24 + 1/24) == 0 # transverse R
# check term definition 2.21.6
variables.ma(2) == - variables.ka(2) == -1/(2 * variables.SlopeRegge) 
variables.diracb((variables.e, variables.k), variables.NS) == np.dot(variables.e, (variables.psib(-1/2) * variables.diracb((0, variables.k), variables.NS))) 
0 == variables.Lb(0) * variables.diracb((variables.e, variables.k), variables.NS) == variables.SlopeRegge * variables.ka(2) * variables.diracb((variables.e, variables.k), variables.NS) 
0 == variables.Gab(variables.m, 1/2) * variables.diracb((variables.e, variables.k), variables.NS) == np.dot((2 * variables.SlopeRegge)**(1/2) * variables.k, (variables.e * variables.diracb((0, variables.k), variables.NS)))
variables.Gab(variables.m, -1/2) * variables.diracb((0, variables.k), variables.NS) == np.dot((2 * variables.SlopeRegge)**(1/2) * variables.k, variables.psib(-1/2) * variables.diracb((0, variables.k), variables.NS))
variables.ka(2) == 0 
np.dot(variables.e, variables.k) == 0
variables.ea(variables.mu) == variables.ea(variables.mu) + variables.lamda * variables.ka(variables.mu) 
variables.diracb((variables.mu, variables.k), variables.R) == variables.diracb((variables.sB, variables.k), variables.R) * variables.ub(variables.sB)
0 == variables.Lb(0) * variables.diracb((variables.u, variables.k), variables.R) == variables.SlopeRegge * variables.ka(2) * variables.diracb((variables.u, variables.k), variables.R) 
variables.Gab(variables.m, 0) * variables.diracb((variables.u, variables.k), variables.R)  == np.dot(variables.SlopeRegge**(1/2) * variables.diracb((variables.sBp, variables.k), variables.R) * variables.k, variables.Rhob(variables.sBp * variables.sB) * variables.ub(variables.sB))
np.dot(variables.k, variables.Rhob(variables.sBp * variables.sB) * variables.ub(variables.sB)) == 0
variables.kb(0) * variables.Rhoa(0) + variables.kb(1) * variables.Rhoa(1) == - variables.kb(1) * variables.Rhoa(0) * (variables.Rhoa(0) * variables.Rhoa(1) - 1) == -2 * variables.kb(1) * variables.Rhoa(0) * (variables.Sb(0) - 1/2)
(variables.Sb(0) - 1/2) * variables.diracb(((variables.sB, 0), variables.k), variables.R) * variables.ub(variables.sB) == 0
variables.B16 == (+1/2, variables.B8) + (-1/2, variables.B8p)
variables.B16p == (+1/2, variables.B8p) + (-1/2, variables.B8)
#reference table on 2.23
variables.SlopeRegge/4 * variables.ma(2) == variables.N - variables.v == variables.Ntilde - variables.vtilde 
variables.dirac((variables.i, variables.sB)) * variables.Rhoab(variables.i, variables.sB * variables.sBp)
variables.Qb(variables.B) == 1/(2 * math.pi * variables.i) * integrate((diff(variables.z) * variables.jb(variables.B) - diff(variables.zl) * variables.jbtilde(variables.B)), 2 * math.pi) # contour integral
variables.jb(variables.B) == variables.c * variables.Tab(variables.m, variables.B) + variables.gamma * variables.Tab(variables.m, variables.F) + 1/2 * (variables.c * variables.Tab(variables.g, variables.B) + variables.gamma * variables.Tab(variables.g, variables.F)) == variables.c * variables.Tab(variables.m, variables.B) + variables.gamma * variables.Tab(variables.m, variables.F) + variables.b * variables.c * diff(variables.c, variables.z) * variables.Beta * variables.gamma + 1/4 * variables.c * (diff(variables.Beta, variables.z)) * variables.gamma - 3/4 * variables.c * variables.Beta * diff(variables.gamma, variables.z) - variables.b * variables.gammaa(2) 
# check 2.24.22 
variables.bracketc(variables.Qb(variables.B), variables.bb(variables.n)) == variables.Lb(variables.n) 
variables.brackets(variables.Qb(variables.B), variables.Betab(variables.r)) == variables.Gb(variables.r) 
variables.Qb(variables.B) == summation(variables.cb(-variables.m) * variables.Lab(variables.m, variables.m) + summation(variables.gammab(-variables.r) * variables.Gab(variables.m, variables.r) - summation(1/2 * (variables.n - variables.m) * variables.point(variables.bb(-variables.m -variables.n) * variables.cb(variables.m) * variables.cb(variables.n)), ((variables.m, variables.n))), (variables.r)), variables.m) + summation(variables.brackets(1/2 * (2 * variables.r - variables.m) * variables.point(variables.Betab(-variables.m-variables.r) * variables.cb(variables.m) * variables.gammab(variables.r)) - variables.point(variables.bb(-variables.m) * variables.gammab(variables.m-variables.r) * variables.gammab(variables.r))), ((variables.m, variables.r))) + variables.aa(variables.g) * variables.cb(0)
variables.bb(0) * variables.dirac(variables.psi) == variables.Lb(0) * variables.dirac(variables.psi) == 0
variables.Betab(0) * variables.dirac(variables.psi) == variables.Gb(0) * variables.dirac(variables.psi) == 0

#2.25.2
variables.alpha == 1- 2 * variables.v 
(variables.alpha, variables.F, variables.alphatilde, variables.Ftilde)
math.exp(math.pi * variables.i * (variables.Fb(1) * variables.alphab(2) - variables.Fb(2)* variables.alphab(1) - variables.Fbtilde(1) * variables.alphabtilde(2) + variables.Fbtilde(2) * variables.alphabtilde(1)))
(variables.Fb(1) * variables.alphab(2) - variables.Fb(2)* variables.alphab(1) - variables.Fbtilde(1) * variables.alphabtilde(2) + variables.Fbtilde(2) * variables.alphabtilde(1)) == 2 * variables.ZB 
(variables.alphab(1) + variables.alphab(2), variables.Fb(1) + variables.Fb(2), variables.alphabtilde(1) + variables.alphabtilde(2), variables.Fbtilde(1) + variables.Fbtilde(2))
#check spectra solving equations 2.26.table
variables.xa(2) == -variables.xa(2) 
variables.psia(2) == - variables.psia(2) 
variables.psiatilde(2) == - variables.psiatilde(2)
# check solutions 2.27.table
#IIA:
variables.brackets(0) + variables.brackets(1) + variables.brackets(2) + variables.brackets(3) + (2) + variables.B8 + variables.B8p + variables.B56 + variables.B56p
#IIB:
variables.brackets(0)**(2) + variables.brackets(2)**2 + variables.bracketsb(4, variables.pos) + (2) + variables.B8p**2 +variables.B56**2
math.exp(math.pi* variables.i * variables.F) == math.exp(math.pi * variables.i * variables.Ftilde) == variables.pos * 1
math.exp(math.pi* variables.i * variables.F) == +1
math.exp(math.pi* variables.i * variables.Ftilde) == (-1)**(variables.alphatilde)
variables.alpha == variables.alphatilde 
math.exp(math.pi * variables.i * variables.F) == math.exp(math.pi * variables.i * variables.Ftilde)
variables.psiab(variables.mu, (-1/2)) * variables.diracb((0, variables.sB, variables.k), variables.NS - variables.R) * variables.ub(variables.mu * variables.sB) 
variables.ub(variables.mu * variables.sB) == variables.ub(variables.mu * variables.sB) + variables.kb(variables.mu) * variables.zetab(variables.sB) 
variables.VOb(variables.sB) * variables.ea(-variables.phitilde) * variables.psiatilde(variables.mu) * variables.ea(np.dot(variables.i * variables.k, variables.X)) 
variables.ea(variables.phi) * variables.psia(variables.mu) * variables.VObtilde(variables.sB) * variables.ea(np.dot(variables.i * variables.k, variables.X))
variables.VOb(variables.sB) * variables.VObtilde(variables.sBp) 
variables.VObtilde(variables.sB) * variables.VOb(variables.sBp) == - variables.VOb(variables.sBp) * variables.VObtilde(variables.sB)
variables.brackets(0) + variables.brackets(2) + (2) + variables.B8p + variables.B56 == variables.B1 + variables.B28 + variables.B35 + variables.B8p + variables.B56
# check spectra, 2.30.table
variables.brackets(0) + variables.brackets(2) + (2) + variables.B8p + variables.B56 + variables.holderb(variables.B8b(variables.v) + variables.B8,variables.S * variables.Of(variables.n))
variables.brackets(0) + variables.brackets(2) + (2) + variables.B8p + variables.B56 + variables.holderb(variables.B8b(variables.v) + variables.B8,variables.S * variables.pf(variables.k))

(variables.Rhoa(0) * diff(0, 1) + variables.Rhoa(1) * diff(1)) * variables.u == 0
variables.B16 == (1, variables.B8) + (-1/2, variables.B8p)
(variables.taub(variables.F1))/(variables.taub(variables.D1)) == variables.g == variables.ea(variables.Phi)
variables.lb(0) == (4 * math.pi**3)**(-1/8) * variables.kappaa(1/4)
(variables.tauab(-1/2, variables.F1))/(variables.ga(-1/4)) == (variables.lb(0)) == (variables.tauab(-1/2, variables.D1))/(variables.ga(1/4))
variables.Phip == - variables.Phi 
variables.Gbp(variables.mu * variables.v) == variables.ea(-variables.Phi) * variables.Gb(variables.mu * variables.v)
variables.Bbp(2) == variables.Cb(2)
variables.Cbp(2) == -variables.Bb(2)
variables.Cbp(4) == variables.Cb(4)
variables.Gb(variables.E * variables.mu * variables.v) == variables.ea(-variables.Phi/2) * variables.Gb(variables.mu * variables.v) == variables.ea(-variables.Phip/2) * variables.Gbp(variables.mu * variables.v)
integrate(variables.Bbp(2), variables.M) == integrate((variables.Bb(2) * variables.d + variables.Cb(2) * variables.c), variables.M)
variables.tauab(2, (variables.p, variables.q)) == variables.lab(-4, 0) * (variables.MOa(-1))**(variables.i * variables.j) * variables.qb(variables.i) * variables.qb(variables.j) == variables.lab(-4, 0) * variables.brackets(variables.ea(variables.Phi) * (variables.p + variables.Cb(0) * variables.q)**2 + variables.ea(-variables.Phi) * variables.qa(2))
variables.Cb(0) == variables.Cb(0) + variables.b 
integrate(variables.da(10) * variables.x * (-variables.G)**(1/2) * variables.ea(-2 * variables.Phi) * (variables.R + 4 * diff(variables.Phi, variables.mu) * diff(variables.Phi, variables.z, variables.mu))) - 1/2 * integrate(variables.ea(2 * variables.alpha * variables.Phi) * Abs(variables.Fb(variables.q))**2)
integrate(variables.Fb(variables.q), variables.Sb(variables.q)) == variables.Q 
variables.d * variables.starf(variables.ea(2 * variables.alpha * variables.Phi) * variables.Fb(variables.q)) == 0
integrate(variables.starf(variables.ea(2 * variables.alpha * variables.Phi) * variables.Fb(variables.q)), variables.Sb(10 - variables.q)) == variables.Qp 
variables.Gb(variables.m * variables.n) == variables.ea(2 * variables.Phi) * variables.deltab(variables.m * variables.n) 
variables.Gb(variables.mu * variables.v) == variables.etab(variables.mu * variables.v) 
variables.Hb(variables.m * variables.n * variables.p) == - variables.epsilonba(variables.m * variables.n * variables.p, variables.q) * diff(variables.Phi, variables.q)
variables.ea(2 * variables.Phi) == variables.ea(2 * variables.Phif(math.inf)) + (variables.Q)/(2 * math.pi**2 * variables.ra(2))
variables.taub(variables.NS, 5) == (2 * math.pi**2 * variables.SlopeRegge)/(variables.kappaa(2)) == 1/((2 * math.pi)**5 * variables.ga(2) * variables.SlopeRegge**3)
variables.ea(2 * variables.Phi) == variables.ea(2 * variables.Phif(math.inf)) + 1/(2 * math.pi**2) * summation((variables.Qb(variables.i))/((variables.x - variables.xbs(variables.i))**2), (variables.i, 1, variables.N))
variables.gab(2, variables.D3) == 2 * math.pi * variables.g 
variables.gab(2, variables.D3) == (4 * math.pi**2)/(variables.gab(2, variables.D3))
1/(4 * math.pi) * integrate(variables.Cb(0) * variables.Trf(WedgeProduct(variables.Fb(2), variables.Fb(2))))
-1/(2 * variables.gab(2, variables.D3)) * integrate(variables.da(4) * variables.x * variables.Trf(Abs(variables.Fb(2))**2)) + variables.theta/(8 * math.pi**2) * integrate(variables.Trf(WedgeProduct(variables.Fb(2), variables.Fb(2))))
variables.B27 == variables.B10 + variables.B16 + variables.B1 
# check duality chains 2.189.2-3 
(variables.n * variables.m > 0)
(variables.N, variables.Ntilde) == (variables.n * variables.m, 0)
(variables.m * variables.n < 0)
(variables.N * variables.Ntilde) == (0, - variables.n * variables.m)
variables.Trf(variables.qa(variables.N)) == 2 ** 8 * np.prod(((1 + variables.qa(variables.k))/(1 - variables.qa(variables.k)))**8, (variables.k, 1, math.inf))
-variables.Beta == - math.exp(variables.brackets(math.pi * variables.i * (variables.sb(1) + variables.sb(2) + variables.sb(3) + variables.sb(4))))
variables.dirac(variables.sb(0), variables.i)
variables.Gb(variables.I * variables.mu * variables.v) == variables.ea(-variables.Phib(variables.h)) * variables.Gb(variables.h * variables.mu * variables.v)
variables.Phib(variables.I) == - variables.Phib(variables.h)
variables.Fbtilde(variables.I, 3) == variables.Hbtilde(variables.h, 3) 
variables.Ab(variables.I, 1) == variables.Ab(variables.h, 1)
variables.Lamdaab(variables.i, 0) == variables.La(-1/2) * quad(diff(variables.xas(1)) * variables.Lamdaaf(variables.i, variables.xas(1)), (0, variables.L))
variables.bracketc(variables.Lamdaab(variables.i, 0), variables.Lamdaab(variables.j, 0)) == variables.deltaa(variables.deltaa(variables.i * variables.j))
variables.taub(variables.D1) == (math.pi**(1/2))/(2**(1/2) * variables.kappa) * (4 * math.pi**2 * variables.SlopeRegge) == (variables.gab(2, variables.YM))/(8 * math.pi * variables.kappaa(2)) # type I
(math.pi**2 * variables.SlopeRegge**2)/(np.cross(2, math.factorial(4) * variables.gab(2, variables.YM))) * (variables.t * variables.Fa(4)) == (variables.gab(2, variables.YM))/(2**10 * math.pi**5 * math.factorial(4) * variables.kappaa(2)) * (variables.t * variables.Fa(4))
1/(2**8 * math.pi * math.factorial(4) * variables.SlopeRegge) * (variables.t * variables.Fa(4)) == (variables.gab(2, variables.YM))/(2**10 * math.pi**5 * math.factorial(4) * variables.kappaa(2)) * (variables.t * variables.Fa(4))
variables.psiab(variables.mu, -1/2) * variables.dirac((0, variables.k), variables.i * variables.j) * variables.lamdab(variables.i * variables.j) 
variables.psiab(variables.m, -1/2) * variables.dirac((0, variables.k),variables.i * variables.j) * variables.lamdabp(variables.i * variables.j)
variables.M * variables.lamda * variables.Ma(-1) == - variables.lamdaa(variables.T) 
variables.M * variables.lamdap * variables.Ma(-1) == variables.lamdapa(variables.T)
variables.lamda == variables.sigmaa(variables.a) 
variables.lamdap == variables.I
variables.hat(variables.omega) * variables.dirac(variables.psi, variables.i * variables.j) == variables.gammab(variables.j * variables.jp) * variables.dirac(variables.omega * variables.psi, variables.jp * variables.ip) * variables.gammaab(-1, variables.ip * variables.i)
variables.V == variables.ea(variables.i * (variables.Hb(3) + variables.Hb(4))/2) 
variables.holderb(variables.psia(6) + variables.i * variables.psia(7), -1/2) * variables.holderb(variables.psia(8) + variables.i * variables.psia(9), -1/2) * variables.dirac(0)
variables.gammaab(variables.T, 9) * variables.gammaab(-1, 9) == variables.omegaab(2, (5, -9)) * variables.gammaab(variables.T, 5) * variables.gammaab(-1, 5)
variables.taub(0) == 1/(variables.g * variables.SlopeRegge**(1/2)) 
variables.n * variables.taub(0) == variables.n/(variables.g * variables.SlopeRegge**(1/2))
variables.Rb(10) == variables.g * variables.SlopeRegge**(1/2)
variables.kappaab(2, 11) == 2 * math.pi * variables.Rb(10) * variables.kappaa(2) == 1/2 * (2 * math.pi)**8 * variables.ga(3) * variables.SlopeRegge**(9/2)
variables.Mb(11) == variables.ga(-1/3) * variables.SlopeRegge**(-1/2) 
variables.g == (variables.Mb(11) * variables.Rb(10))**(3/2) 
variables.SlopeRegge == variables.Mab(-3, 11) * variables.Rab(-1, 10)
variables.d == 9
variables.U == variables.SLf(2, variables.ZB)
variables.d == 8 
variables.U == np.cross(variables.SLf(3, variables.ZB), variables.SLf(2, variables.ZB))
variables.d == 6 
variables.U == variables.S * variables.Of(5, 5, variables.ZB)
variables.taub(variables.D2) == 1/((2 * math.pi)**2 * variables.g  * variables.SlopeRegge**(3/2)) == (variables.Mab(3, 11))/((2 * math.pi)**2)
variables.taub(variables.F1) == 1/(2 * math.pi* variables.SlopeRegge) == 2 * math.pi * variables.Rb(10) * variables.taub(variables.D2)
variables.Sf(variables.brackets(variables.F, variables.lamda, variables.X)) == - variables.taub(2) * integrate(variables.da(3) * variables.x * variables.bracketc(variables.brackets(-sp.Determinant(variables.etab(variables.mu * variables.v) + diff(variables.xa(variables.m), variables.mu) * diff(variables.xa(variables.m), variables.v) + 2* math.pi * variables.SlopeRegge * variables.Fb(variables.mu * variables.v)))**(1/2) + (variables.epsilona(variables.mu * variables.v, variables.rho))/2 * variables.lamda * diff(variables.Fb(variables.v, variables.rho), variables.mu)))
variables.Sf(variables.brackets(variables.lamda, variables.X)) == -variables.taub(2) * integrate(variables.da(3) * variables.x * variables.bracketc(-sp.Determinant(variables.etab(variables.mu * variables.v) + diff(variables.xa(variables.m), variables.mu) * diff(variables.xa(variables.m), variables.v) + (2 * math.pi * variables.SlopeRegge)**(-2) * diff(variables.lamda, variables.mu) * diff(variables.lamda, variables.v)))**(1/2))
variables.delta * variables.psib(variables.M) == variables.Dab(variables.neg, variables.M) * variables.zeta 
variables.delta * variables.psibtilde(variables.M) == variables.Dab(variables.pos, variables.M) * variables.zeta 
variables.S * variables.Of(9, 1) == np.cross(variables.S * variables.Of(5, 1), variables.S * variables.Of(4))
variables.B16 == (variables.B4, variables.B2) + (variables.B4p, variables.B2p)
variables.B16p == (variables.B4, variables.B2p) + (variables.B4p, variables.B2) 
variables.taub(variables.NS, 5) == 1/((2 * math.pi)**2 * variables.ga(2) * variables.SlopeRegge**3) == (variables.tauab(2, variables.D2))/(2 * math.pi) == (variables.Mab(6, 11))/((2 * math.pi)**5)
variables.taub(variables.D2) == variables.taub(variables.M2) 
variables.taub(variables.NS, 5) == variables.taub(variables.M5)
variables.taub(variables.D4) == 2 * math.pi * variables.Rb(10) * variables.taub(variables.M5) 
variables.taub(1) == variables.r * variables.taub(variables.M2) 
variables.Rbp(9) == variables.Rab(-1, 9) 
variables.gp == variables.g * variables.Rab(-1, 9)
1/(variables.ga(2)) * integrate(variables.da(10) * variables.x) == (2 * math.pi * variables.Rb(9))/(variables.ga(2)) * integrate(variables.da(9) * variables.x)
variables.gb(variables.I) * variables.gpa(-1) == variables.ga(-1) * variables.Rb(9) 
variables.Rb(9 * variables.I) == variables.gpa(-1/2) * variables.Rbp(9) == variables.ga(-1/2) * variables.Rab(-1/2, 9) 
variables.gb(variables.Ip) == variables.gb(variables.I) * variables.Rab(-1, 9 * variables.I) == variables.ga(-1/2) * variables.Rab(3/2, 9) 
variables.Rb(9 * variables.Ip) == variables.Rab(-1, 9 * variables.I) == variables.ga(1/2) * variables.Rab(1/2, 9)
variables.Rb(10 * variables.M) == variables.gab(2/3, variables.Ip) == variables.ga(-1/3) * variables.Rb(9) 
variables.Rb(9 * variables.M) == variables.gab(-1/3, variables.Ip) * variables.Rb(9 * variables.Ip) == variables.ga(2/3)
# check duality chains 2.207.6-8

variables.E - variables.n/(variables.Rb(10)) == variables.Of(variables.Rb(10)/variables.n)
variables.H == variables.Rb(10) * variables.Trf(1/2 * variables.pb(variables.i) * variables.pb(variables.i) - (variables.Mab(6, 11))/(16 * math.pi**2) * variables.brackets(variables.xa(variables.i), variables.xa(variables.j))**2 - (variables.Mab(3, 11))/(4 * math.pi) * variables.lamda * variables.Rhoa(0) * variables.Rhoa(variables.i) * variables.brackets(variables.xa(variables.i), variables.lamda))
variables.E == (variables.Rb(10))/(2) * variables.Trf(variables.pb(variables.i) * variables.pb(variables.i)) == (variables.qa(2))/(2 * variables.pb(10))
variables.xa(variables.i) == variables.xab(variables.i, 0) + variables.xas(variables.i) 
variables.xab(variables.i, 0) == variables.Yab(variables.i, 1) * variables.Ib(1) + variables.Yab(variables.i, 2) * variables.Ib(2)
variables.xas(variables.i) == variables.xabs(variables.i, 11) + variables.xabs(variables.i, 22) + variables.xabs(variables.i, 12) + variables.xabs(variables.i, 21)
variables.psif(variables.xbs(11), variables.xbs(22)) == variables.psibf(0, variables.xbs(11)) * variables.psibf(0, variables.xbs(22))
variables.brackets(variables.xab(variables.i, 0), variables.xabs(variables.j, 12)) == (variables.Yab(variables.i, 1) - variables.Yab(variables.i, 2)) * variables.xabs(variables.j, 12)
variables.Lb(variables.eff) == -variables.Vf(variables.r, variables.v) == 4 * math.pi**(5/2) * variables.Rhof(7/2) * variables.SlopeRegge**3 * variables.nb(1) * variables.nb(2) * (variables.va(4))/(variables.ra(7)) == (15 * math.pi**3)/2 * (variables.pb(10) * variables.ppb(10))/(variables.Mab(9, 11) * variables.Rb(10)) * (variables.va(4))/(variables.ra(7))
# check matrcies 2.215.9
variables.Ua(variables.n) == variables.Va(variables.n) == 1 
variables.U * variables.V == variables.alpha * variables.V * variables.U 
variables.xa(variables.i) == summation(variables.xab(variables.i, variables.r * variables.s) * variables.Ua(variables.r) * variables.Va(variables.s), ((variables.r, variables.s), variables.brackets(1 - variables.n/2), variables.brackets(variables.n/2))) 
variables.xa(variables.i) == variables.Xaf(variables.i, (variables.p, variables.q)) == summation(variables.xab(variables.i, variables.r * variables.s) * math.exp(variables.i * variables.p * variables.r + variables.i * variables.q * variables.s), ((variables.r, variables.s), variables.brackets(1 - variables.n/2), variables.brackets(variables.n/2)))
variables.brackets(variables.xa(variables.i), variables.xa(variables.j)) == (2 * math.pi * variables.i)/variables.n * (diff(variables.xa(variables.i), variables.q) * diff(variables.xa(variables.j), variables.p)  - diff(variables.xa(variables.i), variables.p) * diff(variables.xa(variables.j), variables.q)) + variables.Of(variables.na(-2)) == (2 * math.pi * variables.i)/(variables.n) * variables.bracketcb(variables.xa(variables.i), variables.xa(variables.j), (variables.PB)) + variables.Of(variables.na(-2))
variables.Tr == variables.n * integrate((diff(variables.q) * diff(variables.p))/((2 * math.pi)**2))
variables.Rb(10) * integrate(diff(variables.q) * diff(variables.p) * ((variables.n)/(8 * math.pi**2) * variables.Pib(variables.i) * variables.Pib(variables.i) + (variables.Mab(6, 11))/(16 * math.pi**2 * variables.n) * variables.bracketcab((variables.xa(variables.i), variables.xa(variables.j)), 2, variables.PB) - variables.i * (variables.Mab(3, 11))/(8 * math.pi**2) * variables.lamda * variables.Rhoa(0) * variables.Rhoa(variables.i) * variables.bracketc((variables.xa(variables.i), variables.lamda), variables.PB)))
variables.xa(1) == variables.a * variables.q 
variables.xa(2) == variables.b * variables.p 
(variables.Mab(6, 11) * variables.Rb(10) * variables.aa(2) * variables.ba(2))/(2 * variables.n) == (variables.Mab(6, 11) * variables.Aa(2))/(2 * (2 * math.pi)**4 * variables.pb(10)) == (variables.tauab(2, variables.M2) * variables.Aa(2))/(2 * variables.pb(10))
variables.xa(variables.i) == variables.Yab(variables.i, 1) * variables.Ib(1) + variables.Yab(variables.i, 2) * variables.Ib(2)
(variables.xas(0), variables.xas(10)) == (variables.xas(0) - math.pi * variables.Rbtilde(10), variables.xas(10) + math.pi * variables.Rbtilde(10))
(variables.xas(0), variables.xas(10)) == (variables.xas(0) - math.pi * variables.Rbtilde(10), variables.xas(10) + math.pi * variables.Rbtilde(10) + 2 * math.pi * variables.epsilona(2) * variables.Rbtilde(10))
(variables.xaps(0), variables.xaps(10)) == (variables.xaps(0), variables.xaps(10) + 2 * math.pi * variables.epsilon * variables.Rbtilde(10))
variables.xaps(0) + variables.posneg * variables.xaps(10) == variables.epsilona(variables.negpos * 1) * (variables.xas(0) + variables.posneg * variables.xas(10))
variables.Ep, variables.ppb(10) == variables.Of(variables.epsilona(-1))
variables.Ep - variables.ppb(10) == variables.Of(variables.epsilon)
variables.Rb(10) == variables.epsilon * variables.Rbtilde(10)
(variables.Rb(10))/(variables.SlopeRegge**(1/2)) * np.prod((variables.SlopeRegge**(1/2))/(variables.Rb(variables.m)), variables.m) == variables.Rab((3 - variables.k)/2, 10) * (variables.Mab(3, 11))**((1- variables.k)/2) * np.prod(variables.Rab(-1, variables.m), variables.m)
variables.ds2 == variables.Zf(variables.r)**(-1/2) * variables.etab(variables.mu * variables.v) * diff(variables.xas(variables.mu)) * diff(variables.xas(variables.v)) + variables.Zf(variables.r)**(1/2) * diff(variables.xas(variables.m)) * diff(variables.xas(variables.m)) 
variables.ea(2 * variables.Phi) == variables.ga(2) * variables.Zf(variables.r)**((3 - variables.p)/2)
variables.Zf(variables.r) == 1 + (variables.rhoa(7 - variables.p))/(variables.ra(7 - variables.p))
variables.ra(2) == variables.xas(variables.m) * variables.xas(variables.m)
variables.rhoa(7 - variables.p) == variables.SlopeRegge**((7 - variables.p)/2) * variables.g * variables.Q * (4 * math.pi)**((5 - variables.p)/2) * variables.Rho * (7 - variables.p)/2 
variables.ds2 == variables.Zab(-1/2, 1) * variables.Zab(-1/2, 5) * variables.brackets(variables.etab(variables.mu * variables.v) * diff(variables.xas(variables.mu)) * diff(variables.xas(variables.v)) + (variables.Zb(variables.n) - 1) * (diff(variables.t) + diff(variables.xbs(5)))**2) + variables.Zab(1/2, 1) * variables.Zab(1/2, 5) * diff(variables.xas(variables.i)) * diff(variables.xas(variables.i)) + variables.Zab(1/2, 1)* variables.Zab(-1/2, 5) * diff(variables.xas(variables.m)) * diff(variables.xas(variables.m))
variables.ea(-2 * variables.Phi) == (variables.Zb(5))/(variables.Zb(1))
variables.Zb(1) == 1 + (variables.rab(2, 1))/(variables.ra(2)) 
variables.rab(2, 1) == ((2 * math.pi)**4 * variables.g * variables.Qb(1) * variables.SlopeRegge**3)/(variables.Vb(4))
variables.Zb(5) == 1 + (variables.rab(2, 5))/(variables.ra(2)) 
variables.rab(2, 5) == variables.g * variables.Qb(5) * variables.SlopeRegge 
variables.Zb(variables.n) == 1 + (variables.rab(2, variables.n))/(variables.ra(2)) 
variables.rab(2, variables.n) == ((2 * math.pi)**5 * variables.ga(2) * variables.pb(5) * variables.SlopeRegge**4)/(variables.L * variables.Vb(4))
diff(variables.sab(2, variables.E)) == variables.Zab(-3/4, 1) * variables.Zab(-1/4, 5) * variables.brackets(variables.etab(variables.mu * variables.v) * diff(variables.xas(variables.mu)) * diff(variables.xas(variables.v)) + (variables.Zb(variables.n) - 1) * (diff(variables.t) + diff(variables.xbs(5)))**2) + variables.Zab(1/4, 1) * variables.Zab(3/4, 5) * diff(variables.xas(variables.i)) * diff(variables.xas(variables.i)) + variables.Zab(1/4, 1) * variables.Zab(-1/4, 5) * diff(variables.xas(variables.m)) * diff(variables.xas(variables.m))
((variables.rab(2, 1))/(variables.ra(2)))**(1/4) * ((variables.rab(2, 5))/(variables.ra(2)))**(3/4) * variables.ra(2) * diff(variables.omegaa(2)) == variables.rab(1/2, 1) * variables.rab(3/2, 5) * diff(variables.omegaa(2))
variables.A == 2 * math.pi**2 * variables.L * variables.Vb(4) * variables.rb(1) * variables.rb(5) * variables.rb(variables.n) == 2**6 * math.pi**7 * variables.ga(2) * variables.SlopeRegge**4 * (variables.Qb(1) * variables.Qb(5) * variables.nb(5))**(1/2) == variables.kappaa(2) * (variables.Qb(1) * variables.Qb(5) * variables.nb(5))**(1/2)
variables.S == (2 * math.pi * variables.A)/(variables.kappaa(2)) == 2 * math.pi * (variables.Qb(1) * variables.Qb(5) * variables.nb(5))**(1/2)
variables.V == 1/((2 * math.pi * variables.SlopeRegge)**2) * Abs(variables.xb(variables.i * variables.chi) -variables.chi * variables.Yb(variables.i))**2 + (variables.gab(2, 1))/4 * (variables.Dab(variables.A, 1)) * variables.Dab(variables.A, 1) + (variables.gab(2, 5))/(4 * variables.Vb(4)) * variables.Dab(variables.A, 5) * variables.Dab(variables.A, 5) # double check 2.221.9
variables.xa(variables.i) == variables.xas(variables.i) * variables.Ib(variables.Qb(1))
variables.Ya(variables.i) == variables.xas(variables.i) * variables.Ib(variables.Qb(5))
quad(diff(variables.E) * variables.nf(variables.E) * math.exp(-variables.Beta * variables.E), 0, math.inf) == variables.Trf(variables.brackets(math.exp(-variables.Beta * variables.H)))
variables.ds2 == (variables.ra(2))/(variables.rb(1) * variables.rb(5)) * variables.etab(variables.mu * variables.v) * diff(variables.xas(variables.mu)) * diff(variables.xas(variables.v)) + (variables.rb(1) * variables.rb(5))/(variables.ra(2)) * diff(variables.ra(2)) + variables.rb(1) * variables.rb(5) * diff(variables.omegaa(2)) + (variables.rb(1))/(variables.rb(5)) * diff(variables.xas(variables.m)) * diff(variables.xas(variables.m))
np.cross(variables.A * diff(variables.Sb(3)), variables.Sa(3), variables.Ta(4)) 
np.cross(variables.A * diff(variables.Sb(5)), variables.Sa(5))
math.exp(variables.bracketc(math.pi * variables.M * variables.brackets((variables.c + variables.ctilde) * variables.SlopeRegge/3))**(1/2))

#2.293.1
variables.SBb(variables.het) == 1/(2 * variables.kappaab(2, 10)) * integrate(variables.da(10) * variables.x * (-variables.G)**(1/2) * variables.ea(-2 * variables.PhiB) * variables.brackets(variables.R + 4 * diff(variables.PhiB, variables.mu) * diff(variables.PhiB, 1, variables.mu) - 1/2 * Abs(variables.Hbtilde(3))**2 - (variables.SlopeRegge)/4 * variables.Trbf(variables.v, Abs(variables.Fb(2))**2)))
variables.Hbtilde(3) == diff(variables.Bb(2)) - (variables.SlopeRegge)/4 * variables.Trbf(variables.v, WedgeProduct(variables.Ab(1), diff(variables.Ab(1))) - WedgeProduct(2 * variables.i * variables.Ab(1), variables.Ab(1), variables.Ab(1)/3))
variables.Gb(variables.mu * variables.v)
variables.Bb(variables.mu * variables.v) 
variables.PhiB
variables.Gb(variables.i * variables.jl) 
variables.Bb(variables.i * variables.jl)
variables.Aab(variables.a, variables.mu) 
variables.Ab(variables.i, variables.jl * variables.xl)
variables.Ab(variables.il, variables.j * variables.x)
variables.SB == 1/(2 * variables.kappaa(2, 4)) * integrate(variables.da(4) * variables.x * (-variables.G)**(1/2) * variables.brackets(variables.R - 2 * diff(variables.PhiBb(4), variables.mu) * diff(variables.PhiBb(4), 1, variables.mu) - 1/2 * variables.ea(-4 * variables.PhiBb(4)) * Abs(variables.Hb(3))**2 - 1/2 * variables.Ga(variables.i * variables.jl) * variables.Ga(variables.k * variables.ll) * (diff(variables.Gb(variables.i * variables.ll), variables.mu) * diff(variables.Gb(variables.jl * variables.k), 1, variables.mu) + diff(variables.Bb(variables.i * variables.ll), variables.mu) * diff(variables.Bb(variables.jl * variables.k), 1, variables.mu))))
variables.PhiBb(4) == variables.PhiB - 1/4 * sp.Determinant(variables.Gb(variables.m * variables.n))
variables.Gb(variables.mu * variables.v * variables.Einstein) == variables.ea(- 2 * variables.PhiBb(4)) * variables.Gb(variables.mu * variables.v)
variables.Gb(variables.i * variables.jl) == variables.Gb(variables.jl * variables.i) == variables.star(variables.Gb(variables.j * variables.il)) == variables.star(variables.Gb(variables.il * variables.j))
variables.Bb(variables.i * variables.jl) == - variables.Bb(variables.jl * variables.i) == - variables.star(variables.Bb(variables.j * variables.il)) == variables.star(variables.Bb(variables.il * variables.j))
-1/2* integrate(variables.da(4) * variables.x * (-variables.G)**(1/2) * variables.ea(-4 * variables.PhiBb(4)) * Abs(variables.Hb(3))**2) + integrate(variables.a * variables.diff(variables.Hb(3))) == -1/2 * integrate(variables.da(4) * variables.x * (-variables.G)**(1/2) * variables.ea(4 * variables.PhiBb(4)) * diff(variables.a, variables.mu) * diff(variables.a, 1, variables.mu))
1/(2 * variables.kappaab(2, 4)) * integrate(variables.da(4) * variables.x * (-variables.G)**(1/2) * variables.brackets(variables.R - (2 * diff(variables.star(variables.S), variables.mu) * diff(variables.S, 1, variables.mu))/((variables.S + variables.star(variables.S))**2) - 1/2 * variables.Ga(variables.i * variables.jl) * variables.Ga(variables.k * variables.ll) * diff(variables.Tb(variables.i * variables.ll), variables.mu) * diff(variables.Tb(variables.jl * variables.k), 1, variables.mu)))
variables.S == variables.ea(-2 * variables.PhiBb(4)) + variables.i * variables.a 
variables.Tb(variables.i * variables.jl) == variables.Gb(variables.i * variables.jl) + variables.Bb(variables.i * variables.jl)
variables.kappaab(2, 4) * variables.K == -ln(variables.S + variables.star(variables.S)) - ln(sp.Determinant(variables.Tb(variables.i * variables.jl) + variables.star(variables.Tb(variables.i * variables.jl))))
-1/2 * integrate(variables.da(4) * variables.x * (-variables.G)**(1/2) * variables.ea(-4 * variables.PhiBb(4)) * Abs(variables.Hbtilde(3))**2) + integrate(variables.a * variables.brackets(diff(variables.Hbtilde(3)) + (variables.SlopeRegge)/4 * variables.Trbf(variables.v, WedgeProduct(variables.Hb(2), variables.Fb(2)))))
-1/(4 * variables.gab(2, 4)) * integrate(variables.ea(-2 * variables.PhiBb(4)) * variables.Trbf(variables.v, Abs(variables.Fb(2))**2)) + 1/(2 * variables.gab(2, 4)) * integrate(variables.v, WedgeProduct(variables.Fb(2), variables.Fb(2)))
variables.fb(variables.a * variables.b) == (variables.deltab(variables.a * variables.b))/(variables.gab(2, 4)) * variables.S 
variables.kappaab(2, 4) * variables.K == -ln(variables.S + variables.star(variables.S)) - ln(sp.Determinant(variables.brackets(variables.Tb(variables.i * variables.jl) + variables.star(variables.Tb(variables.i * variables.jl)) - variables.SlopeRegge * variables.Trbf(variables.v, variables.star(variables.Ab(variables.i) * variables.Ab(variables.j))))))
variables.W == variables.epsilona(variables.i, variables.j * variables.k) * variables.Trbf(variables.v, variables.Ab(variables.i) * variables.brackets(variables.Ab(variables.j), variables.Ab(variables.k)))
variables.W == variables.epsilona(variables.i, variables.j * variables.k) * variables.epsilona(variables.ll, variables.ml, variables.nl) * variables.da(variables.xl * variables.yl * variables.zl) * variables.Ab(variables.i, variables.ll * variables.xl) * variables.Ab(variables.j, variables.ml * variables.yl) * variables.Ab(variables.k, variables.nl * variables.zl)
variables.Tb(variables.i * variables.jl) == variables.Tb(variables.i) * variables.deltab(variables.i * variables.jl) # no sum on i
variables.kappaab(2, 4) * variables.K == - ln(variables.S + variables.star(variables.S)) - summation(ln(variables.Tb(variables.i) + variables.star(variables.Tb(variables.i))), variables.i) + variables.SlopeRegge * summation((variables.Trbf(variables.v, variables.Ab(variables.i) * variables.star(variables.Ab(variables.i))))/(variables.Tb(variables.i) + variables.star(variables.Tb(variables.i))), variables.i)
variables.Tb(variables.i) == (variables.ab(variables.i) * variables.Tb(variables.i) - variables.i * variables.bb(variables.i))/(variables.i * variables.cb(variables.i) * variables.Tb(variables.i) + variables.db(variables.i))
variables.ab(variables.i) * variables.db(variables.i) - variables.bb(variables.i) * variables.cb(variables.i) == 1
variables.Tb(variables.i) + variables.star(variables.Tb(variables.i)) == (variables.tb(variables.i) + variables.star(variables.Tb(variables.i)))/(Abs(variables.i * variables.cb(variables.i) * variables.Tb(variables.i0 + variables.db(variables.i)))**2)
variables.kappaab(2, 4) * variables.K == variables.kappaab(2, 4) * variables.K + variables.Ref(variables.brackets(summation(ln(variables.i * variables.cb(variables.i) * variables.Tb(variables.i) + variables.db(variables.i)), variables.i)))
variables.Ab(variables.i) == (variables.Ab(variables.i))/(variables.i * variables.cb(variables.i) * variables.Tb(variables.i) + variables.db(variables.i))
variables.W == variables.W/(np.prod(variables.i * variables.cb(variables.i) * variables.Tb(variables.i) + variables.db(variables.i), (variables.i, 2, 4)))
(variables.S * variables.Uf(3, 3))/(np.cross(variables.S * variables.Uf(3), variables.S * variables.Uf(3), variables.S * variables.Uf(3, 3, variables.ZB)))
(variables.S * variables.Uf(1, 1))/(np.cross(variables.Uf(1), variables.P * variables.S * variables.Lf(2, variables.ZB)))
variables.lamdaab(variables.pos * variables.K, -1/6) * variables.alphaab(variables.jl, -1/3) * variables.diracb(0, (variables.NS, variables.NS))
variables.K == 1, 2, 3
variables.Da(variables.a) == variables.star(variables.Mb(variables.K * variables.jl)) * variables.tab(variables.a, variables.K * variables.L) * variables.Mb(variables.L * variables.jl) == variables.Trf(variables.dagger(variables.M) * variables.ta(variables.a) * variables.M)
variables.M * variables.dagger(variables.M) == variables.rhoa(2) * variables.I # ->
variables.M == variables.rho * variables.U 
variables.Cb(variables.alpha) * variables.star(variables.Cb(variables.alpha)) * np.prod((variables.tb(variables.i) + variables.star(variables.tb(variables.i)))**(variables.nab(variables.i, variables.alpha)), (variables.i, 2, 4))
variables.Cb(variables.alpha) == variables.Cb(variables.alpha) * np.prod((variables.i * variables.cb(variables.i) * variables.Tb(variables.i) + variables.db(variables.i))**(variables.nab(variables.i, variables.alpha)), (variables.i, 2, 4)) 
1/(variables.gabf(2, variables.a, variables.mu)) == (variables.S * variables.kb(variables.a))/(variables.gab(2, 4)) + (variables.bb(variables.a))/(16 * math.pi**2) * ln((variables.mab(2, variables.S * variables.U))/(variables.mu**(2))) + 1/(16 * math.pi**2) * variables.Deltabtilde(variables.a)
variables.Betab(variables.a) == (variables.bb(variables.a) * variables.gab(3, variables.a))/(16 * math.pi**2)
variables.Deltabtilde(variables.a) == variables.Deltab(variables.a) + 16 * math.pi**2 * variables.kb(variables.a) * variables.Y 
variables.Deltab(variables.a) == integrate((variables.d2(variables.tau))/(variables.taub(2)) * variables.brackets(variables.BObf(variables.a, (variables.tau, variables.taul)) - variables.bb(variables.a)))
variables.mb(variables.S * variables.U) == (2 * math.exp(variables.brackets((1 - variables.gamma)/2)))/(3**(3/4) * (2 * math.pi * variables.SlopeRegge)**(1/2))
# check conditions 2.299.37
variables.Za(2) == variables.neg * variables.Za(2) 
variables.Za(3) == variables.neg * variables.Za(3) 
variables.Za(4) == variables.pos * variables.Za(4)
# - 
variables.Deltab(variables.a) == variables.cb(variables.a) - summation((variables.bab(variables.i, variables.a) * Abs(variables.Pa(variables.i)))/(Abs(variables.P)) * variables.bracketc(ln(variables.brackets((variables.Tb(variables.i) + variables.star(variables.Tb(variables.i))) * Abs(variables.Tb(variables.i))**4)) + ln(variables.brackets((variables.Ub(variables.i) + variables.star(variables.Ub(variables.i))) * Abs(variables.Ub(variables.i))**4))), variables.i)
ln(variables.brackets((variables.Tb(variables.i) + variables.star(variables.Tb(variables.i))) * Abs(variables.etaf(variables.Tb(variables.i)))**4)) == ln(variables.Tb(variables.i) + variables.star(variables.Tb(variables.i))) + 4 * variables.Ref(variables.brackets(ln(variables.etaf(variables.Tb(variables.i)))))
