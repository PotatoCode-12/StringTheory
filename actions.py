import math
from scipy.integrate import quad
import numpy as np 
import sympy as sp
from sympy import symbols, diff,integrate, Function, Derivative, summation, oo, Abs, ln
from sympy.diffgeom import WedgeProduct
import variables
import constants

ActionNambuGoto = symbols('e')
AppNambuGoto = symbols('e')
def two_path_monotomic():
    variables.xa(variables.mu, variables.t) == variables.xap(variables.mu, variables.taup, variables.tau) #for monotomic functions of time
def NambuGoto():
    ActionNambuGoto == quad(np.diff(variables.t)*np.diff(variables.sigma)*-1/(2*math.pi*variables.SlopeRegge)*(np.linalg.det(-np.array(variables.h(variables.a,variables.b))))**(1/2), variables.WorldSheet, math.inf) 
    AppNambuGoto == Function('AppNambuGoto')(variables.x, variables.GWorldSheetMetric) == -1/(4*math.pi*variables.SlopeRegge)*quad(np.diff(variables.tau)*np.diff(variables.sigma)*(-variables.GWorldSheetMetric)**(1/2)*variables.GWorldSheetMetric**(variables.a*variables.b)*diff(variables.xa(variables.mu), variables.a)*diff(variables.xb(variables.mu), variables.b), variables.WorldSheet, math.inf) # with set GWorldSheetMetric


variables.AppNambuGotoP == AppNambuGoto - variables.lamdas * variables.chi == -quad(np.diff(variables.tau)*np.diff(variables.sigma)*(-variables.GWorldSheetMetric)**(1/2)(1/(4*math.pi*variables.SlopeRegge)*variables.GWorldSheetMetric**(variables.a*variables.b)*diff(variables.xa(variables.mu), variables.a)*diff(variables.xb(variables.mu), variables.b) + variables.lamdas/(4*math.pi)*variables.ScalarRicci), variables.WorldSheet, math.inf) #chi introduction with symetries retained
if variables.xa(variables.pos, variables.tau) == variables.tau :
    variables.ActionWorldLine == 1/2(quad(np.diff(variables.tau)(-2*variables.eta**(-1)*variables.xdiffaneg + variables.eta**(-1)*(2*variables.xdiff**(variables.i))-variables.eta*variables.m**2)))
    variables.pb(variables.mu) == diff(variables.L)/variables.xadiff(variables.mu)
    variables.pb(variables.neg) == -variables.eta**(-1)
    variables.pb(variables.i) == variables.eta**(-1)*variables.xadiff(variables.i)
    variables.H == variables.pb(variables.neg)*variables.xadiff(variables.neg) + variables.pb(variables.i) * variables.xadiff(variables.i) - variables.L == (variables.pa(variables.i)*variables.pa(variables.i) +variables.m**2)/(2*variables.pa(variables.pos))
def actionMasslessScalar():
    variables.S == 1/(4*math.pi*variables.SlopeRegge) * integrate((sp.diff(variables.tau)* sp.diff(variables.sigma))*(diff(variables.xa(variables.mu), variables.D1)*diff(variables.xb(variables.mu), variables.D1) + diff(variables.xa(variables.mu), variables.D2) * diff(variables.xb(variables.mu), variables.D2))) # for D scalar fields in two dimensions

variables.S == 1/(2*math.pi*variables.SlopeRegge) * integrate(variables.d2z*diff(variables.xa(variables.mu,), variables.z)*diff(variables.xb(variables.mu), variables.zl)) # pg 34 contour in complex cordinates

variables.S == 1/(2*math.pi) * integrate(variables.d2z*variables.b*diff(variables.c, variables.zl))
variables.hb(variables.b) == variables.lamda
variables.hb(variables.c) == 1- variables.lamda # 50

variables.S == 1/(2*math.pi) * integrate(variables.d2z * variables.btilde * diff(variables.ctilde)) # for antiholomorphic theory pg 51
variables.j == - variables.point(variables.b * variables.c) #noether current

variables.InvPoincare == 1/2 * integrate(diff(variables.tau) * (variables.eta**(-1) * variables.Gbf(variables.mu*variables.v, variables.X) * diff(variables.xa(variables.mu), variables.tau) * diff(variables.xa(variables.v), variables.tau) - variables.eta* variables.m**2)) # 108.1
variables.Sb(variables.sigma) == 1/(4*variables.SlopeRegge*math.pi)*integrate(variables.d2sigma * variables.g**(1/2) * variables.ga(variables.a*variables.b) * variables.Gbf(variables.mu * variables.v, variables.X) * diff(variables.xa(variables.mu), variables.a) * diff(variables.xa(variables.v), variables.b), variables.M)
variables.Gbf(variables.mu*variables.v, variables.X) == variables.etab(variables.mu*variables.v) + variables.chibf(variables.mu*variables.v, variables.X)
#reference excluded integrand 108.4
variables.chibf(variables.mu*variables.v, variables.X) == -4*math.pi*variables.gb(variables.c)*math.e**(variables.i*variables.k*variables.X) * variables.sb(variables.mu * variables.v)
variables.Sb(variables.sigma) == 1/(4*math.pi*variables.SlopeRegge) * integrate(variables.d2sigma * variables.g**(1/2), variables.M) * variables.brackets((variables.ga(variables.a*variables.b) * variables.Gbf(variables.mu * variables.v, variables.X) + variables.i*variables.epsilona(variables.a*variables.b) * variables.Bbf(variables.mu * variables.v, variables.X))*diff(variables.xa(variables.mu), variables.a) * diff(variables.xa(variables.v), variables.b) + variables.SlopeRegge * variables.R * variables.phif(variables.X))
variables.delta * variables.Bbf(variables.mu * variables.v, variables.X) == diff(variables.zetabf(variables.v, variables.X), variables.mu) - diff(variables.zetabf(variables.mu, variables.X), variables.v)
variables.Hb(variables.omegas * variables.mu * variables.v) == diff(variables.Bb(variables.mu * variables.v), variables.omegas) + diff(variables.Bb(variables.v * variables.omegas), variables.mu) + diff(variables.Bb(variables.omegas * variables.mu), variables.v)
#vertex operators under the action
variables.Gbf(variables.mu*variables.v, variables.X) == variables.etab(variables.mu * variables.v) - 4*math.pi*variables.gb(variables.c) * variables.sb(variables.mu * variables.v) * math.e**(variables.i*variables.k*variables.X)
variables.Bbf(variables.mu * variables.v, variables.X) == -4*math.pi*variables.gb(variables.c) * variables.ab(variables.mu * variables.v) * math.e**(variables.i*variables.k*variables.X)
variables.phif(variables.X) == -4*math.pi*variables.gb(variables.c) * variables.phi * math.e**(variables.i*variables.k*variables.X)
variables.Tab(variables.a, variables.a) == -1/(2*variables.SlopeRegge)  * variables.Betaab(variables.G, variables.mu*variables.v) * variables.ga(variables.a*variables.b) * diff(variables.xa(variables.mu), variables.a) * diff(variables.xa(variables.v), variables.b) - variables.i/(2*variables.SlopeRegge) * variables.Betaab(variables.B, variables.mu*variables.v) * variables.epsilona(variables.a*variables.b) * diff(variables.xa(variables.mu), variables.a) * diff(variables.xa(variables.v), variables.b) - 1/2 * variables.Betaa(variables.phi) * variables.R
variables.Betaab(variables.G, variables.mu* variables.v) == -variables.SlopeRegge/2 * (diff(variables.chib(variables.mu * variables.v)) - diff(diff(variables.chib(variables.mu*variables.omegas), 1, variables.omegas), variables.v) - diff(diff(variables.chib(variables.omegas * variables.v), 1, variables.omegas), variables.mu) + diff(diff(variables.chiab(variables.omegas, variables.omegas), variables.v), variables.mu) + 2*variables.SlopeRegge*diff(diff(variables.phi, variables.v), variables.mu))
variables.Betaab(variables.B, variables.mu * variables.v) == -variables.SlopeRegge/2 * diff(variables.Hb(variables.omegas, variables.mu, variables.v), 1, variables.omegas) # 111.13
variables.Betaa(variables.phi) == (variables.D-26)/6 - variables.SlopeRegge/2 * diff(variables.phi, variables.omegas, 2)
# after renormalizing to two spacetime derivatives with RB == R-bolded
variables.Betaab(variables.G, variables.mu * variables.v) == variables.SlopeRegge * variables.RBb(variables.mu * variables.v) + variables.SlopeRegge * 2 * variables. Deltab(variables.mu) * variables.Deltab(variables.v) * variables.phi - variables.SlopeRegge/4 * variables.Hb(variables.mu * variables.lamda * variables.omegas) * variables.Hba(variables.v, variables.lamda * variables.omegas) + variables.Of(variables.SlopeRegge**2)
variables.Betaab(variables.B, variables.mu * variables.v) == - variables.SlopeRegge/2 * variables.Deltaa(variables.omegas) * variables.Hb(variables.omegas * variables.mu * variables.v) + variables.SlopeRegge * variables.Deltaa(variables.omegas) * variables.phi * variables.Hb(variables.omegas * variables.mu * variables.v) + variables.Of(variables.SlopeRegge**2)
variables.Betaa(variables.phi) == (variables.D-26)/6 - variables.SlopeRegge/2 * variables.Delta**2 * variables.phi + variables.SlopeRegge*variables.Deltab(variables.omegas) * variables.phi * variables.Deltaa(variables.omegas) * variables.phi - variables.SlopeRegge/24 * variables.Hb(variables.mu * variables.v * variables.lamda) * variables.Ha(variables.mu * variables.v * variables.lamda) + variables.Of(variables.SlopeRegge**2)
variables.Betaab(variables.G, variables.mu*variables.v) == variables.Betaab(variables.B, variables.mu * variables.v) == variables.Betaa(variables.phi) == 0 # 112.15
#for points
variables.Gbf(variables.mu * variables.v, variables.X) == variables.etab(variables.mu * variables.v)
variables.Bbf(variables.mu * variables.v, variables.X) == 0
variables.phif(variables.X) == variables.phib(variables.D0)
variables.lamda == variables.phib(variables.D0) #extra condition
#additionally, if
variables.phif(variables.X) == variables.Vb(variables.mu) * variables.xa(variables.mu)
variables.Vb(variables.mu) * variables.Va(variables.mu) == (26-variables.D)/(6*variables.SlopeRegge)
#spacetime action
variables.SB == 1/(2*variables.kappaab(variables.D2, variables.D0)) * integrate(diff(variables.x, 1, variables.D) * (-variables.G)**(1/2) * math.e**(-2*variables.phi)) * variables.brackets(-(2*(variables.D-26))/(3*variables.SlopeRegge) + variables.RB - 1/12 * variables.Hb(variables.mu * variables.v * variables.lamda) * variables.Ha(variables.mu * variables.v * variables.lamda) + 4 * diff(variables.phi, variables.mu) * diff(variables.phi, 1, variables.mu) + variables.Of(variables.SlopeRegge))
variables.delta * variables.SB == -1/(2*variables.kappaab(variables.D2, variables.D0) * variables.SlopeRegge) * integrate(diff(variables.x, 1, variables.D) * (-variables.G)**(1/2) * math.e**(-2*variables.phi)) * variables.brackets(variables.delta * variables.Gb(variables.mu * variables.v) * variables.Betaa(variables.G * variables.mu *variables.v) + variables.delta * variables.Bb(variables.mu * variables.v) * variables.Betaa(variables.B * variables.mu * variables.v) + (2 * variables.delta * variables.phi - 1/2 * (variables.Ga(variables.mu * variables.v) * variables.delta(variables.Gb(variables.mu * variables.v))))*(variables.Betaabp(variables.G * variables.omegas, variables.omegas) - 4 * variables.Betaa(variables.phi)))
variables.Gbftilde(variables.mu * variables.v, variables.x) == math.exp(2 * variables.omegasf(variables.x)) * variables.Gbf(variables.mu * variables.v, variables.x)
variables.RBtilde == math.exp(-2*variables.omegas) * variables.brackets(variables.RB - 2*(variables.D-1) * variables.Deltaa(variables.D2) * variables.omegas - (variables.D - 2) *(variables.D - 1) * diff(variables.omegas, variables.mu) * diff(variables.omegas, 1, variables.mu))
variables.phitilde == variables.phi - variables.phib(variables.D0)
variables.SB == 1/(2 * variables.kappaa(variables.D2)) * integrate(diff(variables.X, 1, variables.D) * (-variables.G)**(1/2)) * variables.brackets(-(2*(variables.D-26))/(variables.SlopeRegge * 3) * math.e**(4*variables.phitilde/(variables.D-2)) + variables.RBtilde - 1/12 * math.e**(-8*variables.phitilde/(variables.D-2)) * variables.Hb(variables.mu * variables.v * variables.lamda) * variables.Hatilde(variables.mu * variables.v * variables.lamda) - 4/(variables.D-2) * diff(variables.phitilde, variables.mu) * diff(variables.phitilde, 1, variables.mutilde) + variables.Of(variables.SlopeRegge))
variables.kappa == (8 * math.pi * variables.Gb(variables.N)) **(1/2) == (8*math.pi)**(1/2)/(variables.Mb(variables.P)) # check constant value 114.26
variables.SB == -1/4 * integrate(diff(variables.x, 1, variables.D) * variables.Trf(variables.Fb(variables.mu * variables.v) * variables.Fa(variables.mu * variables.v)))
variables.Fb(variables.mu * variables.v) == diff(variables.Ab(variables.v), variables.mu) - diff(variables.Ab(variables.mu), variables.v) - variables.i*variables.g*variables.brackets(variables.Ab(variables.mu), variables.Ab(variables.v))
variables.g*variables.Ab(variables.mu) == variables.Abp(variables.mu)
variables.g*variables.Fb(variables.mu * variables.v) == variables.Fbp(variables.mu * variables.v)

#2.85.1
2 * variables.kappaab(2, 11) * variables.SBb(11) == integrate(variables.da(11) * variables.x * (-variables.G)**(1/2) * (variables.R - 1/2 * Abs(variables.Fb(4))**2) -1/6 * integrate(WedgeProduct(variables.Ab(3), variables.Fb(4), variables.Fb(4))))
#check 2.85.2, expanded states
variables.ds2 == variables.Gabf(11, variables.M * variables.N, variables.xas(variables.mu)) * diff(variables.xas(variables.M)) * diff(variables.xas(variables.N)) == variables.Gabf(10, variables.mu * variables.v, variables.xas(variables.mu)) * diff(variables.xas(variables.mu)) * diff(variables.xas(variables.v)) + math.exp(2 * variables.sigmaf(variables.xas(variables.mu))) * variables.brackets(diff(variables.xas(10)) + variables.Abf(variables.v, variables.xas(variables.mu) * diff(variables.xas(variables.v))))**2
variables.SBb(1) == 1/(2 * variables.kappaab(2, 10)) * integrate(variables.da(10) * variables.x * (-variables.G)**(1/2) * (variables.ea(variables.sigma) * variables.R - 1/2 * variables.ea(3 * variables.sigma) * Abs(variables.Fb(2))**2))
variables.SBb(2) == -1/(4 * variables.kappaab(2, 10)) * integrate(variables.da(10) * variables.x * (-variables.G)**(1/2) * (variables.ea(-variables.sigma) * Abs(variables.Fb(3))**2 + variables.ea(variables.sigma) * Abs(variables.Fbtilde(4))**2))
variables.SBb(3) == - 1/(4 * variables.kappaab(2, 10)) * integrate(WedgeProduct(variables.Ab(2), variables.Fb(4), variables.Fb(4))) == - 1/(4 * variables.kappaab(2, 10)) * integrate(WedgeProduct(variables.Ab(3), variables.Fb(3), variables.Fb(4)))
variables.Fbtilde(4) == diff(variables.Ab(3)) - WedgeProduct(variables.Ab(1), variables.Fb(3))
-WedgeProduct(diff(variables.lamdab(0), variables.Fb(3))) == - variables.d * (WedgeProduct(variables.lamdab(0), variables.Fb(3)))
variables.deltap * variables.Ab(3) == WedgeProduct(variables.lamdab(0), variables.Fb(3))
variables.d * variables.Fbtilde(4) == - WedgeProduct(variables.Fb(2), variables.Fb(3))
WedgeProduct(-diff(variables.lamdab(0)), variables.Fb(3)) == -diff(WedgeProduct(variables.lamdab(0), variables.Fb(0)))
variables.deltap * variables.Ab(3) == WedgeProduct(variables.lamdab(0), variables.Fb(3))
diff(variables.Fbtilde(4)) == WedgeProduct(-variables.Fb(2), variables.Fb(3))
variables.Gb(variables.mu * variables.v) == variables.ea(-variables.sigma) * variables.Gbf(variables.new, variables.mu * variables.v)
variables.sigma == 2 * variables.Phi/3
variables.SBb(variables.IIA) == variables.SBb(variables.NS) + variables.SBb(variables.R) + variables.SBb(variables.CS)
variables.SBb(variables.NS) == 1/(2 * variables.kappaab(2, 10)) * integrate(variables.daf(10, variables.x) * (-variables.G)**(1/2) * variables.ea(-2 * variables.Phi) * (variables.R + 4 * diff(variables.Phi, variables.mu) * diff(variables.Phi, 1, variables.mu) - 1/2 * Abs(variables.Hb(3))**2))
variables.SBb(variables.R) == -1/(4 * variables.kappaab(2, 10)) * integrate(variables.daf(10, variables.x) * (-variables.G)**(1/2) * (Abs(variables.Fb(2))**2 + Abs(variables.Fbtilde(4))**2))
variables.SBb(variables.CS) == -1/(4 * variables.kappaab(2, 10)) * integrate(WedgeProduct(variables.Bb(2), variables.Fb(4), variables.Fb(4)))
variables.Cb(1) == variables.ea(-variables.Phi) * variables.Cbp(1)
integrate(variables.daf(10, variables.x) * (-variables.G)**(1/2) * Abs(variables.Fb(2))**2) == integrate(variables.daf(10, variables.x) * (-variables.G)**(1/2) * variables.ea(-2 * variables.Phi) * Abs(variables.Fbp(2))**2)
variables.Fbp(2) == diff(variables.Cbp(1)) - WedgeProduct(diff(variables.Phi), variables.Fbp(2))
diff(variables.Fbp(2)) == WedgeProduct(diff(variables.Phi), variables.Fbp(2))
variables.delta * variables.Cbp(1) == diff(variables.lamdabp(0)) - variables.lamdabp(0) * diff(variables.Phi)
# check ommited operator states 2.88.14 - 16
variables.d * variables.eb(variables.p) == variables.d * variables.Star(variables.eb(variables.p)) == 0
variables.Tb(variables.F) == variables.i * (2/variables.SlopeRegge)**(1/2) * variables.psia(variables.mu) * diff(variables.xa(variables.mu), variables.z) - 2 * variables.i * (variables.SlopeRegge/2)**(1/2) * variables.Phib(0, variables.mu) * diff(variables.psia(variables.mu), variables.z)
(variables.d * variables.eb(variables.p)) * WedgeProduct(variables.d * variables.Phi , variables.eb(variables.p)) == variables.d * variables.Star(variables.eb(variables.p)) - WedgeProduct(variables.d * variables.Phi, variables.Star(variables.eb(variables.p))) == 0
1/(2 * math.pi * variables.SlopeRegge) * integrate(variables.Bb(2), variables.M)
variables.Fbtilde(6) == variables.Star(variables.Fbtilde(4))
variables.Fbtilde(8) == variables.Star(variables.Fb(2))
variables.d * variables.Star(variables.Fb(10)) == 0 
variables.SBbp(variables.IIA) == variables.SBbtilde(variables.IIA) - 1/(4 * variables.kappaab(2, 10)) * integrate(variables.daf(10, variables.x) * (-variables.G)**(1/2) * variables.Ma(2) + 1/(2 * variables.kappaab(2, 10)) * integrate(variables.M * variables.Fb(10)))
variables.Fb(2) == variables.Fb(2) + variables.M * variables.Bb(2) 
variables.Fb(4) == variables.Fb(4) + WedgeProduct(1/2 * variables.M * variables.Bb(2), variables.Bb(2))
variables.Fbtilde(4) == variables.Fbtilde(4) + WedgeProduct(1/2 * variables.M * variables.Bb(2), variables.Bb(2))
variables.SBb(variables.IIA) == variables.SBb(variables.NS) + variables.SBb(variables.R) + variables.SBb(variables.CS)
variables.SBb(variables.NS) == 1/(2 * variables.kappaab(2, 10)) * integrate(variables.daf(10, variables.x) * (-variables.G)**(1/2) * variables.ea(-2 * variables.Phi) * (variables.R + 4 * diff(variables.Phi, variables.mu) * diff(variables.Phi, 1, variables.mu) - 1/2 * Abs(variables.Hb(3))**2))
variables.SBb(variables.R) == -1/(4 * variables.kappaab(2, 10)) * integrate(variables.daf(10, variables.x) * (-variables.G)**(1/2) * (Abs(variables.Fb(1))**2 + Abs(variables.Fbtilde(3))**2 + 1/2 * Abs(variables.Fbtilde(5))**2))
variables.SBb(variables.CS) == -1/(4 * variables.kappaab(2, 10)) * integrate(WedgeProduct(variables.Cb(4), variables.Hb(3), variables.Fb(3)))
variables.Fbtilde(3) == variables.Fb(3) - WedgeProduct(variables.Cb(0), variables.Hb(3))
variables.Fbtilde(5) == variables.Fb(5) - WedgeProduct(1/2 * variables.Cb(2), variables.Hb(3)) + WedgeProduct(1/2 * variables.Bb(2), variables.Fb(3))
variables.d * variables.Star(variables.Fbtilde(5)) == variables.d * variables.Fbtilde(5) == WedgeProduct(variables.Hb(3), variables.Fb(3))
variables.Star(variables.Fbtilde(5)) == variables.Fbtilde(5)
variables.Gb(variables.E * variables.mu * variables.v) == variables.ea(-variables.Phi/2) * variables.Gb(variables.mu * variables.v) 
variables.tau == variables.Cb(0) + variables.i * variables.ea(-variables.Phi)
variables.MOb(variables.i * variables.j) == 1/(variables.Im * variables.tau) * np.array([[Abs(variables.tau)**2, - variables.Re(variables.tau)],
                                                                                        [-variables.Re(variables.tau, 1)]])
variables.Fab(variables.i, 3) == np.array([[variables.Hb(3)],
                                          [variables.Fb(3)]])
variables.SBb(variables.IIB) == 1/(2 * variables.kappaab(2, 10)) * integrate(variables.daf(10, variables.x) * (-variables.Gb(variables.E))**(1/2) * (variables.Rb(variables.E) - (diff(variables.taul, variables.mu) * diff(variables.tau, 1, variables.mu))/(2 * (variables.Im * variables.tau)**2) - np.dot((variables.MOb(variables.i * variables.j))/2 * variables.Fab(variables.i, 3), variables.Fab(variables.j, 3)) - 1/4 * Abs(variables.Fbtilde(5))**2) - (variables.epsilonb(variables.i * variables.j))/(8 * variables.kappaab(2, 10)) * integrate(WedgeProduct(variables.Cb(4), variables.Fab(variables.i, 3), variables.Fab(variables.j, 3))))
variables.taup == (variables.a * variables.tau + variables.b)/(variables.c * variables.tau + variables.d) 
variables.Fabp(variables.i, 3) == variables.Lamdaab(variables.i, variables.j) * variables.Fab(variables.j, 3)
variables.Lamdaab(variables.i, variables.j) == np.array([[variables.d, variables.c],
                                                        [variables.b, variables.a]])
variables.Fbptilde(5) == variables.Fbtilde(5)
variables.Gbp(variables.E * variables.mu * variables.v) == variables.Gb(variables.E * variables.mu * variables.v)
variables.MOp == (variables.Lamdaa(-1))**variables.T * variables.MO * variables.Lamdaa(-1)
variables.SBb(variables.I) == variables.SBb(variables.c) + variables.SBb(variables.o)
variables.SBb(variables.c) == 1/(2 * variables.kappaab(2, 10)) * integrate(variables.daf(10, variables.x) * (-variables.G)**(1/2) * variables.brackets(variables.ea(-2 * variables.Phi) * (variables.R + 4 * diff(variables.Phi, variables.mu) * diff(variables.Phi, 1, variables.mu)) - 1/2 * Abs(variables.Fbtilde(3))**2))
variables.SBb(variables.o) == -1/(2 * variables.gab(2, 10)) * integrate(variables.daf(10, variables.x) * (-variables.G)**(1/2) * variables.ea(-variables.Phi) * variables.Trbf(variables.v, Abs(variables.Fb(2))**2))
variables.Fbtilde(3) == variables.d * variables.Cb(2) - (variables.kappaab(2, 10))/(variables.gab(2, 10)) * variables.omegab(3)
variables.omegab(3) == variables.Trbf(variables.v, WedgeProduct(variables.Ab(1), variables.d * variables.Ab(1)) - WedgeProduct((2 * variables.i)/3 * variables.Ab(1), variables.Ab(1), variables.Ab(1)))
variables.delta * variables.omegab(3) == variables.d * variables.Trbf(variables.v, variables.lamda * variables.d * variables.Ab(1))
variables.delta * variables.Cb(2) == (variables.kappaab(2, 10))/(variables.gab(2, 10)) * variables.Trbf(variables.v, variables.lamda * variables.d * variables.Ab(1))
variables.SBb(variables.het) == 1/(2 * variables.kappaab(2, 10)) * integrate(variables.daf(10, variables.x) * (-variables.G)**(1/2) * variables.ea(-2 * variables.Phi) * variables.brackets(variables.R + 4 * diff(variables.Phi, variables.mu) * diff(variables.Phi, 1, variables.mu) - 1/2 * Abs(variables.Hbtilde(3))**2 - (variables.kappaab(2, 10))/(variables.gab(2, 10)) * variables.Trbf(variables.v, Abs(variables.Fb(2))**2)))
variables.Hbtilde(3) == variables.d * variables.Bb(2) - (variables.kappaab(2, 10))/(variables.gab(2, 10)) * variables.omegab(3)
variables.delta * variables.Bb(2) == (variables.kappaab(2, 10))/(variables.gab(2, 10)) * variables.Trbf(variables.v, variables.lamda * variables.d * variables.Ab(1))
variables.Gb(variables.I * variables.mu * variables.v) == variables.ea(-variables.Phib(variables.h)) * variables.Gb(variables.h * variables.mu * variables.v)
variables.Phib(variables.I) == - variables.Phib(variables.h)
variables.Fbtilde(variables.I * 3) == variables.Hbtilde(variables.h * 3)
variables.Ab(variables.I * 1) == variables.Ab(variables.h * 1)
variables.Sb(variables.int) == integrate(variables.d2z * (variables.jab(variables.a, variables.z) * variables.Aab(variables.a, variables.zl) + variables.jab(variables.a, variables.zl) * variables.Aab(variables.a, variables.z)))
variables.Xf(variables.brackets(variables.A)) == 1/2 * integrate(variables.d2(variables.zb(1)) * variables.d2(variables.zb(2)) *  variables.brackets((variables.hat(variables.kb(variables.L)))/(variables.zab(2, 12)) * variables.Aabf(variables.a, variables.zl, (variables.zb(1), variables.zlb(1))) * variables.Aabf(variables.a, variables.zl, (variables.zb(2), variables.zlb(2))) + (variables.hat(variables.kb(variables.R)))/(variables.zlab(2, 12)) * variables.Aabf(variables.a, variables.z, (variables.zb(1), variables.zlb(1))) * variables.Aabf(variables.a, variables.z, (variables.zb(2), variables.zlb(2)))))
variables.delta * variables.Zf(variables.brackets(variables.A)) == 2 * math.pi * integrate(variables.d2z * variables.lamdaaf(variables.a, (variables.z, variables.zl)) * variables.brackets(variables.hat(variables.kb(variables.L)) * diff(variables.Aabf(variables.a, variables.zl, (variables.z, variables.zl)), variables.z) + variables.hat(variables.kb(variables.R)) * diff(variables.Aabf(variables.a, variables.z, (variables.z, variables.zl)), variables.zl)))
variables.delta * variables.Zf(variables.brackets(variables.A)) == - 2 * math.pi * variables.hat(variables.k) * variables.delta * integrate(variables.d2z * variables.Aabf(variables.a, variables.z, (variables.z, variables.zl)) * variables.Aabf(variables.a, variables.zl, (variables.z, variables.zl)))
variables.Zfp(variables.brackets(variables.A)) == variables.Zf(variables.brackets(variables.A)) + 2 * math.pi * variables.hat(variables.k) * integrate(variables.d2z * variables.Aabf(variables.a, variables.z, (variables.z, variables.zl)) * variables.Aabf(variables.a, variables.zl, (variables.z, variables.zl))) == (variables.hat(variables.k))/2 * integrate(variables.d2(variables.zb(1)) *(variables.d2(variables.zb(2))) * ln(Abs(variables.zab(2, 12))) * variables.Fabf(variables.a, variables.z * variables.zl, (variables.zb(1), variables.zlb(1))) * variables.Fabf(variables.a, variables.z * variables.zl, (variables.zb(2), variables.zlb(2))))
summation(variables.qa(2), variables.L) - summation(variables.qa(2), variables.R) == 0
summation(1, variables.L) - summation(1, variables.R) == 0
summation(variables.q, variables.L) - summation(variables.q, variables.R) == 0
summation(variables.qa(3), variables.L) == 0
summation(variables.q, variables.L) == 0
variables.hat(variables.Ib(variables.d + 2)) == variables.d * variables.hat(variables.Ib(variables.d + 1))
variables.delta * variables.hat(variables.Ib(variables.d + 1)) == variables.d * variables.hat(variables.Ib(variables.d))
variables.delta * ln(variables.Z) == -variables.i/(2 * math.pi)**5 * integrate(variables.hat(variables.Ibf(variables.d, (variables.Fb(2), variables.Rb(2)))))
variables.hat(variables.Ibf(variables.B8, (variables.Fb(2), variables.Rb(2)))) == - (variables.Trf*variables.Fab(6, 2)())/1440 + (variables.Trf(variables.Fab(4, 2)) * variables.trf(variables.Rab(2, 2)))/2304 - (variables.Trf(variables.Fab(2, 2) * variables.trf(variables.Rab(4, 2))))/23040 - (variables.Trf(variables.Fab(2, 2)) * variables.brackets(variables.trf(variables.Rab(2, 2)))**2)/18432 + (variables.n * variables.trf(variables.Rab(6, 2)))/725760 + (variables.n * variables.trf(variables.Rab(4, 2)) * variables.trf(variables.Rab(2, 2)))/552960 + (variables.n * variables.brackets(variables.trf(variables.Rab(2, 2)))**3)/1327104
variables.hat(variables.Ibf(variables.B56, (variables.Fb(2), variables.Rb(2)))) == - 495 * (variables.trf(variables.Rab(6, 2)))/725760 + 225 * (variables.trf(variables.Rab(4, 2)) * variables.trf(variables.Rab(2, 2)))/552960 - 63 * (variables.brackets(variables.trf(variables.Rab(2, 2)))**3)/1327104
variables.hat(variables.Ibf(variables.SD, variables.Rb(2))) == 992 * (variables.trf(variables.Rab(6, 2)))/725760 - 448 * (variables.trf(variables.Rab(4, 2)) * variables.trf(variables.Rab(2, 2)))/552960 + 128 * (variables.brackets(variables.trf(variables.Rab(2, 2)))**3)/1327104
variables.hat(variables.Ibf(variables.IIB, variables.Rb(2))) == - 2 * variables.hat(variables.Ibf(variables.B8, variables.Rb(2))) + 2 * variables.hat(variables.Ibf(variables.B56, variables.Rb(2))) + variables.hat(variables.Ibf(variables.SD, variables.Rb(2))) == 0
variables.SBp == integrate(variables.Bb(2) * variables.Trf(variables.Fab(4, 2)))
variables.delta * variables.SBp == integrate(variables.Trf(variables.lamda * variables.d * variables.Ab(1)) * variables.Trf(variables.Fab(4, 2)))
variables.hat(variables.Ib(variables.d)) == variables.Trf(variables.lamda * variables.d * variables.Ab(1)) * variables.Trf(variables.Fab(4, 2))
variables.hat(variables.Ib(variables.d + 1)) == variables.Trf(variables.Ab(1) * variables.Fb(2)) * variables.Trf(variables.Fab(4, 2))
variables.hat(variables.Ib(variables.d + 2)) == variables.Trf(variables.Fab(2, 2)) * variables.Trf(variables.Fab(4, 2))
variables.SBpp == integrate(variables.Bb(2) * variables.brackets(variables.Trf(variables.Fab(2, 2)))**2)
variables.Trbf(variables.a, variables.ta(2)) == (variables.n - 2) * variables.Trbf(variables.v, variables.ta(2))
variables.Trbf(variables.a, variables.ta(4)) == (variables.n - 8) * variables.Trbf(variables.v, variables.ta(4)) + 3 * variables.Trbf(variables.v, variables.ta(2)) * variables.Trbf(variables.v, variables.ta(2))
variables.Trbf(variables.a, variables.ta(6)) == (variables.n - 32) * variables.Trbf(variables.v, variables.ta(6)) + 15 * variables.Trbf(variables.v, variables.ta(2)) * variables.Trbf(variables.v, variables.ta(4))
variables.Trbf(variables.a, variables.ta(4)) == 1/100 * variables.brackets(variables.Trbf(variables.a, variables.ta(2)))**2
variables.Trbf(variables.a, variables.ta(6)) == 1/7200 * variables.brackets(variables.Trbf(variables.a, variables.ta(2)))**3 
integrate(variables.Bb(2) * variables.Xbf(8, (variables.Fb(2), variables.Rb(2))))
variables.Hbtilde(3) == variables.d * variables.Bb(2) - variables.c * variables.omegabs(3 * variables.Y) - variables.cp * variables.omegabs(3 * variables.L)
variables.omegabs(3 * variables.L) == variables.omegabs(1) * variables.d * variables.omegabs(1) + 2/3 * variables.omegaabs(2, 1)
variables.delta * variables.omegabs(3 * variables.L) == variables.d * variables.trf(variables.Theta * variables.d * variables.omegabs(1))
variables.delta * variables.Ab(1) == variables.d * variables.lamda 
variables.delta * variables.omegabs(1) == variables.d * variables.Theta 
variables.delta * variables.Bb(2) == variables.c * variables.Trf(variables.lamda * variables.d * variables.Ab(1)) + variables.cp * variables.trf(variables.Theta * variables.d * variables.omegabs(1))
variables.brackets(variables.c * variables.Trf(variables.Fab(2, 2)) + variables.cp * variables.Trf(variables.Rab(2, 2))) * variables.Xbf(8, (variables.Fb(2), variables.Rb(2)))
variables.hat(variables.Ib(variables.I)) == variables.hat(variables.Ibf(variables.B56, variables.Rb(2))) + variables.hat(variables.Ibf(variables.B8, variables.Rb(2))) + variables.hat(variables.Tbf(variables.B8, (variables.Fb(2), variables.Rb(2)))) == 1/1440 * variables.bracketc(-variables.Trbf(variables.a, variables.Fab(6, 2)) + 1/48 * variables.Trbf(variables.a, variables.Fab(2, 2)) * variables.Trbf(variables.a, variables.Fab(4, 2)) - (variables.brackets(variables.Trbf(variables.Fab(2, 2)))**3)/14400) + (variables.n - 496) * variables.bracketc((variables.trf(variables.Rab(6, 2)))/725760 + (variables.trf(variables.Rab(4, 2)) * variables.trf(variables.Rab(2, 2)))/552960 + (variables.brackets(variables.trf(variables.Rab(2, 2)))**3)/1327104) + (variables.Yb(4) * variables.Xb(8))/768
variables.Yb(4) == variables.trf(variables.Rab(2, 2)) - 1/30 * variables.Trbf(variables.a, variables.Fab(2, 2)) 
variables.Xb(8) == variables.trf(variables.Rab(4, 2)) + (variables.brackets(variables.trf(variables.Rab(2, 2)))**2)/4 - (variables.Trbf(variables.a, variables.Fab(2, 2)) * variables.trf(variables.Rab(2, 2)))/30 + (variables.Trbf(variables.a, variables.Fab(4, 2)))/3 - (variables.brackets(variables.Trbf(variables.a, variables.Fab(2, 2))**2))/900
(variables.Yb(4) * variables.Xb(8))/768

#2.146.2
variables.AOb(variables.R - variables.R) == - variables.AOb(variables.NS - variables.NS)
2 * variables.i * variables.kappaa(2) * variables.tauab(2, variables.p) * variables.Gbf(9 - variables.p, variables.y)
variables.tauab(2, variables.p) == (math.pi)/(variables.kappaa(2)) * (4 * math.pi * variables.SlopeRegge)**(3-variables.p)
-1/(4 * variables.kappaab(2, 10)) * integrate(variables.da(10) * variables.x * (-variables.G)**(1/2) * Abs(variables.Fb(variables.p+2))**2 + variables.mub(variables.p) *integrate(variables.Cb(variables.p + 1)))
-2 * variables.kappaab(2, 10) * variables.i * variables.muaab(2, variables.p) * variables.Gbf(9-variables.p, variables.y)
variables.muab(2, variables.p) == (math.pi)/(variables.kappaab(2, 10)) * (4 * math.pi**2 * variables.SlopeRegge)**(3-variables.p) == variables.ea(2 * variables.Phib(0)) * variables.tauab(2, variables.p) == variables.Tab(2, variables.p)
integrate(variables.Fb(2), variables.Sb(2)) == variables.mub(variables.m)
math.exp(variables.i * variables.mub(variables.e) * integrate(variables.Ab(1), variables.P)) == math.exp(variables.i * variables.mub(variables.e) * integrate(variables.Fb(2), variables.D)) # contour integral
math.exp(variables.i * variables.mub(variables.e) * integrate(variables.Fb(2),variables.Sb(2))) == math.exp(variables.i * variables.mub(variables.e) * variables.mub(variables.m))
variables.mub(variables.e) * variables.mub(variables.m) == 2 * math.pi * variables.n 
integrate(variables.Fb(variables.p + 2), variables.Sb(variables.p + 2)) == variables.mub(6 - variables.p)/(2 * variables.kappaab(2, 10))
variables.mub(variables.p) * variables.mub(6 - variables.p) == math.pi * variables.n/variables.kappaab(2, 10)
variables.SBb(variables.Db(variables.p)) == - variables.mub(variables.p) * integrate(variables.da(variables.p + 1) * variables.zeta * variables.Trf(variables.bracketc(variables.ea(-variables.Phi) * variables.brackets(-sp.Determinant(variables.Gb(variables.a*variables.b) + variables.Bb(variables.a * variables.b) + 2 * math.pi * variables.SlopeRegge * variables.Fb(variables.a * variables.b)))**(1/2))))
variables.Of(variables.brackets(variables.xa(variables.m), variables.xa(variables.n))**2)
integrate(variables.Cb(2)) == integrate(diff(variables.xas(0) * (diff(variables.xas(1)) * variables.Cb(0, 1) + diff(variables.xas(2)) * variables.Cb(0, 2)))) == integrate(diff(variables.xas(0)) * diff(variables.xas(1)) * (variables.Cb(0, 1) + diff(variables.xa(2), 1) * variables.Cb(0, 2)))
integrate(diff(variables.xas(0)) * diff(variables.xas(1)) * diff(variables.xas(2)) * (variables.Cb(0, 1, 2) + 2 * math.pi * variables.SlopeRegge * variables.Fb(12) * variables.Cb(0)))
variables.i * variables.mub(variables.p) * integrate(variables.Trf(variables.brackets(WedgeProduct(math.exp(2 * math.pi * variables.SlopeRegge * variables.Fb(2) + variables.Bb(2)), summation(variables.Cb(variables.q), variables.q)))))
-variables.i * integrate(variables.da(variables.p + 1) * variables.zeta * variables.Trf(variables.lamdal * variables.Rhoa(variables.a) * variables.Db(variables.a) * variables.lamda))
(variables.taub(variables.F1))/(variables.taub(variables.D1)) == 1/(2 * math.pi * variables.SlopeRegge) * (variables.kappa)/(4 * math.pi**(5/2) * variables.SlopeRegge) == variables.kappa/(8 * math.pi**(7/2) * variables.SlopeRegge**2) # same D1, differetn meaning check
variables.g == (variables.taub(variables.F1))/(variables.taub(variables.D1))
variables.kappaa(2) == 1/2 * (2 * math.pi)**7 * variables.ga(2) * variables.SlopeRegge**4
variables.taub(variables.p) == 1/(variables.gf(2 * math.pi)**variables.p * variables.SlopeRegge**((variables.p + 1)/2)) == (2 * variables.kappaa(2))**(-1/2) * (2 * math.pi)**((7 - 2 * variables.p)/2) * variables.SlopeRegge**((3 - variables.p)/2)
variables.kappaab(2, 10) == 1/2 * (2 * math.pi)**7 * variables.SlopeRegge**4
variables.gab(2, variables.D * variables.p) == 1/((2 * math.pi * variables.SlopeRegge)**2 * variables.taub(variables.p)) == (2 * math.pi)**(variables.p-2) * variables.g  * variables.SlopeRegge**((variables.p - 3)/ 2)
1/(4 * variables.gab(2, variables.D * variables.p)) * variables.Trb(variables.f)
variables.ttilde == np.array([[variables.t, 0],
                             [0, - variables.ta(variables.T)]])
1/(4 * variables.gab(2, variables.D * variables.p)) * variables.Trbf(variables.f, variables.ta(2)) == 1/(4 * variables.gab(2, variables.D * variables.p, variables.S * variables.Of(2 * variables.n))) * variables.Trbf(variables.v, variables.tatilde(2))
# type I
variables.kappaab(2, 10 -variables.k) == (2 * math.pi * variables.R)**(-variables.k) * variables.kappaa(2)
variables.gab(2, 10 - variables.k, variables.Y * variables.M) == (2 * math.pi * variables.R)**(-variables.k) * variables.gab(2,variables.YM)
#(next line)
variables.kappaab(2, 10 - variables.k) == 2 * (2 * math.pi * variables.Rp)**(-variables.k) * variables.kappaa(2)
variables.gab(2, 10 - variables.k, variables.YM) == variables.gab(2, (variables.D * (9 - variables.k), variables.S * variables.Of(32)))
(variables.gab(2, variables.YM))/(variables.kappa) == 2 * (2 * math.pi)**(7/2) * variables.SlopeRegge # type I
variables.SB == -1/((2 * math.pi * variables.SlopeRegge)**2 * variables.gab(2, variables.YM)) * integrate(variables.da(10) * variables.x * variables.Trf(variables.bracketc(variables.brackets(-sp.Determinant(variables.etab(variables.mu * variables.v) + 2 * math.pi * variables.SlopeRegge * variables.Fb(variables.mu * variables.v)))**(1/2))))
((2 * math.pi * variables.SlopeRegge)**2)/(32 * variables.gab(2, variables.YM)) * variables.Trbf(variables.v, (4 * variables.Fb(variables.mu * variables.v) * variables.Fa(variables.v * variables.sigma) * variables.Fb(variables.sigma * variables.rho) * variables.Fa(variables.rho * variables.mu) - variables.Fb(variables.mu * variables.v) * variables.Fa(variables.v * variables.mu) * variables.Fb(variables.sigma * variables.rho) * variables.Fa(variables.rho * variables.sigma)))
# check unexpanded version 2.152.34

# interactions 2.152.1
variables.Qb(variables.alpha) + variables.holderb(variables.Betaa(variables.perp) * variables.Qtilde, variables.alpha)
variables.Betaa(variables.perp) == np.prod(variables.Betaa(variables.m), (variables.m, variables.Sb(variables.D)))
variables.Qb(variables.alpha) + variables.holderb(variables.Betaap(variables.perp) * variables.Qtilde, variables.alpha) == variables.Qb(variables.alpha) + variables.bracketsb(variables.Betaa(variables.perp) * (variables.Betaa(variables.perp - 1) * variables.Betaap(variables.perp)) * variables.Qtilde, variables.alpha) 
variables.Betaap(variables.perp) == np.prod(variables.Betaa(variables.m), (variables.m, variables.Sbp(variables.D)))
# check unexpanded state for Beta 2.153.3
variables.xa(variables.mu, (variables.w, variables.wl)) == variables.xa(variables.mu, variables.w) + variables.Xaftilde(variables.mu, variables.wl)
variables.xa(variables.mu, variables.w) == variables.xas(variables.mu) + (variables.SlopeRegge/2)**(1/2) * variables.brackets(-variables.alphaab(variables.mu, 0) * variables.w + variables.i * summation((variables.alphaab(variables.mu, variables.m))/variables.m * math.exp(variables.i * variables.m * variables.w), (variables.m, variables.Z, variables.m != 0)))
variables.xa(variables.mu, variables.w) == variables.i * (variables.SlopeRegge/2)**(1/2) * summation((variables.alphaab(variables.mu, variables.r))/variables.r * math.exp(variables.i * variables.r * variables.w), (variables.r, variables.ZB + 1/2))
variables.Xaftilde(variables.mu, variables.wl) == variables.Xaf(variables.mu, 2 * math.pi - variables.wl)
variables.Xaftilde(variables.mu, variables.wl) == - variables.Xaf(variables.mu, 2 * math.pi - variables.wl)
(8 - variables.lbND) * (- 1/24 - 1/48) + variables.lbND * (1/48 + 1/24) == -1/2 + (variables.lbND)/(8)
variables.Qb(variables.alpha) + variables.holderb(variables.rhoa(-1) * variables.Betaa(variables.perp) * variables.rho * variables.Qtilde, variables.alpha)
(variables.Betaa(variables.perp))**(-1) * variables.rhoa(-1) * variables.Betaa(variables.perp) * variables.rho == (variables.Betaa(variables.perp))**(-1) * variables.Betaa(variables.perp) * variables.rhoa(2) == variables.rhoa(2)
math.exp(2 * variables.i * summation(variables.sb(variables.a) * variables.phib(variables.a), (variables.a, 1, 4)))
np.diag(variables.brackets(math.exp(variables.i * variables.phib(1)), math.exp(variables.i * variables.phib(2)), math.exp(variables.i * variables.phib(3)), math.exp(variables.i * variables.phib(4))))
variables.sigmaa(1) == 0 #:
diff(variables.Ref(variables.Za(variables.a)), 1) == variables.Imf(variables.Za(variables.a)) == 0
variables.sigmaa(1) == math.pi #:
diff(variables.Ref(variables.brackets(math.exp(- variables.i * variables.phib(variables.a)) * variables.Za(variables.a))), 1) == variables.Imf(variables.brackets(math.exp(-variables.i * variables.phib(variables.a)) * variables.Za(variables.a))) == 0
variables.Zaf(variables.a, (variables.w, variables.wl)) == variables.LOaf(variables.a, variables.w) + variables.LOlaf(variables.a, -variables.wl) == math.exp(-2 * variables.i * variables.phib(variables.a)) * variables.LOaf(variables.a, variables.w + 2 * math.pi) + variables.LOlaf(variables.a, -variables.wl)
variables.LOaf(variables.a, variables.w) == variables.i * (variables.SlopeRegge/2)**(1/2) * summation((variables.alphaab(variables.a, variables.r))/variables.r * math.exp(variables.i * variables.r * variables.w), (variables.r, variables.ZB + variables.vb(variables.a)))
variables.qa(variables.Eb(0)) * np.prod(variables.brackets(1 - variables.qa(variables.m + (variables.phi/math.pi)))**(-1) * variables.brackets(1 - variables.qa(variables.m + 1 - (variables.phi/math.pi)))**(-1), (variables.m, 0, math.inf)) == -variables.i * (math.exp(variables.phia(2) * variables.t /math.pi) * variables.etaf(variables.i * variables.t))/(variables.thetavarb((1, 1), (variables.i * variables.phi * variables.t/math.pi, variables.i * variables.t)))
variables.Eb(0) == 1/24 - 1/2* (variables.phi/math.pi - 1/2)**2 
variables.thetavarb((1, 1), (variables.v, variables.i * variables.t)) == -2 * variables.qa(1/8) * math.sin(math.pi * variables.v) * np.prod((1- variables.qa(variables.m)) * (1 - variables.z * variables.qa(variables.m)) * (1 - variables.za(-1) * variables.qa(variables.m)), (variables.m, 1, math.inf))
variables.thetavarb((1, 1), (-variables.i * variables.v/variables.t, variables.i/variables.t)) == -variables.i * variables.ta(1/2) * math.exp(math.pi * variables.va(2)/variables.t) * variables.thetavarb((1, 1), (variables.v, variables.i * variables.t))
variables.Zabf(variables.alpha, variables.Beta, (variables.phi, variables.i * variables.t)) == (variables.thetavarb((variables.alpha, variables.Beta), (variables.i * variables.phi * variables.t/math.pi, variables.i * variables.t)))/(math.exp(variables.phia(2) * variables.t/math.pi) * variables.etaf(variables.i * variables.t))
1/2 * variables.brackets(np.prod(variables.Zabf(0, 0, (variables.phib(variables.a), variables.i * variables.t)), (variables.a, 1, 4)) - np.prod(variables.Zabf(0, 1, (variables.phib(variables.a), variables.i * variables.t)), (variables.a, 1, 4)) - np.prod(variables.Zabf(1, 0, (variables.phib(variables.a), variables.i * variables.t)), (variables.a, 1, 4)) - np.prod(variables.Zabf(1, 1, (variables.phib(variables.a), variables.i * variables.t)), (variables.a, 1, 4)))
np.prod(variables.Zabf(1, 1, (variables.phipb(variables.a), variables.i * variables.t)), variables.a, 1, 4)
variables.phipb(1) == 1/2 * (variables.phib(1) + variables.phib(2) + variables.phib(3) + variables.phib(4)) 
variables.phipb(2) == 1/2 * (variables.phib(1) + variables.phib(2) - variables.phib(3) - variables.phib(4))
variables.phipb(3) == 1/2 * (variables.phib(1) - variables.phib(2) + variables.phib(3) - variables.phib(4))
variables.phipb(4) == 1/2 * (variables.phib(1) - variables.phib(2) - variables.phib(3) + variables.phib(4))
variables.V == -quad((diff(variables.t))/variables.t * (8 * math.pi**2 * variables.SlopeRegge * variables.t)**(-1/2) * math.exp(-(variables.t * variables.yab(2, 1))/(2 * math.pi * variables.SlopeRegge)) * np.prod((variables.thetavarb((1, 1,), (variables.i * variables.phipb(variables.a) * variables.t/math.pi, variables.i * variables.t)))/(variables.thetavarb((1, 1), (variables.i * variables.phib(variables.a) * variables.t/math.pi, variables.i * variables.t))), (variables.a, 1, 4)), (0, math.inf))
variables.thetavarb((1, 1), (variables.i * variables.phib(variables.a) * variables.t/math.pi, variables.i * variables.t))**(-1) == variables.i * variables.L * variables.etaf(variables.i * variables.t)**(-3) * (8 * math.pi**2 * variables.SlopeRegge * variables.t)**(-1/2)
variables.i * variables.etaf(variables.i * variables.t)**(-3) * math.exp(variables.brackets(-((variables.t * (variables.yab(2, 8) + variables.yab(2, 9)))/(2 * math.pi * variables.SlopeRegge))))
np.prod((variables.thetavarb((1, 1), (variables.i * variables.phipb(variables.a) * variables.t/math.pi, variables.i * variables.t)))/(variables.thetavarb((1, 1), (variables.i * variables.phib(variables.a) * variables.t/math.pi, variables.i * variables.t))), (variables.a, 1, 4)) == np.prod((math.sin(variables.phipb(variables.a)))/(math.sin(variables.phib(variables.a))), (variables.a, 1, 4))
variables.ma(2) == (variables.yab(2, 1))/(4 * math.pi**2 * variables.SlopeRegge**2) - (variables.phib(1))/(2 * math.pi * variables.SlopeRegge)
0 <= variables.phib(1) <= math.pi 

# interactions 2.159.1
variables.xa(3) == variables.Xpa(0) * np.tanh(variables.u)
variables.AO == -variables.i * variables.Vb(variables.p) * quad((diff(variables.t))/variables.t * (8 * math.pi**2 * variables.SlopeRegge * variables.t)**(-variables.p/2) * math.exp(-(variables.t * variables.ya(2))/(2 * math.pi * variables.SlopeRegge)) * (variables.thetavarb((1, 1), (variables.u * variables.t/(2 * math.pi), variables.i * variables.t))**4)/(variables.etaf(variables.i * variables.t)**9 * variables.thetavarb((1, 1), (variables.u * variables.t/math.pi, variables.i * variables.t))), (0, math.inf))
variables.AO == (variables.Vb(variables.p))/(8 * math.pi**2 * variables.SlopeRegge)**(variables.p/2) * quad((diff(variables.t))/variables.t * variables.ta((6 - variables.p)/2) * math.exp(-(variables.t * variables.ya(2))/(2 * math.pi * variables.SlopeRegge)) * (variables.thetavarb((1, 1), (variables.i * variables.u/(2 * math.pi), variables.i/variables.t))**4)/(variables.etaf(variables.i/variables.t)**9 * variables.thetavarb((1, 1), (variables.i * variables.u/math.pi, variables.i/variables.t))), (0, math.inf))
variables.AO == - variables.i * quad(diff(variables.tau) * variables.Vf(variables.rf(variables.tau), variables.v), -math.inf, math.inf)
variables.rf(variables.tau)**2 == variables.ya(2) + variables.va(2) * variables.taua(2)
variables.v == np.tanh(variables.u)
variables.Vf(variables.r, variables.v) == variables.i * (2 * variables.Vb(variables.p))/((8 * math.pi**2 * variables.SlopeRegge)**((variables.p + 1)/2)) * quad(diff(variables.t) * np.cross(variables.ta((5 - variables.p)/2), math.exp(-(variables.t * variables.ra(2))/(2 * math.pi * variables.SlopeRegge)) * ((np.tanh(variables.u)) * variables.thetavarb((1, 1), (variables.i * variables.u/(2 * math.pi), variables.i / variables.t))**4)/(variables.etaf(variables.i/variables.t)**9 * variables.thetavarb((1, 1), (variables.i * variables.u/math.pi, variables.i/variables.t)))), 0, math.inf)
variables.Vf(variables.r, variables.v) == -variables.va(4) * (variables.Vb(variables.p))/((8 * math.pi**2 * variables.SlopeRegge)**((variables.p + 1)/2)) * quad(diff(variables.t) * variables.ta((5 - variables.p)/2) * math.exp(- (variables.t * variables.ra(2))/(2 * math.pi * variables.SlopeRegge)) + variables.Of(variables.va(6)), 0, math.inf) == - (variables.va(4))/(variables.ra(7 - variables.p)) * (variables.Vb(variables.p))/(variables.SlopeRegge**(variables.p - 3)) * 2 **(2 - 2 * variables.p) * math.pi**((5 - 3 * variables.p)/2) * variables.Rhof((7 - variables.p)/2) + variables.Of(variables.va(6))
# check approximations 2.160.8 - 12
variables.L == variables.Trf(variables.bracketc(1/(2 * variables.g * variables.SlopeRegge**(1/2)) * variables.Db(0) * variables.xa(variables.i) * variables.Db(0) * variables.xa(variables.i) + 1/(4 * variables.g * variables.SlopeRegge**(1/2) * (2 * math.pi * variables.SlopeRegge)**2) * variables.brackets(variables.xa(variables.i), variables.xa(variables.j))**2 - variables.i/2 * variables.lamda * variables.Db(0) * variables.lamda + 1/(4 * math.pi * variables.SlopeRegge) * variables.lamda * variables.Rhoa(0) * variables.Rhoa(variables.i) * variables.brackets(variables.xa(variables.i), variables.lamda)))
variables.H == variables.Trf(variables.bracketc((variables.g * variables.SlopeRegge**(1/2))/(2) * variables.pb(variables.i) * variables.pb(variables.i) - 1/(16 * math.pi**2 * variables.g * variables.SlopeRegge**(5/2)) * variables.brackets(variables.xa(variables.i), variables.lamda)**2 - 1/(4 * math.pi * variables.SlopeRegge) * variables.lamda * variables.Rhoa(0) * variables.Rhoa(variables.i) * variables.brackets(variables.xa(variables.i), variables.lamda)))
variables.brackets(variables.pb(variables.i * variables.a * variables.b), variables.xab(variables.j, variables.c * variables.d)) == - variables.i * variables.deltaab(variables.j, variables.l) * variables.deltab(variables.a * variables.d) * variables.deltab(variables.b * variables.c)
variables.xa(variables.i) =variables.ga(1/3) * variables.SlopeRegge**(1/2) * variables.Ya(variables.i) 
variables.H == (variables.ga(1/3))/(variables.SlopeRegge**(1/2)) & variables.Trf(variables.bracketc(1/2 * variables.pb(variables.Y, variables.i) * variables.pb(variables.Y, variables.i) - 1/(16 * math.pi**2) * variables.brackets(variables.Ya(variables.i), variables.Ya(variables.j))**2 - 1/(4 * math.pi) * variables.lamda * variables.Rhoa(0) * variables.Rhoaf(variables.i, variables.brackets(variables.Ya(variables.i), variables.lamda))))
variables.diracb((variables.sb(3), variables.sb(4)), variables.NS)
-math.exp(variables.brackets(math.pi * variables.i * (variables.sb(3) + variables.sb(4)))) == + 1 #goes to
variables.sb(3) == variables.sb(4)
variables.diracb((variables.sb(1), variables.sb(2)), variables.R) 
variables.SB == -1/(4 * variables.gab(2, variables.D9)) * integrate(variables.da(10) * variables.x * variables.Fb(variables.M * variables.n) * variables.Fa(variables.M * variables.N) - 1/(4 * variables.gab(5, variables.D5)) * integrate(variables.da(6) * variables.x * variables.Fpb(variables.M * variables.N) * variables.Fpa(variables.M * variables.N))) - integrate(variables.da(6) * variables.x * variables.brackets(variables.Db(variables.mu) * variables.dagger(variables.chi) * variables.da(variables.mu) * variables.chi + (variables.gab(2, variables.D5))/(2) * summation((variables.dagger(variables.chib(variables.i)) * variables.sigmaab(variables.A, variables.i * variables.j) * variables.chib(variables.j))**2, (variables.A, 1, 3))))
variables.Apb(variables.M) == variables.Apb(variables.mu)
variables.xbp(variables.m)/(2 * math.pi * variables.SlopeRegge)
variables.Ab(variables.i) == (variables.xb(variables.i))/(2 * math.pi * variables.SlopeRegge) 
variables.Apb(variables.i) == (variables.Xpb(variables.i))/(2 * math.pi * variables.SlopeRegge)
((variables.Xb(variables.i) - variables.Xpb(variables.i))/(2 * math.pi * variables.SlopeRegge))**2 * variables.dagger(variables.chi) * variables.chi

# low energy actions 2.360.1
variables.fb(variables.a * variables.b) == (2 * variables.hat(variables.kb(variables.a)) * variables.deltab(variables.a * variables.b))/(variables.gab(2, 4)) * variables.S 
integrate(WedgeProduct(variables.Fb(2), variables.Fb(2)))
variables.S == variables.S * variables.i * variables.epsilon 
variables.S == variables.t * variables.S 
variables.Gb(variables.mu * variables.v * 4 * variables.E) == variables.t * variables.Gb(variables.mu * variables.v * 4 * variables.E)
variables.SB == variables.t * variables.SB 
variables.fb(variables.a * variables.b) == variables.deltab(variables.a * variables.b) * variables.S/(8 * math.pi**2)
(variables.gab(2, variables.YM))/(8 * math.pi**2) == 1/(variables.Ref(variables.bracket(variables.S)))
integrate(variables.da(4) * variables.x * (-variables.Gb(4 * variables.E))**(1/2) * math.exp(variables.kappaa(2) * variables.K) * (variables.Ka(variables.il * variables.j) * variables.star(variables.Wb(variables.cov(variables.i))) * variables.Wb(variables.cov(variables.j)) - 3 * variables.kappaa(2) * variables.star(variables.W) * variables.W))
variables.SBb(variables.L) == variables.ta(1 - variables.L) * variables.SBb(variables.L)
variables.gOb(variables.i * variables.j) == variables.bracket(variables.diracbs((diff(variables.LOb(variables.ws), variables.z))/(diff(variables.Phia(variables.i), variables.z)), (diff(variables.LOb(variables.ws), variables.z))/(diff(variables.Phia(variables.j), variables.z))))
1/2 * variables.gOb(variables.i * variables.j) * diff(variables.Phia(variables.i), variables.mu) * diff(variables.Phia(variables.j), variables.z, variables.mu)

# 2.382.3
variables.K == - 3 * ln(variables.T + variables.star(variables.T))
variables.K == - ln(variables.Imf(summation(variables.star(variables.xa(variables.I)) * diff(variables.Ff(variables.X), variables.I), variables.I)))
variables.Ff(variables.X) == ((variables.xa(1))**3)/(variables.xa(0))
variables.T == (variables.i * variables.xa(1))/(variables.xa(0))
variables.delta * variables.xa(1) == variables.epsilon * variables.xa(0)
variables.delta * variables.F == 3 * variables.epsilon * (variables.xa(1))**2
variables.delta * variables.F == variables.cb(variables.I * variables.J) * variables.xa(variables.I) * variables.xa(variables.J) 
variables.Delta * variables.F == variables.i * variables.lamda * (variables.xa(0))**2
variables.delta * variables.xa(variables.I) == variables.omegaas(variables.I * variables.J) * variables.xa(variables.J)
variables.Ta(variables.A) == (variables.i * variables.xa(variables.A))/(variables.xa(0))
variables.F == (variables.db(variables.A * variables.B * variables.C) * variables.xa(variables.A) * variables.xa(variables.B) * variables.xa(variables.C))/(variables.xa(0)) + variables.i * variables.lamda * (variables.xa(0))**2 
variables.dirac(variables.c, variables.ctilde)
variables.dirac(variables.c, variables.atilde)
variables.dirac(variables.a, variables.ctilde)
variables.dirac(variables.a, variables.atilde)
variables.PhiaB(variables.pos, variables.pos)
variables.PhiaB(variables.pos, variables.neg)
variables.PhiaB(variables.neg, variables.pos)
variables.PhiaB(variables.neg, variables.neg)
variables.h >= 1/2 * (variables.Qb(variables.PhiB) + variables.Qb(variables.PsiB)) == variables.hb(variables.PhiB) + variables.hb(variables.PsiB)
variables.j == variables.psia(variables.i) * variables.psia(variables.il)
variables.jtilde == variables.psiatilde(variables.i) * variables.psiatilde(variables.il)
variables.Q == variables.p 
variables.h == variables.p/2 
variables.Qtilde == - variables.q 
variables.htilde == variables.q/2 
(variables.Gab(variables.pos, 0))**2 == 0 
variables.Tab(variables.top, variables.B) == variables.Tb(variables.B) + 1/2 * diff(variables.j, variables.z)
