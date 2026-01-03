import math
from scipy.integrate import quad
import numpy as np 
from sympy import symbols, diff,integrate, Function, Derivative, summation, oo, Abs
import sympy as sp
import variables
import constants

InvPoincare = symbols('invpoincare')
InvPoincarep = symbols('invpoincare')

def poincare_invariance_conditions():
    InvPoincare == -variables.m*integrate(np.diff(variables.tau)-variables.xadiff(variables.mu)*variables.xbdiff(variables.mu), variables.t)**(1/2) #for invariance under Poincare transformations
    InvPoincare*variables.NambuGotoVar == -variables.m*integrate(np.diff(variables.tau)*variables.ubdiff(variables.mu)*variables.xa(variables.mu)*variables.NambuGotoVar)
    variables.ua(variables.mu) == variables.xadiff(variables.mu)*(integrate(-variables.xa(variables.v)*variables.xb(variables.v)))**(-1/2)
    InvPoincarep  == 1/2*integrate(variables.n**(-1)*diff(variables.xa(variables.mu)*variables.xa(variables.mu))-variables.n*variables.m**2) #using n as a tetrad
    
variables.A == (2-variables.D)/24 # for Weyl invariance condition

diff_resulta = variables.sigma
for _ in range(int(variables.d)):
    diff_resulta = diff(diff_resulta, variables.d)
a = diff_resulta

variables.PIC(diff(variables.phip))* math.exp(-variables.S*variables.PIC(variables.phip)) == variables.PIC(diff(variables.phi)) * math.exp(-variables.S*variables.PIC(variables.phi)) # 41
0 == integrate(variables.PIC(diff(variables.phip)) * math.exp(-variables.S*variables.PIC(variables.phip))) - integrate(variables.PIC(diff(variables.phi))*math.exp(-variables.S*variables.PIC(variables.phi))) == variables.epsilon/(2*math.pi*variables.i)* integrate(a*variables.g**(1/2)*variables.rho(variables.sigma)*variables.bracket(variables.Deltab(variables.a)*variables.ja(variables.a, variables.sigma)))
variables.Deltab(variables.a)*variables.ja(variables.a) == 0

#159.1
variables.da(variables.mu) * variables.tp == Abs(sp.Determinant(diff(variables.tp)/(diff(variables.t)))) * variables.da(variables.mu) * variables.t
np.prod(1/(4*math.pi) * (variables.b, diff(diff(variables.hat(variables.g), variables.k), variables.k)), (variables.k, 1, variables.mu)) == sp.Determinant(diff(variables.t)/(diff(variables.tp))) * np.prod(1/(4 * math.pi) * (variables.b, diff(variables.hat (variables.g), variables.k)), (variables.k, 1, variables.mu))
variables.delta * variables.hat(variables.gbfp(variables.a*variables.b, variables.t, variables.sigma)) == 2 * variables.delta * variables.omegas *  ( variables.t, variables.sigma) * variables.hat(variables.gbf(variables.a * variables.b, variables.t, variables.sigma))
(variables.b, diff(variables.hat(variables.gp), variables.k)) == integrate(variables.d2sigma * (variables.hat(variables.g))**(1/2) * variables.bb(variables.a * variables.b * variables.hat(variables.gap(variables.a * variables.c)) * variables.hat(variables.gap(variables.b * variables.d)) * diff(variables.hat(variables.gbp(variables.c * variables.d)), variables.k))) == integrate(variables.d2sigma * (variables.hat(variables.g))**(1/2) * variables.bb(variables.a *variables.b) * (variables.hat(variables.ga(variables.a *variables.c)) * variables.hat(variables.ga(variables.b * variables.d)) * diff(variables.hat(variables.gb(variables.c *variables.d)), variables.k) + 2 * variables.hat(variables.ga(variables.a *variables.b)) * diff(variables.omegas, variables.k))) == (variables.b, diff(variables.hat(variables.g), variables.k))
variables.delta * (variables.b, diff(variables.hat(variables.g), variables.k)) == -2 * (variables.b, variables.P1 * diff(variables.zeta, variables.k)) == -2 * ((variables.P1a(variables.T) * variables.b), diff(variables.zeta, variables.k)) == 0
variables.deltab(variables.B) * variables.VOb(variables.m) == variables.i * variables.epsilon * diff(variables.ca(variables.a) * variables.VOb(variables.m), variables.a)
variables.deltab(variables.B) * (variables.b, diff(variables.hat(variables.g), variables.k)) == variables.i * variables.epsilon * (variables.T, diff(variables.hat(variables.g), variables.k))
variables.muab(variables.k * variables.a, variables.b) == 1/2 * variables.hat(variables.ga(variables.b * variables.c)) * diff(variables.hat(variables.gb(variables.a *variables.c)), variables.k)
1/(2 * math.pi) * (variables.b, variables.mub(variables.k)) == 1/(2 * math.pi) * integrate(variables.d2z * (variables.bb(variables.z * variables.z) * variables.muba(variables.k * variables.zl, variables.z) + variables.bb(variables.zl * variables.zl) * variables.muba(variables.k*variables.z, variables.zl)))
variables.zbp(variables.m) == variables.zb(variables.m) + variables.delta * variables.ta(variables.k) * variables.vab(variables.z * variables.m, variables.k * variables.m) * (variables.zb(variables.m), variables.zlb(variables.m))
# check approx at 162.10
variables.muba(variables.k * variables.zb(variables.m), variables.zlb(variables.m)) == diff(variables.vab(variables.zlb(variables.m), variables.k *variables.m), variables.zb(variables.m))
variables.muba(variables.k * variables.zlb(variables.m), variables.zb(variables.m)) == diff(variables.vab(variables.zb(variables.m), variables.k *variables.m), variables.zlb(variables.m))
1/(2 * math.pi) * (variables.b, variables.mub(variables.k)) == 1/(2 * math.pi * variables.i) * summation(integrate(variables.d * variables.zb(variables.m) * variables.vab(variables.zb(variables.m) , variables.k *variables.m) * variables.bb(variables.zb(variables.m) * variables.zb(variables.m)) - variables.d * variables.zlb(variables.m) * variables.vab(variables.zlb(variables.m), variables.k *variables.m) * variables.bb(variables.zlb(variables.m) * variables.zlb(variables.m)), variables.Cb(variables.m)), variables.m) # contour integral
variables.d * variables.zb(variables.m)/(variables.d * variables.ta(variables.k)) == variables.vab(variables.zb(variables.m), variables.k *variables.m)
# check conditions for transofrmation as provided 162.14 to .15 with evaluation bounds\
variables.z == variables.zp + variables.zb(variables.v)
integrate((variables.d * variables.zp)/(2 *math.pi * variables.i) * variables.bb(variables.zp * variables.zp) * integrate((variables.d * variables.zlbp(variables.m))/(-2*math.pi * variables.i) * variables.bb(variables.zpl * variables.zpl), variables.C), variables.C) == variables.bb(-1) * variables.bbtilde(-1)
variables.m == variables.mu  + 2 * variables.nb(variables.c) + variables.nb(variables.o) - variables.keta == -3 * variables.chi + 2 * variables.nb(variables.c) + variables.nb(variables.o)
variables.bb(-1) * variables.bbtilde(-1) * variables.ctilde * variables.c * variables.VOb(variables.m) == variables.VOb(variables.m)

#198.1, mobius invariance conditions
variables.zp == (variables.alpha * variables.z + variables.Beta)/(variables.gamma * variables.z + variables.delta)
variables.bracketb(variables.AOb(variables.i, (0, 0)), variables.Sb(2)) == variables.bracketb(variables.AObp(variables.i, (0, 0)), variables.Sb(2)) == variables.gammaa(-variables.hb(variables.i)) * variables.gammaal(-variables.hbtilde(variables.i)) * variables.bracketb(variables.AOb(variables.i, (0, 0)), variables.Sb(2))
variables.bracketb(variables.AOb(variables.i, (variables.zb(1), variables.zlb(1))) * variables.AOb(variables.j, (variables.zb(2), variables.zlb(2))), variables.Sb(2)) == variables.zab(-variables.hb(variables.i) - variables.hb(variables.j), 12) * variables.zlab(-variables.hbtilde(variables.i)-variables.hbtilde(variables.j), 12) * variables.bracketb(variables.AOb(variables.i, (1, 1)) * variables.AOb(variables.j, (0, 0)), variables.Sb(2))
if variables.hb(variables.p) != variables.hb(variables.q):
    if variables.hbtilde(variables.p) != variables.hbtilde(variables.q):
        variables.bracketb(variables.OOb(variables.p, (variables.zb(1), variables.zlb(1))) * variables.OOb(variables.zb(2), variables.zlb(2)), variables.Sb(2)) == 0
variables.bracketb(np.prod(variables.OOb(variables.pb(variables.i), (variables.zb(variables.i), variables.zlb(variables.i))), (variables.i, 1, 3)), variables.Sb(2)) =variables.Cb(variables.pb(1) * variables.pb(2) * variables.pb(3)) * np.prod(variables.zab(variables.h- 2 * (variables.hb(variables.i) + variables.hb(variables.j)), variables.i * variables.j) * variables.zlab(variables.htilde - 2 * (variables.hbtilde(variables.i) + variables.hbtilde(variables.j)), variables.i * variables.j), ((variables.i, variables.j), 1, 3, variables.i<variables.j))
variables.bracketb(np.prod(variables.OOb(variables.pb(variables.i), (variables.zb(variables.i), variables.zlb(variables.i))), (variables.i, 1, 4)), variables.Sb(2)) == np.cross(variables.Cb(variables.pb(1) * variables.pb(2) * variables.pb(3) * variables.pb(4)) * (variables.zb(variables.c), variables.zlb(variables.c)) * variables.holdera((variables.zb(12), variables.zb(34)), variables.h) * variables.holdera((variables.zlb(12), variables.zlb(34)), variables.htilde), np.prod(variables.zab(-variables.hb(variables.i) - variables.hb(variables.j), variables.i * variables.j) * variables.zlab(-variables.hbtilde(variables.i), variables.i * variables.j), ((variables.i, variables.j), 1, 4, variables.i<variables.j)))

# 2.31.1
variables.ZBb(variables.Tb(2)) == variables.Vb(10) * integrate((variables.d2(variables.tau))/(4 * variables.taub(2)) * integrate((variables.da(10) * variables.k)/(2 * math.pi)**(10) * summation((-1)**(variables.FBb(variables.i)) * variables.qa(variables.SlopeRegge * (variables.ka(2) + variables.mab(2, variables.i))/4) * variables.qla(variables.SlopeRegge * (variables.ka(2) + variables.mabtilde(2, variables.i))/4), (variables.i, variables.hilberta(variables.perp)))), variables.F)
variables.ma(2) == 4 * variables.hilberta(variables.perp)/variables.SlopeRegge 
variables.matilde(2) == 4 * (variables.Hatilde(variables.perp))/variables.SlopeRegge 
variables.Zbf(variables.X, variables.tau) == (4 * math.pi**(2) * variables.SlopeRegge * variables.taub(2))**(-1/2) * (variables.q * variables.ql)**(-1/24)  * np.prod((summation(variables.qa(variables.n * variables.Nb(variables.n)) * variables.qla(variables.n * variables.Nbtilde(variables.n)), ((variables.Nb(variables.n), variables.Nbtilde(variables.n)), 1, math.inf))), (variables.n, 1, math.inf)) == (4 * math.pi**2 * variables.SlopeRegge * variables.taub(2))**(-1/2) * Abs(variables.etaf(variables.taub))**(-2) 
variables.psif(variables.w+ 2* math.pi) == math.exp(variables.brackets(math.pi * variables.i * (1 - variables.alpha))) * variables.psif(variables.w) 
# check rising operators 2.32.5 
variables.Trbf(variables.alpha, variables.qa(variables.H)) == variables.qa((3 * variables.alphaa(2) - 1)/24) * np.prod(variables.brackets(1 + variables.qa(variables.m-(1-variables.alpha)/2)) * variables.brackets(1 + variables.qa(variables.m - (1 + variables.alpha)/2)), (variables.m, 1, math.inf))
variables.Zabf(variables.alpha, variables.Beta, variables.tau) == variables.Trbf(variables.alpha, variables.brackets(variables.qa(variables.H) * math.exp(math.pi * variables.i * variables.Beta * variables.Q))) == np.cross(variables.qa((3 * variables.alphaa(2) - 1)/24) * math.exp(math.pi * variables.i * variables.alpha * variables.Beta/2), np.prod(variables.brackets(1 + math.exp(math.pi * variables.i * variables.Beta) * variables.qa(variables.m-(1 - variables.alpha)/2)) * variables.brackets(1 + math.exp(-math.pi * variables.i * variables.Beta) * variables.qa(variables.m-(1 + variables.alpha)/2)), (variables.m, 1, math.inf))) == 1/(variables.etaf(variables.tau)) * variables.jacobithetavar(((variables.alpha/2), (variables.Beta/2)), ((0, variables.tau)))
variables.Zabf(0, 0 ,variables.tau) == variables.Trbf(variables.NS, variables.brackets(variables.qa(variables.H)))
variables.Zabf(0, 1, variables.tau) == variables.Trbf(variables.NS, variables.brackets(math.exp(math.pi * variables.i * variables.F) * variables.qa(variables.H)))
variables.Zabf(1, 0, variables.tau) == variables.Trbf(variables.R,variables.brackets(variables.qa(variables.H)))
variables.Zabf(1,1, variables.tau) == variables.Trbf(variables.R, variables.brackets(math.exp(math.pi * variables.i * variables.F ) * variables.qa(variables.H)))
variables.Zabf(variables.posneg, variables.psi, variables.tau) == 1/2 * variables.brackets(variables.Zabf(0, 0, variables.tau)**4 - variables.Zabf(0, 1, variables.tau)**4 - variables.Zabf(1, 0, variables.tau)**4  + variables.negpos * variables.Zabf(1, 1, variables.tau)**4)
variables.Zb(variables.Tb(2)) == variables.i * variables.Vb(10) * integrate((variables.d2(variables.tau))/(16 * math.pi**2 * variables.SlopeRegge * variables.tauab(2, 2)) * variables.Zab(8, variables.X) * variables.Zabf(variables.pos, variables.psi, variables.tau) * variables.star(variables.Zabf(variables.posneg, variables.psi, variables.tau)), variables.F)
variables.psif(variables.w + 2 * math.pi) == - math.exp(-math.pi * variables.i * variables.alpha) * variables.psif(variables.w)
variables.psif(variables.w + 2 * math.pi * variables.tau) == -math.exp(-math.pi * variables.i * variables.Beta) * variables.psif(variables.w) 
variables.psif(variables.brackets(variables.w + 2 * math.pi * (variables.tau + 1))) == math.exp(variables.brackets(-math.pi * variables.i * (variables.alpha * variables.Beta))) * variables.psif(variables.w)
variables.psifp(variables.wp + 2 * math.pi) == - math.exp(-math.pi * variables.i * variables.Beta) * variables.psifp(variables.wp) 
variables.psifp(variables.wp - 2 * math.pi/variables.tau) == - math.exp(math.pi * variables.i * variables.alpha) * variables.psifp(variables.wp)
variables.Zabf(variables.alpha, variables.Beta, variables.tau) == variables.Zabf(variables.Beta, -variables.alpha, (-1/variables.tau)) == math.exp(variables.brackets(-math.pi * variables.i * (3 * variables.alphaa(2) - 1)/12)) * variables.Zabf(variables.alpha, variables.alpha + variables.Beta - 1, variables.tau + 1)
1/2 * variables.brackets(Abs(variables.Zabf(0, 0, variables.tau))**variables.N + Abs(variables.Zabf(0, 1, variables.tau))**variables.N + Abs(variables.Zabf(1, 0, variables.tau))**variables.N + variables.negpos * Abs(variables.Zabf(1, 1, variables.tau))**variables.N)
variables.Zabf(0, 0, variables.tau)**4 - variables.Zabf(0, 1, variables.tau)**4 - variables.Zabf(1, 0, variables.tau)**4 == 0
variables.alpha == variables.alphatilde 
math.exp(math.pi * variables.i * variables.F) == math.exp(math.pi * variables.i * variables.Ftilde) 
(variables.kb(variables.R), variables.kb(variables.L)) == (variables.nb(1), variables.nb(2))
# or
(variables.kb(variables.R), variables.kb(variables.L)) == (variables.nb(1) + 1/2, variables.nb(2) + 1/2)
diff(variables.H, variables.z) * diff(variables.H, variables.zl) == - variables.psil * variables.psi * variables.psiltilde * variables.psitilde
variables.kb(variables.R) == variables.m/(3**1/2) 
variables.kb(variables.L) == variables.n/(3**(1/2))
variables.m - variables.n == 3 * variables.ZB 
math.exp(variables.brackets(variables.posneg * variables.i * 3**(1/2) * variables.Hf(variables.z)))
math.exp(variables.brackets(variables.posneg * variables.i * 3**(1/2) * variables.Hftilde(variables.zl)))
(variables.H, variables.Htilde) == (variables.H, variables.Htilde) + (2 * math.pi)/(np.cross(2, 3**(1/2))) * (1, -1)
# 2.37.1, divergences of type I theory
variables.Zb(variables.Cb(2)) == variables.Zb(variables.Cb(2), 0) + variables.Zb(variables.Cb(2), 1)
variables.Zb(variables.Cb(2), 0) == variables.i * variables.Vb(10) * variables.na(2) * quad((diff(variables.t))/(8 * variables.t) * (8 * math.pi**2 * variables.SlopeRegge * variables.t)**(-5) * variables.etaf(variables.i * variables.t)**(-8) * variables.brackets(variables.Zabf(0, 0, variables.i * variables.t)**4 - variables.Zabf(1, 0, variables.i * variables.t)**4), 0, math.inf)
variables.Zb(variables.Cb(2), 1) == variables.i * variables.Vb(10) * variables.na(2) * quad((diff(variables.t))/(8 * variables.t) * (8 * math.pi**2 * variables.SlopeRegge * variables.t)**(-5) * variables.etaf(variables.i * variables.t)**(-8) * variables.brackets(-variables.Zabf(0, 1, variables.i * variables.t)**4 - variables.Zabf(0, 1, variables.i * variables.t)**4), 0, math.inf)
variables.etaf(variables.i * variables.t) == variables.t**(-1/2) * variables.etaf(variables.i/variables.t)
variables.Zabf(variables.alpha, variables.Beta, variables.i * variables.t) == variables.Zabf(variables.Beta, -variables.alpha, variables.i/variables.t)
variables.Zb(variables.Cb(2), 0) == variables.i * (variables.Vb(10) * variables.na(2))/(8 * math.pi * ( 8 * math.pi**2 * variables.SlopeRegge)**5) * quad(diff(variables.s) * variables.etaf(variables.i * variables.s/math.pi)**(-8) * variables.brackets(variables.Zabf(0, 0, variables.i * variables.s/math.pi)**4 - variables.Zabf(0, 1, variables.i * variables.s/math.pi)**4), 0, math.inf) == variables.i * (variables.Vb(10) * variables.na(2))/(8 * math.pi * (8 * math.pi * variables.SlopeRegge)**(5)) * quad(diff(variables.s) * variables.brackets(16 + variables.Of(math.exp(-2 * variables.s))), 0, math.inf)
variables.mub(10) * integrate(variables.Cb(10)) 
np.dot((1 + math.exp(math.pi * variables.i * variables.F))/2, (1 + math.exp(math.pi * variables.i * variables.Ftilde))/2) 
variables.psif(variables.w + 2 * math.pi * variables.i * variables.t) == - variables.R * variables.psif(variables.w) * variables.Ra(-1) == - math.exp(math.pi * variables.i * variables.Beta) * variables.psiftilde(variables.wl)
variables.psiftilde(variables.wl + 2 * math.pi * variables.i * variables.t) == - variables.R * variables.psiftilde(variables.wl) * variables.Ra(-1) == - math.exp(math.pi * variables.i * variables.Betatilde) * variables.psif(variables.w)
variables.psif(variables.w + 4 * math.pi * variables.i * variables.t) == math.exp(variables.brackets(math.pi * variables.i * (variables.Beta + variables.Betatilde))) * variables.psif(variables.w) 
variables.qa(-1/3) * np.prod((1 + variables.qa(2 * variables.m - 1))**(8), (variables.m, 1, math.inf)) == variables.Zabf(0, 0, 2 * variables.i * variables.t)**4 
-12 * variables.qa(2/3) * np.prod((1 + variables.qa(2 * variables.m))**8, (variables.m, 1, math.inf)) == - variables.Zabf(1, 0, 2 * variables.i * variables.t)**4 
variables.Zb(variables.Kb(2), 0) == variables.i * variables.Vb(10) * quad((diff(variables.t))/(8 * variables.t) * (4 * math.pi**2 * variables.SlopeRegge * variables.t)**(-5) * variables.etaf(2 * variables.i * variables.t)**(-8) * variables.brackets(variables.Zabf(0, 0, 2 * variables.i * variables.t)**4 - variables.Zabf(1, 0, 2 * variables.i * variables.t)**4), 0, math.inf) == variables.i * (2**10 * variables.Vb(10))/(8 * math.pi * (8 * math.pi**2 * variables.SlopeRegge)**5) * quad(diff(variables.s) * variables.etaf(variables.i * variables.s/math.pi)**(-8) * variables.brackets(variables.Zabf(0, 0, variables.i * variables.s/math.pi)**4 - variables.Zabf(0, 1, variables.i * variables.s/math.pi)**4), 0, math.inf) == variables.i * (2 **10 * variables.Vb(10))/(8 * math.pi * (8 * math.pi **2 * variables.SlopeRegge)**5) * quad(diff(variables.s) * variables.brackets(16 + variables.Of(math.exp(-2 * variables.s))), 0, math.inf) 
variables.omega * variables.psiaf(variables.mu, variables.w) * variables.omegaa(-1) == variables.psiaftilde(variables.mu, math.pi - variables.wl) == variables.psiaf(variables.mu, variables.w - math.pi)
variables.omega * variables.psiab(variables.mu, variables.r) * variables.omegaa(-1) == math.exp(-math.pi * variables.i * variables.r) * variables.psiab(variables.mu, variables.r) 
variables.omegaa(2) == math.exp(math.pi * variables.i * variables.F)
(1 + variables.omega + variables.omegaa(2) + variables.omegaa(3))/4 
variables.psiaf(variables.mu, variables.w + 4 * math.pi * variables.i * variables.t) == -math.exp(math.pi * variables.i * variables.Beta) * variables.psiaf(variables.mu, variables.w - 2 * math.pi * variables.i * variables.t - math.pi) == variables.psiaf(variables.mu, variables.w - 2 * math.pi)
-16 * variables.qa(1/2) * np.prod(variables.brackets(1 + (-1)**variables.m * variables.qa(variables.m))**8 - (1 - 1)**4 * variables.qa(1/2) * np.prod(variables.brackets(1 - (-1)**variables.m * variables.qa(variables.m))**8, (variables.m, 1, math.inf)), (variables.m, 1, math.inf)) == variables.Zabf(0, 1, 2 * variables.i * variables.t)**4 * variables.Zabf(1, 0, 2 * variables.i * variables.t)**4
variables.Zb(variables.Mb(2), 1) == variables.i * variables.n * variables.Vb(10) * quad((diff(variables.t))/(8 * variables.t) * (8 * math.pi**2 * variables.SlopeRegge * variables.t)**(-5) * (variables.Zabf(0, 1, 2 * variables.i * variables.t)**4 * variables.Zabf(1, 0, 2 * variables.i * variables.t)**4)/(variables.etaf(2 * variables.i * variables.t)*88 * variables.Zabf(0, 0, 2 * variables.i * variables.t)**4), 0, math.inf) == variables.posneg * 2 * variables.i * variables.n * (2**5  * variables.Vb(10))/(8 * math.pi * (8 * math.pi**2 * variables.SlopeRegge)**5) * quad(diff(variables.s) * (variables.Zabf(0, 1, 2 * variables.i * variables.t)**4 * variables.Zabf(1, 0, 2 * variables.i * variables.t)**4)/(variables.etaf(2 * variables.i * variables.t)*88 * variables.Zabf(0, 0, 2 * variables.i * variables.t)**4), 0, math.inf) == variables.posneg * 2 * variables.i * variables.n * (2**5 * variables.Vb(10))/(8 * math.pi * (8 * math.pi**2 * variables.SlopeRegge)**5) * quad(diff(variables.s) * variables.brackets(16 + variables.Of(math.exp(-2 * variables.s))), 0, math.inf) 
variables.Zb(1) == - variables.i * (variables.n + variables.negpos * 32)**2 * (variables.Vb(10))/(8 * math.pi * (8 * math.pi**2 * variables.SlopeRegge)**5) * quad(diff(variables.s) * variables.brackets(16 + variables.Of(math.exp(-2 * variables.s))), 0, math.inf) 
