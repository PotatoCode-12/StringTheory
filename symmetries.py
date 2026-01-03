import math
from scipy.integrate import quad
import numpy as np 
from sympy import symbols, diff,integrate, Function, Derivative, summation, oo, Abs, ln
from sympy.diffgeom import WedgeProduct
import sympy as sp
import variables
import constants

def isometries_1():
    if (variables.xap(variables.mu, variables.tau, variables.sigma)) == variables.lamda(variables.mu, variables.v)*variables.xav(variables.v, variables.tau, variables.sigma)+variables.a(variables.mu): #for translations a(mu)and lamda,isosymmetry of flat group
        variables.WorldSheetMetricP(variables.tau, variables.sigma) == variables.WorldSheetMetricE(variables.a*variables.b,variables.tau, variables.sigma) #for D-dimensional poincare invariance
def isometries_2():
    if variables.xap(variables.mu, variables.taup, variables.sigmap) == variables.xa(variables.mu, variables.tau, variables.sigma):
        diff(variables.sigmap**variables.c)/(diff(variables.sigma**variables.a))*diff(variables.sigmap**variables.d)/(diff(variables.sigma**variables.b))*variables.NWorldSheetMetric(variables.c*variables.d, variables.taup, variables.sigmap) == variables.WorldSheetMetric(variables.a*variables.b,variables.tau, variables.sigma)#diff invariance

if variables.xap(variables.mu, variables.tau, variables.sigma) == variables.xa(variables.mu, variables.tau, variables.sigma):
    variables.WorldSheetMetricP(variables.a*variables.b, variables.tau, variables.sigma) == math.e**(2*variables.omega(variables.tau, variables.sigma))*variables.WorldSheetMetric(variables.a*variables.b, variables.tau, variables.sigma)

variables.jf(variables.z) == variables.i*variables.vf(variables.z)*variables.Tf(variables.z)
variables.jftilde(variables.zl) == variables.i*variables.vf(variables.z) * variables.trans*variables.Tftilde(variables.zl) # check conditions pg 44 and similarities

#2.46.2, world-sheet supersymmetries
variables.cb(variables.h) == (-1)**(2 * variables.h + 1) * variables.brackets(3 * (2 * variables.h - 1)**2 - 1)
variables.cb(2) == - 26 
variables.cb(3/2) == +11 
variables.cb(1) == - 2 
variables.cb(1/2) == -1
variables.cb(0) == -2
variables.Tab(variables.posneg, variables.F) == 2**(-1/2) * (variables.Tb(variables.F1 + variables.posneg * variables.i * variables.Tb(variables.F2)))
# check algebra summing to without OPE calculations 2.47.4
variables.S == 1/(2 * math.pi) * integrate(variables.d2z * (diff(variables.Zl, variables.z) * diff(variables.Z, variables.zl) + variables.psil * diff(variables.psi, variables.zl) + variables.psiltilde * diff(variables.psitilde, variables.z)))
variables.Tb(variables.B) == - diff(variables.Zl, variables.z) * diff(variables.Z, variables.z) - 1/2 * (variables.psil * diff(variables.psi, variables.z) + variables.psi * diff(variables.psil, variables.z)) 
variables.j == - variables.psil * variables.psi 
variables.Tab(variables.pos, variables.F) == 2**(1/2) * variables.i * variables.psi * diff(variables.Zl, variables.z) 
variables.Tab(variables.neg, variables.F) == 2**(1/2) * variables.i * variables.psil * diff(variables.Z, variables.z)
variables.xa(variables.mu, (variables.z, variables.zl))
variables.psiaftilde(variables.mu, variables.zl)
#check values of mu 2.49.1
variables.lamdaaf(variables.A, variables.z) 
# check values of A ''
variables.S == 1/(4 * math.pi) * integrate(variables.d2z * (2/(variables.SlopeRegge) * diff(variables.xa(variables.mu), variables.z) * diff(variables.xb(variables.mu), variables.zl) + variables.lamdaa(variables.A) * diff(variables.lamdaa(variables.A), variables.zl) + variables.psiatilde(variables.mu) * diff(variables.psibtilde(variables.mu), variables.z)))
# check operator product approximations ''
variables.Tb(variables.B) == -1/(variables.SlopeRegge) * diff(variables.xa(variables.mu), variables.z) * diff(variables.xb(variables.mu), variables.z) - 1/2 * variables.lamdaa(variables.A) * diff(variables.lamdaa(variables.A), variables.z)
variables.Tbtilde(variables.B) == -1/(variables.SlopeRegge) * diff(variables.xa(variables.mu), variables.zl) * diff(variables.xb(variables.mu), variables.zl) - 1/2 * variables.psiatilde(variables.mu) * diff(variables.psibtilde(variables.mu), variables.zl)
variables.Tbtilde(variables.F) == variables.i * (2/variables.SlopeRegge)**(1/2) * variables.psiatilde(variables.mu) * diff(variables.xb(variables.mu), variables.zl)
variables.lamdaaf(variables.A, variables.w + 2 * math.pi) == variables.Oa(variables.A * variables.B) * variables.lamdaaf(variables.B, variables.w)
np.dot(variables.sB, variables.sBp) + 1/2 == variables.ZB 
math.exp(math.pi * variables.i * variables.Ftilde) == 1 
variables.lamdaaf(variables.A, variables.w + 2 * math.pi) == variables.posneg * variables.lamdaaf(variables.A, variables.w)
math.exp(math.pi * variables.i * variables.F) == 1
variables.Lamdaa(variables.posneg * variables.K) == 2**(1/2) * (variables.lamdaa(2 * variables.K - 1) + variables.posneg * variables.i * variables.lamdaa(2 * variables.K))
#check K states 2.51.11
variables.F == summation(variables.qb(variables.K), (variables.K, 1, 16))
variables.Zbf(16, variables.tau) == 1/2 * variables.brackets(variables.Zabf(0, 0, variables.tau)**16 + variables.Zabf(0, 1, variables.tau)**16 + variables.Zabf(1, 0, variables.tau)**16 + variables.Zabf(1, 1, variables.tau)**16)
# NS:
-1
# R: 
+1
variables.lamdaab(variables.A, -1/2) * variables.diracb(0, variables.NS) 
variables.alphaab(variables.i, -1) * variables.diracb(0, variables.NS) 
variables.lamdaab(variables.A, -1/2) * variables.lamdaab(variables.B, -1/2) * variables.diracb(0, variables.NS) 
np.cross((variables.B8b(variables.v), variables.B1), (variables.B8b(variables.v) + variables.B8)) == (variables.B1, variables.B1) + (variables.B28, variables.B1) + (variables.B35, variables.B1) + (variables.B56, variables.B1) + (variables.B8p, variables.B1)
np.cross((variables.B1, variables.B496), (variables.B8b(variables.v) + variables.B8)) == (variables.B8b(variables.v), variables.B496) + (variables.B8, variables.B496)
# check different functions for different A bounds, 2.53.19
variables.lamdaaf(variables.A, variables.w + 2 * math.pi) == variables.eta * variables.lamdaaf(variables.A, variables.w)
variables.lamdaaf(variables.A, variables.w + 2 * math.pi) == variables.etap * variables.lamdaaf(variables.A, variables.w)
math.exp(math.pi * variables.i * variables.Fb(1))
math.exp(math.pi * variables.i * variables.Fbp(1))
math.exp(math.pi * variables.i * variables.Fb(1)) == math.exp(math.pi * variables.i * variables.Fbp(1)) == math.exp(math.pi * variables.i * variables.Ftilde) == 1
variables.Zbf(8, variables.tau)**2 == 1/4 * variables.brackets(variables.Zabf(0, 0, variables.tau)**8 + variables.Zabf(0, 1, variables.tau)**8 + variables.Zabf(1, 0, variables.tau)**8 + variables.Zabf(1, 1, variables.tau)**8)**2
variables.alphaab(variables.i, -1) * variables.diracb(0, (variables.NS, variables.NSp)) 
variables.lamdaab(variables.A, -1/2) * variables.lamdaab(variables.B, -1/2) * variables.diracb(0, (variables.NS, variables.NSp))
# ordering constant == 0
(variables.B8b(variables.v), variables.B1, variables.B1) + (variables.B1, variables.B120, variables.B1) + (variables.B1, variables.B1, variables.B120) + (variables.B1, variables.B128, variables.B1) + (variables.B1, variables.B1, variables.B128)
math.exp(variables.brackets(variables.i * summation(variables.qb(variables.K) * variables.Haf(variables.K, variables.z), (variables.K, 1, 16))))
# check different charge states 2.55.27
variables.qb(variables.K) == +1/2 
variables.qb(variables.K) == 0
summation(variables.qb(variables.K), (variables.K, 1, 16)) == 2 * variables.ZB 
(variables.B1, variables.B1, variables.B1) + (variables.B28, variables.B1, variables.B1) + (variables.B35, variables.B1, variables.B1) + (variables.B56, variables.B1, variables.B1) + (variables.B8p, variables.B1, variables.B1) + (variables.B8b(variables.v), variables.B248, variables.B1) + (variables.B8, variables.B248, variables.B1) + (variables.B8b(variables.v), variables.B1, variables.B248) + (variables.B8, variables.B1, variables.B248)
#2.55.1
1/2 * variables.brackets(variables.Zabf(0,0, variables.tau)**16 * variables.star(variables.Zabf(0, 0, variables.tau))**4 - variables.Zabf(0, 1, variables.tau)**16 * variables.star(variables.Zabf(0, 1, variables.tau))**4 - variables.Zabf(1, 0, variables.tau)**16 * variables.star(variables.Zabf(1, 0, variables.tau))**4 - variables.Zabf(1, 1, variables.tau)**16 * variables.star(variables.Zabf(1, 1, variables.tau))**4)
math.exp(variables.brackets(math.pi * variables.i * (variables.F + variables.Ftilde))) == 1
variables.lamdaab(variables.A, -1/2) * variables.diracb(0, (variables.NS, variables.NS)) 
variables.ma(2) == - 2/variables.SlopeRegge 
math.exp(math.pi * variables.i * variables.F) == math.exp(math.pi * variables.i * variables.Ftilde) == -1
variables.alphaab(variables.i, -1) * variables.psiabtilde(variables.j, -1/2) * variables.diracb(0, (variables.NS, variables.NS))
variables.lamdaab(variables.A, -1/2) * variables.lamdaab(variables.B, -1/2) * variables.psiabtilde(variables.j, -1/2) * variables.diracb(0, (variables.NS, variables.NS))
np.dot((1 + math.exp(math.pi * variables.i * (variables.F + variables.Ftilde)))/(2), (1 + math.exp(math.pi * variables.i * variables.Ftilde))/2) == np.dot((1 + math.exp(math.pi * variables.i * variables.F))/2, (1 + math.exp(math.pi * variables.i * variables.Ftilde))/2) 
# zero point energy, 2.57.6
-1 + variables.k/16
variables.Z == 1/(variables.order(variables.H)) * summation(variables.Zb(variables.hb(1), variables.hb(2)), ((variables.hb(1), variables.hb(2)), variables.H)) # check advanced conditions 2.57.7 & 8
variables.Z == 1/(variables.order(variables.H)) * summation(variables.epsilonf(variables.hb(1), variables.hb(2)) * variables.Zb(variables.hb(1), variables.hb(2)), ((variables.hb(1), variables.hb(2)), variables.H))
variables.epsilonf(variables.hb(1), variables.hb(2)) == variables.epsilonf(variables.hb(2), variables.hb(1))**(-1)
variables.epsilonf(variables.hb(1), variables.hb(2)) * variables.epsilonf(variables.hb(1), variables.hb(3)) == variables.epsilonf(variables.hb(1), variables.hb(2), variables.hb(3))
variables.epsilonf(variables.h, variables.h) == 1 
variables.hat(variables.hb(2)) * variables.diracb(variables.psi, variables.hb(1)) == variables.epsilonf(variables.hb(1), variables.hb(2))**(-1) * variables.diracb(variables.psi, variables.hb(1))
variables.hat(variables.h) == variables.epsilonf(variables.hb(1), variables.h) * variables.hat(variables.h)
(variables.hb(1), variables.hb(2)) == (math.exp(variables.brackets(math.pi * variables.i * (variables.kb(1) * variables.Fb(1) + variables.lb(1) * variables.Ftilde))), math.exp(variables.brackets(math.pi * variables.i * (variables.kb(2) * variables.Fb(1) + variables.lb(2) * variables.Ftilde))))
variables.epsilonf(variables.hb(1), variables.hb(2)) == (-1)**(variables.kb(1) * variables.lb(2) + variables.kb(2) * variables.lb(1))
math.exp(math.pi * variables.i * variables.Fb(1)) == math.exp(math.pi * variables.i * variables.Fbp(1)) == math.exp(math.pi * variables.i * variables.Ftilde) == 1
math.exp(variables.brackets(math.pi * variables.i * (variables.Fb(1) + variables.alphabp(1) + variables.alphatilde))) == math.exp(variables.brackets(math.pi * variables.i * (variables.Fbp(1) + variables.alphab(1) + variables.alphatilde))) == math.exp(variables.brackets(math.pi * variables.i * (variables.Ftilde + variables.alphab(1) + variables.alphabp(1)))) == 1
# check spectrum tables, 2.58 & 9.tables

#2.282.4
variables.Qb(variables.alpha) == variables.Dfb(variables.Phi, variables.alpha * variables.Beta) * variables.Qb(variables.Beta)
variables.Qb(variables.s) == math.exp(np.dot(2 * math.pi * variables.i * variables.sB, variables.Phi)) * variables.Qb(variables.sB)
variables.Phib(2) + variables.Phib(3) + variables.Phib(4) == 0
variables.S * variables.Of(9, 1) == np.cross(variables.S * variables.Of(3, 1), variables.S * variables.Of(6)) == np.cross(variables.S * variables.Of(3, 1), variables.S * variables.Uf(3))
variables.B16 == (variables.B2, variables.B4) + (variables.B2l, variables.B4l) == (variables.B2, variables.B3) + (variables.B2, variables.B1) + (variables.B2l, variables.B3l) + (variables.B2l, variables.B1)
variables.Phib(2) + variables.Phib(3) == variables.Phib(4) == 0
variables.P == variables.S * variables.Uf(2) == variables.S * variables.Uf(3) == variables.S * variables.Of(6)

#2.302.1
variables.Gb(variables.M * variables.N) == np.array([[variables.ff(variables.y) * variables.etab(variables.mu * variables.v), 0],
                                                    [0, variables.Gbf(variables.m * variables.n, variables.y)]])
variables.delta * variables.psib(variables.mu) == variables.Deltab(variables.mu) * variables.epsilon 
variables.delta * variables.psib(variables.m) == diff(1, variables.m) + 1/4 * variables.omegab(variables.neg, variables.m * variables.n * variables.p) * variables.Rhoa(variables.n * variables.p) * variables.epsilon
variables.delta * variables.chi == variables.Rhoa(variables.m) * diff(variables.PhiB, variables.m) - 1/12 * variables.Rhoa(variables.m * variables.n * variables.p) * variables.Hbtilde(variables.m * variables.n * variables.p) * variables.epsilon
variables.delta * variables.lamda == variables.Fb(variables.m * variables.n) * variables.Rhoa(variables.m * variables.n) * variables.epsilon
variables.omegaab(variables.posneg, (variables.M, variables.N, variables.P)) == variables.omegabs(variables.M, variables.N, variables.P) + variables.posneg * 1/2 * variables.Hb(variables.M, variables.N, variables.P)
variables.B16 == (variables.B2, variables.B4) + (variables.B2l, variables.B4l)
variables.epsilonf(variables.y) == variables.epsilonbf(variables.alpha * variables.Beta, variables.y) + variables.star(variables.epsilonbf(variables.alpha * variables.Beta, variables.y))
variables.epsilonb(variables.alpha * variables.Beta) == variables.ub(variables.alpha) * variables.zetabf(variables.Beta, variables.y)
variables.Hbtilde(variables.m * variables.n * variables.p) == 0
diff(variables.PhiB, variables.m) == 0
variables.Gb(variables.mu * variables.v) == variables.etab(variables.mu * variables.v)
variables.Deltab(variables.m) * variables.zeta == 0
variables.brackets(variables.Deltab(variables.m), variables.Deltab(variables.n)) * variables.zeta == 1/4 * variables.Rb(variables.m * variables.n * variables.p * variables.q) * variables.Rhoa(variables.p * variables.q) * variables.zeta == 0
variables.Fb(variables.i * variables.j) == variables.Fb(variables.il * variables.jl) == 0
variables.Ga(variables.i * variables.jl) * variables.Fb(variables.i * variables.jl) == 0
diff(variables.Hbtilde(3)) == (variables.SlopeRegge)/4 * variables.brackets(variables.trf(WedgeProduct(variables.Rb(2), variables.Rb(2))) - variables.Trbf(variables.v, WedgeProduct(variables.Fb(2), variables.Fb(2))))
variables.Rb(variables.m * variables.n) == 0
variables.Deltaa(variables.m) * variables.Fb(variables.m * variables.n) == 0
# 2.305.2, Calabi-Yau manifolds
variables.chif(variables.K) == summation((-1)**(variables.p) * variables.bb(variables.p), (variables.p, 0, variables.d))
variables.Deltab(variables.d) == variables.starf(variables.d) * variables.d + variables.d * variables.star(variables.d) == (variables.d + variables.starf(variables.star(variables.d)))**2
variables.bb(variables.p) == variables.b(variables.d - variables.p)
integrate(WedgeProduct(variables.omegabs(variables.p), variables.alphab(variables.d - variables.p)), variables.K) == integrate(variables.alphab(variables.d - variables.p), (variables.Nf(variables.omegas)))
variables.Gb(variables.i * variables.j) == variables.Gb(variables.il * variables.jl) == 0
diff(1, variables.z) == diff(variables.za(variables.i)) * diff(1, variables.i)
variables.diff(1, variables.zl) == diff(variables.zla(variables.il)) * diff(1, variables.il)
diff(1, variables.z, 2) == diff(1, variables.zl, 2) == 0
variables.Deltab(diff(1, variables.z)) == diff(variables.dagger(diff(1, variables.z)), variables.z) + variables.dagger(diff(diff(1, variables.z), variables.z)) 
variables.Deltab(diff(1, variables.zl)) == diff(variables.dagger(diff(1, variables.zl)), variables.zl) + variables.dagger(diff(diff(1, variables.zl), variables.zl))
variables.Jb(1, 1) == variables.i * variables.Gb(variables.i * variables.jl) * diff(variables.za(variables.i)) * diff(variables.zla(variables.jl))
diff(variables.Jb(1, 1)) == 0
variables.Gb(variables.i * variables.jl) == (diff(1, variables.z))/(diff(variables.za(variables.i), variables.z)) * (diff(1, variables.z))/(diff(variables.zla(variables.jl), variables.z)) * variables.Kf(variables.z, variables.zl)
variables.Kfp(variables.z, variables.zl) == variables.Kf(variables.z, variables.zl) + variables.ff(variables.z) + variables.star(variables.ff(variables.z))
variables.Deltab(variables.d) == 2 * variables.Deltab(diff(1, variables.zl)) == 2 * variables.Deltab(diff(1, variables.z))
variables.Habf((variables.p, variables.q), diff(1, variables.zl), variables.K) == variables.Habf((variables.p, variables.q), diff(1, variables.z), variables.K) == variables.Haf((variables.p, variables.q), variables.K)
variables.bb(variables.k) == summation(variables.ha(variables.p, variables.k-variables.p), (variables.p, 0, variables.k))
variables.ha(variables.p, variables.q) == variables.ha(variables.q, variables.p)
variables.ha(variables.n - variables.p, variables.n - variables.q) == variables.ha(variables.p, variables.q)
variables.Jb(1, 1) == summation(variables.va(variables.A) *variables.omegabs(1, (1, variables.A)), variables.A)
variables.ROb(1, 1) == variables.Rb(variables.i * variables.jl) * diff(variables.za(variables.i)) * diff(variables.zl, variables.jl)
variables.ha(variables.p, 0) == variables.ha(3 - variables.p, 0)
variables.bb(1) == variables.ha(1, 0) == variables.ha(0, 1) == 0
# 2.309.29, Hodge diamond, temp ommit
variables.chi == 2 * (variables.ha(1, 1) - variables.ha(2, 1))
# 2.312.1, massless spectrum
variables.Deltab(variables.M) * variables.Deltaa(variables.M) == diff(1, variables.mu) * diff(1, variables.z, variables.mu) + variables.Deltab(variables.m) * variables.Deltaa(variables.m)
variables.Rhob(variables.M) * variables.Deltaa(variables.M) == variables.Rhob(variables.mu) * diff(1, variables.z, variables.mu) + variables.Rhob(variables.m) * variables.Deltaa(variables.m)
np.cross(variables.Eb(8), variables.Eb(8)) == np.cross(variables.S * variables.Uf(3), variables.Eb(6), variables.Eb(6))
(variables.B1, variables.B78, variables.B1) + (variables.B1, variables.B1, variables.B248)
(variables.B3, variables.B27, variables.B1)
(variables.B3l, variables.B27l, variables.B1)
(variables.B8, variables.B1, variables.B1)
variables.ab(variables.i * variables.ll * variables.ml * variables.x) == variables.ab(variables.i, variables.j * variables.x) * variables.Ga(variables.j * variables.kl) * variables.omegab(variables.kl * variables.ll * variables.ml)
variables.B16 == (variables.B2, variables.B1) + (variables.B2, variables.B3) + (variables.B2l, variables.B1l) + (variables.B2l, variables.B3l)
Abs(variables.ha(2, 1) - variables.ha(1, 1)) == (Abs(variables.chi))/2
variables.gb(variables.i * variables.ll * variables.ml) == variables.gb(variables.i, variables.j) * variables.Ga(variables.j * variables.kl) * variables.omegab(variables.kl * variables.ll * variables.ml)
# 2.315.1, new thing idk
variables.phif(variables.x, variables.y) == summation(variables.Phibf(variables.m, variables.x) * variables.fbf(variables.m, variables.y), variables.m)
variables.LObf(4, variables.Phi) == integrate(variables.da(6) * variables.y * variables.LObf(10, variables.phi))
variables.phi == variables.phib(1) + variables.phib(variables.h)
variables.m * variables.phiab(2, variables.h) + variables.g * variables.phib(variables.h) * variables.phiab(2, 1)
-(variables.ga(2))/(4 * variables.m) * variables.phiab(4, 1)
variables.abf((variables.i, variables.jl * variables.xl), (variables.x, variables.y)) == summation(variables.Phiabf(variables.A, variables.xl, variables.x) * variables.omegabfs(variables.A * variables.i * variables.jl, variables.y), variables.A)
variables.lamdabf((variables.i, variables.jl * variables.xl), (variables.x, variables.y)) == summation(variables.lamdaabf(variables.A, variables.xl, variables.x) * variables.omegabfs(variables.A * variables.il * variables.jl, variables.y), variables.A)
integrate(variables.da(6) * variables.y * variables.Trbf(variables.v, variables.lamdal * variables.Rhoa(variables.m) * variables.brackets(variables.Ab(variables.m),variables.lamda)))
variables.da(variables.xl * variables.yl * variables.zl) * variables.lamdalab(variables.A, variables.xl) * variables.lamdaab(variables.B, variables.yl) * variables.Phiab(variables.C, variables.zl) * integrate(WedgeProduct(variables.omegabs(1, (1, variables.A)), variables.omegabs(1, (1, variables.B)), variables.omegabs(1, (1, variables.C))), variables.K)
variables.Wf(variables.Phi) == variables.da(variables.xl * variables.yl * variables.zl) * variables.Phiab(variables.A, variables.xl) * variables.Phiab(variables.B, variables.yl) * variables.Phiab(variables.C, variables.zl) * integrate(WedgeProduct(variables.omegabs(1, (1, variables.A)), variables.omegabs(1, (1, variables.B)), variables.omegabs(1, (1, variables.C))), variables.K)
# check 2.317.10, how to implement #
integrate(WedgeProduct(variables.omegabs(1, (1, variables.A)), variables.omegabs(1, (1, variables.B)), variables.omegabs(1, (1, variables.C))), variables.K)
integrate(variables.omegabs(1, (1, variables.B)), variables.Na(variables.A)) == variables.deltaab(variables.A, variables.B)
# - 
(variables.gb(variables.i * variables.j) + variables.bb(variables.i * variables.j)) * (variables.x, variables.y) == summation(variables.Ta(variables.A, variables.x) * variables.omegabfs(variables.A * variables.i * variables.jl, variables.y), variables.A)
variables.Gb(variables.A * variables.Bl) == 1/(variables.V) * integrate(variables.da(6) * variables.y * (sp.Determinant(variables.G))**(1/2) * variables.Ga(variables.i * variables.kl) * variables.Ga(variables.l * variables.jl) * variables.omegabs(variables.A * variables.i * variables.jl) * variables.star(variables.omegabs(variables.B * variables.k * variables.l)))
variables.Jb(1, 1) + variables.i * variables.Bb(1, 1) == variables.Ta(variables.A) * variables.omegabs(1, (1, variables.A))
variables.Ta(variables.A) == variables.va(variables.A) + variables.i * variables.ba(variables.A)
variables.Gb(variables.A * variables.Bl) == - (diff(1, variables.z, 2))/(diff(variables.Ta(variables.A), variables.z) * diff(variables.star(variables.Ta(variables.B)), variables.z)) * ln(variables.Wf(variables.v))
variables.Wf(variables.c) == integrate(WedgeProduct(variables.Jb(1, 1), variables.Jb(1, 1), variables.Jb(1, 1)), variables.K) # check # thing
variables.Kbf(1, (variables.T, variables.star(variables.T))) == - ln(variables.Wf(variables.v))
variables.Gbp(variables.A * variables.Bl) == math.exp(variables.brackets(variables.kappaab(2, 4) * (variables.Kb(2) - variables.Kb(1))/3)) * variables.Gb(variables.A * variables.Bl)
variables.abf((variables.i, variables.j * variables.x), (variables.x, variables.y)) == 1/2 * summation(variables.chiabf(variables.a, variables.x, variables.x) * variables.omegabfs(variables.a * variables.i * variables.jl * variables.ll, variables.y) * variables.omegaabf(variables.kl * variables.ll, variables.j, variables.y), variables.a)
variables.Gb(variables.a * variables.bl) == -(integrate(WedgeProduct(variables.omegabs(1, (2, variables.a)), variables.star(variables.omegabs(1, (2, variables.b)))), variables.K))/(integrate(WedgeProduct(variables.omegab(3, 0), variables.star(variables.omegab(3, 0))), variables.K)) == - (diff(1, variables.z))/(diff(variables.xa(variables.a), variables.z)) * (diff(1, variables.z))/(diff(variables.xa(variables.bl), variables.z)) * variables.Kbf(2, (variables.X, variables.star(variables.X)))
variables.Kbf(2, (variables.X, variables.star(variables.X))) == ln(variables.i * integrate(WedgeProduct(variables.omegab(3, 0), variables.star(variables.omegab(3, 0))), variables.K))
variables.bracketc(variables.Aa(variables.I), variables.Bb(variables.J))
# check 2.319.23
variables.Za(variables.I) == integrate(variables.omegab(3, 0), variables.Aa(variables.I))
variables.gOb(variables.I) == (diff(variables.gO, variables.z))/(diff(variables.Za(variables.I), variables.z)) 
variables.gOf(variables.lamda * variables.Z) == variables.lamdaa(2) * variables.gOf(variables.Z)
variables.Kbf(2, (variables.Z, variables.star(variables.Z))) == ln(variables.Im(variables.star(variables.Za(variables.I)) * diff(variables.gOf(variables.Z), variables.I)))
variables.Wf(variables.Z, variables.chi) == (variables.chia(variables.a) * variables.chia(variables.b) * variables.chia(variables.c))/(math.factorial(3)) * (diff(variables.gOf(variables.Z), variables.z, 3))/(diff(variables.Za(variables.a), variables.z) * diff(variables.Za(variables.b), variables.z) * diff(variables.Za(variables.c), variables.z))
variables.Gbp(variables.a * variables.bl) == math.exp(variables.brackets(variables.kappaab(2, 4) * (variables.Kb(1) - variables.Kb(2))/3)) * variables.Gb(variables.a * variables.bl)
np.array([[variables.Apa(variables.I)],
          [variables.Bpb(variables.J)]]) == variables.S * np.array([[variables.Aa(variables.I)],
                                                                   [variables.Bb(variables.J)]])
np.array([[variables.Zpa(variables.I)],
          [variables.gOpb(variables.J)]]) == variables.S * np.array([[variables.Za(variables.I)],
                                                                    [variables.gOb(variables.J)]])
# new thing, 2.322.1
1/(2 * math.pi * variables.SlopeRegge) * integrate(variables.Bb(1, 1), variables.M) == (variables.nb(variables.A) * variables.ba(variables.A))/(2 * math.pi * variables.SlopeRegge)
variables.Ta(variables.A) == variables.Ta(variables.A) + variables.i * variables.epsilona(variables.A)
variables.Ta(variables.A) == variables.ca(variables.A) * variables.T
-3 * ln(variables.T + variables.star(variables.T))
1/(2 * math.pi * variables.SlopeRegge) * integrate(variables.d2z * variables.Gb(variables.i * variables.jl) * (diff(variables.Za(variables.i), variables.z) * diff(variables.Za(variables.jl), variables.zl) + diff(variables.Za(variables.il), variables.z) * diff(variables.Za(variables.j), variables.zl)))
1/(2 * math.pi * variables.SlopeRegge) * integrate(variables.Jb(1, 1), variables.m) == (variables.nb(variables.A) * variables.va(variables.A))/(2 * math.pi * variables.SlopeRegge)
math.exp(-variables.nb(variables.A) * variables.Ta(variables.A)/(2 * math.pi * variables.SlopeRegge))

# 2.327.1
variables.Q == 1/(2 * math.pi * variables.i) * integrate((diff(variables.z) * variables.jb(variables.z) - diff(variables.zl) * variables.jb(variables.zl)), 2 * math.pi) # contour
variables.jb(variables.z) * diff(variables.xa(variables.mu), variables.zl) * variables.ea(np.dot(variables.i * variables.k, variables.X)) 
diff(variables.xa(variables.mu), variables.z) * variables.jb(variables.zl) * variables.ea(np.dot(variables.i * variables.k, variables.X))
variables.Q == 1/(2 * math.pi * variables.i) * integrate((diff(variables.z) * diff(variables.theta) * variables.J - diff(variables.zl) * diff(variables.thetal) * variables.Jtilde), 2 * math.pi) # contour
variables.Rbf(8, 2 * math.pi * variables.Rb(9)) == variables.SlopeRegge/(variables.Rbf(8, 0))
variables.PhiBf(2 * math.pi * variables.Rb(9)) == - variables.PhiBf(0)
variables.gf(2 * math.pi * variables.Rb(9)) == 1/(variables.gf(0))
variables.Nb(variables.pos * 1/2, variables.r) - variables.Nb(variables.neg * 1/2, variables.r) == variables.Trbf((variables.r, variables.kerf(variables.Gb(0))), variables.brackets(variables.i * math.exp(math.pi * variables.i * variables.Fbtilde(variables.K))))
variables.Nb(variables.pos * 1/2, variables.r) - variables.Nb(variables.neg * 1/2, variables.r) == variables.Trbf(variables.r, variables.brackets(variables.i * math.exp(math.pi * variables.i * variables.Fbtilde(variables.K))))
variables.bracket((variables.alpha, variables.out), (variables.Beta, variables.inn)) == variables.bracket(variables.VObl(variables.alpha) * variables.VOb(variables.Beta))
variables.bracket((variables.alpha, variables.out), (variables.Beta, variables.inn)) == variables.bracket(np.dot(variables.theta, variables.VObl(variables.alpha) * variables.theta, variables.VOb(variables.Beta))) == variables.bracket((variables.theta * variables.Betal, variables.out), (variables.theta * variables.alphal, variables.inn))
variables.bracket((np.dot(variables.C * variables.P * variables.T, variables.Beta), variables.out), (np.dot(variables.C * variables.P * variables.T, variables.alpha), variables.inn)) == variables.bracket((variables.alpha, variables.out), (variables.Beta, variables.inn))
variables.SBb(variables.theta) == (variables.theta)/(8 * math.pi**2) * integrate(WedgeProduct(variables.Fb(2), variables.Fb(2)))
Abs(variables.theta) < 10**(-9)
1/(2 * variables.gab(2, 4)) * integrate(WedgeProduct(variables.a * variables.Fab(variables.a, 2), variables.Fab(variables.a, 2)))
variables.a == variables.a + variables.epsilon 
math.exp(-variables.nb(variables.A) * variables.Ta(variables.A)/(2 * math.pi * variables.SlopeRegge)) == math.exp(variables.brackets(-variables.nb(variables.A) * (variables.va(variables.A) + variables.i * variables.ba(variables.A))/(2 * math.pi * variables.SlopeRegge)))
# 2.335.1, gauge symmetries
variables.gab(2, 10) == (4 * variables.kappaab(2, 10))/(variables.SlopeRegge)
variables.gab(2, 4) == (variables.gab(2, 10))/(variables.V)
variables.kappaab(2, 4) == (variables.kappaab(2, 10))/(variables.V)
variables.gab(2, 4) == (4 * variables.kappaab(2, 4))/(variables.SlopeRegge)
variables.gab(2, variables.YM) == (4 * variables.kappaa(2))/(variables.SlopeRegge)
variables.gab(2, variables.YM) == (2 * variables.kappaa(2))/(variables.hat(variables.k) * variables.SlopeRegge)
-1/(4 * variables.gab(2, variables.YM)) * variables.Fab(variables.a, variables.mu * variables.v) * variables.Fa(variables.a * variables.mu * variables.v)
(variables.gab(2, variables.YM))/(variables.kappa) == 2 * (2 * math.pi)**(7/2) * variables.SlopeRegge # check condition 2.336.7
(variables.gab(2, variables.YM))/(variables.kappa) == (2 * (2 * math.pi)**(7/2) * variables.SlopeRegge)/(variables.Va(1/2)) # condition
variables.Tb(variables.B) == variables.Tab(variables.s, variables.B) + variables.Tbp(variables.B)
np.dot(variables.jab(variables.a, -1), 1) == variables.ja(variables.a)
diff(variables.xa(variables.mu), variables.z) * variables.psiatilde(variables.a) * variables.ea(np.dot(variables.i * variables.k, variables.X))
np.dot(variables.Gbtilde(-1/2), variables.psiatilde(variables.a)) == variables.jatilde(variables.a)
np.dot(variables.Gbtilde(1/2), variables.jatilde(variables.a)) == np.dot(2 * variables.Lbtilde(0), variables.psiatilde(variables.a)) == variables.psiatilde(variables.a)
variables.jatilde(variables.a) == variables.jabtilde(variables.a, variables.psi) + variables.japtilde(variables.a)
variables.jabtilde(variables.a, variables.psi) == -(variables.i)/(2 * variables.k) * variables.fa(variables.a * variables.b * variables.c) * variables.psiatilde(variables.b) * variables.psiatilde(variables.c)
variables.Tbtilde(variables.F) == variables.Tabtilde(variables.s, variables.F) + variables.Tbpptilde(variables.F)
variables.Tabtilde(variables.s, variables.F) == -(variables.i)/(6 * variables.ka(2)) * variables.fa(variables.a * variables.b * variables.c) * variables.psiatilde(variables.a) * variables.psiatilde(variables.b) * variables.psiatilde(variables.c) + 1/(variables.k) * variables.psiatilde(variables.a) * variables.japtilde(variables.a)
variables.Tbtilde(variables.B) == variables.Tabtilde(variables.psi, variables.B) + variables.Tbptilde(variables.B) + variables.Tbpptilde(variables.B)
variables.Tabtilde(variables.psi, variables.B) == -1/(2 * variables.k) * variables.psiatilde(variables.a) * diff(variables.psiatilde(variables.a), variables.zl)
variables.Tbptilde(variables.B) == 1/(2 * (variables.kp + variables.hf(variables.g))) * variables.point(variables.jptilde * variables.jptilde)
variables.catilde(variables.psi) == (variables.dimf(variables.g))/2 
variables.cptilde == (variables.kp * variables.dimf(variables.g))/(variables.kp + variables.hf(variables.g))
variables.cpptilde == variables.ctilde - variables.catilde - variables.cptilde
variables.catilde(variables.psi) + variables.cptilde == ((3 * variables.kp + variables.hf(variables.g)) * variables.dimf(variables.g))/(2 * (variables.kp + variables.hf(variables.g)))
(variables.dimf(variables.g))/2 <= variables.catilde(variables.psi) + variables.cptilde <= (3 * variables.dimf(variables.g))/2
variables.Tabtilde(variables.s, variables.F) == variables.i * variables.psitilde * diff(variables.H, variables.zl)
variables.Tabtilde(variables.psi, variables.B) + variables.Tbptilde(variables.B) == -1/2 * variables.psitilde * diff(variables.psitilde, variables.zl) - 1/2 * diff(variables.H, variables.zl) * diff(*variables.H, variables.zl)
variables.Sb(variables.alpha) * variables.VOb(variables.K) * variables.ea(variables.i * variables.kb(variables.mu) * variables.xa(variables.mu))
variables.Gabtilde(2, 0) == variables.Lbtilde(0) - (variables.ctilde)/24 >= 0
variables.Labtilde(variables.i, 0) == (variables.catilde(variables.i))/(24)
variables.Labtilde(variables.psi, 0) + variables.Lbptilde(0) - (variables.catilde(variables.psi) + variables.cptilde)/24 >= (variables.hf(variables.g) * variables.dimf(variables.g))/(24 * (variables.kp + variables.hf(variables.g))) > 0
variables.psia(variables.a) * variables.psiatilde(variables.mu) * variables.ea(np.dot(variables.i * variables.k, variables.X))
variables.psia(variables.mu) * variables.psiatilde(variables.a) * variables.ea(np.dot(variables.i * variables.k, variables.X))
variables.ctilde >= 8/2 + variables.catilde(variables.S * variables.Uf(3), 1) + 3/2 + variables.catilde(variables.S * variables.Uf(2), 1) + variables.catilde(variables.Uf(1)) == 10
# 2.343.1
(variables.mb(variables.s))/(variables.mb(variables.grav)) == variables.gb(variables.YM) * (variables.hat(variables.k)/2)**(1/2)
variables.Y/2 == variables.diag(-1/3, -1/3, -1/3, 1/2, 1/2)
variables.gb(3) == variables.gb(2) == variables.gb(1) == variables.gb(variables.S * variables.Uf(5))
1/2 * variables.gp * variables.Y == variables.gb(variables.Uf(1)) * variables.ta(variables.Uf(1)) #changes to
variables.gp == (3/5)**(1/2) * variables.gb(1)
(5/3)**(1/2) * variables.gp == variables.gb(2) == variables.gb(3)
variables.mu * (diff(1, variables.z))/(diff(variables.mu, variables.z)) * variables.gb(variables.i) == (variables.bb(variables.i))/(16 * math.pi**2) * variables.gab(3, variables.i)
variables.alphaabf(-1, variables.i, variables.mu) == variables.alphaabf(-1, variables.i, variables.mb(variables.GUT)) + (variables.bb(variables.i))/(4 * math.pi) * ln((variables.mab(2, variables.GUT))/(variables.mua(2)))
# check 2.346.8, how to implement summation bounds
(math.sin(variables.thetabf(variables.w, variables.mb(variables.Z))))**2 == 0.212 + variables.posneg * 0.003
(math.sin(variables.thetabf(variables.w, variables.mb(variables.Z))))**2 == 0.234 + variables.posneg * 0.003
(math.sin(variables.thetabf(variables.w, variables.mb(variables.Z))))**2 == 0.2313 + variables.posneg * 0.0003
variables.mb(variables.GUT) == 10**(16.1 + variables.posneg * 0.3) 
variables.mb(variables.SU) == np.cross(variables.ka(1/2) * variables.gb(variables.YM), 5.27 * 10**(17)) == 3.8 * 10**17 # GeV
np.cross(variables.S* variables.Ufp(5), variables.Uf(1)) == variables.S * variables.Of(10)
np.cross(variables.S * variables.Uf(4), variables.S * variables.Ufb(2, variables.L), variables.S * variables.Ufb(2, variables.R)) == 2 * variables.Of(10)
np.cross(variables.S * variables.Ufb(3, variables.C), variables.S * variables.Ufb(3, variables.L), variables.S * variables.Ufb(3, variables.R)) == variables.Eb(6)
variables.Deltab(variables.a) == variables.cb(variables.a) - summation((variables.bab(variables.i, variables.b) * Abs(variables.Ga(variables.i)))/(Abs(variables.G)) * ln(variables.brackets((variables.Tb(variables.i) + variables.star(variables.Tb(variables.i))) * Abs(variables.etaf(variables.Tb(variables.i)))**4 * (variables.Ub(variables.i) + variables.star(variables.Ub(variables.i))) * Abs(variables.etaf(variables.Ub(variables.i)))**4)), variables.i)
np.cross(variables.K, (variables.Sb(1))/(variables.ZBb(2)))
# mehr auf unification 2.351.1
(variables.alphab(2))/(variables.alphab(3)) == (variables.hat(variables.kb(3)))/(variables.hat(variables.kb(2))) == (variables.kb(3))/(variables.kb(2))
diff(variables.H, variables.z) * variables.psitilde * variables.ea(np.dot(variables.i * variables.k, variables.X))
variables.Qp == variables.Qb(variables.EM) + variables.T/3
variables.Qp == variables.Qb(variables.EM) + 1/3 * variables.T == 1/2 * variables.Y + variables.Ib(3) + 1/3 * variables.T 
np.diag([-1/3, -1/3, -1/3, 1/2, 1/2]) + np.diag([0, 0, 0, 1/2, 1/2]) + variables.diag([1/3, 1/3, -2/3, 0, 0]) == np.diag([0, 0, -1, 1, 0])
variables.jp == variables.lamdaa(variables.neg * 6) * variables.lamdaa(variables.pos * 6) - variables.lamdaa(variables.neg * 7) * variables.lamdaa(variables.pos * 7) == variables.i * diff(variables.Hb(7) - variables.Hb(6), variables.z)
variables.lamdaaf(variables.pos * variables.K, variables.sigmab(1) + 2 * math.pi) == math.exp(2 * math.pi * variables.i * variables.vb(variables.K)) * variables.lamdaaf(variables.pos * variables.K, variables.sigmab(1))
math.exp(variables.brackets(variables.i * summation((1/2 - variables.vb(variables.K)) * variables.Hb(variables.K), variables.K)))
variables.Qp == variables.vb(6) - variables.vb(7)
variables.jab(3, variables.S * variables.Uf(3)) == variables.i/2 * diff(variables.Ha(4) - variables.Ha(5), variables.z)
variables.jab(8, variables.S * variables.Uf(3)) == (variables.i)/(np.cross(2, 3**(1/2))) * diff(variables.Ha(4) + variables.Ha(5) - 2 * variables.Ha(6), variables.z)
variables.jab(3, variables.S * variables.Uf(2)) == variables.i/2 * diff(variables.Ha(7) - variables.Ha(8), variables.z)
variables.jb(variables.Y/2) == variables.i/6 * diff(variables.brackets(-2 * (variables.Ha(4) + variables.Ha(5) + variables.Ha(6)) + 3 * (variables.Ha(7) + variables.Ha(8))), variables.z)
variables.jp == variables.jb(variables.Y/2) + variables.jab(3, variables.S * variables.Uf(2)) + 2/(3**(1/2)) * variables.jab(8, variables.S * variables.Uf(3)) == variables.i * diff(variables.Ha(7) - variables.Ha(6), variables.z)
math.exp(variables.brackets(variables.i * (variables.Ha(6) - variables.Ha(7))))
variables.mub(1) * variables.Hb(1) * variables.L + variables.etab(1) * variables.Ua(variables.c) * variables.Da(variables.c) * variables.Da(variables.c) + variables.etab(2) * variables.Q * variables.L * variables.Da(variables.c) + variables.etab(3) * variables.L * variables.L * variables.Ea(variables.c) + (variables.lamdab(1))/(variables.M) * variables.Q * variables.Q * variables.Q * variables.L + (variables.lamdab(2))/(variables.M) * variables.Ua(variables.c) * variables.Ua(variables.c) * variables.Da(variables.c) * variables.Ea(variables.c) + (variables.lamdab(3))/(variables.M) * variables.L * variables.L * variables.Hb(2) * variables.Hb(2)
# spacetime supersymmetry 2.356.1
variables.Jbtilde(variables.alpha) == variables.ea(-variables.Phitilde/2) * variables.Sbtilde(variables.alpha) * variables.Sigmatilde
variables.Jbtilde(diff(variables.alpha, variables.tau)) == variables.ea(-variables.Phitilde/2) * variables.Sbtilde(diff(variables.alpha)) * variables.Sigmaltilde
# 2.357.2 - 4 out
variables.Sigmaftilde(variables.zl) * variables.Sigmaftilde(0) == variables.Of(variables.zla(3/4))
variables.bracket(variables.Sigmaftilde(variables.zlb(1)) * variables.Sigmalftilde(variables.zlb(2)) * variables.Sigmaftilde(variables.zlb(3)) * variables.Sigmalftilde(variables.zlb(4))) == ((variables.zlb(13) * variables.zlb(24))/(variables.zlb(12) * variables.zlb(14) * variables.zlb(23) * variables.zlb(34)))**(3/4) * variables.ff(variables.zlb(1), variables.zlb(2), variables.zlb(3), variables.zlb(4))
variables.brackets(variables.jftilde(variables.zlb(2)) * variables.Sigmaftilde(variables.zlb(3)) * variables.Sigmalftilde(variables.zlb(4))) == (3 * variables.zlab(1/4, 34))/(2 * variables.zlb(23) ** variables.zlb(24))
variables.jftilde(variables.zl) == 3**(1/2) * variables.i * diff(variables.Hftilde(variables.zl), variables.zl)
variables.Tbtilde(variables.B) == - 1/2 * diff(variables.Htilde, variables.zl) * diff(variables.Htilde, variables.zl) + variables.Tbptilde(variables.B)
variables.Sigmatilde == math.exp(3**(1/2) * variables.i * variables.Htilde/2) * variables.Sigmaptilde
variables.Sigmatilde == math.exp(3**(1/2) * variables.i * variables.Htilde/2)
variables.Sigmaltilde == math.exp(-3**(1/3) * variables.i * variables.Htilde/2)
variables.Tbftilde(variables.F, variables.zl) * variables.Sigmaftilde(0) == variables.Of(variables.zla(-1/2))
variables.Tbftilde(variables.zl) * variables.Sigmalftilde(0) == variables.Of(variables.zla(-1/2))
variables.Tbtilde(variables.F) == variables.Tabtilde(variables.pos, variables.F) + variables.posneg * variables.Tabtilde(variables.neg, variables.F)
variables.Tabtilde(variables.pos, variables.F) == math.exp(variables.i * variables.Htilde/(3**(1/2))) 
variables.Tabtilde(variables.neg, variables.F) == -variables.i * variables.Htilde/(3**(1/2))
variables.Jbtilde(1/2, 1/2) == math.exp(variables.brackets(1/2 * (-variables.Phitilde + variables.i * variables.Hbtilde(0) + variables.i * variables.Hbtilde(1) + 3 **(1/2) * variables.i * variables.Htilde)))
variables.Jbtilde(variables.GSO) == 1/2 * diff(-variables.Phitilde + variables.i * variables.Hbtilde(0) + variables.i * variables.Hbtilde(1) + 3**(1/2) * variables.i * variables.Htilde, variables.zl)

# 2.363.1, breaking
(variables.Phib(2), variables.Phib(3), variables.Phib(4)) == (0, 0, 1/2 * math.pi)
(diff(diff(variables.K, variables.tl), variables.t))**(-1) * Abs(diff(variables.W, variables.t) + variables.kappaa(2) * diff(variables.K, variables.t) * variables.W)**2 == 3 * variables.kappaa(2) * Abs(variables.W)**2
variables.Wf(variables.Phi) == diff(variables.Wf(variables.Phi), variables.i) == variables.Daf(variables.Phi, variables.star(variables.Phi)) == 0
variables.V == variables.Ref(variables.brackets(variables.S/(8 * math.pi**2) + variables.fbf(1, variables.T))) * (variables.Da(2))/2
variables.D == 1/(variables.Ref(variables.brackets((variables.S/(8 * math.pi**2)) + variables.fbf(1, variables.T)))) * (2 * variables.zeta - variables.i * variables.kappaa(2) * variables.Kb(variables.comma(variables.i)) * (variables.delta * variables.Phia(variables.i))/(variables.delta * variables.lamda))
variables.Sa(-variables.L - 1)
variables.delta * variables.S == variables.i * variables.q * variables.delta * variables.lamda 
variables.V == (variables.qa(2))/(variables.S + variables.star(variables.S))**3
variables.delta * 1/(4 * math.pi**2) * integrate(WedgeProduct(variables.Imf(variables.S) * variables.Fab(variables.a, 2), variables.Fab(variables.a, 2))) == (variables.q * variables.delta * variables.lamda)/(4 * math.pi**2) * integrate(WedgeProduct(variables.Fab(variables.a, 2), variables.Fab(variables.a, 2)))
variables.D == variables.q/(variables.S + variables.star(variables.S)) + summation(variables.qb(variables.i) * variables.star(variables.Phia(variables.i)) * variables.Phia(variables.i), (variables.Phia(variables.i) != variables.S))
# ommitted exampple section 2.366.(18.8)

#2.402.1
variables.holderb((variables.ha(1, 1), variables.ha(2, 1)), variables.MO) == variables.holderb((variables.ha(2, 1), variables.ha(1, 1)), variables.WO)
math.exp(variables.bracketc(2 * math.pi * variables.i * variables.brackets(variables.r * (variables.Qb(2) - variables.Qb(3)) + variables.s * (variables.Qb(3) - variables.Qb(4)) + variables.t * (variables.Qb(4) - variables.Qb(5)))))
variables.X == variables.X + 2 * math.pi* (variables.SlopeRegge/variables.n)**(1/2) 
variables.Rp == (variables.SlopeRegge/variables.n)**(1/2)
(variables.zb(1), variables.zb(2), variables.zb(3), variables.zb(4), variables.zb(5)) == (variables.zb(1), variables.alphaa(variables.r) * variables.zb(2), variables.alphaa(variables.s - variables.r) * variables.zb(3), variables.alphaa(variables.t - variables.s) * variables.zb(4), variables.alphaa(-variables.t) * variables.zb(5))
variables.WO == variables.MO/variables.Rho 
variables.Gf(variables.z) == variables.zab(5, 1) + variables.zab(5, 2) + variables.zab(5, 3) + variables.zab(5, 4) + variables.zab(5, 5) - 5 * variables.psi * variables.zb(1) * variables.zb(2) * variables.zb(3) * variables.zb(4) * variables.zb(5)
variables.Za(variables.I) == integrate(variables.omegab(3, 0), variables.Aa(variables.i))
variables.gObf(variables.I, variables.Z) == integrate(variables.omegab(3, 0), variables.Bb(variables.I))
variables.F == (variables.xa(0))**2 * variables.brackets(5 * variables.i/6 * variables.Ta(3) - 25 * variables.i /(2 * math.pi**3) * variables.zetaf(3) + summation(variables.Cb(variables.k) * math.exp(-2 * math.pi * variables.k * variables.T), (variables.k, 1, math.inf)))
variables.nb(variables.k) == 2875, 609250, 317206375, 242467530000
variables.Ref(variables.Ta(variables.A)) == integrate(variables.Jb(1, 1), variables.Na(variables.A)) == integrate(variables.d2(variables.w) * variables.gb(variables.i * variables.jl) * (diff(variables.xa(variables.i), variables.z))/(diff(variables.w, variables.z)) * (diff(variables.xa(variables.jl), variables.z))/(diff(variables.wl, variables.z)), variables.Na(variables.A)) > 0
np.dot(variables.star(variables.Phi), variables.Phi) - np.dot(variables.star(variables.rho), variables.rho) - variables.r == 0 
(variables.Phib(variables.i), variables.rhob(variables.i)) == (variables.ea(variables.i * variables.lamda) * variables.Phib(variables.i), variables.ea(-variables.i * variables.lamda) * variables.rhob(variables.i))
# new thing 2.409.1
variables.zab(5, variables.i) == variables.psi * variables.zb(1) * variables.zb(2) * variables.zb(3) * variables.zb(4) * variables.zb(5) 
variables.i == 1, 2, 3, 4, 5
variables.psia(5) == 1
summation(variables.wab(2, variables.i), variables.i) == 0
summation(Abs(variables.wab(2, variables.i), variables.i)) == 2 * variables.rhoa(2)
np.dot(variables.x, variables.x) == variables.rhoa(2) 
np.dot(variables.y, variables.y) == variables.rhoa(2)
np.dot(variables.x, variables.y) == 0
np.cross(variables.Sa(3), variables.Sa(2), variables.RaB(variables.pos))
summation(variables.wab(2, variables.i), variables.i) == variables.psi - 1
np.dot(variables.x, variables.x) == Abs(variables.psi - 1) 
variables.y == 0
# check 2.410.9
variables.gOb(1) == variables.gOb(1) + variables.Za(1)
variables.cbf(variables.mu * variables.n * variables.p * variables.q, (variables.x, variables.y)) == variables.cabf(1, variables.mu, variables.x) * variables.omegabfs((1, variables.n * variables.p * variables.q), variables.y)
integrate(variables.cb(4), variables.D) == integrate(variables.cab(1, 1), variables.P)
-1/(32 * math.pi**2) * ln((variables.Lamdaa(2))/(variables.Ma(2))) * variables.Fb(variables.mu * variables.v) * variables.Fa(variables.mu * variables.v)
1/(8 * math.pi) * variables.Ref(variables.i * diff(variables.gOb(1), 1))
variables.dagger(variables.PhiBb(variables.alpha)) * variables.sigmaab(variables.A, variables.alpha * variables.Beta) * variables.PhiBb(variables.Beta) == 0
variables.A == 1, 2, 3
variables.za(1) * variables.Hbf(1, variables.z) + variables.za(2) * variables.Hbf(2, variables.z) == 0
variables.za(1) == variables.za(2) == variables.Hbf(1, variables.z) == variables.Hbf(2, variables.z) == 0
variables.qab(variables.I, variables.i) == variables.deltaab(variables.I, variables.i)
variables.i == 1, 2, ..., 15 
variables.qab(variables.I, 16) == - 1
variables.dagger(variables.PhiBb(variables.i * variables.alpha)) * variables.sigmaab(variables.A, variables.alpha * variables.Beta) * variables.PhiBb(variables.i * variables.Beta) - variables.dagger(variables.PhiBb(16 * variables.alpha)) * variables.sigmaab(variables.A, variables.alpha * variables.Beta) * variables.PhiBb(16 * variables.Beta) == 0
variables.A == 1, 2, 3
variables.i == 1, 2, ..., 15
variables.PhiBb(variables.i * variables.alpha) == variables.PhiBb(16 * variables.alpha) 
variables.i == 1, 2, ..., 15 
(variables.ha(1, 1), variables.ha(2, 1)) == (2, 86)
variables.B16 == (variables.B4, variables.B2) + (variables.B4p, variables.B2p) 
variables.B16p == (variables.B4, variables.B2p) + (variables.B4p, variables.B2)
# check hodge diamon 2.416.2
variables.star(variables.omegabs(2)) == variables.posneg * variables.omegabs(2)
variables.gb(variables.i * variables.j) == variables.omegaab(variables.kl, (variables.i, variables.j)) * variables.omegabs((variables.i, variables.j)* variables.kl) # check 2.416.4
variables.Hb(variables.mu * variables.v * variables.sigma * variables.p * variables.q) == variables.Hb(variables.mu * variables.v * variables.sigma) * variables.omegabs(variables.p * variables.q)
variables.star(10) == variables.star(4) * variables.star(6)
(variables.S * variables.Of(20, 4, variables.RB))/(np.cross(variables.S * variables.Of(20, variables.RB), variables.S * variables.Of(4, variables.RB)))
variables.Tb(variables.F) == variables.i * variables.psib(variables.m) * diff(variables.xa(variables.m), variables.z) == variables.i * variables.psib(variables.r) * variables.eab(variables.r, variables.m) * diff(variables.xa(variables.m), variables.z)
(variables.B56, variables.B1)**(10) + (variables.B1, variables.B1)**(65)
(variables.B28, variables.B2)**(10) + (variables.B1, variables.B1)**(65)
variables.nb(variables.H) + 29 * variables.nb(variables.T) - variables.nb(variables.V) == 273 
variables.F == variables.starf(variables.F)
integrate(variables.trf(WedgeProduct(variables.Rb(2), variables.Rb(2))), variables.K3) == integrate(variables.Trbf(variables.v, WedgeProduct(variables.Fb(2), variables.Fb(2))), variables.K3)
variables.nb(1) + variables.nb(2) + variables.nb(5) == 24 
# final thing 2.422.1
variables.Gb(variables.I * variables.mu * variables.v) == variables.gab(-1, variables.h) * variables.Gb(variables.h * variables.mu * variables.v)
variables.gb(variables.I) == variables.gab(-1, variables.h)
variables.Rb(variables.m * variables.I) == variables.gab(-1/2, variables.h) * variables.Rb(variables.m * variables.h)
variables.gp == variables.Vab(-1, 1) * variables.gb(variables.I) == variables.Vab(-1, variables.h) * variables.gab((variables.k - 2)/2, variables.h)
variables.Rbp(variables.m) == variables.Rab(-1, variables.m * variables.I) == variables.gab(1/2, variables.h) * variables.Rab(-1, variables.m * variables.h)
(variables.Ta(variables.k))/(variables.ZBb(2)) 
variables.ZBb(2) == variables.bracketc(1, variables.omega * variables.hat(variables.Beta))
variables.Rb(10, variables.M) == variables.gpa(2/3) == variables.gab(1/3, variables.h) * variables.Vab(-2/3, variables.h)
variables.Rb(variables.m * variables.M) == variables.gpa(-1/3) * variables.Rbp(variables.m) == variables.gab(1/3, variables.h) * variables.Rab(-1, variables.m * variables.h) * variables.Vab(1/3, variables.h)
variables.gpp == variables.gpa(-1) == variables.gab(-1, variables.h) * variables.Vb(variables.h)
variables.Rppb == variables.gpa(-1/2) * variables.Rpb(variables.m) == variables.Rab(-1, variables.m * variables.h) * variables.Vab(1/2, variables.h)
variables.omega * variables.hat(variables.Beta) == variables.Betaab(-1, variables.R) * variables.omega * variables.Betab(variables.R) == variables.omega * variables.Betaab(-1, variables.L) * variables.Betab(variables.R)
variables.hat(variables.Beta) == variables.Betaab(-1, variables.L) * variables.Betab(variables.R) == variables.Betaab(-2, variables.L) * variables.Betab(variables.L) * variables.Betab(variables.R) == math.exp(-2 * math.pi * variables.i * variables.Jb(variables.L)) * variables.Beta == math.exp(math.pi * variables.i * variables.n * variables.FBb(variables.L)) * variables.Beta 
variables.Gb(variables.mu * variables.v) * variables.pos 
variables.Bb(variables.mu * variables.v) * variables.neg 
variables.PhiB * variables.pos
variables.C * variables.neg 
variables.Cb(variables.mu * variables.v) * variables.pos 
variables.Cb(variables.mu * variables.v, variables.rho * variables.sigma) * variables.neg 
variables.Gb(variables.mu * variables.v) * variables.pos 
variables.Bb(variables.mu * variables.v) * variables.pos 
variables.Phi * variables.pos 
variables.C * variables.neg 
variables.Cb(variables.mu * variables.v) * variables.neg 
variables.Cb(variables.mu * variables.v, variables.rho * variables.sigma) * variables.neg 
(variables.Ta(variables.k))/(variables.ZBb(2)) 
variables.ZBb(2) == variables.bracketc(1, math.exp(math.pi * variables.i * variables.FBb(variables.L)) * variables.Beta) 
variables.gb(variables.A) == variables.gpp * variables.Rppba(9, -1) == variables.gab(-1, variables.h) * variables.Rb(9 * variables.h) * variables.Vab(1/2, variables.h)
variables.Rb(9, variables.A) == variables.Rppba(9, -1) == variables.Vab(-1/2, variables.h) * variables.Rb(9, variables.h)
variables.Rb(variables.m, variables.A) == variables.Rppb(variables.m) == variables.Vab(1/2, variables.h) * variables.Rab(-1, variables.m, variables.h)
variables.m == 6, 7, 8
(variables.Ta(4))/(variables.ZBb(2)) 
variables.ZBb(2) == variables.bracketc(1, variables.Beta)
variables.ea(-2 * variables.PhiBb(6)) == variables.V * variables.ea(-2 * variables.PhiB)
variables.PhiBb(6) == - variables.PhiBb(6) 
variables.Gb(variables.mu * variables.v) == variables.ea(-2 * variables.PhiBb(6)) * variables.Gb(variables.mu * variables.v) 
variables.Hbtilde(3) == variables.ea(-2 * variables.PhiBb(6)) * variables.starf(6, variables.Hbtilde(3)) 
variables.Fab(variables.a, 2) == variables.Fab(variables.a, 2)
variables.SBb(variables.het) == 1/(2 * variables.kappaab(2, 6)) * integrate(variables.da(6) * variables.x * (-variables.Gb(6))**(1/2) * variables.ea(-2 * variables.PhiBb(6)) * (variables.R + 4 * diff(variables.PhiBb(6), variables.mu) * diff(variables.PhiBb(6), variables.z, variables.mu) - 1/2 * Abs(variables.Hbtilde(3))**2 - (variables.kappaab(2, 6))/(2 * variables.gab(2, 6)) * Abs(variables.Fb(2))**2))
variables.SBb(variables.IIA) == 1/(2 * variables.kappaab(2, 6)) * integrate(variables.da(6) * variables.x * (-variables.Gb(6))**(1/2) * (variables.ea(-2 * variables.PhiBb(6)) * variables.R + 4 * variables.ea(-2 * variables.PhiBb(6)) * diff(variables.PhiBb(6), variables.mu) * diff(variables.PhiBb(6), variables.z, variables.mu) - 1/2 * Abs(variables.Hbtilde(3))**2 - (variables.kappaab(2, 6))/(2 * variables.gab(2, 6)) * variables.ea(-2 * variables.PhiBb(6)) * Abs(variables.Fb(2))**2))
(variables.Ta(4))/(variables.ZBb(2)) == variables.K3 
np.cross((variables.S * variables.Of(19, 3, variables.RB))/(np.cross(variables.S * variables.Of(19, variables.RB), variables.S * variables.Of(3, variables.RB))), variables.RBa(variables.pos))
variables.ea(-2 * variables.PhiBb(4)) == variables.Rb(4) * variables.Rb(5) * variables.ea(-2 * variables.PhiBb(6))
variables.ea(variables.PhiBb(4)) == (variables.Rb(4) * variables.Rb(5))**(-1/2) 
variables.starf(4, variables.Htilde) == variables.ea(-2 * variables.PhiBb(4)) * diff(variables.Bb(45))
variables.starf(4, variables.Htilde) == variables.ea(-2 * variables.PhiBb(4)) * diff(variables.a)
variables.S == variables.i * variables.star(variables.rho)
np.cross((variables.S * variables.Uf(1, 1))/(np.cross(variables.Uf(1), variables.S * variables.Uf(1, 1, variables.ZB))), (variables.Of(22, 6, variables.RB))/(np.cross(variables.Of(22, variables.RB), variables.Of(6, variables.RB), variables.Of(22, 6, variables.ZB))))
variables.dirac(1, -1) 
variables.dirac(1, -1/2)**2 
variables.dirac(1, 0)
