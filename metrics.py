import math
from scipy.integrate import quad
import numpy as np 
from sympy import symbols, diff,integrate, Function, Derivative, summation, oo, ln, Abs
import variables
import constants

variables.aa(variables.mu) * variables.bb(variables.mu) == -variables.aapos*variables.baneg - variables.aaneg * variables.bapos + variables.a**variables.i * variables.b**variables.i
variables.abneg == -variables.aapos
abpos = -variables.aaneg
variables.ab(variables.i) == variables.aa(variables.i)
if np.array([variables.GWorldSheetMetrica(variables.tau, variables.tau), variables.GWorldSheetMetrica(variables.tau, variables.sigma)], [variables.GWorldSheetMetrica(variables.sigma, variables.tau), variables.GWorldSheetMetrica(variables.sigma, variables.sigma)]) == np.array([-variables.GWorldSheetMetricb(variables.sigma, variables.sigma)(variables.tau), variables.GWorldSheetMetricb(variables.tau, variables.sigma)(variables.tau, variables.sigma)], [variables.GWorldSheetMetricb(variables.tau, variables.sigma)(variables.tau, variables.sigma), variables.GWorldSheetMetricb(variables.tau, variables.sigma)**(-1)(variables.tau)(1-(variables.GWorldSheetMetricb(variables.tau, variables.sigma))**2(variables.tau, variables.sigma))]):
    L=-1/(4*math.pi*variables.SlopeRegge)*quad(np.diff(variables.sigma), 0 ,variables.l) * (variables.GWorldSheetMetricb(variables.sigma, variables.sigma)((2*diff(variables.xas(variables.neg), variables.tau)) - diff(variables.xa(variables.i), variables.tau)*diff(variables.xa(variables.i), variables.tau)) - 2*variables.GWorldSheetMetricb(variables.sigma, variables.tau)*(diff(variables.ya(variables.neg), variables.tau) - diff(variables.xa(variables.i), variables.tau) * diff(variables.xa(variables.i), variables.sigma)) + (variables.GWorldSheetMetricb*(variables.sigma, variables.sigma)**(-1)(1-variables.GWorldSheetMetricb(variables.tau, variables.sigma)**2)*diff(variables.xa(variables.i), variables.sigma)*diff(variables.xa(variables.i), variables.sigma)))
    variables.xas(variables.neg, variables.tau) == 1/variables.l * quad(np.diff(variables.sigma)*variables.xa(variables.neg, variables.tau, variables.sigma), 0, variables.l)
    variables.ya(variables.neg, variables.tau, variables.sigma) == variables.xa(variables.neg, variables.tau, variables. sigma) - variables.xas(variables.neg, variables.tau) #transformation
if variables.sigma == (0, variables.l):
    variables.GWorldSheetMetricb(variables.tau, variables.sigma)*(diff(variables.xa(variables.mu), variables.tau))-variables.GWorldSheetMetric(variables.tau, variables.tau) *(diff(variables.xa(variables.mu), variables.sigma)) == 0 # open string boundary condition
    if variables.mu == variables.pos:
        variables.GWorldSheetMetricb(variables.tau, variables.sigma) # open string boundary condition
    if variables.mu == variables.i:
        diff(variables.xa(variables.i), variables.sigma) == 0 # boundry condition
    L == -variables.l/(2*math.pi*variables.SlopeRegge)*variables.GWorldSheetMetricb(variables.sigma, variables.sigma) * diff(variables.xas(variables.neg), variables.tau) + 1/(4*math.pi*variables.SlopeRegge)*quad(np.diff(variables.sigma)*(variables.GWorldSheetMetric(variables.sigma, variables.sigma)*diff(variables.xa(variables.i), variables.tau)*diff(variables.xa(variables.i), variables.tau) - variables.GWorldSheetMetric(variables.sigma, variables.sigma)**(-1)*diff(variables.xa(variables.i), variables.sigma)*diff(variables.xa(variables.i), variables.sigma)))
    variables.pb(variables.neg) == -variables.pa(variables.pos) == diff(variables.L)/(diff(diff(variables.xas(variables.neg), variables.tau))) == -variables.l/(2*math.pi*variables.SlopeRegge)*variables.GWorldSheetMetric(variables.sigma, variables.sigma)

integrate(variables.brackets(diff(variables.X) * diff(variables.g)) * math.exp(-variables.S)) # definition for spacetime 82.1
variables.S == variables.Sb(variables.X) + variables.lamdab(variables.chi) 
variables.Sb(variables.X) == 1/(4*math.pi*variables.SlopeRegge) * quad(variables.d2sigma*variables.g**(1/2) * variables.ga(variables.a, variables.b) * diff(variables.xa(variables.mu), variables.a) * diff(variables.xb(variables.mu), variables.b), variables.M)
variables.chi == 1/(4*math.pi) * quad(variables.d2sigma * variables.g**(1/2) * variables.R + 1/(2*math.pi) * integrate(diff(variables.s) * variables.k, 0, 2*math.pi), variables.M) # gepdisoc integration for second

#redefinition 85.8
variables.dsp2 == math.exp(2*variables.omegas) * Abs(diff(variables.f, variables.z))**(-2) * diff(variables.zp) * diff(variables.zpl)
variables.omegas == ln(Abs(diff(variables.f, variables.z)))

#92.8
variables.Tab(variables.a, variables.a) == variables.ab(variables.D1) *variables.R

#231.1, periodic metric:
variables.xas(4) == variables.xas(4) + 2 * math.pi * variables.R 
variables.d * variables.sa(2) == variables.Gab(variables.D, variables.M * variables.N) * variables.d * variables.xas(variables.M) * variables.d * variables.xas(variables.N) == variables.Gb(variables.mu * variables.v) * variables.d * variables.xas(variables.mu) * variables.d * variables.xas(variables.v) + variables.Gbf(variables.d * variables.d, variables.d * variables.xas(variables.d) + variables.Ab(variables.mu) *  variables.d * variables.xas(variables.mu))**2
variables.xaps(variables.d) == variables.xas(variables.d) + variables.lamdaf(variables.xas(variables.mu))
variables.Abp(variables.mu) == variables.Ab(variables.mu) - diff(variables.lamda, variables.mu)
variables.phif(variables.xas(variables.M)) == summation(variables.phibf(variables.n, variables.xas(variables.mu)) * math.exp(variables.i * variables.n * variables.xas(variables.d)/variables.R), (variables.n, -math.inf, math.inf))
diff(diff(variables.phibf(variables.n, variables.xas(variables.mu)), 1, variables.mu), variables.mu) == variables.na(2)/(variables.Ra(2)) * variables.phibf(variables.n, variables.xas(variables.mu))
-variables.pa(variables.mu) * variables.pb(variables.mu) == variables.na(2)/(variables.Ra(2))
# ricci scalar for metric 232.8
variables.RB == variables.RBb(variables.d) - 2 * variables.e**(-variables.sigma) * variables.Deltaa(2) * variables.ea(variables.sigma) - 1/4 * variables.ea(2 * variables.sigma) * variables.Fb(variables.mu * variables.v) * variables.Fa(variables.mu * variables.v)
# graviton dilation action
variables.Sb(1) == 1/(2 * variables.ketaab(2, 0)) * integrate(variables.da(variables.D) * variables.x(-variables.Gb(variables.D)**(1/2)) * variables.ea(-2 * variables.Phi) * (variables.Rb + (4 * variables.Deltab(variables.mu) * variables.Phi * variables.Deltaa(variables.mu) * variables.Phi))) == np.cross(((math.pi * variables.R)/(variables.ketaab(2, 0)) * integrate(variables.da(variables.d) * variables.x * (-variables.Gb(variables.d))**(1/2) * variables.ea(-2 * variables.Phi + variables.sigma))),(variables.RBb(variables.d) - 4 * diff(variables.Phi, variables.mu) * diff(variables.sigma, 1, variables.mu) + 4 * diff(variables.Phi, variables.mu) * diff(variables.Phi, 1, variables.mu) - 1/4 * variables.ea(2 * variables.sigma) * variables.Fb(variables.mu * variables.v) * variables.Fa(variables.mu * variables.v))) == np.cross((math.pi * variables.R/(variables.ketaab(2, 0)) * integrate(variables.da(variables.d) * variables.x * (-variables.G)**(1/2) * variables.ea(-2 * variables.Phib(variables.d)))),(variables.RBb(variables.d) - diff(variables.sigma, variables.mu) * diff(variables.sigma, 1, variables.mu) + 4 * diff(variables.Phib(variables.d), variables.mu) * diff(variables.Phib(variables.d), 1, variables.mu) - 1/4 * variables.ea(2 * variables.sigma) * variables.Fb(variables.mu * variables.v) * variables.Fa(variables.mu * variables.v)))
diff(1, variables.mu) + variables.i * variables.pb(variables.d) * variables.Ab(variables.mu) == diff(1, variables.mu) + variables.i * variables.n * variables.Abtilde(variables.mu)
variables.gab(2, variables.d) == (variables.ketaab(2, 0) * variables.ea(2 * variables.Phib(variables.d)))/(math.pi * variables.Ra(3) * variables.ea(2 * variables.sigma)) == (2 * variables.ketaab(2, variables.d))/(variables.rhoa(2))
1/(variables.ketaab(2, variables.d)) == (2 * math.pi * variables.rho)/(variables.ketaa(2))
variables.Sb(2) == -1/(24 * variables.ketaab(2, 0)) * integrate(variables.da(variables.D) * variables.x * (-variables.Gb(variables.D))**(1/2) * variables.ea(-2 * variables.Phi) * variables.Hb(variables.M * variables.N * variables.L) * variables.Ha(variables.M * variables.N * variables.L)) == -(math.pi * variables.R)/(12 * variables.ketaab(2, 0)) * integrate(variables.da(variables.d) * variables.x * (-variables.Gb(variables.d))**(1/2) * variables.ea(-variables.Phib(variables.d)) * (variables.Hbtilde(variables.mu * variables.v * variables.lamda) * variables.Hatilde(variables.mu * variables.v * variables.lamda) + 3 * variables.ea(-2 * variables.sigma) * variables.Hb(variables.d * variables.mu * variables.v) * variables.Hab(variables.mu * variables.v, variables.d)))
# check for Hbtilde definition 234.14
variables.Bbp(variables.v * variables.lamda) == variables.Bb(variables.v * variables.lamda) - variables.lamda * variables.Hb(variables.d * variables.v * variables.lamda)
#for CFT 235.1, field =
variables.X == variables.X + 2 * math.pi * variables.R
variables.k =variables.n/variables.R 
variables.N == variables.ZB 
variables.Xf(variables.sigma + 2 * math.pi) == variables.Xf(variables.sigma) + 2 * math.pi * variables.R * variables.w 
variables.w == variables.ZB 
diff(variables.Xf(variables.z), variables.z) == -variables.i * (variables.SlopeRegge/2)**(1/2) * summation((variables.alphab(variables.m))/(variables.za(variables.m + 1)), (variables.m, -math.inf, math.inf))
diff(variables.Xf(variables.zl), variables.zl) == -variables.i * (variables.SlopeRegge/2)**(1/2) * summation((variables.alphabtilde(variables.m))/(variables.zla(variables.m + 1)), (variables.m, -math.inf, math.inf))
2 * math.pi * variables.R * variables.w == integrate(diff(variables.z) * diff(variables.X, variables.z) + diff(variables.zl) * diff(variables.X, variables.zl), math.pi * 2) == 2 * math.pi * (variables.SlopeRegge/2)**(1/2) * (variables.alphab(0) - variables.alphabtilde(0)) #  contour integral
variables.p == 1/(2 * math.pi * variables.SlopeRegge) * integrate(diff(variables.z) * diff(variables.X, variables.z) - diff(variables.zl) * diff(variables.X, variables.zl), 2 * math.pi) == (2 * variables.SlopeRegge)**(-1/2) * (variables.alphab(0) + variables.alphabtilde(0)) # contour integral
# for periodic dimension 236.7
variables.pb(variables.L) == (2/variables.SlopeRegge)**(1/2) * variables.alphab(0) == variables.n/variables.R + (variables.w * variables.R)/(variables.SlopeRegge)
variables.pb(variables.R) == (2/variables.SlopeRegge)**(1/2) * variables.alphabtilde(0) == variables.n/variables.r - (variables.w * variables.R)/(variables.SlopeRegge)
variables.Lb(0) == (variables.SlopeRegge * variables.pab(2, variables.l))/(4) + summation(variables.alphab(-variables.n) * variables.alphab(variables.n), (variables.n, 1, math.inf))
variables.Lbtilde(0) == (variables.SlopeRegge * variables.pab(2, variables.R))/4 + summation(variables.alphabtilde(-variables.n) * variables.alphabtilde(variables.n), (variables.n, 1, math.inf))
#partition function 237.9
(variables.q * variables.ql)**(-1/24) * variables.Trf(variables.qa(variables.Lb(0)) * variables.qla(variables.Lbtilde(0))) == Abs(variables.etaf(variables.tau))**(-2) * summation(variables.qa(variables.SlopeRegge * variables.pab(2, variables.L/4)) * variables.qla(variables.SlopeRegge * variables.pab(2, variables.R)/4), ((variables.n, variables.w), -math.inf, math.inf)) == Abs(variables.etaf(variables.tau))**(-2) * summation(math.exp(variables.brackets(-math.pi * variables.taub(2) * ((variables.SlopeRegge * variables.na(2))/(variables.Ra(2)) + (variables.wa(2) * variables.Ra(2))/(variables.SlopeRegge)) + 2 * math.pi * variables.i * variables.taub(1) * variables.n * variables.w)), ((variables.n, variables.w), -math.inf, math.inf))
summation(math.exp(-math.pi * variables.a * variables.na(2) + 2 * math.pi * variables.i * variables.b * variables.n), (variables.n, -math.inf, math.inf)) == variables.aa(-1/2) * summation(variables.brackets((math.pi * (variables.m - variables.b)**2)/variables.a), (variables.m, -math.inf, math.inf))
# new partition function
2 * math.pi * variables.R * variables.Zbf(variables.X, variables.tau) * summation(math.exp(-(math.pi * variables.Ra(2) * Abs(variables.m - variables.w * variables.tau)**2)/(variables.SlopeRegge * variables.taub(2))), ((variables.m, variables.w), -math.inf, math.inf))
variables.Xf(variables.sigmaa(1) + 2 * math.pi, variables.sigmaa(2)) == variables.Xf(variables.sigmaa(1), variables.sigmaa(2)) + 2 * math.pi * variables.w * variables.R 
variables.Xf(variables.sigmaa(1) + 2 * math.pi * variables.taub(1), variables.sigmaa(2) + 2 * math.pi * variables.taub(2)) == variables.Xf(variables.sigmaa(1), variables.sigmaa(2)) + 2 * math.pi * variables.m * variables.R 
variables.xb(variables.c * variables.l) == variables.sigmaa(1) * variables.w * variables.R + variables.sigmaa(2) * (variables.m - variables.w * variables.taub(1)) * variables.R/(variables.taub(2))
variables.brackets(variables.xbs(variables.L), variables.pb(variables.L)) == variables.brackets(variables.xbs(variables.R), variables.pb(variables.R)) == variables.i 
variables.Xf(variables.z, variables.zl) == variables.Xbf(variables.L, variables.z) + variables.Xbf(variables.R, variables.zl)
variables.Xbf(variables.L, variables.z) == variables.xbs(variables.L) - variables.i * (variables.SlopeRegge)/(2) * variables.pb(variables.L) * ln(variables.z) + variables.i * (variables.SlopeRegge/2)**(1/2) * (summation((variables.alphab(variables.m)/(variables.m * variables.za(variables.m), variables.m, -math.inf, -1))) + summation((variables.alphab(variables.m)/(variables.m * variables.za(variables.m), variables.m, 1, math.inf))))
variables.Xbf(variables.R, variables.zl) == variables.xbs(variables.R) - variables.i * variables.SlopeRegge/2 * variables.pb(variables.R) * ln(variables.zl) + variables.i * (variables.SlopeRegge/2)**(1/2) * (summation((variables.alphabtilde(variables.m)/(variables.m * variables.zla(variables.m))), (variables.m, -math.inf -1)) + summation((variables.alphabtilde(variables.m)/(variables.m * variables.zla(variables.m))), (variables.m, 1, math.inf)))
# check antiholoimorphic approximations 238.17
variables.VObf(variables.kb(variables.L) * variables.kb(variables.R), (variables.z, variables.zl)) == variables.point(variables.ea(variables.i * variables.kb(variables.L) * variables.Xbf(variables.L, variables.z) + variables.i * variables.kb(variables.R) * variables.Xbf(variables.R, variables.zl)))
# OPE approx. 238.19
math.exp(math.pi * variables.i * variables.SlopeRegge * (variables.kb(variables.L) * variables.kbp(variables.L) - variables.kb(variables.R) * variables.kbp(variables.R))) == math.exp(2 * math.pi * variables.i * (variables.n * variables.wp + variables.w * variables.nprime)) == 1
variables.brackets(variables.Xbf(variables.L, variables.zb(1)), variables.Xbf(variables.L, variables.zb(2))) == (math.pi * variables.i * variables.SlopeRegge)/2 * variables.sign(variables.sigmaab(1, 1) - variables.sigmaab(1, 2))
variables.VObf(variables.kb(variables.L) * variables.kb(variables.R), (variables.z, variables.zl)) == math.exp(variables.brackets(math.pi * variables.i * (variables.kb(variables.L) - variables.kb(variables.R)) * (variables.pb(variables.L) + variables.pb(variables.R) * variables.SlopeRegge/4))) * variables.point(variables.ea(variables.i * variables.kb(variables.L) * variables.Xbf(variables.L, variables.z) + variables.i * variables.kb(variables.R) * variables.Xbf(variables.R, variables.zl)))
math.exp(variables.bracketc(math.pi * variables.i * variables.brackets((variables.kb(variables.L) - variables.kb(variables.R)) * (variables.kbp(variables.L) + variables.kbp(variables.R)) - (variables.kbp(variables.L) - variables.kbp(variables.R)) * (variables.kb(variables.L) + variables.kb(variables.R))) * variables.SlopeRegge/4)) == math.exp(variables.brackets(math.pi * variables.i * (variables.n * variables.wp - variables.w * variables.nprime)))
# check replacements for expansion 239.24/5

