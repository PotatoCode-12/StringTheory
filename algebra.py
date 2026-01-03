import math
from scipy.integrate import quad
import numpy as np 
from sympy import symbols, diff,integrate, Function, Derivative, summation, oo, Abs
import sympy as sp
import variables
import constants

m = symbols('n') # for summation
- math.inf < variables.sigmaa(variables.D2) < math.inf #check conditions pg 52
variables.w == variables.simgaa(variables.D1) + variables.i*variables.sigmaa(variables.D2)
variables.z == math.exp(-variables.i*variables.w) == math.exp(-variables.i*variables.sigmaa(variables.D1) + variables.sigmaa(variables.D2))
variables.Tbf(variables.z, variables.z, variables.z) == summation((variables.Lb(m))/(variables.za(m+2)), m, -math.inf, math.inf) # laurent expansions
variables.Tbftilde(variables.zl, variables.zl, variables.zl) == summation((variables.Lbtilde(m))/(variables.zla(m+2)))
variables.Lb(m) == quad((diff(variables.z))/(2*math.pi*variables.i*variables.z) * variables.za(variables.m+2) * variables.Tbf(variables.z, variables.z, variables.z), 0, 2*math.pi) # change z systematically to fit countour bounds
variables.Tbf(variables.w, variables.w, variables.w) == -summation(math.exp(variables.i*m*variables.sigmaa(variables.D2) - m*variables.sigmaa(variables.D2)) * variables.Tb(m), m, -math.inf, math.inf) # fourier transformation
variables.Tbf(variables.wl, variables.wl, variables.wl) =-summation(math.exp(-variables.i*m*variables.sigmaa(variables.D2) - m*variables.sigmaa(variables.D2)) * variables.Tbtilde(m), m, -math.inf, math.inf)
variables.Tb(m) == variables.Lb(m) - variables.deltab(m, 0) * variables.c/24
variables.Tbtilde(m) == variables.Lbtilde(m) - variables.deltab(m, 0) * variables.ctilde/24
variables.Tb(variables.w, variables.w) == (diff(variables.z, variables.w))**(variables.D2)*variables.Tb(variables.z, variables.z) + variables.c/24 # additive shift
variables.H == quad((diff(variables.sigmaa(variables.D1)))/(math.pi*2) * variables.Tb(variables.D2, variables.D2), 0, math.pi*2) == variables.Lb(0) + variables.Lbtilde(0) - (variables.c +variables.ctilde)/24 # 54
variables.Qb(variables.i, variables.C) == quad((diff(variables.z))/(math.pi*2*variables.i) * variables.jb(variables.i), 0, math.pi*2) # contour integral 54
variables.Qm1*variables.Qm2 - variables.Qm2*variables.Qm1 == np.array([variables.Qm1, variables.Qm2]) # matric element
variables.Qftilde(variables.C) == -quad((diff(variables.zl))/(2*math.pi*variables.i) * variables.jtilde, 0, 2*math.pi) #contour
variables.brackets(variables.Lb(m), variables.Lb(variables.n)) == (m - variables.n) * variables.Lb(m + variables.n) + variables.c/12 * (m**3-m) * variables.deltab(m, -variables.n) # Virasaro algebra
variables.brackets(variables.Lb(variables.D0), variables.Lb(variables.n)) == -variables.n*variables.Lb(variables.n) # generator condition
variables.Lb(variables.D0) * variables.Lb(variables.n) * variables.dirac(variables.psi) == variables.Lb(variables.n) * (variables.Lb(variables.D0) - variables.n) * variables.dirac(variables.psi) == (variables.h - variables.n) * variables.Lb(variables.n) * variables.dirac(variables.psi) # check eigenstates value condition 56.21
variables.brackets(variables.Lb(variables.D0), variables.Lb(variables.D1)) == -variables.Lb(variables.D1) # generator closed algebra
variables.brackets(variables.Lb(variables.D0, -variables.D1)) == variables.Lb(-variables.D1)
variables.brackets(variables.Lb(variables.D1), variables.Lb(-variables.D1)) == 2*variables.Lb(variables.D0)
variables.OO(variables.z) == summation((variables.OOb(m))/(variables.za(m == variables.h)), (m, -math.inf, math.inf)) # weight (h, 0)
variables.brackets(variables.Lb(m), variables.OOb(variables.n)) == variables.brackets((variables.h - 1) * m - variables.n) * variables.OOb(m + variables.n)
variables.Tb(variables.a, variables.b) * variables.na(variables.a) * variables.ta(variables.b) == 0 # tensor conditions
variables.Tbf(variables.z, variables.z, variables.z) == variables.Tbf(variables.zl, variables.zl, variables.zpl)
variables.i*variables.z < 0 # redefine for all imaginary conditions pg 57
variables.Lb(m) == 1/(2*math.pi*variables.i) * quad(diff(variables.z) * variables.za(m+1) * variables.tb(variables.z, variables.z) - diff(variables.zl) * variables.zla(m+1) * variables.Tb(variables.zl, variables.zl), variables.C) == 1/(2*math.pi*variables.i) * quad(diff(variables.z) * variables.za(m + 1) * variables.Tbf(variables.z, variables.z, variables.z), 0, 2*math.pi) # contour
variables.dirac(variables.Lb(m), variables.Lb(variables.n)) == (m - variables.n) * variables.Lb(m + variables.n) + variables.c/12 *(m**3-m) * variables.deltab(m, -variables.n)

# under contitions pg 61
variables.Lb(variables.m) == 1/2*summation(variables.ANO(variables.alphaab(variables.mu, variables.m-variables.n) * variables.alphaab(variables.mu, variables.n)) + variables.i*(variables.SlopeRegge/2)**(1/2) * (variables.m + 1) * variables.Va(variables.mu) * variables.alphab(variables.mu, variables.m), (variables.n, -math.inf, math.inf))
variables.Lb(variables.m) == summation((variables.m*variables.lamda - variables.n) * variables.ANO(variables.bb(variables.n)  * variables.cb(variables.m-variables.n)) + variables.deltab(variables.m, 0) * variables.aa(variables.g), (variables.n, -math.inf, math.inf)) # virasaro generators
2*variables.Lb(0, variables.dirac(variables.down)) == (variables.Lb(variables.D1)*variables.Lb(-variables.D1) - variables.Lb(-variables.D1) * variables.Lb(variables.D1)) * variables.dirac(variables.down) == (variables.lamda*variables.bb(variables.D0) * variables.cb(variables.D1)) * variables.brackets((1-variables.lamda) * variables.bb(-variables.D1) * variables.cb(variables.D0)) * variables.dirac(variables.down) == variables.lamda(1-variables.lamda)*variables.dirac(variables.down) # ordering constants #61.20-21
variables.Lb(variables.m) == summation((variables.m*variables.lamda - variables.n)* variables.ANO(variables.bb(variables.n) * variables.cb(variables.m-variables.n)) + (variables.lamda*(1-variables.lamda))/2 * variables.deltab(variables.m, 0), (variables.n, -math.inf, math.inf))
# ghost current
variables.Na(variables.g) == -1/(2*math.pi*variables.i) * quad(variables.diff(variables.w) * variables.jb(variables.w), 0, 2*math.pi) == summation(variables.cb(-variables.n) * variables.bb(variables.n) - variables.bb(-variables.n)*variables.cb(variables.n), (variables.n, 1, math.inf)) + variables.cb(variables.D0) * variables.bb(variables.D0) - 1/2
variables.brackets(variables.Na(variables.g), variables.bb(variables.m)) == -variables.bb(variables.m)
variables.brackets(variables.Na(variables.g), variables.cb(variables.m)) == variables.cb(variables.m)
#with ghost number +-1/2
variables.Na(variables.g)*variables.dirac(variables.down) == -1/2*variables.dirac(variables.down)
variables.Na(variables.g)*variables.dirac(variables.up) == 1/2*variables.dirac(variables.up)

#viasaro 72.14
2*variables.hb(variables.OO) * variables.dirac([variables.OO], [variables.OO]) == 2*variables.dirac([variables.OO], [variables.Lb(variables.D0)], [variables.OO]) == variables.dirac([variables.OO], [variables.brackets(variables.Lb(variables.D1), variables.Lb(-variables.D1))], [variables.OO]) == np.linalg.norm(variables.Lb(-variables.D1) * variables.dirac(variables.OO))**2 >= 0
variables.Lb(-variables.D1) * variables.OO == variables.Lbtilde(-variables.D1) * variables.OO == variables.Lbtilde(-variables.D1) * variables.OO == 0

#2.1.1
variables.pb(variables.mu) * variables.pa(variables.mu) + variables.ma(2) == 0
variables.Lb(0) * variables.dirac(variables.psi) == 0
variables.i * variables.pb(variables.mu) * variables.Rhoa(variables.mu) + variables.m == 0
variables.bracketc(variables.Rhoa(variables.mu), variables.Rhoa(variables.v)) == 2 * variables.etaa(variables.mu *variables.v)
variables.S == 1/(4 * math.pi) * integrate(variables.d2z * (2/variables.SlopeRegge * diff(variables.xa(variables.mu), variables.z) * diff(variables.xb(variables.mu), variables.zl) + variables.psia(variables.mu) * diff(variables.psib(variables.mu), variables.zl) + variables.psiatilde(variables.mu) * diff(variables.psibtilde(variables.mu), variables.z)))
#excluded 2.2.6-7
variables.Tbf(variables.F, variables.z) == variables.i * (2/variables.SlopeRegge)**(1/2) * variables.psiaf(variables.mu, variables.z) * diff(variables.xb(variables.mu, variables.z), variables.z)
variables.Tbftilde(variables.F, variables.zl) == variables.i * (2/variables.SlopeRegge)**(1/2) * variables.psiaftilde(variables.mu, variables.zl) * diff(variables.xb(variables.mu, variables.zl), variables.zl)
variables.jaf(variables.eta, variables.z) * variables.etaf(variables.z) * variables.Tbf(variables.F, variables.z) 
variables.jaftilde(variables.eta, variables.zl) == variables.star(variables.etaf(variables.z)) * variables.Tbftilde(variables.F, variables.zl)
variables.epsilona(-1) * (2/variables.SlopeRegge)**(1/2) * variables.delta * variables.xa(variables.mu, (variables.z, variables.zl)) == + variables.etaf(variables.z) * variables.psiaf(variables.mu, variables.z) + variables.star(variables.etaf(variables.z)) * variables.psiaftilde(variables.mu, variables.zl) 
variables.epsilona(-1) * (variables.SlopeRegge/2)**(1/2) * variables.delta * variables.psiaf(variables.mu, variables.z) == - variables.etaf(variables.z) * diff(variables.xa(variables.mu, variables.z), variables.z)
variables.epsilona(-1) * (variables.SlopeRegge/2)**(1/2) * variables.delta * variables.psiaftilde(variables.mu, variables.zl) == - variables.star(variables.etaf(variables.z)) * diff(variables.xa(variables.mu, variables.zl), variables.zl)
variables.deltab(variables.etab(1)) * variables.deltab(variables.etab(2)) - variables.deltab(variables.etab(2)) * variables.deltab(variables.etab(1)) == variables.deltab(variables.v) 
variables.Tb(variables.B) == - 1/variables.SlopeRegge * diff(variables.xa(variables.mu), variables.z) * diff(variables.xb(variables.mu), variables.z) - 1/2 * variables.psia(variables.mu) * diff(variables.psib(variables.mu), variables.z)
#check approximation 2.3.13 & 4.15
variables.c == (1 + 1/2) * variables.D == 3/2 * variables.D 
variables.hb(variables.b) == variables.lamda 
variables.hb(variables.c) == 1 - variables.lamda 
variables.hb(variables.Beta) == variables.lamda - 1/2 
variables.hb(variables.gamma) == 3/2 - variables.lamda 
variables.Sb(variables.B * variables.C) == 1/(2 * math.pi) * integrate(variables.d2z * (variables.b * diff(variables.c, variables.zl) + variables.Beta * diff(variables.gamma, variables.zl)))
variables.Tb(variables.B) == (diff(variables.b, variables.z)) * variables.c - variables.lamda * diff(variables.b * variables.c, variables.z) + (diff(variables.Beta, variables.z)) * variables.gamma - 1/2 * (2 * variables.lamda - 1) * diff(variables.Beta * variables.gamma, variables.z)
variables.Tb(variables.F) == - 1/2 * (diff(variables.Beta, variables.z)) * variables.c + (2 * variables.lamda - 1)/(2) * diff(variables.b * variables.c, variables.z) - 2 * variables.b * variables.gamma 
variables.brackets(-3 * (2 * variables.lamda -1)**2 + 1) + variables.brackets(3 * (2 * variables.lamda -2)**2 - 1) == 9 - 12 * variables.lamda 
0 == 3/2 * variables.D - 15
variables.D == 10
#pour variables.lamda == 2
variables.Tb(variables.B) == -(diff(variables.b, variables.z)) * variables.c - 2 * variables.b * diff(variables.c, variables.z) - 1/2 * (diff(variables.Beta, variables.z)) * variables.gamma - 3/2 * variables.Beta * diff(variables.gamma, variables.z) 
variables.Tb(variables.F) == (diff(variables.Beta, variables.z)) * variables.c + 3/2 * variables.Beta * diff(variables.c, variables.z) - 2 * variables.b * variables.gamma
variables.Tbf(variables.B, variables.z) == - 1/variables.SlopeRegge * diff(variables.xa(variables.mu), variables.z) * diff(variables.xb(variables.mu), variables.z) + variables.Vb(variables.mu) * diff(variables.xa(variables.mu), variables.z, 2) - 1/2 * variables.psia(variables.mu) * diff(variables.psiab(variables.mu), variables.z)
variables.Tbf(variables.F, variables.z) == variables.i * (2/variables.SlopeRegge)**(1/2) * variables.psia(variables.mu) * diff(variables.xb(variables.mu), variables.z) - variables.i *(2 * variables.SlopeRegge)**(1/2) * variables.Vb(variables.mu) * diff(variables.psia(variables.mu), variables.z)
variables.c == 3/2 * variables.D + 7 * variables.SlopeRegge * variables.Va(variables.mu) * variables.Vb(variables.mu)

#2.5.1
1/(4 * math.pi) * integrate(variables.d2(variables.w) * (variables.psia(variables.mu) * diff(variables.psib(variables.mu), variables.wl) + variables.psiatilde(variables.mu) * diff(variables.psibtilde(variables.mu), variables.w)))
variables.psiaf(variables.mu, variables.w + 2 * math.pi) == + variables.psiaf(variables.mu(variables.w))
variables.psiaf(variables.mu, variables.w + 2 * math.pi) == - variables.psiaf(variables.mu(variables.w))
variables.psiaf(variables.w + 2 * math.pi) == math.exp(2 * math.pi * variables.i * variables.v) * variables.psiaf(variables.mu, variables.w)
variables.psiaftilde(variables.mu, variables.wl + 2 * math.pi) == math.exp(-2 * math.pi * variables.i * variables.vtilde) * variables.psiaftilde(variables.mu, variables.wl)
variables.Tbf(variables.F, variables.w + 2 * math.pi) == math.exp(2 * math.pi * variables.i * variables.v) * variables.Tbf(variables.F, variables.w)
variables.Tbftilde(variables.F, variables.wl + 2 * math.pi) == math.exp(- 2 * math.pi * variables.i * variables.vtilde) * variables.Tbftilde(variables.F, variables.wl)
variables.psiaf(variables.mu, variables.w) == variables.i**(-1/2) * summation(variables.psiab(variables.mu, variables.r) * math.exp(variables.i * variables.r * variables.w), (variables.r, variables.ZB + variables.v))
variables.psiaftilde(variables.mu, variables.wl) == variables.i**(1/2) * summation(variables.psiabtilde(variables.mu, variables.r) * math.exp(-variables.i * variables.r * variables.wl), (variables.r, variables.ZB + variables.vtilde))
variables.psiabf(variables.mu, variables.za(1/2), variables.z) == (diff(variables.w, variables.z))**(1/2) * variables.psiabf(variables.mu, variables.wa(1/2), variables.w) == variables.ia(1/2) * variables.za(-1/2) * variables.psiabf(variables.mu, variables.wa(1/2), variables.w)
variables.psiaf(variables.mu, variables.z) == summation((variables.psiab(variables.mu, variables.r))/(variables.za(variables.r + 1/2)), (variables.r, variables.ZB + variables.v))
variables.psiaftilde(variables.mu, variables.zl) == summation((variables.psiabtilde(variables.mu, variables.r))/(variables.zla(variables.r + 1/2)), (variables.r, variables.ZB + variables.vtilde))
diff(variables.xa(variables.mu, variables.z), variables.z) == - variables.i * (variables.SlopeRegge/2)**(1/2) * summation((variables.alphab(variables.mu, variables.m))/(variables.za(variables.m + 1)), (variables.m, -math.inf, math.inf)) 
diff(variables.xa(variables.mu, variables.zl), variables.zl) == -variables.i * (variables.SlopeRegge/2)**(1/2) * summation((variables.alphaabtilde(variables.mu, variables.m))/(variables.za(variables.m + 1)), (variables.m, -math.inf, math.inf)) 
variables.bracketc(variables.psiab(variables.mu, variables.r), variables.psiab(variables.v, variables.s)) == variables.bracketc(variables.psiabtilde(variables.mu, variables.r), variables.psiabtilde(variables.v, variables.s)) == variables.etaa(variables.mu * variables.v) * variables.deltab(variables.r, -variables.s)
variables.brackets(variables.alphaab(variables.mu, variables.m), variables.alphaab(variables.v, variables.n)) == variables.brackets(variables.alphaabtilde(variables.mu, variables.m), variables.alphaabtilde(variables.v, variables.n)) == variables.m * variables.etaa(variables.mu * variables.v) * variables.deltab(variables.m, -variables.n)
variables.Tbf(variables.F, variables.z) == summation((variables.Gb(variables.r))/(variables.za(variables.r + 3/2)), (variables.r, variables.ZB + variables.v))
variables.Tbftilde(variables.F, variables.zl) == summation((variables.Gbtilde(variables.r))/(variables.za(variables.m + 3/2)), (variables.r, variables.ZB + variables.vtilde))
variables.Tbf(variables.B, variables.z) == summation((variables.Lb(variables.m))/(variables.za(variables.m + 2)), (variables.m, -math.inf, math.inf)) 
variables.Tbftilde(variables.B, variables.zl) == summation((variables.Lbtilde(variables.m))/(variables.zla(variables.m + 2)), (variables.m, -math.inf, math.inf))
variables.brackets(variables.Lb(variables.m), variables.Lb(variables.n)) == (variables.m - variables.n) * variables.Lb(variables.m + variables.n) + variables.c/12 * (variables.ma(3) - variables.m) * variables.deltab(variables.m, -variables.n)
variables.bracketc(variables.Gb(variables.r), variables.Gb(variables.s)) == 2 * variables.Lb(variables.r + variables.s) + variables.c/12 * (4 * variables.ra(2) - 1) * variables.deltab(variables.r, -variables.s)
variables.brackets(variables.Lb(variables.m), variables.Gb(variables.r)) == (variables.m - 2 * variables.r)/2 * variables.Gb(variables.m + variables.r)
variables.Lb(variables.m) == 1/2 * summation(variables.point(variables.alphaab(variables.mu, variables.m-variables.n) * variables.alphab(variables.mu * variables.n)), (variables.n, variables.ZB)) + 1/4 * summation((2 * variables.r - variables.m) * variables.point(variables.psiab(variables.mu, variables.m-variables.r) * variables.psib(variables.mu * variables.r)), (variables.r, variables.ZB + variables.v)) + variables.aa(variables.m) * variables.deltab(variables.m, 0)
variables.Gb(variables.r) == summation(variables.alphaab(variables.mu, variables.n) * variables.psib(variables.mu * variables.r - variables.n), (variables.n, variables.ZB))
variables.aa(variables.m) == 1/16 * variables.D 
variables.aa(variables.m) == 0
variables.psiaf(variables.mu, (0, variables.sigmaa(2))) == math.exp(2 * math.pi * variables.i * variables.v) * variables.psiaftilde(variables.mu, (0, variables.sigmaa(2))) 
variables.psiaf(variables.mu, (math.pi, variables.sigmaa(2))) == math.exp(2 * math.pi * variables.i * variables.vp) * variables.psiaftilde(variables.mu, (variables.math.pi, variables.sigmaa(2)))
variables.psiaf(variables.mu, (variables.sigmaa(1), variables.sigmaa(2))) == variables.psiaftilde(variables.mu, (2 * math.pi - variables.sigmaa(1), variables.sigmaa(2)))
variables.psiab(variables.mu, variables.r) * variables.diracb(0, variables.NS) == 0
variables.r > 0
variables.Rhoa(variables.mu) == 2**(1/2) * variables.psiab(variables.mu, 0)
#check 2.8.18
math.exp(math.pi * variables.i * variables.F)
variables.Sigmaa(variables.mu * variables.lamda) == - variables.i/2 * summation(variables.brackets(variables.psiab(variables.mu, variables.r), variables.psiab(variables.lamda, - variables.r)), (variables.r, variables.ZB + variables.v))
variables.Sb(variables.a) == variables.ia(variables.deltab(variables.a, 0)) *variables.Sigmaa(2 * variables.a, 2 * variables.a + 1)
variables.F == summation(variables.Sb(variables.a), (variables.a, 0, 4))
variables.Sb(1) * (variables.psiab(2, variables.r) + variables.posneg * variables.psiab(3, variables.r)) == (variables.psiab(2, variables.r) + variables.posneg * variables.i * variables.psiab(3, variables.r)) * (variables.Sb(1) + variables.posneg * 1)
math.exp(math.pi * variables.i * variables.f) * variables.diracb(0, variables.NS) == - variables.diracb(0, variables.NS)
math.exp(math.pi * variables.i * variables.F) * variables.diracb(variables.sB, variables.R) == variables.diracb(variables.sBp, variables.R) * variables.Rhob(variables.sBp * variables.sB)
# ommitted 2.9.25

# lie algebra 2.29.1
variables.brackets(variables.Ta(variables.a), variables.Ta(variables.b)) == variables.i * variables.fab(variables.a * variables.b, variables.c) * variables.Ta(variables.c)
variables.brackets(variables.Ta(variables.a), variables.brackets(variables.Ta(variables.b), variables.Ta(variables.c))) + variables.brackets(variables.Ta(variables.b), variables.brackets(variables.Ta(variables.c), variables.Ta(variables.a))) + variables.brackets(variables.Ta(variables.c), variables.brackets(variables.Ta(variables.a), variables.Ta(variables.b))) == 0
math.exp(variables.i * variables.thetab(variables.a) * variables.Ta(variables.a))
(variables.Ta(variables.a), variables.Ta(variables.b)) == variables.da(variables.a * variables.b)
(variables.brackets(variables.T, variables.Tp), variables.Tpp) + (variables.Tp, variables.brackets(variables.T, variables.Tpp)) == 0
variables.Trf(variables.tab(variables.a, variables.r) * variables.tab(variables.b, variables.r)) == variables.Tb(variables.r) * variables.da(variables.a * variables.b) 
variables.tab(variables.a, variables.r) * variables.tab(variables.b, variables.r) * variables.db(variables.a * variables.b) == variables.Qb(variables.r)
variables.M * variables.T * variables.Ma(-1) == - variables.Ta(variables.T)
variables.M == variables.i * np.array([[0, variables.Ib(variables.k)],
                                      [-variables.Ib(variables.k), 0]])
variables.M * variables.U * variables.Ma(-1) == (variables.Ua(variables.T))**(-1)
variables.brackets(variables.Ha(variables.i), variables.Ea(variables.alpha)) == variables.SlopeRegge * variables.Ea(variables.alpha)
#check conditions for which 2.61.12
variables.brackets(variables.Ea(variables.alpha), variables.Ea(variables.Beta)) == variables.epsilonf(variables.alpha, variables.Beta) * variables.Ea(variables.alpha + variables.Beta)
variables.brackets(variables.Ea(variables.alpha), variables.Ea(variables.Beta)) == np.dot(2 * variables.alpha, variables.H/variables.alphaa(2))
variables.brackets(variables.Ea(variables.alpha), variables.Ea(variables.Beta)) == 0
# check combined vectors
np.array([[0, variables.i],
          [-variables.i, 0]])
(variables.posneg * 1, 0**(variables.k-1))
(+1, +1, 0**(variables.k-2)) 
(+1, -1, 0**(variables.k-2))
(-1, -1, 0**(variables.k-2))
(variables.posneg * 1, 0**(variables.k-1)) 
(variables.posneg * 2, 0**(variables.k-1))
(+1, -1, 0**(variables.n-2))
(+1/2, +1/2, +1/2, +1/2, +1/2, +1/2, +1/2, +1/2)
summation(variables.alphaa(variables.i), (variables.i, 1, 8)) == 2 * variables.ZB 
summation((variables.alphaa(variables.i))**2, (1, variables.i, 8))
-summation(variables.fab(variables.a * variables.c, variables.d) * variables.fab(variables.b * variables.d, variables.c), ((variables.c, variables.d))) == variables.hf(variables.g) * variables.psia(2) * variables.da(variables.a * variables.b)
variables.Eb(8) == np.cross(variables.SU(3), variables.Eb(6))
variables.B248 == (variables.B8, variables.B1) + (variables.B1, variables.B78) + (variables.B3, variables.B27) + (variables.B3l, variables.B27l)
variables.B27 == (variables.holderb((variables.B3, variables.B2), 1)) + variables.holderb((variables.B3l, variables.B1, -4)) + variables.holderb((variables.B1, variables.B1), 6) + variables.holderb((variables.B3l, variables.B1), 2) + variables.holderb((variables.B1, variables.B2), -3) + variables.brackets(variables.B1b(0)) + variables.brackets(variables.holderb((variables.B3l, variables.B1), 2) + variables.holderb((variables.B3, variables.B1), -2)) + variables.brackets(variables.holderb((variables.B1, variables.B2), -3) + variables.holderb((variables.B1, variables.B2), 3)) + variables.brackets(variables.B1b(0))
# current algebras 2.66.1 ommited, aprrox
variables.jaf(variables.a, variables.z) == summation((variables.jab(variables.a, variables.m))/((variables.za(variables.m + 1))), (variables.m, -math.inf, math.inf))
variables.brackets(variables.jab(variables.a, variables.m), variables.jab(variables.b, variables.n)) == variables.m * variables.ka(variables.a * variables.b) * variables.deltab(variables.m, -variables.n) + variables.i * variables. cab(variables.a * variables.b, variables.c) * variables.jab(variables.c, variables.m + variables.n)
variables.brackets(variables.jab(variables.a, 0), variables.jab(variables.b, variables.n)) == variables.i * variables.cab(variables.a * variables.b, variables.c) * variables.jab(variables.c, 0)
variables.cab(variables.a * variables.b, variables.c) == variables.fab(variables.a * variables.b, variables.c) 
variables.fab(variables.b * variables.c, variables.d) * variables.ka(variables.a * variables.d) + variables.fab(variables.b * variables.a, variables.d) * variables.ka(variables.d * variables.c) == 0
variables.ka(variables.a * variables.b) == variables.hat(variables.k) * variables.da(variables.a * variables.b)
variables.brackets(variables.Lb(variables.m), variables.jab(variables.a, variables.n)) == -variables.n * variables.jab(variables.a, variables.m + variables.n)
variables.hat(variables.k) * variables.da(variables.a * variables.a) == variables.bracket(1, variables.brackets(variables.jab(variables.a, 1), variables.jab(variables.a, -1)), 1) == np.abs(variables.dirac(variables.jab(variables.a, -1), 1))**2
variables.ja(3) == (np.dot(variables.alpha, variables.H))/(variables.alphaa(2)) 
variables.Ja(variables.posneg) == variables.Ea(variables.posneg * variables.alpha)
variables.brackets(variables.Ja(3), variables.Ja(variables.posneg)) == variables.posneg * variables.Ja(variables.posneg)
variables.brackets(variables.Ja(variables.pos), variables.Ja(variables.neg)) == 2 * variables.Ja(3)
(np.dot(variables.alpha, variables.Hb(0)))/(variables.alphaa(2))
variables.Eab(variables.alpha, 0) 
variables.Eab(-variables.alpha, 0)
(np.dot(variables.alpha, variables.Hb(0)) + variables.hat(variables.k))/(variables.alphaa(2)) 
variables.Eab(variables.alpha, 1) 
variables.Eab(-variables.alpha, -1)
variables.k == (2 * variables.hat(variables.k))/(variables.psia(2))
#2.68.14, OPE current approx.
variables.ja(variables.a) == variables.i * diff(variables.Ha(variables.a), variables.z)
variables.i * variables.lamdaa(variables.A) * variables.lamdaa(variables.B)
1/2 * variables.lamdaa(variables.A) * variables.lamdaa(variables.B) * variables.tab(variables.a, (variables.r, variables.A * variables.B))
variables.point(variables.j * variables.jf(variables.zb(1))) == sp.limit((variables.jaf(variables.a, variables.zb(1)) * variables.jaf(variables.a, variables.zb(2)) - (variables.hat(variables.k) - variables.dimf(variables.g))/(variables.zab(2, 12))), variables.zb(2), variables.zb(1))
# check product expansion 2.69.19 & 70.24, with holomorphic terms, and antisymetry approximation .20 & 22
variables.Tabf(variables.s, variables.B, variables.z) == 1/((variables.k + variables.hf(variables.g)) * variables.psia(2)) * variables.point(variables.j * variables.jf(variables.z))
variables.ca(variables.g, variables.k) == (variables.k * variables.dimf(variables.g))/(variables.k + variables.hf(variables.g))
variables.Lab(variables.s, 0) == 1/((variables.k + variables.hf(variables.g)) * variables.psia(2)) * (variables.jab(variables.a, 0) * variables.jab(variables.a, 0) + 2 * summation(variables.jab(variables.a, -variables.n) * variables.jab(variables.a, variables.n), (variables.n, -math.inf, math.inf)))
variables.Lab(variables.s, variables.m) == 1/((variables.k + variables.hf(variables.g)) * variables.psia(2)) * summation(variables.jab(variables.a, variables.n) * variables.jab(variables.a, variables.m-variables.n), (variables.n, -math.inf, math.inf))
variables.m != 0
variables.Tab(variables.s, variables.B) == 1/2 * variables.point(variables.j * variables.j)
variables.Tbp(variables.B) == variables.Tb(variables.B) - variables.Tab(variables.s, variables.B)
#check ope aprrox. 2.71.29 & 30
variables.cp == variables.c - variables.ca(variables.g, variables.k)
variables.ca(variables.g, variables.k) <= variables.c 
variables.ca(variables.g, variables.k) == variables.c 
variables.ca(variables.g, variables.k) == (variables.k * variables.dimf(variables.g) * variables.rank(variables.g))/(variables.dimf(variables.g) + (variables.k - 1) * variables.rank(variables.g))
variables.ca(variables.g, 1) == variables.rank(variables.g)
variables.ca(variables.g, variables.k) == (3 * variables.k)/(2 + variables.k)
variables.rank(variables.g) <= variables.ca(variables.g, variables.k) <= variables.dimf(variables.g)
variables.jab(variables.a, 0) * variables.dirac(variables.r, variables.i) == variables.dirac(variables.r, variables.j) * variables.tab(variables.a, (variables.r, variables.j * variables.i))
variables.Lab(variables.s, 0) * variables.dirac(variables.r, variables.i) == 1/((variables.k + variables.hf(variables.g)) * variables.psia(2)) * variables.dirac(variables.r, variables.k) * variables.tab(variables.a, (variables.r, variables.k * variables.j)) * variables.tab(variables.a, (variables.r, variables.j * variables.i)) == (variables.Qb(variables.r))/((variables.k + variables.hf(variables.g)) * variables.psia(2)) * variables.dirac(variables.r, variables.i)
variables.hb(variables.r) == (variables.Qb(variables.r))/((variables.k + variables.hf(variables.g)) * variables.psia(2)) == (variables.Qb(variables.r))/(2 * variables.hat(variables.k) + variables.Qb(variables.g))
variables.hb(variables.r) == (variables.Qb(variables.r))/((variables.k + variables.hf(variables.g)) * variables.psia(2)) == (variables.Qb(variables.r))/(2 * variables.hat(variables.k) + variables.Qb(variables.g)) 
variables.hb(variables.j) == (variables.j * (variables.j + 1))/(variables.k + 2)
variables.bracket(variables.dirac(variables.r, variables.lamda), variables.brackets(variables.Eab(variables.alpha, 1), variables.Eab(-variables.alpha, -1)), variables.dirac(variables.r, variables.lamda)) == 2 * variables.bracket(variables.dirac(variables.r, variables.lamda), (np.dot(variables.alpha, variables.Hb(0)) + variables.hat(variables.k)), variables.dirac(variables.r, variables.lamda))/variables.alphaa(2) == 2 * (np.dot(variables.alpha, variables.lamda) + variables.hat(variables.k))/variables.alphaa(2)
variables.hat(variables.k) >= Abs(np.dot(variables.alpha, variables.lamda))
variables.k >= (2 * Abs(np.dot(variables.psi, variables.lamda)))/(variables.psia(2)) == 2 * Abs(variables.Ja(3))
diff(variables.xa(variables.mu), variables.tau) * variables.lamdaa(variables.a) * variables.ea(np.dot(variables.i * variables.k, variables.X)) 
variables.lamdaaf(variables.a, variables.yb(1)) * variables.lamdaaf(variables.b, variables.yb(2)) == variables.brackets(variables.thetaf(variables.yb(1) - variables.yb(2)) * variables.dab(variables.a * variables.b, variables.c) + variables.thetaf(variables.yb(2) - variables.yb(1)) * variables.dab(variables.b * variables.a, variables.c)) * variables.lamdaaf(variables.c, variables.yb(2))
variables.lb(variables.L, variables.R) == (variables.SlopeRegge/2)**(1/2) * variables.kb(variables.L, variables.R)
(variables.lab(variables.m, variables.L), variables.lab(variables.n, variables.R)) 
variables.d <= variables.m <= 25
variables.d <= variables.n <= 9
variables.lf(variables.lp) == np.dot(variables.lb(variables.L), variables.lbp(variables.L)) - np.dot(variables.lb(variables.R), variables.lbp(variables.R))
variables.lf(variables.l) == 2 * variables.ZB 
variables.Rho == variables.star(variables.Rho)
# check lattice summations, 2.74.5-6
variables.Rhob(variables.r) == variables.star(variables.Rhob(variables.g))
variables.Rhob(variables.w) == variables.star(variables.Rhob(variables.g))
#check sublattices 2.75.9
variables.Rho == summation(np.cross(variables.Rhob(variables.r), variables.Rhobtilde(variables.r)), (variables.r))
variables.Rho == variables.Lamda * variables.Rhob(0)
variables.Lamda == variables.Of(26 - variables.d, 10 - variables.d, variables.RB)
variables.Lamdab(1) * variables.Lamda * variables.Lamdab(2) * variables.Rhob(0) == variables.Lamda * variables.Rhob(0)
variables.Lamdab(1) == np.cross(variables.Of(26 - variables.d, variables.RB), variables.Of(10 - variables.d, variables.RB))
variables.Lamdab(2) == variables.Of(26 - variables.d, 10 - variables.d, variables.ZB)
(variables.Of(26 - variables.d, 10 - variables.d, variables.RB))/(np.cross(variables.Of(26 - variables.d, variables.RB), variables.Of(10 - variables.d, variables.RB), variables.Of(26 - variables.d, 10 - variables.d, variables.ZB)))
variables.lab(2, variables.L) == 2
variables.lb(variables.R) == 0
1/2 * variables.brackets((36 - 2 * variables.d) * (35 - 2 * variables.d) - (26 - variables.d) * (25 - variables.d) - (10 - variables.d) * (9 - variables.d)) == (26 - variables.d) * (10 - variables.d)
variables.kb(variables.L * variables.m) == (variables.nb(variables.m))/variables.R + (variables.wa(variables.n) * variables.R)/(variables.SlopeRegge) * (variables.Gb(variables.m * variables.n) + variables.Bb(variables.m * variables.n)) - variables.qa(variables.I) * variables.Aab(variables.I, variables.m) - (variables.wa(variables.n) * variables.R)/2 * variables.Aab(variables.I, variables.n) * variables.Aab(variables.I, variables.m)
variables.kab(variables.I, variables.L) == (variables.qa(variables.I) + variables.wa(variables.m) * variables.R * variables.Aab(variables.I, variables.m)) * (2/variables.SlopeRegge)**(1/2)
variables.kb(variables.R * variables.m) == (variables.nb(variables.m))/variables.R + (variables.wa(variables.n) * variables.R)/(variables.SlopeRegge) * (-variables.Gb(variables.m * variables.n) + variables.Bb(variables.m * variables.n)) - variables.qa(variables.I) * variables.Aab(variables.I, variables.m) - (variables.wa(variables.n) * variables.R)/2 * variables.Aab(variables.I, variables.n) * variables.Aab(variables.I, variables.m)
variables.R * variables.Aab(variables.I, 9) == variables.diag((1/2)**8, 0**8)
variables.Rp * variables.Aab(variables.I, 9) == variables.diag(1, 0**7, 1, 0**7)
variables.kb(variables.L, variables.R) == variables.ntilde/variables.R + variables.posneg * (2 * variables.m * variables.R)/variables.SlopeRegge
variables.kbp(variables.L, variables.R) == variables.nptilde/variables.Rp + variables.posneg * (2 * variables.mp * variables.Rp)/variables.SlopeRegge
variables.B8b(variables.v) == +1, 0**6, -1
variables.B8 == +1/2**4, -1/2**4
np.cross(variables.B8b(variables.v), variables.B8b(variables.v)) == +2, +1**12, 0**38, -1**12, -2
np.cross(variables.B8, variables.B8b(variables.v)) == (3/2)**4,  1/2**28, -1/2**28, -(3/2)**4
variables.bracketc(variables.Qb(variables.alpha), variables.dagger(variables.Qb(variables.Beta))) == 2 * variables.Pb(variables.mu) * variables.holderb(variables.Rhoa(variables.mu), variables.Rhoa(0), variables.alpha * variables.Beta) + 2 * variables.Pb(variables.R * variables.m) * variables.holderb(variables.Rhoab(variables.m) * variables.Rhoa(0), variables.alpha * variables.Beta)
2 * variables.holderb(variables.M + variables.kb(variables.R * variables.m) * variables.Rhoa(variables.m) * variables.Rhoa(0), variables.alpha * variables.Beta)
2 * (variables.M + variables.posneg * Abs(variables.kb(variables.R)))
variables.Ma(2) == variables.kab(2, variables.R) + 4 * (variables.Ntilde - 1/2)/variables.SlopeRegge # NS
variables.Ma(2) == variables.kab(2, variables.R) + 4 * variables.Ntilde/variables.SlopeRegge # R
variables.Ma(2) == variables.kab(2, variables.L) + 4 * (variables.N - 1)/variables.SlopeRegge 
variables.N == 1 + variables.SlopeRegge * (variables.kab(2, variables.R) - variables.kab(2, variables.L))/4 == 1 - variables.nb(variables.m) * variables.wa(variables.m) - variables.qa(variables.I) * variables.qa(variables.I)/2 
variables.bracketc(variables.Qb(variables.alpha), variables.dagger(variables.Qb(variables.Beta))) == 2 * variables.Pb(variables.M) * variables.holderb(variables.Rhoa(variables.M) * variables.Rhoa(0), variables.alpha * variables.Beta) - 2 * (variables.Delta * variables.xb(variables.m))/(2 * math.pi * variables.SlopeRegge) * variables.holderb(variables.Rhoa(variables.m) * variables.Rhoa(0))
1/(2 * math.pi * variables.SlopeRegge) * integrate(variables.B, variables.M) == 1/2 * integrate(variables.da(10) * variables.x * variables.jaf(variables.M * variables.N, variables.x) * variables.Bbf(variables.M * variables.N, variables.x))
variables.jaf(variables.M * variables.N, variables.x) == 1/(2 * math.pi * variables.SlopeRegge) * integrate(variables.d2sigma * (diff(variables.xa(variables.M), 1) * diff(variables.xa(variables.N), 2) - diff(variables.xa(variables.N), 1) * diff(variables.xa(variables.M), 2)) * variables.deltaaf(10, variables.x - variables.Xf(variables.sigma)), variables.M)
variables.Qa(variables.M) == integrate(variables.da(9) * variables.x * variables.ja(variables.M * 0)) == 1/(2 * math.pi * variables.SlopeRegge) * integrate(diff(variables.xa(variables.M)))
variables.bracketc(variables.Qb(variables.alpha), variables.dagger(variables.Qb(variables.Beta))) =2 * (variables.Pb(variables.M) - variables.Qb(variables.M)) * variables.holderb(variables.Rhoa(variables.M) * variables.Rhoa(0), variables.alpha * variables.Beta)

# 2.229.1 
variables.Lb(0) * variables.dirac(variables.h) == variables.h * variables.dirac(variables.h)
variables.Lb(variables.m) * variables.dirac(variables.h) == 0
variables.m > 0 
variables.brackets(variables.Lb(variables.m), variables.Lb(variables.n)) == (variables.m - variables.n) * variables.Lb(variables.m + variables.n) + (variables.c)/12 * (variables.ma(3) - variables.m) * variables.deltab(variables.m, -variables.n) 
variables.dirac(variables.phi) == variables.dirac(variables.phi) - variables.dirac(variables.i) * variables.diracbs(variables.i, variables.phi)
variables.bracket(variables.phi, variables.Lb(-variables.n) * variables.Lb(variables.n), variables.phi) > 0 
variables.MOa(1) == variables.bracket(variables.h, variables.Lb(1) * variables.Lb(-1), variables.h) == 2 * variables.h 
variables.MOa(2) == np.array([[variables.dirac(variables.h, variables.Lab(2, 1))],
                             [variables.dirac(variables.h, variables.Lb(2))]]) * np.array([[variables.Lab(2, -1) * variables.dirac(variables.h), variables.Lb(-2) * variables.dirac(variables.h)]])
variables.MOa(2) == (sp.Matrix([[8 * variables.ha(2) + 4 * variables.h, 6 * variables.h],
                               [6 * variables.h, 4 * variables.h + variables.c/2]]))
sp.Determinant(variables.MOa(2)) == 32 * variables.h * (variables.h - variables.hb(variables.pos)) * (variables.h - variables.hb(variables.neg)) 
16 * variables.hb(variables.posneg) == (5 - variables.c) + variables.posneg * variables.brackets((1 - variables.c) * (25 - variables.c))**(1/2) 
variables.MOabf(variables.N, (variables.bracketc(variables.k), variables.bracketc(variables.kp)), (variables.c, variables.h))
summation(variables.kb(variables.i), variables.i) == variables.N 
sp.Determinant(variables.brackets(variables.MOaf(variables.N, (variables.c, variables.h)))) == variables.Kb(variables.N) * np.prod(variables.holdera((variables.h - variables.hb((variables.r, variables.s))), (variables.P * (variables.N - variables.r * variables.s))), (variables.r * variables.s, 1, variables.N))
variables.hb((variables.r, variables.s)) == (variables.c - 1)/24 + 1/4 * (variables.r * variables.alphab(variables.pos) + variables.s & variables.alphab(variables.neg))**2 
variables.alphab(variables.posneg) == (24)**(-1/2) * variables.brackets((1 - variables.c)**(1/2) + variables.posneg * (25  - variables.c) **(1/2))
np.prod(1/(1 - variables.qa(variables.n)), (variables.n, 1, math.inf)) == summation(variables.Pf(variables.k) * variables.qa(variables.k), (variables.k, 0, math.inf))
variables.h == variables.hb(variables.r, variables.s) == (variables.brackets(variables.r * (variables.m + 1) - variables.s * variables.ma)**2 - 1)/(4 * variables.m * (variables.m + 1))
variables.hb(variables.r, variables.s) == (25 - (3 * variables.r - 2 * variables.s)**2)/24 
variables.h == variables.hb(variables.r, variables.s) + variables.r * variables.s == (25 - (3 * variables.r - 2 * variables.s)**2)/24 
# conformal bootstrap 2.233.2
variables.Tf(variables.z) * variables.OObf(1, variables.zb(1)) == summation(np.dot((variables.z - variables.zb(1))**(variables.k - 2) * variables.Lb(-variables.k), variables.OObf(1, variables.zb(1))), (variables.k, -math.inf, math.inf))
variables.OObf(variables.m, (variables.z, variables.zl)) * variables.OObf(variables.n, (0, 0)) == summation(np.cross(variables.za(-variables.hb(variables.m) - variables.hb(variables.n) + variables.hb(variables.i) + variables.N) * variables.zla(-variables.hbtilde(variables.m) - variables.hbtilde(variables.n) + variables.hbtilde(variables.i) + variables.Ntilde), np.dot(variables.cab(variables.i * variables.bracketc(variables.k, variables.ktilde), variables.m * variables.n) * variables.Lb(-variables.bracketc(variables.k)) * variables.Lbtilde(-variables.bracketc(variables.ktilde)), variables.OObf(variables.i, (0, 0)))))
# reference evaluation abilties on 2.234.7
variables.cab(variables.i * variables.bracketc(variables.k, variables.ktilde), variables.m * variables.n) == variables.Betaab(variables.i * variables.bracketc(variables.k), variables.m * variables.n) * variables.Betaabtilde(variables.i * variables.bracketc(variables.ktilde), variables.m * variables.n) * variables.cab(variables.i, variables.m * variables.n)
variables.bracketb(variables.OObf(variables.j, (math.inf, math.inf)) * variables.OObf(variables.l, (1, 1)) * variables.OObf(variables.m, (variables.z, variables.zl)) * variables.OObf(variables.n, (0, 0)), variables.Sb(2)) == summation(variables.cab(variables.i, variables.j * variables.l) * variables.cb(variables.i * variables.m * variables.n) * variables.FOabf(variables.j * variables.l, variables.m * variables.n, (variables.i, variables.z)) * variables.FOabftilde(variables.j * variables.l, variables.m * variables.n, (variables.i, variables.zl)), variables.i)
variables.FOabf(variables.j * variables.l, variables.m * variables.n, (variables.i, variables.z)) == summation(variables.za(-variables.hb(variables.m) - variables.hb(variables.n) + variables.hb(variables.i) + variables.N) * variables.Betaab(variables.i * variables.bracketc(variables.k), variables.j * variables.l) * variables.MOb((variables.bracketc(variables.k), variables.bracketc(variables.kp))) * variables.Betaab(variables.i * variables.bracketc(variables.kp), variables.m * variables.n), ((variables.bracketc(variables.k), variables.bracketc(variables.kp))))
variables.Zf(variables.tau) == summation(variables.qa(-variables.c/24 + variables.hb(variables.i) + variables.N) * variables.qla(-variables.ctilde/24 + variables.hbtilde(variables.i) + variables.Ntilde), ((variables.i, variables.bracketc(variables.k, variables.kp)))) == summation(variables.chibf((variables.c, variables.hb(variables.i)), variables.q) * variables.chibf((variables.ctilde, variables.hlb(variables.i)), variables.ql), variables.i) 
variables.chibf((variables.c, variables.h), variables.q) == variables.qa(-variables.c/24 + variables.h) * summation(variables.qa(variables.N), variables.bracketc(variables.k))
variables.chibf((variables.c, variables.h), variables.q) == variables.qa(-variables.c/24 + variables.h) * np.prod((1)/(1 - variables.qa(variables.n)), (variables.n, 1, math.inf)) 
# checl approximations 2.236.14-15
variables.Zf(variables.i * variables.lO) <= variables.NO * variables.lO * math.exp(math.pi/(6 * variables.lO))
# minimal models 2.236.1
variables.h == variables.hb(1, 1) == (variables.c - 1)/24 + ((variables.alphab(variables.pos) + 2 * variables.alphab(variables.neg))**2)/4 
variables.NOb(1, 2) == np.dot(variables.brackets(variables.Lb(- 2) - 3/(2 * (2  * variables.hb(1, 2) + 1)) * variables.Lab(2, -1)), variables.OOb(1, 2)) == 0
0 == variables.bracketb(variables.NObf((1,2 ), variables.zb(1)) * np.prod(variables.OObf(variables.i, variables.zb(variables.i)), (variables.i, 2, variables.n)), variables.Sb(2)) == variables.brackets(variables.LOb(-2) - 3 /(2 * (2 * variables.hb(1, 2) + 1)) * variables.LOab(2, -1)) * variables.AOb(variables.n) == variables.brackets(summation((variables.hb(variables.i))/((variables.zb(variables.i) - variables.zb(1))**2), (variables.i, 1, variables.n)) - summation(1/(variables.zb(variables.i) - variables.zb(1)) * (diff(1, variables.z))/(diff(variables.zb(variables.i), variables.z)) - 3/(2 * (2 * variables.hb(1, 2) + 1)) * (diff(1, variables.z, 2))/(diff(variables.zab(2, 1), variables.z)), (variables.i, 2, variables.n))) * variables.AOb(variables.n) 
variables.AOb(variables.n) == variables.bracketb(variables.OObf((1, 2), (variables.zb(1), variables.zlb(1))) * np.prod(variables.OObf(variables.i, (variables.zb(variables.i), variables.zb(variables.i))), (variables.i, 2, variables.n)), variables.Sb(2))
variables.brackets(summation((variables.hb(variables.i))/((variables.zb(variables.i) - variables.zb(1))**2), (variables.i, 2, 4)) - summation((variables.hb(1,2) - variables.hb(2) - variables.hb(3) - variables.hb(4) + 2 * (variables.hb(variables.i) + variables.hb(variables.j)))/((variables.zb(variables.i) - variables.zb(1)) * (variables.zb(variables.j) - variables.zb(1))), (2 <= variables.i < variables.j <= 4)) + summation(1/(variables.zb(variables.i) - variables.zb(1)) * (diff(1, variables.z))/(diff(variables.zb(1), variables.z)) - 3/(2 * (2 * variables.hb(1, 2) + 1)) * (diff(1, variables.z, 2))/(diff(variables.zab(2, 1), variables.z)), (variables.i, 2, 4))) * variables.AOb(4) == 0
variables.bracketb(variables.OObf((variables.pair(variables.r, variables.s), variables.pair(variables.rtilde, variables.stilde)), (variables.zb(1), variables.zlb(1))) * np.prod(variables.OObf(variables.i, (variables.zb(variables.i), variables.zlb(variables.i))), (variables.i, 2, 4)), variables.Sb(2)) == summation(summation(variables.ab(variables.i * variables.j) * variables.fbf(variables.i, variables.z) * variables.fbftilde(variables.j, variables.zl), (variables.j, 1, variables.rtilde * variables.stilde)), (variables.i, 1, variables.r * variables.s))
3/(2 * (2 * variables.hb(1, 2) + 1)) * variables.keta * (variables.keta - 1) + variables.keta - variables.hb(variables.i) == 0 
variables.hb(variables.posneg) == variables.hb(1, 2) + variables.hb(variables.i) + variables.kappab(variables.posneg)
variables.hb(variables.i) == (variables.c - 1)/24 + (variables.gammaa(2))/4
variables.hb(variables.posneg) == (variables.c - 1)/24 + ((variables.gamma + variables.posneg * variables.alphab(variables.neg))**2)/4 
variables.OOb(1, 2) * variables.OOb(variables.gamma) == variables.brackets(variables.OOb(variables.gamma + variables.alphab(variables.neg))) + variables.brackets(variables.OOb(variables.gamma - variables.alphab(variables.neg)))
variables.OOb(2, 1) * variables.OOb(variables.gamma) == variables.brackets(variables.OOb(variables.gamma + variables.alphab(variables.pos))) + variables.brackets(variables.gamma - variables.alphab(variables.pos)) 
variables.OOb(1, 2) * variables.OOb(variables.r, variables.s) == variables.brackets(variables.r, variables.s) == variables.brackets(variables.OObs(variables.r, variables.s + 1)) + variables.brackets(variables.OOb(variables.r, variables.s - 1)) 
variables.OOb(2, 1) * variables.OOb(variables.r, variables.s) == variables.brackets(variables.OOb(variables.r + 1, variables.s)) + variables.brackets(variables.OOb(variables.r - 1, variables.s))
variables.c == 1 - 6 * ((variables.p - variables.q)**2)/(variables.p * variables.q) 
variables.hb(variables.r, variables.s) == ((variables.r * variables.q - variables.s * variables.p)**2 - (variables.p - variables.q)**2)/(4 * variables.p * variables.q)
variables.hb(variables.p - variables.r, variables.q - variables.s) == variables.hb(variables.r, variables.s)
1 <= variables.r <= variables.p - 1 
1 <= variables.s <= variables.q - 1
variables.p == variables.m 
variables.q == variables.m + 1
variables.hb(1, 1) == 0 
variables.hb(2, 1) == 1/2 
variables.hb(1, 2) == 1/16
variables.OOb(variables.rb(1), variables.sb(1)) * variables.OOb(variables.rb(2), variables.sb(2)) == summation(variables.brackets(variables.OOb(variables.r, variables.s)))
# check states fpr r & s, 2.240.19b/c 
2 * math.pi * (variables.hb(variables.ip) - variables.hb(variables.i) - variables.h) == 2 * math.pi * variables.Qb(variables.i) 
variables.OOb(variables.p - 1, 1) * variables.OOb(variables.p - 1, 1) == variables.brackets(variables.OOb(1, 1))
variables.Qb(variables.r, variables.s) == (variables.p * (1 - variables.s) + variables.q * (1 - variables.r))/(2) % 1
variables.c == 1 - 24 * variables.alphaab(2, 0)
variables.T == -1/2 * diff(variables.phi, variables.z) * diff(variables.phi, variables.z) + 2**(1/20) * variables.i * variables.alphab(0) * diff(variables.phi, variables.z, 2)
variables.Vb(variables.alpha) == math.exp(2**(1/2) * variables.i * variables.SlopeRegge * variables.phi)
variables.alpha == variables.alphab(0) - variables.gamma/2 
(variables.c - 1)/24 + (variables.gammaa(2))/4
variables.bracket(variables.Vb(variables.alphab(1)) * variables.Vb(variables.alphab(2)) * variables.Vb(variables.alphab(3)) * variables.Vb(variables.alphab(4)))
summation(variables.alphab(variables.i), variables.i) == 2 * variables.alphab(0)
variables.Jb(variables.posneg) == math.exp(2 **(1/2) * variables.i * variables.alphab(variables.posneg) * variables.phi)
variables.Qb(variables.posneg) == integrate(variables.Jb(variables.posneg, variables.z, 2 * math.pi))# contour 2.242.32
variables.nb(variables.pos) == 1/2 * summation(variables.rb(variables.i) - 2, (variables.i))
variables.nb(variables.neg) == 1/2 * summation(variables.sb(variables.i) - 2, variables.i)
# current alegbras 
variables.Lb(variables.m) * variables.dirac(variables.r, variables.i) == variables.jab(variables.a, variables.m) * variables.dirac(variables.r, variables.i) 
variables.m >0 
variables.jab(variables.a, 0) * variables.dirac(variables.r, variables.i) == variables.dirac(variables.r, variables.j) * variables.tab(variables.a, (variables.r, variables.j * variables.i)) 
variables.Tf(variables.z) == 1/((variables.k + variables.hf(variables.g)) * variables.psia(2)) * variables.point(variables.j * variables.jf(variables.z))
variables.ca(variables.g, variables.k) == (variables.k * variables.dim(variables.g))/(variables.k + variables.hf(variables.g))
variables.hb(variables.r) == (variables.Qb(variables.r))/((variables.k + variables.hf(variables.g)) * variables.psia(2))
variables.JOab(variables.a, -variables.m) == -summation((variables.ta(variables.af(variables.i)))/((variables.zb(variables.i) - variables.zb(1))**variables.m), (variables.i, 2, variables.n))
variables.Lb(variables.m) == 1/((variables.k + variables.hf(variables.g)) * variables.psia(2)) * summation(variables.jab(variables.a, variables.n) * variables.jab(variables.a, variables.m - variables.n), (variables.n, - math.inf, math.inf))
variables.LOb(-1) - 2/((variables.k + variables.hf(variables.g)) * variables.psia(2)) * summation(variables.ta(variables.af(variables.i)) * variables.JOab(variables.a, -1), variables.a)
variables.jab(variables.pos, -1) 
variables.jab(3, 0) - variables.k/2 
variables.jab(variables.neg, 1) 
variables.dirac(variables.j, variables.m)
(variables.jab(variables.pos, -1))**(variables.k - 2 * variables.m + 1) * variables.dirac(variables.j, variables.m) == 0
0 == summation(variables.bracketsb((variables.ta(variables.pos * 2))**(variables.lb(2)), (variables.mb(2), variables.nb(2))) * variables.bracketsb((variables.ta(variables.pos * 3))**variables.lb(3), (variables.mb(3), variables.nb(3))) * variables.bracketb(variables.OOb(variables.jb(1), variables.mb(1)) * variables.OOb(variables.jb(2), variables.mb(2)) * variables.OOb(variables.jb(3), variables.mb(3)), variables.Sb(2)), ((variables.mb(2), variables.mb(3))))
variables.lb(2) + variables.lb(3) >= variables.k - 2 * variables.mb(1) + 1 
variables.bracketb(variables.OOb(variables.jb(1), variables.jb(1)) * variables.OOb(variables.jb(2), variables.mb(2)) * variables.OOb(variables.jb(3), variables.mb(3)), variables.Sb(2)) == 0
variables.jb(1) + variables.jb(2) +variables.jb(3) > variables.k 
np.cross(variables.brackets(variables.jb(1)), variables.brackets(variables.k/2)) == variables.brackets(variables.k/2 - variables.jb(1))
variables.Zf(variables.tau) == summation(variables.nb(variables.r * variables.rtilde) * variables.chibf(variables.r, variables.q) * variables.star(variables.chibf(variables.rtilde, variables.q)), ((variables.r, variables.rtilde)))
variables.chibf(variables.r, variables.qp) == summation(variables.Sb(variables.r * variables.rp) * variables.chibf(variables.rp, variables.q), variables.rp)
variables.dagger(variables.S) * variables.n * variables.S == variables.n 
variables.Sb(variables.j * variables.jp) == (2/(variables.k + 2))**(1/2) * math.sin((math.pi * (2 * variables.j + 1) * (2 * variables.jp + 1))/(variables.k + 1))
variables.nb(variables.j * variables.jtilde) == variables.deltab(variables.j * variables.jtilde) 
# check evaluated conditions 2.246.24
variables.S == 1/(2 * math.pi * variables.SlopeRegge) * integrate(variables.d2z * (variables.Gb(variables.m * variables.n) + variables.Bb(variables.m * variables.n) * diff(variables.xa(variables.m), variables.z) * diff(variables.xa(variables.n), variables.zl)))
variables.Hb(variables.m * variables.n * variables.p) == (variables.q)/(variables.ra(3)) * variables.epsilonb(variables.m * variables.n * variables.p)
variables.Rb(variables.m * variables.n) == 2/(variables.ra(2)) * variables.Gb(variables.m * variables.n)
variables.Betaab(variables.G< variables.m * variables.n) == variables.SlopeRegge * variables.Gb(variables.m * variables.n) * (2/(variables.ra(2)) - (variables.qa(2))/(2 * variables.ra(6))) 
variables.Betaa(variables.Phi) == 1/2 - (variables.SlopeRegge * variables.qa(2))/(4 * variables.ra(6))
variables.ra(2) == Abs(variables.q)/2 + variables.Of(variables.SlopeRegge)
variables.c == 6 * variables.Betaa(variables.Phi) == 3 - (6 * variables.SlopeRegge)/(variables.ra(2)) + variables.Of((variables.SlopeRegge**2)/(variables.ra(4)))
variables.c == 3 - 6/(variables.k + 2)
math.exp(variables.i/(2 * math.pi * variables.SlopeRegge) * integrate(variables.B, variables.M)) == math.exp(variables.i/(2 * math.pi * variables.SlopeRegge) * integrate(variables.H, variables.N))
1 == math.exp(variables.i/(2 * math.pi * variables.SlopeRegge) * integrate(variables.H, variables.Sb(3))) == math.exp((math.pi * variables.i * variables.q)/(variables.SlopeRegge))
variables.q == 2 * variables.SlopeRegge * variables.n 
variables.ra(2) == variables.SlopeRegge * Abs(variables.n)
variables.c == 3 - 7/(Abs(variables.n)) + variables.Of(1/variables.na(2)) 
variables.g == variables.xa(4) + variables.i * variables.xas(variables.i) * variables.sigmaa(variables.i) 
summation((variables.xas(variables.i))**2, (variables.i, 1, 4)) == 1
variables.S == (Abs(variables.n))/(4 * math.pi) * integrate(variables.d2z * variables.Trf(diff(variables.ga(-1), variables.z) * diff(variables.g, variables.zl)), variables.M) + (variables.i * variables.n)/(12 * math.pi) * integrate(variables.Trf(variables.omegaas(3)), variables.N)
diff(variables.omegaas(3)) == 0
(variables.n)/(12 * math.pi) * integrate(variables.Trf(variables.chi), variables.M) 
variables.delta * variables.S == (Abs(variables.n))/(2* math.pi) * integrate(variables.d2z * variables.Trf(variables.brackets(diff(variables.g, variables.zl) * variables.ga(-1) * diff(variables.ga(-1) * variables.delta * variables.g, variables.z)))) == (Abs(variables.n))/(2 * math.pi) * integrate(variables.d2z * variables.Trf(variables.brackets(variables.ga(-1) * diff(variables.g, variables.z) * diff(variables.delta * variables.g * variables.ga(-1), variables.zl))))
variables.delta * variables.gf((variables.z, variables.zl)) == variables.i * variables.epsilonb(variables.L) * variables.gf((variables.z, variables.zl)) - variables.i * variables.gf(variables.z, variables.zl) * variables.epsilonb(variables.R)
variables.delta * variables.gf((variables.z, variables.zl)) == variables.i * variables.epsilonbf(variables.L, variables.z) * variables.gf((variables.z, variables.zl)) - variables.i * variables.gf((variables.z, variables.zl)) * variables.epsilonbf(variables.R, variables.zl)
Abs(variables.n) * variables.Trf(variables.epsilonb(variables.R) * variables.ga(-1) * diff(variables.g)) 
Abs(variables.n) * variables.Trf(variables.epsilonb(variables.L) * diff(variables.g, variables.zl) * variables.ga(-1))
variables.LO == 1/(4 * math.pi) * diff(variables.phia(variables.a), variables.z) * diff(variables.phia(variables.a), variables.zl) + variables.Of(variables.phia(3))
variables.jab(variables.a, variables.R) == Abs(variables.n)**(1/2) * diff(variables.phia(variables.a), variables.z) + variables.Of(variables.phia(2))
variables.jab(variables.a, variables.L) == Abs(variables.n)**(1/2) * diff(variables.phia(variables.a), variables.zl) + variables.Of(variables.phia(2))
variables.Dabf(variables.r, variables.i * variables.j, variables.g) == variables.OOabf(variables.r, variables.i, variables.z) * variables.OOabftilde(variables.r, variables.j, variables.zl)
variables.Ta(variables.G) == variables.Ta(variables.H) + variables.Ta(variables.G/variables.H) 
variables.ca(variables.G/variables.H) == variables.ca(variables.G/variables.H) == variables.ca(variables.G) - variables.ca(variables.H)
variables.G == variables.SumDirect(variables.S * variables.Ufb(2, variables.k), variables.S * variables.Ufb(2, 1)) 
variables.ca(variables.G) == 4 - 6/(variables.k + 2)
variables.H == variables.S * variables.Ufb(2, variables.k + 1) 
variables.ca(variables.H) == 3 - 6/(variables.k + 3)
variables.ca(variables.G/variables.H) == 1 - 6/((variables.k + 2) * (variables.k + 3))
variables.chiabf(variables.G, variables.r, variables.q) == summation(variables.nab(variables.r, variables.rp * variables.rpp) * variables.chiabf(variables.H, variables.rp, variables.q) * variables.chiabf(variables.G/variables.H, variables.rpp, variables.q), ((variables.rp, variables.rpp)))
variables.Sb((variables.r * variables.s, variables.rp * variables.sp)) == variables.brackets(8/((variables.q + 1)* (variables.q + 1)) * (-1)**((variables.r + variables.s) * (variables.rp + variables.sp)))**(1/2) * math.sin((math.pi * variables.r * variables.rp)/variables.p) * math.sin((math.pi * variables.s * variables.sp)/ variables.q)
(variables.S * variables.Ufb(2, variables.k))/(variables.Uf(1))
variables.c == 2 - 6/(variables.k + 2)
variables.ja(3) == variables.i * (variables.k/2)**(1/2) * diff(variables.H, variables.z)
variables.Ta(variables.H) == -1/2 * diff(variables.H, variables.z) * diff(variables.H, variables.z)
variables.ja(variables.pos) == math.exp(variables.brackets(variables.i * variables.Hf(2/variables.k)**(1/2))) * variables.psib(1) 
variables.ja(variables.neg) == math.exp(variables.brackets(-variables.i * variables.Hf(2/variables.k)**(1/2))) * variables.dagger(variables.psib(1))
variables.point((variables.ja(variables.pos))**(variables.l)) == np.dot((variables.jab(variables.pos, -1))**variables.l, variables.l) == math.exp(variables.i * variables.l * variables.Hf(2/variables.k)**(1/2)) * variables.psib(variables.l)
variables.OOb((variables.j, variables.m)) == math.exp(variables.brackets(variables.i * variables.m * variables.Hf(2/(variables.k))**(1/2))) * variables.psiab(variables.j, variables.m)
variables.point(variables.jab(variables.a, (1)) * variables.jab(variables.a, (1)))
variables.point(variables.jab(variables.a, (1)) * variables.jab(variables.a, (2)))
variables.point(variables.jab(variables.a, (2)) * variables.jab(variables.a, (2)))
variables.brackets(variables.jabf(variables.b, (1), variables.z) + variables.jabf(variables.b, (2), variables.z)) * variables.point(variables.jab(variables.a, (variables.i) * variables.jabf(variables.a, (variables.j), 0))) == summation(1/(variables.za(variables.k + 1)) * variables.brackets(variables.jab(variables.b, (variables.k, 1)) + variables.jab(variables.b, (variables.k, 2))) * np.dot(variables.jab(variables.a, (-1, variables.i)) * variables.jab(variables.a, (-1, variables.j)), 1), (variables.k, 0, math.inf))
variables.G == variables.SumDirect(variables.S * variables.Ufb(variables.n, variables.kb(1)), variables.S * variables.Ufb(variables.n, variables.kb(2)))
variables.ca(variables.G) == (variables.na(2) - 1) * variables.brackets((variables.kb(1))/(variables.kb(1) + variables.n) + (variables.kb(2))/(variables.kb(2) + variables.n))
variables.H == variables.S * variables.Ufb(variables.n, variables.kb(1) + variables.kb(2))
variables.ca(variables.H) == (variables.na(2) - 1) * (variables.kb(1) + variables.kb(2))/(variables.kb(1) + variables.kb(2) + variables.n)
variables.da(variables.a * variables.b * variables.c) == variables.Trf(variables.ta(variables.a) * variables.bracketc(variables.ta(variables.b), variables.ta(variables.c)))
variables.c == 2 - (24)/((variables.k + 3) * (variables.k + 4))
#2.254.1
sp.Determinant(variables.holderb(variables.MOa(variables.N)), variables.R, variables.NS) == (variables.h - (variables.epsilon * variables.hat(variables.c))/16) * variables.Kb(variables.N) * np.prod((variables.h - variables.hb(variables.r, variables.s))**(variables.Pb(variables.R, variables.NS) * (variables.N - variables.r * variables.s)/2), (1 <= variables.r * variables.s <= 2 * variables.N))
variables.hb(variables.r, variables.s) == (variables.hat(variables.c) - 1 + variables.epsilon)/(16) + 1/4 * (variables.r * variables.hat(variables.alphab(variables.pos)) + variables.s * variables.hat(variables.alphab(variables.neg)))**2
variables.hat(variables.alphab(variables.posneg)) == 1/4 * variables.brackets((1 - variables.hat(variables.c))**(1/2) + variables.posneg * (9 - variables.hat(variables.c))**(1/2))
np.prod((1 + variables.qa(variables.n - 1))/(1 - variables.qa(variables.n)), (variables.n, 1, math.inf)) == summation(variables.Pbf(variables.R, variables.k) * variables.qa(variables.k), (variables.k, 0, math.inf))
np.prod((1 + variables.qa(variables.n - 1/2))/(1 - variables.qa(variables.n)), (variables.n, 1, math.inf)) == summation(variables.Pbf(variables.NS, variables.k) * variables.qa(variables.k), (variables.k, 0, math.inf))
variables.hat(variables.c) >= 1 
variables.h >= variables.epsilon * (variables.hat(variables.c))/16
variables.h == variables.hb(variables.r, variables.s) == (variables.brackets(variables.r * (variables.m + 2) - variables.s * variables.m)**2 - 4)/(8 * variables.m * (variables.m + 2)) + (variables.epsilon)/16
variables.G == variables.SumDirect(variables.S * variables.Ufb(2, variables.k), variables.S * variables.Ufb(2, 2)) 
variables.H == variables.S * variables.Ufb(2, variables.k + 2)
#2.255.1
variables.Zf(variables.tau) == summation(variables.nb(variables.i * variables.jtilde) * variables.chibf(variables.i, variables.q) * variables.star(variables.chibf(variables.jtilde, variables.q)), ((variables.i, variables.jtilde)))
variables.Nb(variables.i, variables.j * variables.k * variables.l) == variables.Nab(variables.r, (variables.i, variables.j)) * variables.Nb(variables.r * variables.k * variables.l)
variables.Nab(variables.r, variables.i * variables.j) * variables.Nb(variables.r * variables.k * variables.l) == variables.Nab(variables.r, variables.i * variables.k) * variables.Nb(variables.r, variables.j * variables.l) == variables.Nab(variables.r, variables.i * variables.l) * variables.Nb(variables.r, variables.j * variables.k)
variables.Nb(variables.i, variables.j * variables.k * variables.l) * (variables.hb(variables.i) + variables.hb(variables.j) + variables.hb(variables.k) + variables.hb(variables.l)) - summation((variables.Nab(variables.r, (variables.i, variables.j)) * variables.Nb(variables.r * variables.k * variables.l) + variables.NNab(variables.r, variables.i * variables.k) * variables.Nb(variables.r, variables.j * variables.l) + variables.Nab(variables.r, variables.i * variables.l) * variables.Nb(variables.r, variables.j * variables.k)) * variables.hb(variables.r), (variables.r)) == variables.ZB 
summation(variables.Nab(variables.r, variables.i * variables.i) * variables.Nb(variables.r * variables.i * variables.i) * (4 * variables.hb(variables.i) - 3 * variables.hb(variables.r)), variables.r) == variables.ZB 
variables.Nab(0, (0, 0)) == variables.Nab(0, (1/2, 1/2)) == variables.Nab(1, (1/2, 1/2)) == variables.Nab(0, (1, 1)) == variables.Nab(1, (1, 1)) == variables.Nab(0, (3/2, 3/2)) == 1
8 * variables.hb(1/2) - 3 * variables.hb(1) 
5 * variables.hb(1) 
4 * variables.hb(3/2) 
variables.Sa(4) == (variables.S * variables.T)**3 == 1
1 == variables.brackets((sp.Determinant(variables.S))**4)**(-3) * variables.brackets((sp.Determinant(variables.S) * sp.Determinant(variables.T))**(3))**(4) == (sp.Determinant(variables.T))**(12)
# T: (2.258.12)
variables.chibf(variables.i, variables.q) == math.exp(variables.brackets(2 * math.pi * variables.i * (variables.hb(variables.i) - variables.c/24))) * variables.chibf(variables.i, variables.q)
(variables.NO * variables.C)/2 - 12 * summation(variables.hb(variables.i), variables.i) == variables.ZB 
variables.Nab(variables.i, variables.j * variables.k) == summation((variables.Sab(variables.r, variables.j) * variables.Sab(variables.r, variables.k) * variables.dagger(variables.Sb(variables.r))**variables.i)/(variables.Sab(variables.r, 0)), (variables.r)) 
variables.Tp == variables.Lb(variables.a * variables.b) * variables.point(variables.ja(variables.a) * variables.ja(variables.b))
variables.Lb(variables.a * variables.b) == 2 * variables.Lb(variables.a * variables.c) * variables.ka(variables.c * variables.d) * variables.Lb(variables.d * variables.b) - variables.Lb(variables.c * variables.d) * variables.Lb(variables.e * variables.f) * variables.fab(variables.c * variables.e, variables.a) * variables.fab(variables.d * variables.f, variables.b) - variables.Lb(variables.c * variables.d) * variables.fab(variables.c * variables.e, variables.f) * variables.fab(variables.d * variables.f, (variables.a, variables.b)) * variables.Lb(variables.b, variables.e) # check 2.258.16 
variables.c == 2* variables.ka(variables.a * variables.b) * variables.Lb(variables.a * variables.b)

# conformal algebra check <333 2.376.2
variables.Tbf(variables.B, variables.z) == summation((variables.Lb(variables.n))/(variables.za(variables.n + 2)), (variables.n, variables.ZB))
variables.jf == summation((variables.Jb(variables.n))/(variables.za(variables.n + 1)), (variables.n, variables.ZB))
variables.Tabf(variables.pos, variables.F, variables.z) == summation((variables.Gab(variables.pos, variables.r))/(variables.za(variables.r + 3/2)), (variables.r, variables.ZB + variables.v))
variables.Tabf(variables.neg, variables.F, variables.z) == summation((variables.Gab(variables.neg, variables.r))/(variables.za(variables.r + 3/2)), (variables.r, variables.ZB - variables.v))
variables.brackets(variables.Lb(variables.m), variables.Gab(variables.posneg, variables.r)) == (variables.m/2 - variables.r) * variables.Gab(variables.posneg, variables.m + variables.r)
variables.brackets(variables.Lb(variables.m), variables.Jb(variables.n)) == - variables.n * variables.Jb(variables.m + variables.n)
variables.bracketc(variables.Gab(variables.pos, variables.r), variables.Gab(variables.neg, variables.s)) == 2 * variables.Lb(variables.r + variables.s) + (variables.r - variables.s) * variables.jb(variables.r + variables.s) + (variables.c)/3 * (variables.ra(2) - 1/4) * variables.deltab(variables.r, - variables.s)
variables.bracketc(variables.Gab(variables.pos, variables.r), variables.Gab(variables.pos, variables.s)) == variables.bracketc(variables.Gab(variables.neg, variables.r), variables.Gab(variables.neg, variables.s)) == 0
variables.brackets(variables.Jb(variables.n), variables.Gab(variables.posneg, variables.r)) == variables.posneg * variables.Gab(variables.posneg, variables.r + variables.n)
variables.brackets(variables.Jb(variables.m), variables.Jb(variables.n)) == variables.c/3 * variables.m * variables.deltab(variables.m, -variables.n)
variables.Tbtilde(variables.F) == variables.Tabtilde(variables.pos, variables.F) + variables.Tabtilde(variables.neg, variables.F)
math.exp(variables.brackets(variables.l * variables.Phitilde + variables.i * variables.sb(0) * variables.Hatilde(0) + variables.i * variables.sb(1) * variables.Hatilde(1) + variables.i * variables.Qtilde * (variables.Htilde/(3**(1/2)))))
variables.l + variables.sb(0) + variables.sb(1) + variables.Qtilde == 2 * variables.ZB 
variables.uO * math.exp(variables.i * variables.Htilde/(3**(1/2)))
variables.uOl * math.exp(-variables.i * variables.Htilde/(3**(1/2)))
variables.ja(variables.a) * math.exp(3**(1/2) * variables.i * variables.Htilde /2)
variables.uO * math.exp(variables.brackets(variables.i * variables.Htilde/(np.cross(2, 3**(1/2)))))
variables.uOl * math.exp(variables.brackets(variables.i * variables.Htilde/(np.cross(2, 3**(1/2)))))
variables.ja(variables.a) * math.exp(-3**(1/2) * variables.i * variables.Htilde/2)
variables.bracketc(variables.Gabtilde(variables.pos, 1/2), variables.Gabtilde(variables.neg, -1/2)) == 2 * variables.Lbtilde(0) + variables.Jbtilde(0)
variables.bracketc(variables.Gabtilde(variables.pos, -1/2), variables.Gabtilde(variables.neg, 1/2)) == 2 * variables.Lbtilde(0) + variables.Jbtilde(0)
variables.jtilde == variables.i * (variables.ctilde/3)**(1/2) * diff(variables.Htilde, variables.zl)
2  * variables.htilde >= Abs(variables.Qtilde) 
variables.Gabtilde(variables.posneg, variables.r) * variables.dirac(variables.c) == 0
variables.r > 0
variables.Lbtilde(variables.n) * variables.dirac(variables.c) == variables.Jbtilde(variables.n) * variables.dirac(variables.c) == 0 
variables.n > 0
variables.Gabtilde(variables.pos, -1/2) * variables.dirac(variables.c) == 0
(3 * variables.Qatilde(2))/(2 * variables.ctilde) <= (Abs(variables.Qtilde))/(2) # becomes
Abs(variables.Qtilde) <= (variables.ctilde)/3
np.dot(variables.Jbtilde(1), variables.VOa(0)) == np.dot(variables.Jbtilde(1) * variables.Gb(-1/2), variables.VOa(-1)) == np.dot((variables.Gbtilde(-1/2) * variables.Jbtilde(1) + variables.Gabtilde(variables.pos, 1/2) - variables.Gabtilde(variables.neg, 1/2)), variables.VOa(-1)) == 0
math.exp(variables.brackets(variables.i * (3/variables.ctilde)**(1/2) * variables.Qtilde * variables.Htilde)) == math.exp(variables.brackets(variables.i * (3/variables.ctilde)**(1/2) * variables.Qtilde * variables.Htilde - variables.i * variables.eta * (variables.ctilde/3)**(1/2) * variables.Htilde)) 
variables.v == variables.v + variables.eta 
variables.Gabtilde(variables.posneg, variables.n) * variables.dirac(variables.psi) == variables.Lbtilde(variables.n) * variables.dirac(variables.psi) == variables.Jbtilde(variables.n) * variables.dirac(variables.psi) == 0
variables.n >= 0 
