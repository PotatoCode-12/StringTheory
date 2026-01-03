import math
from scipy.integrate import quad
import numpy as np 
from sympy import symbols, diff,integrate, Function, Derivative, summation, oo, ln, Abs
import sympy as sp
import variables
import constants

#orbifolds 256.1
variables.xa(25) == -variables.xa(25)
variables.xa(variables.m) == -variables.xa(variables.m)
26 - variables.k <= variables.m <= 25
#group:
variables.ta(variables.m)
variables.xa(25) == variables.xa(25) + 2 * math.pi * variables.R * variables.m 
#group:
variables.ta(variables.m) * variables.r 
variables.xa(25)  == 2 * math.pi * variables.R * variables.m -variables.xa(25) 
#256.4
variables.xa(25, variables.sigmaa(1) + 2 * math.pi) == -variables.xa(25, variables.sigmaa(1))
variables.dirac((variables.N, variables.Ntilde), (variables.ka(variables.mu), variables.n, variables.w)) == (-1)**(summation)(variables.Nab(25, variables.m) + variables.Nabtilde(25, variables.m), (variables.m, 1, math.inf)) * variables.dirac((variables.N, variables.Ntilde), (variables.ka(variables.mu), -variables.n, -variables.w))
variables.xa(25, (variables.z, variables.zl)) == variables.i * (variables.SlopeRegge/2)**(1/2) * variables.summation((1)/(variables.m + 1/2) * (((variables.alphaab(25, variables.m + 1/2))/(variables.za(variables.m + 1/2))) + ((variables.alphaabtilde(25, variables.m+1/2))/(variables.zla(variables.m + 1/2)))), (variables.m, -math.inf, math.inf))
variables.xa(25, (variables.sigmaa(1) + 2 * math.pi)) == 2 * math.pi * variables.R - variables.xa(25, variables.sigmaa(1))
variables.ma(2) == 4/(variables.SlopeRegge) * (variables.N - 15/16) 
variables.N == variables.Ntilde
(variables.q * variables.ql)**(-1/24) * variables.Trbf(variables.U, (1 + variables.r)/2 * variables.qa(variables.Lb(0)) * variables.qla(variables.Lbtilde(0)))
# check partition function, 258.10
(variables.q * variables.ql)**(1/48) * variables.Trbf(variables.T, (1 + variables.r)/2 * variables.qla(variables.Lb(0) * variables.qla(variables.Lbtilde(0)))) == (variables.q * variables.qla)**(1/48) * variables.brackets(np.prod(Abs(1 - variables.qa(variables.m - 1/2))**(-2), (variables.m, 1, math.inf)) + np.prod(Abs(1 + variables.qa(variables.m - 1/2))**(-2), (variables.m, 1, math.inf)))
variables.Zbf(variables.orb, (variables.R, variables.tau)) == 1/2 * variables.Zbf(variables.tor, (variables.R, variables.tau)) + Abs((variables.etaf(variables.tau))/(variables.thetavarb((1, 0), 0, variables.tau))) + Abs((variables.etaf(variables.tau))/(variables.thetavarb((0, 1), 0, variables.tau))) + Abs((variables.etaf(variables.tau))/(variables.thetavarb((0, 0), 0, variables.tau)))
# path integral 259.13
variables.xa(25, (variables.sigmaa(1) + 2 * math.pi, variables.sigmaa(2)))= (-1)**(variables.a + 1) * variables.xa(25, (variables.sigmaa(1), variables.sigmaa(2)))
variables.xa(25, (variables.sigmaa(1) + 2 * math.pi * variables.taub(1), variables.sigmaa(2) + 2 * math.pi * variables.taub(2))) == (-1)**(variables.b + 1) * variables.xa(25, (variables.sigmaa(1), variables.sigmaa(2)))
#twisting
variables.phi * (variables.sigmaa(1) + 2 * math.pi) == np.cross(variables.h, variables.phif(variables.sigmaa(1)))
variables.Pb(variables.H) == 1/(variables.order(variables.H)) * summation(variables.hat(variables.hb(2)), (variables.hb(2), variables.H))
variables.Z == 1/(variables.order(variables.H)) * summation(variables.Zb(variables.hb(1), variables.hb(2)), ((variables.hb(1), variables.hb(2)), variables.H))
variables.phifp(variables.sigmaa(1)) == np.cross(variables.hb(2), variables.phif(variables.sigmaa(1)))
variables.phifp(variables.sigmaa(1) + 2 * math.pi) == np.cross(variables.hbp(1), variables.phifp(variables.sigmaa(1)))
variables.hbp(1) == variables.hb(2) * variables.hb(1) * variables.hab(-1, 2)
# CFTs 261.20
variables.rp 
variables.xa(25) == variables.xa(25 + math.pi * variables.SlopeRegge**(1/2))
variables.Zbf(variables.orb, (variables.SlopeRegge**(1/2), variables.tau)) == variables.Zbf(variables.tor, (2 * variables.SlopeRegge**(1/2), variables.tau))
variables.xa(25) == variables.xa(25) + (2 * math.pi * variables.SlopeRegge**(1/2))/variables.k 
# check surviving massless scalars 263.23

#2.275.1
variables.xa(variables.m) == variables.thetaa(variables.m * variables.n) * variables.xa(variables.n) + variables.va(variables.m) 
variables.psiatilde(variables.m) == variables.thetaa(variables.m * variables.n) * variables.psiatilde(variables.n) 
variables.lamdaa(variables.A) == variables.gamma(variables.A * variables.B) * variables.lamdaa(variables.A)
variables.K == (variables.Ra(6))/(variables.S)
variables.Ta(6) == (variables.Ra(6))/(variables.Lamda)
np.dot((variables.theta, variables.w), (1, variables.v), (variables.theta, variables.w)**(-1)) == (1, variables.theta * variables.v)
variables.Pl == variables.S/variables.Lamda
variables.K == (variables.Ta(6))/(variables.Pl)
variables.gammaf(variables.thetab(1), variables.vb(1)) * variables.gammaf(np.dot((variables.thetab(1), variables.vb(1)), (variables.thetab(2), variables.vb(2))))
#realized phi issue
variables.phif(variables.sigmaa(1) + 2 * math.pi) == np.dot(variables.h, variables.phif(variables.sigmaa(1)))
variables.theta == math.exp(variables.brackets(2 * math.pi * variables.i * (variables.Phib(2) * variables.Jb(45) + variables.Phib(3) * variables.Jb(67) + variables.Phib(4) * variables.Jb(89))))
variables.Za(variables.i) == 2**(-1/2) * (variables.xa(2 * variables.i) + variables.i * variables.xa(2 * variables.i + 1))
variables.i == 2, 3, 4
variables.Za(variables.il) == variables.Zla(variables.i) == 2**(-1/2) * (variables.xa(2 * variables.i) - variables.i * variables.xa(2 * variables.i + 1))
variables.Zaf(variables.i, variables.sigma + 2 * math.pi) == math.exp(2 * math.pi * variables.i * variables.Phib(variables.i)) * variables.Zaf(variables.i, variables.sigma)
variables.psiaf(variables.i, variables.sigma + 2 * math.pi) == math.exp(variables.brackets(2 * math.pi * variables.i * (variables.Phib(variables.i) + variables.v))) * variables.psiaftilde(variables.i, variables.sigma)
variables.n + variables.Phib(variables.i) 
variables.n - variables.Phib(variables.i) 
variables.n - variables.Phib(variables.i) 
variables.n + variables.Phib(variables.i) 
variables.n - variables.Phib(variables.i)
variables.n - variables.Phib(variables.i) + 1/2 
variables.n + variables.Phib(variables.i) 
variables.n + variables.Phib(variables.i) + 1/2
variables.lamdaa(variables.posneg * variables.K) == math.exp(variables.posneg * 2 * math.pi * variables.i * variables.Betab(variables.K)) * variables.lamdaa(variables.posneg * variables.K)
variables.Phib(variables.i) == (variables.rb(variables.i))/(variables.N) 
variables.Betab(variables.K) == (variables.Sb(variables.K))/(variables.N) 
1/2 * summation(variables.Phib(variables.i), (variables.i, 2, 4)) 
1/2 * summation(variables.Betab(variables.K), (variables.K, 1, 8))
1/2 * summation(variables.Betab(variables.K), (variables.K< 9, 16))
summation(variables.rb(variables.i), (variables.i, 2, 4)) == summation(variables.Sb(variables.K), (variables.K, 1, 8)) == summation(variables.sb(variables.K), (variables.K, 9, 16)) == 0 % 2
1/24 * 1/8 * (2 * variables.theta - 1)**2
variables.Lb(0) - variables.Lbtilde(0) == - summation((variables.Na(variables.i) + variables.Natilde(variables.i) + variables.Nabtilde(variables.i, variables.psi)) * variables.Phib(variables.i), (variables.i, 2, 4)) - summation(variables.Na(variables.K) * variables.Betab(variables.K), (variables.K, 1, 16)) - 1/2 * summation(variables.Phib(variables.i) * (1 - variables.Phib(variables.i)), (variables.i, 2, 4)) + 1/2 * summation(variables.Betab(variables.K) * (1 - variables.Betab(variables.K)) % 1, (variables.K, 1, 16))
- 1/2 * summation(variables.Phib(variables.i) * (1 - variables.Phib(variables.i)), (variables.i, 2, 4)) + 1/2 * summation(variables.Betab(variables.K) * (1 - variables.Betab(variables.K)), (variables.K, 1, 16)) == (variables.m)/(variables.N)
summation((variables.Na(variables.i) + variables.Natilde(variables.i) + variables.Nabtilde(variables.i, variables.psi)) * variables.Phib(variables.i), (variables.i, 2, 4)) + summation(variables.Na(variables.K) * variables.Betab(variables.K), (variables.K, 1, 16)) == (variables.m)/(variables.N) % 1 
variables.delta == 1/2 * (-variables.Nbtilde(variables.psi) - 1 + summation(variables.Phib(variables.i), (variables.i, 2, 4)))
variables.delta == variables.k/variables.N 
summation(variables.rab(2, variables.i), (variables.i, 2, 4)) - summation(variables.sab(2, variables.K), (variables.K, 1, 16)) == 0 % (2 * variables.N)
variables.Tbtilde(variables.F) == variables.i * summation(variables.chibtilde(variables.I) * variables.chibtilde(variables.J) * variables.chibtilde(variables.K) * variables.cb(variables.I * variables.J *variables.K), ((variables.I, variables.J, variables.K), 1, 18))
variables.cb(variables.I * variables.J * variables.M) * variables.cb(variables.K * variables.L * variables.M) + variables.cb(variables.J * variables.K * variables.M) * variables.cb(variables.I * variables.L * variables.M) + variables.cb(variables.K * variables.I *variables.M) * variables.cb(variables.J * variables.L * variables.M) == 0
18 * variables.cb(variables.I * variables.K * variables.L) * variables.cb(variables.J * variables.K * variables.L) == variables.deltab(variables.I * variables.J)
variables.ea(np.dot(variables.i * variables.k, variables.xb(variables.R))) 
variables.ka(2) == 6/(variables.SlopeRegge)
variables.ea(np.dot(variables.i * variables.l, variables.xb(variables.R))) * diff(variables.xb(variables.R), variables.zl)
variables.la(2) == 2/(variables.SlopeRegge)
