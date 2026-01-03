import math
from scipy.integrate import quad
import numpy as np 
from sympy import symbols, diff,integrate, Function, Derivative, summation, oo
import variables
import actions
import constants

diff_resulta = variables.xa(variables.mu)
for _ in range(int(variables.a)):
    diff_resulta = diff(diff_resulta, variables.mu)
a = diff_resulta

diff_resultb = variables.xb(variables.mu)
for _ in range(int(variables.b)):
    diff_resultb = diff(diff_resultb, variables.mu)
b = diff_resultb

diff_resultc = variables.xb(variables.mu)
for _ in range(int(variables.c)):
    diff_resultc = diff(diff_resultc, variables.mu)
c = diff_resultc
def TensionRegge():
    T=1/(2*math.pi*variables.SlopeRegge)
(variables.TF(variables.tau,variables.sigma))**(variables.a*variables.b) =-4*math.pi*(-variables.GWorldSheetMetric)**(-1/2)*variables.NambuGotoVar/(variables.NambuGotoVar*variables.WorldSheetMetric(variables.a*variables.b))*actions.AppNambuGoto(0) == -1/variables.SlopeRegge*a*b-1/2*variables.GWorldSheetMetric**(variables.a*variables.b)*diff(variables.xa(variables.mu),variables.c)*c #check value of AppNambuGoto(0)

variables.Tb(variables.z, variables.zl) == 0 # energy momentum tensor state pg 43
diff(variables.Tb(variables.z, variables.z), variables.zl) == diff(variables.Tb(variables.zl, variables.zl), variables.z) == 0
variables.Tf(variables.z) == variables.Tbf(variables.z, variables.z, variables.z)
variables.Tftilde(variables.zl) == variables.Tbf(variables.zl, variables.zl, variables.zl)
variables.Tf(variables.z) =-1/(variables.SlopeRegge)*variables.point(diff(variables.xa(variables.mu), variables.z)*diff(variables.xb(variables.mu), variables.z))
variables.Tftilde(variables.zl) == -1/(variables.SlopeRegge)*variables.point(diff(variables.xa(variables.mu), variables.zl)*diff(variables.xb(variables.mu), variables.zl))

variables.Tf(variables.z) == 1/(-variables.SlopeRegge) * variables.point(diff(variables.xa(variables.mu), variables.z) * diff(variables.xb(variables.mu), variables.z)) + variables.Vb(variables.mu) * diff(variables.xa(variables.mu), variables.z, 2)
variables.Tftilde(variables.zl) == -1/variables.SlopeRegge * variables.point(diff(variables.xa(variables.mu), variables.zl) * diff(variables.xb(variables.mu), variables.zl)) + variables.Vb(variables.mu) * diff(variables.xa(variables.mu), variables.zl, 2) # for energy mometum tensor pg 49
constants.c == variables.ctilde == variables.D + 6*variables.SlopeRegge*variables.Vb(variables.mu)*variables.Va(variables.mu) # central charge
variables.Tf(variables.z) == variables.point(diff(variables.b, variables.z) * variables.c) - variables.lamda*diff(variables.point(variables.b*variables.c), variables.z) # noethers theorom energy-momentum tensor
variables.Tftilde(variables.zl) == 0
variables.c == -3(2*variables.lamda-1)**2 +1 # for TT OPE pg 51
variables.ctilde == 0

variables.hb(variables.Beta) == variables.lamda # check third CFT family
variables.hb(variables.gamma) == 1 - variables.lamda 
variables.S == 1/(2*math.pi) * integrate(variables.d2z * variables.Beta * diff(variables.gamma, variables.zl))
diff(variables.gammaf(variables.z), variables.zl) == diff(variables.Betaf(variables.z), variables.zl) == 0
variables.T == variables.point(diff(variables.Beta, variables.z) * variables.gamma) - variables.lamda*diff(variables.point(variables.Beta*variables.gamma))
variables.Ttilde == 0
variables.c == 3(2*variables.lamda-1)**2 -1
variables.ctilde == 0
