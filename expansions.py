import math
from scipy.integrate import quad
import numpy as np 
from sympy import symbols, diff,integrate, Function, Derivative, summation, oo, ln
import variables
import constants

diff(variables.xa(variables.mu, variables.z), variables.z) == -variables.i*(variables.SlopeRegge/2)**(1/2) * summation((variables.alphaab(variables.mu, variables.m))/(variables.za(variables.m+1)), (variables.m, -math.inf, math.inf)) # Laurent expansions for free scalar pg 58.1
diff(variables.xa(variables.mu, variables.zl), variables.zl) == -variables.i*(variables.SlopeRegge/2)**(1/2) * summation((variables.alphaabtilde(variables.mu, variables.m))/(variables.zla(variables.m+1)))
variables.alphaab(variables.mu, variables.m) == (2/variables.SlopeRegge)**(1/2)*integrate(diff(variables.z)/(2*math.pi) * variables.za(variables.m)*diff(variables.xa(variables.mu, variables.z), variables.z), 0, 2*math.pi) #contour boundsless
variables.alphaabtilde(variables.mu, variables.m) == -(2/variables.SlopeRegge)**(1/2) * integrate((diff(variables.zl))/(2*math.pi) * variables.zla(variables.m)) * diff(variables.xa(variables.mu, variables.zl), variables.zl, 0, 2*math.pi) # 58.2 a,b

#laurent 61.16
variables.bf(variables.z) == summation((variables.bb(variables.m))/(variables.za(variables.m+variables.lamda)), (variables.m, -math.inf, math.inf))
variables.cf(variables.z) == summation((variables.cb(variables.m))/(variables.za(variables.m+1-variables.lamda)), (variables.m, -math.inf, math.inf))
variables.bracketc(variables.bb(variables.m), variables.cb(variables.n)) == variables.deltab(variables.m, -variables.n) # anticommunicators

#for unit circle
variables.xb(variables.b, variables.theta) == summation(variables.xb(variables.n) * math.e**(variables.i*variables.n*variables.theta), variables.n, -math.inf, math.inf) # 67.19
variables.xb(variables.n)*variables.trans == variables.xb(-variables.n)
