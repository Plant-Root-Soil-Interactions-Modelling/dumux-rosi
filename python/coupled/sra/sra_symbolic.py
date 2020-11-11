import sys; sys.path.append("../../modules/"); sys.path.append("../../../../CPlantBox/");  sys.path.append("../../../build-cmake/cpp/python_binding/")

import van_genuchten as vg
from sympy import *

init_printing()

wilting_point = 15000

sand = [0.045, 0.43, 0.15, 3, 1000]
loam = [0.08, 0.43, 0.04, 1.6, 50]
clay = [0.1, 0.4, 0.01, 1.1, 10]
sp = vg.Parameters(sand)

h = symbols('h')  # soil matric potential h:=Abs(h)
hh = symbols('hh')  # soil matric potential as integration bound

theta_R = symbols('theta_R')
# theta_R = sp.theta_R
theta_S = symbols('theta_S')
# theta_S = sp.theta_S
alpha = symbols('alpha')
# alpha = sp.alpha
n = symbols('n')
# n = sp.n
m = 1 - 1 / n  # symbols('m') # self.m = 1. - 1. / self.n
ksat = symbols('K')
# ksat = sp.Ksat

water_content = theta_R + (theta_S - theta_R) / pow(1 + pow(alpha * h, n), m)
effective_saturation = (water_content - theta_R) / (theta_S - theta_R)
hydraulic_conductivity = ksat * (effective_saturation ** Rational(1, 2)) * ((1 - pow(1 - pow(effective_saturation, 1 / m), m)) ** 2)
phi = Integral(hydraulic_conductivity, (h, hh, oo))

print("water_content             ", water_content)
print("water_content             ", simplify(water_content))
print("effective_saturation      ", effective_saturation)
print("effective_saturation      ", simplify(effective_saturation))
print("hydraulic_conductivity    ", hydraulic_conductivity)
print("hydraulic_conductivity    ", simplify(hydraulic_conductivity))
print("phi                       ", phi)
# print("matric_flux_potential     ", simplify(matric_flux_potential))  # generally too hard?
# print("matric_flux_potential-200 ", matric_flux_potential.evalf(subs = {hh:-200.})) # generally too hard?

# q_out = symbols('q_out')
q_out = Integer(0)
r_in = symbols('r_in')
r_out = symbols('r_out')
rho = r_out / r_in
r = symbols('r')  # r = r_in # -> mfp_ = 0 -> h = -15000
phi = symbols("phi")
f = (phi + q_out * r_out * ln(1 / rho)) * ((r ** 2 / r_in ** 2 - 1 + 2 * rho ** 2 * ln(r_in / r)) / (rho ** 2 - 1 + 2 * rho ** 2 * ln(1 / rho))) + q_out * r_out * ln(r / r_in)

dfdr = simplify(diff(f, r))

print("f                       ", f)
print("f                       ", simplify(f))

print("df/dr                   ", dfdr)
print("df/dr                   ", simplify(dfdr.subs(r, r_in)))

# print("matric_potential_mfp   ", matric_potential_mfp)
# print("i try...")
# mfp = 6.6e-7
# matric_potential_mfp = solveset(Eq(matric_flux_potential - mfp, 0), hh, domain = S.Reals)
#
#

# print(latex(water_content))

