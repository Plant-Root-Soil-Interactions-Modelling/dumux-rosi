import sys; sys.path.append("../../modules/"); sys.path.append("../../../../CPlantBox/");  sys.path.append("../../../../CPlantBox/src/python_modules")

import van_genuchten as vg
from sympy import *

init_printing()

wilting_point = 15000

sand = [0.045, 0.43, 0.15, 3, 1000]
loam = [0.08, 0.43, 0.04, 1.6, 50]
clay = [0.1, 0.4, 0.01, 1.1, 10]
sp = vg.Parameters(loam)

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
phi_integral = Integral(hydraulic_conductivity, (h, hh, oo))

print("water_content             ", water_content)
print("water_content             ", simplify(water_content))
print("effective_saturation      ", effective_saturation)
print("effective_saturation      ", simplify(effective_saturation))
print("hydraulic_conductivity    ", hydraulic_conductivity)
print("hydraulic_conductivity    ", simplify(hydraulic_conductivity))
print("phi                       ", phi_integral)
# print("matric_flux_potential     ", simplify(matric_flux_potential))  # generally too hard?
# print("matric_flux_potential-200 ", matric_flux_potential.evalf(subs = {hh:-200.})) # generally too hard?

q_out = symbols('q_out')
# q_out = Integer(0)
r_in = symbols('r_in')
r_out = symbols('r_out')
rho = r_out / r_in
r = symbols('r')  
phi = symbols("phi")

f = (phi + q_out * r_out * ln(1 / rho)) * ((r ** 2 / r_in ** 2 - 1 + 2 * rho ** 2 * ln(r_in / r)) / (rho ** 2 - 1 + 2 * rho ** 2 * ln(1 / rho))) + q_out * r_out * ln(r / r_in)

dfdr = simplify(diff(f, r))
print()
print("f                       ", f)
print("f                       ", simplify(f))

print()
print("df/dr                   ", simplify(dfdr))
der = (phi + q_out * r_out * ln(1 / rho)) / (rho ** 2 - 1 + 2 * rho ** 2 * ln(1 / rho)) * (2 * r / r_in ** 2 - 2 * rho ** 2 / r) + q_out * r_out / r
print("der                     ", simplify(der))
dfdr_r_in = simplify(dfdr.subs(r, r_in))
print("df/dr, r=r_in          ", simplify(dfdr_r_in))
dfdr0 = simplify(dfdr_r_in.subs(q_out, Integer(0)))
print("df/dr, q_out=0          ", simplify(dfdr0))

print()
h_ = 500  # = -500 cm  
phi_integral_ = phi_integral.subs({hh: h_, theta_R: sp.theta_R, theta_S: sp.theta_S, alpha: sp.alpha, n: sp.n, ksat: sp.Ksat})
# print(phi_integral_)
phi_ = phi_integral_.evalf()
print("phi(h_)                 ", phi_)

print()
dfdr_ = dfdr0.evalf(subs={phi:phi_, r_in: 0.02, r_out: 0.6})
print("df/dr                   ", dfdr_)
k_ = hydraulic_conductivity.subs({h: 15000, theta_R: sp.theta_R, theta_S: sp.theta_S, alpha: sp.alpha, n: sp.n, ksat: sp.Ksat})
print("k_                      ", k_)
print("dh/dr                   ", dfdr_ / k_)  # isn't this value by k_ in order to get the matric flux potential?

# print(latex(water_content))
