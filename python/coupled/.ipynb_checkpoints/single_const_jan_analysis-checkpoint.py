import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src");

from functional.xylem_flux import XylemFluxPython  # Python hybrid solver
import plantbox as pb
import functional.van_genuchten as vg

import numpy as np
from scipy import optimize
from scipy import integrate
from scipy.optimize import fsolve
import matplotlib.pyplot as plt


def soil_root_inerface(rx, sx, r, sp):
    assert rx.shape == sx.shape

    hintmin = -1.e5
    hintmax = -2.
    outer_r = r.rs.segOuterRadii()  #

    k_soilfun = lambda hsoil, hint: (vg.fast_mfp[sp](hsoil) - vg.fast_mfp[sp](hint)) / (hsoil - hint)

    rsx = rx * 0
    for i in range(0, len(rx)):

        s = r.rs.segments[i]
        z = 0.5 * (r.rs.nodes[s.x].z + r.rs.nodes[s.y].z)  # segment mid point

        if (sx[i] - z) < hintmax:

            if (sx[i] - z) > hintmin:

                a = r.rs.radii[i]
                kr = r.kr_f(0., 0, 0, 0)  #  kr_f = [](double age, int type, int orgtype, int numleaf)

                rho = outer_r[i] / a  # Eqn [5]
                rho2 = rho * rho
                # b = 2 * (rho2 - 1) / (1 + 2 * rho2 * (np.log(rho) - 0.5))  # Eqn [4]
                b = 2 * (rho2 - 1) / (1 - 0.53 * 0.53 * rho2 + 2 * rho2 * (np.log(rho) + np.log(0.53)))  # Eqn [7]
                print(b, a, outer_r[i], rho)

                b = 0.649559423771715

                fun = lambda x: (a * kr * rx[i] + b * sx[i] * k_soilfun(sx[i], x)) / (b * k_soilfun(sx[i], x) + a * kr) - x
                rsx[i] = fsolve(fun, rx[i])

            else:
                rsx[i] = rx[i]
        else:
            rsx[i] = sx[i]

    return rsx

""" 
    script, that should mimic matlab script single_const_krkx.m 
"""

""" soil parameter """

alpha = 0.04
n = 1.6
Ks = 50

# alpha = 0.018
# n = 1.8
# Ks = 28.46

soil = [0.08, .43, alpha, n, Ks]
sp = vg.Parameters(soil)  # soil parameters to use with vg functions
vg.create_mfp_lookup(sp, -1.e5, 1000)

psitop = -15000
psibottom = -50

hup = -np.logspace(0, 5, 101);
hup = hup[0:100]  # why?
fluxmp = np.array([vg.matric_flux_potential(h, sp) for h in hup])
print(fluxmp[0], fluxmp[60], fluxmp[90], fluxmp[99], hup[90])  # real low values do not agree with matlab

""" root parameter """
r_root = 0.05  # [cm]
kr = 1.81e-04  # units? []
kx = 0.1730  # units? []
collar_p = -16000  # cm

""" soil (per segment) """
zbottom = -50
z = np.linspace(0, -50, 101)  # (-0.25, -49.75, 100)  # mid values of segments
z = z[1:]
psis1 = np.ones(z.shape) * psibottom + (z - np.ones(z.shape) * zbottom) * (psitop - psibottom) / (-zbottom)

""" single root (straight 0 to -50, 100 segments) """
n, segs = [], []
for i in range(0, 101):  # nodes
    n.append(pb.Vector3d(0, 0, -.5 * i))
for i in range(0, len(n) - 1):  # connect nodes
    segs.append(pb.Vector2i(i, i + 1))
rs = pb.MappedSegments(n, segs, [r_root] * len(segs))  # a single root, with radius r_root
r = XylemFluxPython(rs)
r.setKr([kr])
r.setKx([kx])

psixlin = r.solve_dirichlet(0., collar_p, 0, psis1, False)
sink = np.array(r.segFluxes(0., psixlin, psis1, True, False))
suf = r.get_suf(0.)

# # print(psixlin)
# fig, (ax1) = plt.subplots(1, 1)
# ax1.plot(psixlin[1:], z)
# ax1.legend()
# plt.show()

lroot = np.ones(z.shape) * 0.5
rroot = np.ones(z.shape) * r_root
b = np.ones(z.shape) * 0.649559423771715
k_root = np.ones(z.shape) * kr

# _, psiinterface1 = krsoilrootfunction2(psixlin[1:], psis1, lroot, rroot, b, k_root, z, sp)  # Hx, Hsoil, lroot, rroot, B, k_root, zs, sp
psiinterface1 = soil_root_inerface(psixlin[1:], psis1, r, sp)

# plotme(psixlin[1:], psis1, lroot, rroot, b, k_root, z, sp)

psiinterface = psiinterface1.copy()
for i in range(0, 9):
    psixlin = r.solve_dirichlet(0., collar_p, 0, psiinterface, False)
    # _, psiinterface = krsoilrootfunction2(psixlin[1:], psis1, lroot, rroot, b, k_root, z, sp)  # Hx, Hsoil, lroot, rroot, B, k_root, zs, sp
    psiinterface1 = soil_root_inerface(psixlin[1:], psis1, r, sp)

# plotme(psixlin[1:], psis1, lroot, rroot, b, k_root, z, sp)

psixlin = r.solve_dirichlet(0., collar_p, 0, psiinterface, False)
sink2 = np.array(r.segFluxes(0., psixlin, psiinterface, True, False))

fig, (ax1) = plt.subplots(1, 1)
ax1.set_title("Figure 4")
ax1.plot(psis1, psiinterface1, label = "1st iteration")
ax1.plot(psis1, psiinterface, label = "10th iteration")
ax1.legend()
ax1.set_xlabel(r"$\psi$ soil")
ax1.set_ylabel(r"$\psi$ interface")
plt.show()

fig, (ax1) = plt.subplots(1, 1)
ax1.set_title("Figure 5")
ax1.plot(z, -4 * sink, label = "classic sink")  # (-4) instead of upsacling 0.5 cm line to 2 cm layer
ax1.plot(z, -4 * sink2, label = "nonlinear exact")
ax1.legend()
ax1.set_ylabel(r"depth (cm)")
ax1.set_xlabel(r"sink")
plt.show()

print("fin")

