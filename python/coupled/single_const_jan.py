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


def krsoilrootfunction2(Hx, Hsoil, lroot, rroot, B, k_root, zs, sp):

    # fig, (ax1) = plt.subplots(1, 1)

    hintmin = -8.912409381337460e+04
    hintmax = -2

    k_soilfun = lambda hsoil, hint: (vg.fast_mfp[sp](hsoil) - vg.fast_mfp[sp](hint)) / (hsoil - hint)

    ktemp = Hx * 0
    hin = Hx * 0
    print(len(Hx))
    for i in range(0, len(Hx)):

        if (Hsoil[i] - zs[i]) < hintmax:

            if (Hsoil[i] - zs[i]) > hintmin:

                fun = lambda x: (rroot[i] * k_root[i] * Hx[i] + B[i] * Hsoil[i] * k_soilfun(Hsoil[i], x)) / (B[i] * k_soilfun(Hsoil[i], x) + rroot[i] * k_root[i]) - x
                # fun = lambda x: rroot[i] * k_root[i] * (x - Hx[i]) + B[i] * k_soilfun(Hsoil[i], x) * (x - Hsoil[i])

                ktemp[i] = 2 * np.pi * lroot[i] * rroot[i] * k_root[i] * B[i] * k_soilfun(Hsoil[i] - zs[i], hin[i] - zs[i]) / (B[i] * k_soilfun(Hsoil[i] - zs[i], hin[i] - zs[i]) + rroot[i] * k_root[i])

#                 if i % 10 == 0:
#                     h_ = np.linspace(-16000, 0, 1000)
#                     h2_ = [fun(h) for h in h_]
#                     ax1.plot(h_, h2_, label=str(i) + ", root " + str(Hx[i]))

#                 r = optimize.least_squares(fun, Hx[i])
#                 hin[i] = r.x

                hin[i] = fsolve(fun, Hx[i])

            else:
                hin[i] = Hx[i]
                ktemp[i] = 2 * np.pi * lroot[i] * rroot[i] * k_root[i] * B[i] * k_soilfun(hintmin, Hx[i] - zs[i]) / (B[i] * k_soilfun(hintmin, Hx[i] - zs[i]) + rroot[i] * k_root[i])
        else:
            hin[i] = Hsoil[i]
            ktemp[i] = 2 * np.pi * lroot[i] * rroot[i] * k_root[i]

#     ax1.legend()
#     plt.show()
    return ktemp, hin


def plotme(Hx, Hsoil, lroot, rroot, B, k_root, zs, sp):
    fig, (ax1) = plt.subplots(1, 1)

    hintmin = -8.912409381337460e+04
    hintmax = -2

    k_soilfun = lambda hsoil, hint: (vg.fast_mfp[sp](hsoil) - vg.fast_mfp[sp](hint)) / (hsoil - hint)

    ktemp = Hx * 0
    hin = Hx * 0
    for i in range(0, len(Hx)):

        if (Hsoil[i] - zs[i]) < hintmax:

            if (Hsoil[i] - zs[i]) > hintmin:
                fun = lambda x: (rroot[i] * k_root[i] * Hx[i] + B[i] * Hsoil[i] * k_soilfun(Hsoil[i], x)) / (B[i] * k_soilfun(Hsoil[i], x) + rroot[i] * k_root[i]) - x
                # fun = lambda x: rroot[i] * k_root[i] * (x - Hx[i]) + B[i] * k_soilfun(Hsoil[i], x) * (x - Hsoil[i])
                if i % 10 == 0:
                    h_ = np.linspace(-16000, 0, 1000)
                    h2_ = [fun(h) for h in h_]
                    ax1.plot([-16000, 0], [0., 0.])
#                 r = optimize.least_squares(fun, Hx[i])
#                 x = r.x
                    x = fsolve(fun, Hx[i])
                    ax1.plot(h_, h2_, label = "bulk {:g}, root {:g}, interface {:g} (cm)".format(Hsoil[i], Hx[i], x[0]))
                    ax1.plot(x, fun(x), "*")

    plt.xlabel("bulk soil [cm]")
    plt.ylabel("fun(bulk soil) [cm]")
    plt.legend()
    plt.show()


def soil_root_inerface(rx, sx, r, sp):

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

                rho = outer_r[i] / a
                b = 2 * (rho * rho - 1) / (2 * rho * rho * (np.log(rho) - 0.5) + 1);
                # print(rho, a, outer_r)
                # b = outer_r[i] - a
                b = 0.649559423771715

                # fun = lambda x: (a * kr * rx[i] + b * sx[i] * k_soilfun(sx[i], x)) / (b * k_soilfun(sx[i], x) + a * kr) - x
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
kr = 1.81e-04  # units? [1/day]
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

psiinterface1 = soil_root_inerface(psixlin[1:], psis1, r, sp)
# _, psiinterface1 = krsoilrootfunction2(psixlin[1:], psis1, lroot, rroot, b, k_root, z, sp)  # Hx, Hsoil, lroot, rroot, B, k_root, zs, sp

# psiinterface1 = soil_root_inerface(psixlin[1:], psis1, r, sp)

# plotme(psixlin[1:], psis1, lroot, rroot, b, k_root, z, sp)

fig, (ax1) = plt.subplots(1, 1)
ax1.set_title("Figure 4")
ax1.plot(psis1, psiinterface1, label = "1st iteration")

psiinterface = psiinterface1.copy()
for i in range(0, 5):
    psixlin = r.solve_dirichlet(0., collar_p, 0, psiinterface, False)
    # _, psiinterface = krsoilrootfunction2(psixlin[1:], psis1, lroot, rroot, b, k_root, z, sp)  # Hx, Hsoil, lroot, rroot, B, k_root, zs, sp
    psiinterface = soil_root_inerface(psixlin[1:], psis1, r, sp)
    ax1.plot(psis1, psiinterface, label = str(i))

# plotme(psixlin[1:], psis1, lroot, rroot, b, k_root, z, sp)

psixlin = r.solve_dirichlet(0., collar_p, 0, psiinterface, False)
sink2 = np.array(r.segFluxes(0., psixlin, psiinterface, True, False))

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

